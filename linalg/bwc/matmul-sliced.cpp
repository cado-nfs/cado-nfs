/* Manage the in-memory data for the matrix */
/* It's in C++ because the STL is handy, but that's really all there is
 * to it... */
#include "cado.h" // IWYU pragma: keep
#include <cstdio>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <cstdarg>         // for va_list, va_end, va_start

#include <utility>          // for pair, make_pair
#include <vector>
#include <algorithm>    // sort
#include <iostream>     // cout

#include "matmul.hpp"         // for matmul_ptr, matmul_public_s
#include "arith-hard.hpp"
#include "matmul-common.hpp"
#include "matmul_facade.hpp"
#include "macros.h"
#include "params.h"

using namespace std;

// assembly is now disabled for this code because the semantics of the
// asm code have changed.

// #define L1_CACHE_SIZE   32768
// take only 3/4 of the L1 cache.
#define L1_CACHE_SIZE   24576
// for 8-bytes arith_hard::elt values, this gives 3072 items.

/* Here is how the matrix is stored in memory, in the
 * "matmul_sliced_data_s" type.  The matrix is cut in slices, where each
 * slice is a set of contiguous rows. The size of a slice is tuned so as
 * to fill the L1 cache as per the macro above, with some adjustment to
 * handle non-exact divisibility: some slices will have packbase rows and
 * some others packbase+1 rows.  Within a slice, entries are stored as a
 * list of pairs
 *   (column index, row index),
 * sorted according to column index. Then, only the difference between
 * two consecutive column indices is actually stored, so that it will fit
 * in 16 bits. Also the row index fits in 16 bits, so that a slice is
 * actually a list of 16-bit unsigned integers. 
 * All the slices are stored in a big array of uint16_t, in the data
 * field. Within this data field, the information is organized like so:
 *   data[0] : total number of slices
 *   data[1] : number of rows in slice number 1
 *   data[2..3]: number of entries in slice number 1, say k_1
 *   data[4..(4+2*k_1)]: list of entries of slice number 1
 *   data[(4+2k_1+1)]: number of rows in slice number 2
 *   etc...
 * Additionally, for each slice, a slice_info structure is stored, giving
 * statistics about the slice, and the offset in the data array where it
 * is stored.
 */

/* This extension is used to distinguish between several possible
 * implementations of the product. The upper word correspond to the
 * implementation itself, the lower one to the n-th binary incompatible
 * change (make sure to bump it) */
#define MM_EXTENSION   "-sliced"

#define MM_MAGIC_FAMILY        0xa000UL
#define MM_MAGIC_VERSION       0x1007UL
#define MM_MAGIC (MM_MAGIC_FAMILY << 16 | MM_MAGIC_VERSION)

/* see matmul-basic.c */
#define MM_DIR0_PREFERS_TRANSP_MULT   1

struct slice_info {
    unsigned int nrows;
    unsigned int ncoeffs;
    double spent_time;
    unsigned int npasses;
    double dj_avg;
    double dj_sdev;
    double di_avg;
    double di_sdev;
    int32_t di_max;
    size_t data_offset;
};

using data_t = vector<uint16_t>;

template<typename Arith>
struct matmul_sliced : public matmul_interface {
    /* now our private fields */
    Arith * xab;
    data_t data;
    vector<slice_info> dslices_info;
    unsigned int npack;
    void push(uint32_t x)
    {
        ASSERT_ALWAYS(x >> 16 == 0);
        data.push_back(x);
    }
    void push32(uint64_t x)
    {
        ASSERT_ALWAYS(x >> 32 == 0);
        data.push_back(x & ((1U << 16) - 1));
        x >>= 16;
        data.push_back(x & ((1U << 16) - 1));
    }
    static uint32_t read32(data_t::const_iterator & q) {
        uint32_t res;
        res = *q++;
        res |= ((uint32_t) *q++) << 16;
        return res;
    }
    static uint32_t read32(const uint16_t * & q) {
        uint32_t res;
        res = *q++;
        res |= ((uint32_t) *q++) << 16;
        return res;
    }

    void build_cache(matrix_u32 &&) override;
    int reload_cache_private() override;
    void save_cache_private() override;
    void mul(void *, const void *, int) override;
    void report(double scale) override;

    matmul_sliced(matmul_public &&, arith_concrete_base *, cxx_param_list &, int);
    ~matmul_sliced() override = default;

    matmul_sliced(matmul_sliced const &) = delete;
    matmul_sliced& operator=(matmul_sliced const &) = delete;
    matmul_sliced(matmul_sliced &&) noexcept = default;
    matmul_sliced& operator=(matmul_sliced &&) noexcept = default;

    private:
    static unsigned int npack_initial(cxx_param_list & pl) {
        unsigned int npack = L1_CACHE_SIZE;
        param_list_parse(pl, "l1_cache_size", npack);
        npack /= sizeof(typename Arith::elt);
        return npack;
    }
};

template<typename Arith>
matmul_sliced<Arith>::matmul_sliced(matmul_public && P, arith_concrete_base * pxx, cxx_param_list & pl, int optimized_direction)
    : matmul_interface(std::move(P))
    , xab((Arith *) pxx) // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
    , npack(npack_initial(pl))
{
    int const suggest = optimized_direction ^ MM_DIR0_PREFERS_TRANSP_MULT;
    store_transposed = suggest;
    param_list_parse(pl, "mm_store_transposed", store_transposed);
    if (store_transposed != suggest) {
        fprintf(stderr, "Warning, mm_store_transposed"
                " overrides suggested matrix storage ordering\n");
    }           
}


template<typename Arith>
void matmul_sliced<Arith>::build_cache(matrix_u32 && m)
{
    ASSERT_ALWAYS(!m.p.empty());

    uint32_t const i0 = 0;
    uint32_t const i1 = dim[ store_transposed];

    unsigned int const nslices = iceildiv(i1-i0, npack);
    unsigned int const nslices_index = data.size();
    push(0);
    push(0);        // placeholder for alignment

    ASSERT_ALWAYS(nslices);
    unsigned int const packbase = (i1-i0) / nslices;

    unsigned int const nslices1 = (i1-i0) % nslices; /* 1+packbase-sized slices */
    unsigned int const nslices0 = nslices - nslices1;  /* packbase-sized slices */
    unsigned int current;
    unsigned int next = i0;

    uint32_t * ptr = m.p.data();

    for(unsigned int s = 0 ; s < nslices ; s++) {
        current = next;
        unsigned int const npack = packbase + (s < nslices1);
        next = current + npack;

        slice_info si {};
        memset(&si,0,sizeof(si));
        si.nrows = npack;
        si.data_offset = data.size();

        /*
           std::cout << "Packing " << npack << " rows from " << current
           << " to " << next << "\n";
           */
        using L_t = std::vector<std::pair<uint32_t, uint32_t> >;
        L_t L;
        for(unsigned int di = 0 ; di < npack ; di++) {
            for(unsigned int j = 0 ; j < *ptr ; j++) {
                L.emplace_back(ptr[1+j],di);
            }
            ncoeffs += *ptr;
            ptr += 1 + *ptr;
        }
        /* L is the list of (column index, row index) of all
         * coefficients in the current horizontal slice */
        std::ranges::sort(L);

        push32(npack);
        push32(L.size());

        uint32_t j = 0;
        uint32_t i = 0;
        int32_t di_max = 0;
        double sumdj=0;
        double sumdj2=0;
        double sumdi=0;
        double sumdi2=0;
        uint32_t weight=0;
        for(auto const & ji : L) {
            uint32_t const dj = ji.first - j; j += dj;
            int32_t di = ji.second - i; i += di;
            weight++;

            if (di<0) di = -di;
            if (di > di_max) di_max = di;

            double const ddj = dj;
            double const ddi = di;
            sumdj+=ddj;
            sumdj2+=ddj*ddj;
            sumdi+=ddi;
            sumdi2+=ddi*ddi;

            /* If di exceeds +/- 128, then we can emulate the
             * wide displacement by something like:
             * (di_max,0)(0,0)(di_max,0)(0,0)...(di_remainder,dj)
             *
             * It's also really likely that one would gain by
             * moving _both_ up and down in the strip, so that
             * the negative offsets remain significant
             * (otherwise, we'd rather do a cmov based on the
             * sign, or something).
             */

            push(dj); push(i);
        }
        data[nslices_index]++;
        double const dj2_avg = sumdj2 / weight;
        double const dj_avg = sumdj / weight;
        double const dj_sdev = sqrt(dj2_avg-dj_avg*dj_avg);
        double const di2_avg = sumdi2 / weight;
        double const di_avg = sumdi / weight;
        double const di_sdev = sqrt(di2_avg-di_avg*di_avg);
        si.dj_avg = dj_avg;
        si.dj_sdev = dj_sdev;
        si.di_avg = di_avg;
        si.di_sdev = di_sdev;
        si.ncoeffs = weight;

#ifdef  SPARSE_STRIPS
        if (si.dj_avg > 0.5) {
            std::cerr << "hstrip #" << s
                << " is too sparse ; switching to vstrips\n";
            // data.erase(data.begin()+si.data_offset, data.end());
            // break;
        }
#endif
        dslices_info.push_back(si);
    }

    /* There's an other option for the dense strips. For a given
     * set of source values (fitting in L2), store delta_j's for
     * a given destination value. Could actually be eight bits.
     * When the direction changes, increase the destination
     * pointer
     */
#ifdef  SPARSE_STRIPS
    /* We haven't finished the set of horizontal strips. Treat
     * them as vertical sparse strips then. */
#endif

    /* TODO gather info directly from the dslices_info vector
     * now that it has so much data */
    if (nslices1) {
        std::cout
            << "// " << nslices1
            << " sub-slices of " << (packbase+1)
            << " " << rowcol[store_transposed] << "s\n";
    }
    std::cout
        << "// " << nslices0
        << " sub-slices of " << packbase << " "
        << rowcol[store_transposed] << "s\n";
#ifdef  SLICE_STATS
    std::cout << "info per sub-slice:";
    using si_t = std::vector<slice_info>::const_iterator;
    for(si_t s = dslices_info.begin() ; s != dslices_info.end() ; s++) {
        std::cout
            << (s-dslices_info.begin())
            << " " << s->dj_avg << " " << s->dj_sdev
            << " " << s->di_avg << " " << s->di_sdev
            << " " << s->di_max
            << "\n";
    }
    // std::cout << "\n";
#endif
    std::cout << std::flush;
}

template<typename Arith>
int matmul_sliced<Arith>::reload_cache_private()
{
    auto f = matmul_common_reload_cache_fopen(sizeof(typename Arith::elt), *this, MM_MAGIC);
    if (!f) return 0;

    size_t n;
    MATMUL_COMMON_READ_ONE32(n, f.get());

    resize_and_check_meaningful(data, n, f.get());
    MATMUL_COMMON_READ_MANY16(data.data(), n, f.get());

    return 1;
}

template<typename Arith>
void matmul_sliced<Arith>::save_cache_private()
{
    auto f = matmul_common_save_cache_fopen(sizeof(typename Arith::elt), *this, MM_MAGIC);
    if (!f) return;

    size_t const n = data.size();
    MATMUL_COMMON_WRITE_ONE32(n, f.get());
    MATMUL_COMMON_WRITE_MANY16(data.data(), n, f.get());
}

template<typename Arith>
void matmul_sliced<Arith>::mul(void * xdst, void const * xsrc, int d)
{
    ASM_COMMENT("multiplication code");
    const uint16_t * q = data.data();

    uint16_t const nhstrips = *q++;
    q++;        // alignment.
    uint32_t i = 0;
    Arith * x = xab;
    auto const * src = (typename Arith::elt const *) xsrc;
    auto * dst = (typename Arith::elt *) xdst;

    if (d == !store_transposed) {
#ifdef  SLICE_STATS
        std::vector<slice_info>::iterator sit = dslices_info.begin();
        double tick = oncpu_ticks();
#endif
        for(uint16_t s = 0 ; s < nhstrips ; s++) {
            uint32_t const nrows_packed = read32(q);
            typename Arith::elt * where = x->vec_subvec(dst, i);
            x->vec_set_zero(where, nrows_packed);
            ASM_COMMENT("critical loop");
            /* The external function must have the same semantics as this
             * code block */
            uint32_t const ncoeffs_slice = read32(q);
            uint32_t j = 0;
            for(uint32_t c = 0 ; c < ncoeffs_slice ; c++) {
                j += *q++;
                uint32_t const di = *q++;
                x->add(x->vec_item(where, di), x->vec_item(src, j));
            }
            ASM_COMMENT("end of critical loop");
            i += nrows_packed;
#ifdef  SLICE_STATS
            if (!dslices_info.empty()) {
                double ntick = oncpu_ticks();
                sit->spent_time += ntick - tick;
                tick = ntick;
                sit->npasses++;
                sit++;
            }
#endif
        }
    } else {
        /* d == 0 */
        /* BEWARE, it's a priori sub-optimal ! In practice, the
         * difference isn't so striking though. */
        if (iteration[d] == 10) {
            fprintf(stderr, "Warning: Doing many iterations with bad code\n");
        }

        x->vec_set_zero(dst, dim[!d]);
        for(uint16_t s = 0 ; s < nhstrips ; s++) {
            uint32_t j = 0;
            uint32_t const nrows_packed = read32(q);
            ASM_COMMENT("critical loop");
            uint32_t const ncoeffs_slice = read32(q);
            for(uint32_t c = 0 ; c < ncoeffs_slice ; c++) {
                j += *q++;
                uint32_t const di = *q++;
                x->add(x->vec_item(dst, j), x->vec_item(src, i+di));
            }
            i += nrows_packed;
            ASM_COMMENT("end of critical loop");
        }
    }
    ASM_COMMENT("end of multiplication code");
    iteration[d]++;
}

template<typename Arith>
void matmul_sliced<Arith>::report(double scale MAYBE_UNUSED) {
#ifdef  SLICE_STATS
    if (dslices_info.empty())
        return;
    std::ofstream o("dj.stats");
    o << "// Report of timing per slice\n";
    o << "// <snum> <nrows> <ncoeffs> <dj> <npasses> <spent> <M0>\n";
    for(unsigned int s = 0 ; s < dslices_info.size() ; s++) {
        const slice_info& t(dslices_info[s]);
        o << s
            << " " << t.nrows
            << " " << t.ncoeffs
            << " " << t.dj_avg
            << " " << t.npasses
            << " " << t.spent_time
            << " " << (t.spent_time/t.npasses/t.ncoeffs * 1.0e9)
            << " " << t.dj_sdev
            << " " << t.di_avg
            << " " << t.di_sdev
            << " " << t.di_max
            << "\n";
    }
    o << std::flush;
#endif
}

// NOLINTNEXTLINE(misc-use-internal-linkage)
matmul_interface * CADO_CONCATENATE4(new_matmul_, ARITH_LAYER, _, MM_IMPL)(
        matmul_public && P,
        arith_generic * arith,
        cxx_param_list & pl,
        int optimized_direction)
{
    return new matmul_sliced<arith_hard>(std::move(P), arith->concrete(), pl, optimized_direction);
}

/* vim: set sw=4: */
