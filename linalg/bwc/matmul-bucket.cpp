/* Manage the in-memory data for the matrix */

/* It's in C++ because the STL is handy, but that's really all there is
 * to it ; a conversion to C would not be extremely difficult */

#include "cado.h" // IWYU pragma: keep

#include <cinttypes>
#include <climits>
#include <cmath>
#include <cstdarg>         // for va_list, va_end, va_start
#include <cstddef>      // size_t, NULL
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

// IWYU pragma: no_include <memory>
#include <algorithm>    // sort
#include <deque>
#include <limits>       // underlying_type_{min,max}
#include <list>
#include <sstream> // ostringstream // IWYU pragma: keep
#include <string>           // for basic_string
#include <type_traits>      // for __strip_reference_wrapper<>::__type
#include <utility>          // for pair, make_pair, swap
#include <vector>

#include "fmt/format.h"     // for fmt::format, fmt::print

#include "matmul.hpp"       // for matmul_ptr, matmul_public_s, MATMUL_AUX_Z...
#include "macros.h"
#include "verbose.h"    // CADO_VERBOSE_PRINT_BWC_CACHE_BUILD
#include "timing.h"     // wct_seconds
#include "arith-hard.hpp"

/* Make sure that the assembly function is only called if it matches
 * correctly the abase header !! */
#if defined(HAVE_GAS_SYNTAX_ASSEMBLY_SOURCES) && defined(ARITH_LAYER_b64)
// disabling one particular assembly code is done by simply disabling the
// header file (and optionally removing the assembly file from the link
// list is the CMakeLists.txt file, but in reality it's not needed. This
// file being C++, the names are mangled if no header file marks them as
// having C linkage. So commenting out is enough).
#include "matmul-sub-small1.h" // IWYU pragma: keep
#include "matmul-sub-small2.h" // IWYU pragma: keep
#include "matmul-sub-large-fbd.h" // IWYU pragma: keep
#include "matmul-sub-large-fbi.h" // IWYU pragma: keep
#include "matmul-sub-vsc-dispatch.h" // IWYU pragma: keep
#include "matmul-sub-vsc-combine.h" // IWYU pragma: keep
// #define ENABLE_ASM
#endif

#include "matmul-common.hpp"

#include "matmul_facade.hpp"
#include "portability.h" // strdup // IWYU pragma: keep
#include "params.h"

using namespace std;

/* FIXME: this code was modernized in several passes over the years.
 * The current situation is a bit hybrid. When the code is compiled,
 * arith_hard is defined as something specific. But ideally, nothing here
 * but the final function should depend on the arith_hard type, and they
 * should be templated on the typename Arith instead. While this
 * conversion is mostly done for other backends, it's not complete for
 * this one, and there are quite a few places where more work is needed.
 */

/* {{{ Documentation
 *
 * Parameters referred to by this documentation text may be tuned via
 * #define's below.
 *
 * This implementation builds upon the ``sliced'' variant. 
 * The format used with the -sliced implementation is still used, but
 * there is also another format used for very dense strips. We refer to
 * matmul-sliced.cpp for the first format, which is called ``small1''
 * here (and there also).
 *
 * The ``small2'' format is specialised for denser strips.  Within an
 * horizontal strip of less than 4096 rows, when two non-zero columns are
 * never more than 16 positions apart, then the information (dj,i) is
 * stored in 16 bits instead of 32. This reduces the memory footprint
 * quite considerably, since many coefficients are in dense strips.
 * Speed-wise, small2 strips provide a nice speed-up sometimes (not
 * always).
 *
 * When we reach sparser areas, we use another method. Dense (small2 and
 * small1) slices having been discarded, the rest is split in ``large''
 * slices of LSL_NBUCKETS_MAX*256 rows at most (in fact we arrange for the
 * slices to have equal size, so it's a bit less).  For each such large
 * slice, we use scratch space equal to the L2 cache size for performing
 * matrix multiplication. While reading column vector data, we copy
 * coefficient to (up to) LSL_NBUCKETS_MAX ``buckets'', each representing
 * 256 possible row indices. The dispatching information, telling in
 * which bucket a given coefficient has to be copied, is an 8-bit integer
 * stored in a table called ``main''. Each entry in table ``main'' also
 * has a second 8-bit integer indicating the offset to the next useful
 * column entry (because we're considering very sparse blocks, it is not
 * at all uncommon to see average values ~ 10 here, even when
 * LSL_NBUCKETS_MAX==256, so that 65536 rows are considered
 * simultaneously).
 *
 * The number of columns considered while filling the scratch space is
 * limited by the scratch space size. Once it's filled, we ``apply the
 * buckets'', so as to be able to reuse our scratch space buffer for
 * another batch of columns. Thus a ``large slice'' is typically split
 * into several ``vertical blocks''. Extremely sparse large slices will
 * have only one vertical block.
 *
 * ``Applying the buckets'' means taking a list of column coefficients
 * which have been directed to the bucket because it is known that they
 * affect at least one of the 256 corresponding row indices on output. So
 * in addition to the N coefficients stored in the bucket, we have a
 * list of N 8-bit integers indicating which row index is to be modified.
 *
 * [the section below is superseded by the next one]
 * That's it for the basic idea. Now on top of that, we use an extra
 * indirection level for very sparse blocks. Several ``large blocks'' are
 * considered together as part of a ``huge block''. Call this number
 * ``nlarge''. Source coefficients are first dispatched in nlarge ``big
 * buckets'' of appropriate size. The number of columns considered is
 * guided by the fact that we want none of these big buckets to exceed
 * the L2 cache size.  Afterwards, we do a second bispatching from each
 * of these big buckets (in turn) to (up to) LSL_NBUCKETS_MAX small
 * buckets in the same manner as above (except that we no longer need the
 * information on the link to the next useful column.
 *
 * [documentation of the new format]
 * After ``large'' blocks, we consider ``deferred'' blocks, also referred
 * to as (vertical) staircase-shaped blocks. The idea is the following.
 * The portion of the matrix under study is split into several vertical
 * blocks of at most 2^16 entries. In the other direction, we make splits
 * according to the average row density. A first batch of row has some
 * indicative density X, the next one X/2, and so on (the ratio 2 being
 * in fact VSC_BLOCKS_FLUSHFREQ_RATIO) until the sparsest block at the
 * bottom.  Such horizontal blocks need not have the same size. 
 *
 * Vertical strips are processed in order. Buffers corresponding to the
 * sub-matrices in the different sub-matrices created by the horizontal
 * splitting are filled. We thus populate scratch space with coefficients
 * from the source vector. After each such processing of a vertical
 * strip, we consider whether we have to perform flushing w.r.t. the
 * horizontal strips. The buffers corresponding to the densest horizontal
 * strip are flushed after each vstrip processing. Those which are twice
 * sparser are flushed twice less often. Thus the next vstrip occupies
 * some more buffer space. This pattern goes until the sparsest block,
 * which is probably flushed only once. The total amount of buffer space
 * needed is proportional to the number of coefficients in an area having
 * a staircase shape. Assume we have N vstrips of width W. Assume
 * horizontal blocks have a number of rows R0, R1, R2, ... corresponding
 * to densities X, X/2, X/4, ...  (the ratio 2 here is in fact given by
 * VSC_BLOCKS_FLUSHFREQ_RATIO). The formula
 * for the tbuf space is:
 *   max weight of the N matrices of size R0*W in h. strip #0
 * + max weight of the N/2 matrices of size R1*2W in h. strip #1
 * + max weight of the N/4 matrices of size R2*4W in h. strip #2
 * ...
 * It should be noted that the numbers added up in this formula are
 * expected to be of the same magnitude asymptotically. For real-life
 * examples though, it's pretty common to observe somewhat degenerate
 * cases.
 *
 * }}} */


/* {{{ Documentation for parameters */
// take only a portion of the L1 cache.
// #define L1_CACHE_SIZE   192000
// for 8-bytes abt values, this gives 3500 items.
// note that allowing more than 4096 items here (or anything
// SMALL_SLICES_I_BITS dictates) effectively disables small2 slices,
// which pack (i,dj) values in 16 bits.

// tunables for small2 format.
#define SMALL_SLICES_I_BITS     12
#define SMALL_SLICES_DJ_BITS    (16 - SMALL_SLICES_I_BITS)

// make sure that it's something each core can use (i.e. divide by two
// the L2 cache size for a dual-core w/ shared L2).
// #define L2_CACHE_SIZE   1600000

/* To what should we compare the average dj value for determining when we
 * should switch from dense (small) to sparse (large) blocks */
// #define DJ_CUTOFF1   8.0

/* Same, but for stopping large blocks and go to the next type
 * (currently, deferred blocks). Note thatthe time for large blocks grows
 * faster with multi-cores, probably because of the strong dependency on
 * TLB availability.
 */
// #define DJ_CUTOFF2   24.0

/* How many buckets in large slices, and in large slices which are
 * contained in huge slices. */
// #define LSL_NBUCKETS_MAX      256

/* How many large slices MINIMUM must go in a huge slice ? Below this
 * value, we keep going producing large slices. */
// note that this part of the code is currently inactive
#define HUGE_MPLEX_MIN        2 /* default so that the code compiles */

/* Here, I've seen 97 slices be handled best with 2 mplexed rounds. So
 * let's fix 64. */
// note that this part of the code is currently inactive
#define HUGE_MPLEX_MAX        64        /* default so that the code compiles */

/* read batches of VSC_BLOCKS_ROW_BATCH rows, in order to have something
 * which is less sensible to variations incurred by splitting the matrix
 * */
// #define VSC_BLOCKS_ROW_BATCH 512

/* blocks smaller than this amount will be merged with previous / next
 * block */
// #define VSC_BLOCKS_TOO_SMALL_CUTOFF     8192

/* The number of flush frequencies (number of staircase steps) is
 * inversely proportional to this. The larger this value, the taller the
 * steps, and in return, the larger the temp buffers. The effect is
 * generally mild but noticeable.
 */
// #define VSC_BLOCKS_FLUSHFREQ_RATIO 3
/* }}} */

/* {{{ Examples of parameter sets */
/* These defaults work well for 13.4M rows/cols 4x4 submatrix of
 * snfs247.small on a Xeon L5640 (Westmere), for a single core. Note that
 * we define L1_CACHE_SIZE to something insane, because in fact we're
 * utilizing L2.
 */
#if 0
#define L1_CACHE_SIZE   192000
#define L2_CACHE_SIZE   1600000
#define DJ_CUTOFF1   8.0
#define DJ_CUTOFF2   24.0
#define LSL_NBUCKETS_MAX      256
#define VSC_BLOCKS_ROW_BATCH 512
#define VSC_BLOCKS_TOO_SMALL_CUTOFF     8192
#define VSC_BLOCKS_FLUSHFREQ_RATIO 3
#endif

/* This parameter set is successful for 4.6M rows x 2.8 cols, 12x20
 * submatrix of the same snfs247.sparse, for 4 simultaneous cores of a
 * Xeon X3440. In effect, we're disabling large slices here, and use
 * taller steps for vsc.
 */
#if 1
#define L1_CACHE_SIZE   262144
#define L2_CACHE_SIZE   800000
#define DJ_CUTOFF1   8.0
#define DJ_CUTOFF2   24.0 /* This is so small that large slices don't show up */
#define LSL_NBUCKETS_MAX      32
#define VSC_BLOCKS_ROW_BATCH 256
#define VSC_BLOCKS_TOO_SMALL_CUTOFF     8192
#define VSC_BLOCKS_FLUSHFREQ_RATIO 8.0
#endif

/* These used to be the cado-nfs defaults for a long time. They are
 * visibly inferior to the previous parameters on the L5640. */
#if 0
#define L1_CACHE_SIZE   28000
#define L2_CACHE_SIZE   1600000
#define DJ_CUTOFF1   4.5
#define DJ_CUTOFF2   8.0
#define LSL_NBUCKETS_MAX      256
#define VSC_BLOCKS_ROW_BATCH 256
#define VSC_BLOCKS_TOO_SMALL_CUTOFF     8192
#define VSC_BLOCKS_FLUSHFREQ_RATIO 1.7
#endif
/* }}} */

#define xxxDEBUG_BUCKETS

/* This extension is used to distinguish between several possible
 * implementations of the matrix-times-vector product.
 */
#define MM_EXTENSION   "-bucket"

/* MM_MAGIC is stored in files, so as to check compatibility. The upper
 * word (MM_MAGIC_FAMILY) correspond to the implementation itself, the
 * lower one (MM_MAGIC_VERSION) to the n-th binary incompatible change
 * (make sure to bump it) */
#define MM_MAGIC_FAMILY        0xa003UL
#define MM_MAGIC_VERSION       0x1015UL
#define MM_MAGIC (MM_MAGIC_FAMILY << 16 | MM_MAGIC_VERSION)

/* see matmul-basic.c */
#define MM_DIR0_PREFERS_TRANSP_MULT   1

template<typename T>
static inline T * ptrbegin(vector<T>& v) { return v.data(); }
template<typename T>
static inline T const * ptrbegin(vector<T> const & v) { return v.data(); }
template<typename T>
static inline T * ptrend(vector<T>& v) { return v.size() + ptrbegin(v); }
template<typename T>
static inline T const * ptrend(vector<T> const & v) { return v.size() + ptrbegin(v); }

#if 0
static unsigned int idiotic_sum(void * p, unsigned int nbytes)
{
    uint8_t * q = (uint8_t *) p;
    unsigned int res = 0x12345678;
    for( ; nbytes-- ; q++) {
        unsigned int z = *q;
        z ^= z << 24;
        res = (((res << 8) ^ (res * z)) - (z ^ 0xabababab)) ^ res >> 13;
    }
    return res;
}
#endif

#define SLICE_TYPE_NONE           0
#define SLICE_TYPE_SMALL1 1
#define SLICE_TYPE_SMALL1_VBLOCK         2
#define SLICE_TYPE_SMALL2         3

#define SLICE_TYPE_LARGE_ENVELOPE 4

/* When coefficients go through several passes, they are considered as
 * one slice for each pass. Note that the ordering of some slices is
 * mandated by the design, e.g. slices of type LARGE_ASB must follow
 * immediately their parent LARGE_FBI
 */

#if 0
#define SLICE_TYPE_LARGE_FBI      5
#define SLICE_TYPE_LARGE_ASB      6
#endif

#define SLICE_TYPE_HUGE_ENVELOPE  7

#if 0
#define SLICE_TYPE_HUGE_FBI       8
#define SLICE_TYPE_HUGE_FBD       9
#define SLICE_TYPE_HUGE_ASB       10
#endif

#define SLICE_TYPE_DEFER_ENVELOPE       11
#define SLICE_TYPE_DEFER_COLUMN         12
#define SLICE_TYPE_DEFER_DIS            13
#define SLICE_TYPE_DEFER_ROW            14
#define SLICE_TYPE_DEFER_CMB        15

#define SLICE_TYPE_MAX  16


static inline const char* slice_name(int s) { 
    switch(s) {
        case SLICE_TYPE_SMALL1_VBLOCK: return "small1v";
        case SLICE_TYPE_SMALL1: return "small1";
        case SLICE_TYPE_SMALL2: return "small2";
        case SLICE_TYPE_LARGE_ENVELOPE: return "large";
        case SLICE_TYPE_HUGE_ENVELOPE: return "huge";
        case SLICE_TYPE_DEFER_ENVELOPE: return "defer";
        case SLICE_TYPE_DEFER_COLUMN: return "d-col";
        case SLICE_TYPE_DEFER_DIS: return "d-dis";
        case SLICE_TYPE_DEFER_ROW: return "d-row";
        case SLICE_TYPE_DEFER_CMB: return "d-cmb";
        default: return "other";
    }
}

struct slice_header_t {
    uint16_t t;
    uint16_t nchildren;  // w.r.t the tree-like shape of the block structure.
    uint8_t pad[4];      // make the structure 256-bit wide.
    uint32_t i0, i1;
    uint32_t j0, j1;
    uint64_t ncoeffs;
};

/* Ok, there's only one field for the moment. But there might be more
 * eventually */
struct slice_runtime_stats {
    double t = 0;
};

template<typename Arith>
struct matmul_bucket_methods {
    int small1;
    int small2;
    int large;
    int huge;
    int vsc;
    explicit matmul_bucket_methods(const char * desc = nullptr) {
        small1=small2=large=vsc=huge=0;
        if (!desc) {
            /* default configuration */
            small1=small2=large=vsc=1;
            return;
        } 
        const char * p = desc;
        for( ; ; ) {
            const char * q = strchr(p, ',');
            int const n = q ? q-p : strlen(p);
            if (strncmp(p, "small1", n) == 0) {
                ASSERT_ALWAYS(!small1);
                small1=1;
            } else if (strncmp(p, "small2", n) == 0) {
                ASSERT_ALWAYS(!small2);
                small2=1;
            } else if (strncmp(p, "large", n) == 0) {
                ASSERT_ALWAYS(!large);
                large=1;
            } else if (strncmp(p, "huge", n) == 0) {
                ASSERT_ALWAYS(!huge && !vsc);
                huge=1;
            } else if (strncmp(p, "vsc", n) == 0) {
                ASSERT_ALWAYS(!huge && !vsc);
                vsc=1;
            } else {
                fmt::print(stderr, "Parse error: {}\n", p);
            }
            if (!q) break;
            p = q + 1;
        }
        ASSERT_ALWAYS(small1||!small2);
        ASSERT_ALWAYS(small1||small2||large||vsc||huge);
    }
    bool operator<(matmul_bucket_methods const& o) {
        if (vsc != o.vsc) return vsc < o.vsc;
        if (huge != o.huge) return huge < o.huge;
        if (large != o.large) return large < o.large;
        if (small1 != o.small1) return small1 < o.small1;
        if (small2 != o.small2) return small2 < o.small2;
        return false;
    }
    bool something_beyond(const char * s) const {
        matmul_bucket_methods const o(s);
        if (o.huge || o.vsc) return false;
        if (o.large) return huge || vsc;
        if (o.small1 || o.small2) return large || huge || vsc;
        return true;
    }
};

template<typename Arith>
struct matmul_bucket : public matmul_interface {
    using elt = typename Arith::elt;
    /* now our private fields */
    Arith * xab;
    size_t npack = 0;
    size_t scratch1size = 0;
    size_t scratch2size = 0;
    size_t scratch3size = 0;
    elt * scratch1 = nullptr;
    elt * scratch2 = nullptr;
    elt * scratch3 = nullptr;
    vector<uint16_t> t16;       /* For small (dense) slices */
    vector<uint8_t> t8;         /* For large (sparse) slices */
    vector<unsigned int> auxiliary;   /* Various descriptors -- fairly small */
    /* headers are the first thing found in memory */
    vector<slice_header_t> headers;
    slice_runtime_stats main_timing;
    vector<slice_runtime_stats> slice_timings;
    matmul_bucket_methods<Arith> methods;

    void finish_init();

    void build_cache(matrix_u32 &&) override;
    int reload_cache_private() override;
    void save_cache_private() override;
    void mul(void *, const void *, int) override;
    void report(double scale) override;
    void aux(int op, ...) override;

    matmul_bucket(matmul_public &&, arith_concrete_base *, cxx_param_list &, int);
    ~matmul_bucket() override;

    matmul_bucket(matmul_bucket const &) = delete;
    matmul_bucket& operator=(matmul_bucket const &) = delete;
    matmul_bucket(matmul_bucket &&) noexcept = default;
    matmul_bucket& operator=(matmul_bucket &&) noexcept = default;

    private:
    static unsigned int npack_initial(cxx_param_list & pl) {
        unsigned int npack = L1_CACHE_SIZE;
        param_list_parse(pl, "l1_cache_size", npack);
        npack /= sizeof(elt);
        return npack;
    }

    static size_t scratch1size_initial(cxx_param_list & pl) {
        size_t scratch1size = L2_CACHE_SIZE/2;
        param_list_parse(pl, "l2_cache_size", scratch1size);
        scratch1size /= sizeof(elt);
        return scratch1size;
    }
};

template<typename Arith>
matmul_bucket<Arith>::~matmul_bucket()
{
    arith_hard * ab = xab;
    /* yes, it's absolutely intentional! */
    // NOLINTBEGIN(readability-static-accessed-through-instance)
    if (scratch1) ab->free(scratch1);
    if (scratch2) ab->free(scratch2);
    if (scratch3) ab->free(scratch3);
    // NOLINTEND(readability-static-accessed-through-instance)
}

template<typename Arith>
matmul_bucket<Arith>::matmul_bucket(matmul_public && P, arith_concrete_base * pxx, cxx_param_list & pl, int optimized_direction)
    : matmul_interface(std::move(P))
    , xab((arith_hard *) pxx) // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
    , npack(npack_initial(pl))
    , scratch1size(scratch1size_initial(pl))
    , scratch2size(scratch1size * HUGE_MPLEX_MAX) // (1 << 24)/sizeof(elt));
{
    int const suggest = optimized_direction ^ MM_DIR0_PREFERS_TRANSP_MULT;
    store_transposed = suggest;

    param_list_parse(pl, "mm_store_transposed", store_transposed);

    if (store_transposed != suggest) {
        fmt::print(stderr, "Warning, mm_store_transposed"
                " overrides suggested matrix storage ordering\n");
    }   

    methods = matmul_bucket_methods<Arith>(param_list_lookup_string(pl, "matmul_bucket_methods"));
}

/* This moves an element at the tail of a list with no copy, transferring
 * the ownership to the container argument */
template<typename T>
static void transfer(list<T> * ctr, T * elem)
{
    ctr->push_back(T());
    swap(ctr->back(), *elem);
}

template<typename T>
static T underlying_type_max(T const&) { return numeric_limits<T>::max(); }
template<typename T>
static T underlying_type_min(T const&) { return numeric_limits<T>::min(); }
template<typename T, typename U>
static T safe_assign(T & a, U const& b)
{
    ASSERT_ALWAYS(b <= numeric_limits<T>::max());
    ASSERT_ALWAYS(b >= numeric_limits<T>::min());
    return a = b;
}

struct small_slice {
    uint32_t i0 = 0;
    uint32_t i1 = 0;
    unsigned int ncoeffs = 0;
    double dj_avg = 0;
    uint32_t dj_max = 0;
    int is_small2 = 0;

    using Lv_t = std::vector<std::pair<uint16_t, uint16_t> >;
    using Lvci_t = Lv_t::const_iterator;
    using Lvi_t = Lv_t::iterator;

    using Lu_t = std::deque<Lv_t>;
    using Luci_t = Lu_t::const_iterator;
    using Lui_t = Lu_t::iterator;

    Lu_t Lu;
};

struct large_slice {
    slice_header_t hdr[1];
    double dj_avg;
    uint32_t dj_max;

    struct vblock {
        uint32_t j0;
        uint32_t j1;
        unsigned int n;             // number of coeffs.
        unsigned int pad;           // (half the) number of padding coeffs.
        vector<uint8_t> t8c;        // read: contribution to t8.
        vector<unsigned int> auxc;  // read: contribution to auxiliary.
    };

    list<vblock> vbl;

    struct raw {
        vector<uint8_t> ind[LSL_NBUCKETS_MAX];
        vector<uint8_t> main;
        vector<uint32_t> col_sizes;
        vector<uint32_t> pad_sizes;
    };
};

struct huge_slice {
    slice_header_t hdr[1];
    double dj_avg;
    uint32_t dj_max;
    unsigned int nlarge;
    list<large_slice::vblock> vbl;

    struct raw {
        struct subslice {
            vector<uint8_t> ind[LSL_NBUCKETS_MAX];
            vector<uint8_t> main;
        };

        vector<uint8_t> super;
        vector<subslice> subs;
        vector<uint32_t> col_sizes;
        vector<uint32_t> pad_sizes;
    };
};

struct vsc_slice {
    /* This header is mostly a placeholder, and is not here to reflect an
     * additional pass on the data. */
    slice_header_t hdr[1];

    struct step {
        unsigned int defer;
        uint32_t nrows;
        unsigned int density_upper_bound;
        unsigned long tbuf_space;
    };

    vector<step> steps;

    struct middle_slice {
        slice_header_t hdr[1];

        struct sub_slice {
            slice_header_t hdr[1];
            vector<uint16_t> x;
            vector<uint8_t> c;
        };

        vector<sub_slice> sub;
    };

    vector<middle_slice> dispatch;
    vector<uint32_t> transpose_help;
    unsigned long tbuf_space;
};


template<typename Arith>
struct builder {
    matmul_bucket<Arith> * mm;
    std::vector<uint32_t> data;
    // rowhead is a pointer which is naturally excpected to *move* within
    // the data[0] array.
    const uint32_t * rowhead;
    uint32_t nrows_t;
    uint32_t ncols_t;
    const char * rowname;
    const char * colname;
    /* Statistics on the prime parameter affecting performance of
     * subroutines. All are computed as (total sum, #samples) */
    builder(matmul_bucket<Arith> * mm, matrix_u32 && m);

    void do_all_small_slices(uint32_t * p_i0, uint32_t imax, unsigned int npack);
    void do_all_large_slices(uint32_t * p_i0, unsigned int imax, unsigned int scratch1size);
    int prepare_vsc_slices(struct vsc_slice * V, uint32_t i0, uint32_t imax);
    void do_all_huge_slices(uint32_t * p_i0, unsigned int imax, unsigned int scratch2size);
    static void push_vsc_slices(matmul_bucket<Arith> * mm, vsc_slice * V);

    private:
    int do_small_slice(small_slice * S, uint32_t i0, uint32_t i1);

    void split_large_slice_in_vblocks(large_slice * L, large_slice::raw * R, unsigned int scratch1size);
    uint32_t * do_partial_transpose(vector<uint32_t> & cs, uint32_t i0, uint32_t i1);
    int do_large_slice(large_slice * L, uint32_t i0, uint32_t i1, uint32_t imax, unsigned int scratch1size);

    void split_huge_slice_in_vblocks(huge_slice * H, huge_slice::raw * R, unsigned int scratch2size);
    int do_huge_slice(huge_slice * H, uint32_t i0, uint32_t i1, unsigned int scratch2size);

    vector<unsigned long> rowblock_weights(vsc_slice * V);
    void compute_staircase(struct vsc_slice * V);
    void vsc_fill_buffers(vsc_slice * V);


    static void push_small_slice(matmul_bucket<Arith> * mm, small_slice * S);
    static void push_large_slice(matmul_bucket<Arith> * mm, large_slice * L);
    static void push_huge_slice(matmul_bucket<Arith> * mm, huge_slice * H);
};

template<typename Arith>
builder<Arith>::builder(matmul_bucket<Arith> * mm, matrix_u32 && m)
    : mm(mm)
    , data(std::move(m.p))
    , rowhead(data.data())
    , nrows_t(mm->dim[ mm->store_transposed])
    , ncols_t(mm->dim[!mm->store_transposed])
    , rowname(rowcol[ mm->store_transposed])
    , colname(rowcol[!mm->store_transposed])
{
}

/************************************/
/* Elementary (!) buliding routines */

/* there's room for a slice type where data is stored row-major. It would
 * make it possible to avoid most stores, for a nice performance gain. We
 * might even consider accumulating four rows at a time to four different
 * registers, so that the loads are compensated. The in-memory format
 * must account for that.
 */

/* {{{ small slices */

template<typename Arith>
int builder<Arith>::do_small_slice(small_slice * S, uint32_t i0, uint32_t i1)
{
    S->i0 = i0;
    S->i1 = i1;
    S->ncoeffs = 0;
    ASSERT_ALWAYS(i1-i0 <= (1 << 16) );

    /* Make enough vstrips */
    S->Lu.assign(iceildiv(ncols_t, 1UL << 16), small_slice::Lv_t());
    const uint32_t * ptr0 = rowhead;
    /* We're doing a new slice */
    for(uint32_t i = i0 ; i < i1 ; i++) {
        for(unsigned int j = 0 ; j < *rowhead ; j++) {
            uint32_t const jj = rowhead[1+j];
            S->Lu[jj >> 16].emplace_back(jj % (1 << 16), i-i0);
            S->ncoeffs++;
        }
        rowhead += 1 + *rowhead;
    }
    /* L is the list of (column index, row index) of all
     * coefficients in the current horizontal slice */

    /* Convert all j indices to differences */
    S->dj_max = 0;
    S->dj_avg = ncols_t / (double) (S->ncoeffs + !S->ncoeffs);

    uint32_t lu_offset = 0;
    uint32_t j = 0;
    for(auto & lu : S->Lu) {
        std::ranges::sort(lu);
        for(auto const & lv : lu) {
            uint32_t const dj = (lu_offset + lv.first) - j;
            j = lu_offset + lv.first;
            if (dj > S->dj_max) { S->dj_max = dj; }
        }
        lu_offset += 1 << 16;
    }

    /* There are two possible reasons for this slice to be discarded.
     * Either the max dj is too large -- larger than 16 bits -- or the
     * average dj is too large -- larger than some guessed value.
     */
    /* Now that we're splitting in vblocks anyway, it's not important to
     * constrain dj below 2^16 for small1 */
    int keep1 = 1; // S->dj_max < (1UL << 16);
    int keep2 = S->dj_max < (1 << SMALL_SLICES_DJ_BITS);

    if ((i1-i0) >> SMALL_SLICES_I_BITS) keep2=0;

    if (verbose_enabled(CADO_VERBOSE_PRINT_BWC_CACHE_BUILD)) {
        // if (!keep1) fmt::print(" [cannot be small1, beyond impl limits]");
        if (!keep2) fmt::print(" [cannot be small2, beyond impl limits]");
    }

    if (!mm->methods.small2) keep2=0;

    if (mm->methods.something_beyond("small1,small2")) {
        /* Then we have more than just small slices. So we're enforcing
         * our separation criterion based on DJ */
        keep1 = keep1 && S->dj_avg < DJ_CUTOFF1;
        keep2 = keep2 && S->dj_avg < DJ_CUTOFF1;
    }

    S->is_small2 = keep2;

    if (!keep1 && !keep2) {
        rowhead = ptr0;
    }
    return keep1 || keep2;
}

/* }}} */

/* {{{ large slices */

template<typename Arith>
void builder<Arith>::split_large_slice_in_vblocks(large_slice * L, large_slice::raw * R, unsigned int scratch1size)
{
    /* Now split into vslices */
    uint8_t * mp = ptrbegin(R->main);
    uint8_t * ip[LSL_NBUCKETS_MAX];
    for(int k = 0 ; k < LSL_NBUCKETS_MAX ; k++) {
        ip[k] = ptrbegin(R->ind[k]);
    }
    unsigned int vblocknum = 0;
    double vbl_ncols_variance = 0;
    for(uint32_t j = 0 ; j < ncols_t ; vblocknum++) {
        uint32_t const j0 = j;
        size_t n = 0;
        size_t np = 0;
        for( ; j < ncols_t ; j++) {
            size_t const dn = R->col_sizes[j];
            size_t const dnp = R->pad_sizes[j];
            if (n + dn + 2 * (np + dnp) > scratch1size) {
                /* If dn doesn't fit, then we're in really, really big
                 * trouble ! */
                ASSERT_ALWAYS(n != 0);
                break;
            }
            n += dn;
            np += dnp;
        }
        uint32_t const j1 = j;

        large_slice::vblock V;
        V.j0 = j0;
        V.j1 = j1;
        /* First the main block, then all the bucket blocks. These come
         * in order. */
        /* We have nf coeffs, counting padding. This makes 2nf bytes in
         * main, and nf bytes in the bucket blocks */
        V.n = n;
        V.pad = np;
        size_t nf = n + 2 * np;
        V.auxc.reserve(258);
        V.auxc.push_back(V.j1 - V.j0);
        /* XXX Notice that the reader process will consider as a j0 index
         * the farthest possible index from the previous vblock, which
         * means that we _must_ start with a zero in the main block here
         * no matter what -- even if with regard to the way data is set
         * up at the upper level (huge unique vblock), we would have had
         * a positive integer. This digression does not apply to the
         * first vblock.
         *
         * In the case where we're starting with a padding coefficient,
         * then we might actually end up needing to cancel out _several_
         * coefficients, since the padding coeffs need to be cancelled as
         * well.
         */
        V.t8c.resize(3 * nf);
        uint8_t * q = ptrbegin(V.t8c);
        memcpy(q, mp, 2 * nf * sizeof(uint8_t));
        if (vblocknum) {
            uint8_t * pmp = q;
            for(size_t k = 0 ; k + 2 <= nf ; k+=2, pmp+=4) {
                if (pmp[0] != 255 || pmp[1] || pmp[2] || pmp[3])
                    break;
                pmp[0] = 0;
            }
            if (pmp-q) {
                /* It's ok to cancel padding coefficients, anyway it
                 * should be very rare. However, it would be bad to
                 * start cancelling a huge number of these, because then we
                 * would be better off adjusting the pointers (peeling
                 * off the phony computation at the beggining of mf[]).
                 * Problem is that in such a case, we would also have to
                 * adjust the indirection pointer as well, which would be
                 * messy. */
                verbose_printf(CADO_VERBOSE_PRINT_BWC_CACHE_BUILD,
                        "*** cancelled %td padding coefficient\n", (pmp-q)/4);
            }
            pmp[0] = 0;
        }
        unsigned int ind_sizes[LSL_NBUCKETS_MAX] = {0,};
        q += 2*nf;
        for(size_t k = 0 ; k < nf ; k++) {
            ind_sizes[mp[2 * k + 1]]++;
        }
        mp += 2 * nf;
        V.auxc.push_back(nf);
        for(int k = 0 ; k < LSL_NBUCKETS_MAX ; k++) {
            ASSERT(ip[k]+ind_sizes[k] <= ptrend(R->ind[k]));
            memcpy(q, ip[k], ind_sizes[k] * sizeof(uint8_t));
            q += ind_sizes[k];
            ip[k] += ind_sizes[k];
            V.auxc.push_back(ind_sizes[k]);

            // asb_avg[0] += ind_sizes[k];
            // asb_avg[1] += MIN(L->hdr->i1 - L->hdr->i0 - k * 256, 256);
        }
        ASSERT(q == ptrend(V.t8c));

        vbl_ncols_variance += ((double)(j1-j0)) * ((double)(j1-j0));
#if 0
        if (verbose_enabled(CADO_VERBOSE_PRINT_BWC_CACHE_BUILD)) {
            fmt::print(" vbl{}", vblocknum);
            fmt::print(": {} {}..{} ", colname, j0, j1);
            fmt::print("; w {}", V.n);
            fmt::print("; avg dj {:.1f}", (j1 - j0) / (double) V.n);
            if (V.pad) fmt::print("; pad 6*{}", V.pad);
            fmt::print("\n");
        }
#endif
        transfer(&(L->vbl), &V);
    }
    double vbl_ncols_mean = ncols_t;
    vbl_ncols_mean /= vblocknum;
    vbl_ncols_variance /= vblocknum;
    double const vbl_ncols_sdev = sqrt(vbl_ncols_variance - vbl_ncols_mean * vbl_ncols_mean);
    verbose_printf(CADO_VERBOSE_PRINT_BWC_CACHE_BUILD,
            " %u vblocks, sdev/avg = %.2f\n", vblocknum, vbl_ncols_sdev / vbl_ncols_mean);

}

/* Do the partial transpose. Since we're handling sparse areas, it's
 * quite a bit faster than doing the full-blown transpose (even though
 * transposing may be needed anyway because the input matrix is in the
 * wrong order, it saves us hassles to only view the matrix in one given
 * way.
 *
 * Note though that there _is_ some loss ; when testing, we no longer
 * have access to the pre-computed transposed file, since we recompute
 * the transposed parts each time.
 */

template<typename Arith>
uint32_t * builder<Arith>::do_partial_transpose(vector<uint32_t> & cs, uint32_t i0, uint32_t i1) /*{{{*/
{
    uint32_t * cols;
    uint32_t * qptr;

    const uint32_t * ptr = rowhead;
    for(uint32_t i = i0 ; i < i1 ; i++) {
        uint32_t w = *ptr++;
        for( ; w-- ; ) {
            uint32_t const j = *ptr++;
            cs[j]++;
        }
    }
    cols = new uint32_t[ptr - rowhead - (i1 - i0)];
    qptr = cols;

    vector<uint32_t *> colptrs;
    for(uint32_t j = 0 ; j < ncols_t ; j++) {
        colptrs.push_back(qptr);
        qptr += cs[j];
    }
    ASSERT(qptr-cols == ptr - rowhead - (ptrdiff_t) (i1 - i0));
    ptr = rowhead;
    for(uint32_t i = i0 ; i < i1 ; i++) {
        uint32_t w = *ptr++;
        for( ; w-- ; ) {
            uint32_t const j = *ptr++;
            *colptrs[j]++ = i - i0;
        }
    }
    colptrs.clear();
    rowhead = ptr;
    return cols;
}/*}}}*/

template<typename Arith>
int builder<Arith>::do_large_slice(large_slice * L, uint32_t i0, uint32_t i1, uint32_t imax, unsigned int scratch1size)
{
    memset(L->hdr, 0, sizeof(slice_header_t));
    L->hdr->t = SLICE_TYPE_LARGE_ENVELOPE;
    L->hdr->j0 = 0;
    L->hdr->j1 = ncols_t;
    L->hdr->i0 = i0;
    L->hdr->i1 = i1;
    L->hdr->ncoeffs = 0;

    L->dj_max = 0;

    large_slice::raw R[1];

    R->col_sizes.assign(ncols_t, 0);
    R->pad_sizes.assign(ncols_t, 0);

    const uint32_t * ptr0 = rowhead;
    uint32_t * cols = do_partial_transpose(R->col_sizes, i0, i1);
    uint32_t * qptr = cols;

    uint32_t last_j = 0;

    /* First we create a huge unique vblock, and later on decide on
     * how to split it. */
    for(uint32_t j = 0 ; j < ncols_t ; j++) {
        uint32_t len  = R->col_sizes[j];

        if (!len) continue;

        uint32_t diff = j - last_j;
        if (diff > L->dj_max) {
            L->dj_max = diff;
        }
        for( ; diff > 255 ; ) {
            R->main.push_back(255);
            R->main.push_back(0);
            R->ind[0].push_back(0);
            R->main.push_back(0);
            R->main.push_back(0);
            R->ind[0].push_back(0);
            diff -= 255;
            R->pad_sizes[j] += 1;
        }

        /* Only the first entry in the column has a non-zero dj. So place
         * 0 instead inside the loop, and fix it in the end. */
        R->main.push_back(diff);
        last_j = j;
        for( ; len-- ; ) {
            uint32_t i = *qptr++;
            uint8_t const w = i / 256;
            R->main.push_back(w);
            R->main.push_back(0);
            R->ind[w].push_back((uint8_t) i);
            L->hdr->ncoeffs ++;
        }
        R->main.pop_back();
    }
    ASSERT(qptr-cols == rowhead - ptr0 - (ptrdiff_t) (i1 - i0));
    qptr = nullptr;
    delete[] cols;
    cols = nullptr;

    /* The average dj is quite easy */
    L->dj_avg = ncols_t / (double) L->hdr->ncoeffs;

    verbose_printf(CADO_VERBOSE_PRINT_BWC_CACHE_BUILD,
            " w=%" PRIu64 ", avg dj=%.1f, max dj=%u, bucket hit=1/%.1f",
            L->hdr->ncoeffs, L->dj_avg, L->dj_max, LSL_NBUCKETS_MAX * L->dj_avg);

    if (mm->methods.something_beyond("large")) {
        /* Then we may enforce our separation criterion */
        if (L->dj_avg > DJ_CUTOFF2) {
            verbose_printf(CADO_VERBOSE_PRINT_BWC_CACHE_BUILD,
                    "-> too sparse");
            if (imax - i0 < HUGE_MPLEX_MIN * LSL_NBUCKETS_MAX * 256) {
                verbose_printf(CADO_VERBOSE_PRINT_BWC_CACHE_BUILD,
                        "; kept because of short tail\n");
            } else {
                /* We won't keep this slice */
                verbose_printf(CADO_VERBOSE_PRINT_BWC_CACHE_BUILD, "\n");
                rowhead = ptr0;
                return 0;
            }
        }
    }

    verbose_printf(CADO_VERBOSE_PRINT_BWC_CACHE_BUILD,
            "\n");

    split_large_slice_in_vblocks(L, R, scratch1size);

    safe_assign(L->hdr->nchildren, L->vbl.size());

    return 1;
}

/* }}} */

/* {{{ huge slices */

template<typename Arith>
void builder<Arith>::split_huge_slice_in_vblocks(huge_slice * H, huge_slice::raw * R, unsigned int scratch2size)/*{{{*/
{
    struct ptrblock {
        uint8_t * ip[LSL_NBUCKETS_MAX];
        uint8_t * mp;
    };

    /* Now split into vslices */
    // unsigned int lsize = iceildiv(H->hdr->i1 - H->hdr->i0, H->nlarge);
    uint8_t * sp = ptrbegin(R->super);
    vector<ptrblock> ptrs(R->subs.size());
    for(unsigned int i = 0 ; i < R->subs.size() ; i++) {
        ptrs[i].mp = ptrbegin(R->subs[i].main);
        for(int k = 0 ; k < LSL_NBUCKETS_MAX ; k++) {
            ptrs[i].ip[k] = ptrbegin(R->subs[i].ind[k]);
        }
    }
    double vbl_ncols_variance = 0;
    unsigned int vblocknum = 0;
    for(uint32_t j = 0 ; j < ncols_t ; vblocknum++) {
        uint32_t const j0 = j;
#if 0
        /* First way to advance to the next j -- compare the total number
         * of coeffs to the scratch2size argument */
        size_t n = 0;
        size_t np = 0;
        for( ; j < ncols_t ; j++) {
            size_t dn = R->col_sizes[j];
            size_t dnp = R->pad_sizes[j];
            if (n + dn + 2 * (np + dnp) > scratch2size) {
                /* If dn doesn't fit, then we're in really, really big
                 * trouble ! */
                ASSERT_ALWAYS(n != 0);
                break;
            }
            n += dn;
            np += dnp;
        }
        uint32_t j1 = j;
#else
        /* However, there's another way. We might as well arrange so that
         * no mega bucket exceeds the L2 cache. */
        uint32_t j1 = j;
        uint8_t * sp0 = sp;
        uint8_t * spc = sp0;
        /* XXX Notice that the reader process will consider as a j0 index
         * the farthest possible index from the previous vblock, which
         * means that we _must_ start with a zero in the main block here
         * no matter what -- even if with regard to the way data is set
         * up at the upper level (huge unique vblock), we would have had
         * a positive integer. This digression does not apply to the
         * first vblock. */
        if (vblocknum) {
            sp0[0]=0;
        }
        {
            vector<unsigned int> Lsizes(H->nlarge, 0);
            for( ; j < ncols_t && sp != ptrend(R->super) ; ) {
                unsigned int const dj = *sp;
                j += dj;
                if (dj) {
                    /* beginning of a column. */
                    spc = sp;
                    j1 = j;
                }
                sp++;
                unsigned int const k = *sp++;
                Lsizes[k]++;
                if (Lsizes[k] > scratch2size / HUGE_MPLEX_MAX) {
                    Lsizes[k]--;
                    sp = spc;
                    break;
                }
            }
            if (sp == ptrend(R->super)) {
                spc = sp;
                j1 = ncols_t;
            }
        }
        /* TODO: abandon the n/np distinction. It's useful for debugging
         * at best, and it would not hurt to have it displayed only at the
         * upper level.
         */
        size_t n = 0;
        size_t np = 0;
        for(j = j0 ; j < j1 ; j++) {
            n += R->col_sizes[j];
            np += R->pad_sizes[j];
        }
        ASSERT_ALWAYS((n+2*np)*2 == (size_t) (spc - sp0));
        sp = sp0;
#endif

        /* That's the same type, although it contains more stuff. */
        large_slice::vblock V;
        V.j0 = j0;
        V.j1 = j1;
        V.auxc.reserve(1 + H->nlarge * (1 + LSL_NBUCKETS_MAX));
        V.auxc.push_back(V.j1 - V.j0);
        /* First the super block, and then for each contained large
         * slice, the main block, then all the bucket blocks. These come
         * in order. */
        /* We have nf coeffs, counting padding. This makes 2nf bytes in
         * super, nf in mains, and nf in the bucket blocks */
        V.n = n;
        V.pad = np;
        size_t nf = n + 2 * np;
        V.t8c.resize(4 * nf);
        uint8_t * q = ptrbegin(V.t8c);
        memcpy(q, sp, 2 * nf * sizeof(uint8_t));
        /* TODO: If the solution above is to be kept, then of course we'd
         * rather reuse this thing above, since it has already been
         * computed.
         */
        vector<unsigned int> L_sizes(H->nlarge, 0);
        q += 2*nf;
        for(size_t k = 0 ; k < nf ; k++) {
            ASSERT(sp[2 * k + 1] < H->nlarge);
            L_sizes[sp[2 * k + 1]]++;
        }
        sp += 2 * nf;
        V.auxc.push_back(nf);
        for(unsigned int l = 0 ; l < H->nlarge ; l++) {
            V.auxc.push_back(L_sizes[l]);
        }
        for(unsigned int l = 0 ; l < H->nlarge ; l++) {
            unsigned int ind_sizes[LSL_NBUCKETS_MAX] = {0,};
            memcpy(q, ptrs[l].mp, L_sizes[l] * sizeof(uint8_t));
            q += L_sizes[l];
            for(unsigned int k = 0 ; k < L_sizes[l] ; k++) {
                ind_sizes[ptrs[l].mp[k]]++;
            }
            ptrs[l].mp += L_sizes[l];

            // unsigned int i_size = MIN(H->hdr->i1 - H->hdr->i0 - l * lsize, lsize);

            for(int k = 0 ; k < LSL_NBUCKETS_MAX ; k++) {
                ASSERT(ptrs[l].ip[k]-(ptrbegin(R->subs[l].ind[k]))+(ptrdiff_t) ind_sizes[k] <= (ptrdiff_t) R->subs[l].ind[k].size());
                memcpy(q, ptrs[l].ip[k], ind_sizes[k] * sizeof(uint8_t));
                q += ind_sizes[k];
                ptrs[l].ip[k] += ind_sizes[k];
                V.auxc.push_back(ind_sizes[k]);
                // asb_avg[0] += ind_sizes[k];
                // asb_avg[1] += MIN(i_size - k * 256, 256);
            }
        }
        ASSERT(q == ptrend(V.t8c));

        vbl_ncols_variance += ((double)(j1-j0)) * ((double)(j1-j0));

        /* TODO: have some sort of verbosity level */
#if 0
        if (verbose_enabled(CADO_VERBOSE_PRINT_BWC_CACHE_BUILD)) {
            fmt::print(" vbl{}", vblocknum);
            fmt::print(": {} {}..{} ", colname, j0, j1);
            fmt::print("; w {}", V.n);
            fmt::print("; avg dj {:.1f}", (j1 - j0) / (double) V.n);
            if (V.pad) fmt::print("; pad 6*{}", V.pad);
            fmt::print("\n");
        }
#endif

        transfer(&(H->vbl), &V);
    }
    if (!vblocknum) {
        verbose_printf(CADO_VERBOSE_PRINT_BWC_CACHE_BUILD,
                " 0 vblocks\n");
        return;
    }
    double vbl_ncols_mean = ncols_t;
    vbl_ncols_mean /= vblocknum;
    vbl_ncols_variance /= vblocknum;
    double const vbl_ncols_sdev = sqrt(vbl_ncols_variance - vbl_ncols_mean * vbl_ncols_mean);
    verbose_printf(CADO_VERBOSE_PRINT_BWC_CACHE_BUILD,
            " %u vblocks, sdev/avg = %.2f\n", vblocknum, vbl_ncols_sdev / vbl_ncols_mean);
}/*}}}*/

template<typename Arith>
int builder<Arith>::do_huge_slice(huge_slice * H, uint32_t i0, uint32_t i1, unsigned int scratch2size)
{
    memset(H->hdr, 0, sizeof(slice_header_t));
    H->hdr->t = SLICE_TYPE_HUGE_ENVELOPE;
    H->hdr->j0 = 0;
    H->hdr->j1 = ncols_t;
    H->hdr->i0 = i0;
    H->hdr->i1 = i1;
    H->hdr->ncoeffs = 0;

    H->dj_max = 0;

    huge_slice::raw R[1];

    R->col_sizes.assign(ncols_t, 0);
    R->pad_sizes.assign(ncols_t, 0);

    const uint32_t * ptr0 = rowhead;
    uint32_t * cols = do_partial_transpose(R->col_sizes, i0, i1);
    uint32_t * qptr = cols;

    /* How many large slices in this huge slice ? */

    H->nlarge = iceildiv(i1 - i0, LSL_NBUCKETS_MAX * 256);
    ASSERT(H->nlarge >= HUGE_MPLEX_MIN);
    ASSERT(H->nlarge <= HUGE_MPLEX_MAX);
    unsigned int const lsize = iceildiv(i1 - i0, H->nlarge);
    ASSERT(lsize <= LSL_NBUCKETS_MAX * 256);
    ASSERT(lsize * H->nlarge >= i1 - i0);
    R->subs.assign(H->nlarge, huge_slice::raw::subslice());

    verbose_printf(CADO_VERBOSE_PRINT_BWC_CACHE_BUILD,
            " (%u*%u) ", H->nlarge, lsize);
    uint32_t last_j = 0;

    /* First we create a huge unique vblock, and later on decide on
     * how to split it. */
    uint32_t next_dot = ncols_t / H->nlarge;
    for(uint32_t j = 0 ; j < ncols_t ; j++) {
        uint32_t len  = R->col_sizes[j];

        if (verbose_enabled(CADO_VERBOSE_PRINT_BWC_CACHE_BUILD)) {
            if (j > next_dot) {
                printf(".");
                fflush(stdout);
                next_dot += ncols_t / H->nlarge;
            }
        }

        if (!len)
            continue;

        uint32_t diff = j - last_j;
        if (diff > H->dj_max) {
            H->dj_max = diff;
        }
        for( ; diff > 255 ; ) {
            R->super.push_back(255);    // dj
            R->super.push_back(0);      // k
            R->subs[0].main.push_back(0);
            R->subs[0].ind[0].push_back(0);
            R->super.push_back(0);      // dj
            R->super.push_back(0);      // k
            R->subs[0].main.push_back(0);
            R->subs[0].ind[0].push_back(0);
            diff -= 255;
            R->pad_sizes[j] += 1;
        }

        /* Only the first entry in the column has a non-zero dj. So place
         * 0 instead inside the loop, and fix it in the end. */
        R->super.push_back(diff);
        last_j = j;
        for( ; len-- ; ) {
            uint32_t const i = *qptr++;
            uint8_t const w = i / lsize;
            uint8_t const wq = (i % lsize) / 256;
            uint8_t const wr = i % lsize;
            R->super.push_back(w);
            R->super.push_back(0);
            R->subs[w].main.push_back(wq);
            R->subs[w].ind[wq].push_back(wr);
            H->hdr->ncoeffs ++;
        }
        R->super.pop_back();
    }
    verbose_printf(CADO_VERBOSE_PRINT_BWC_CACHE_BUILD, "\n");
    ASSERT(qptr-cols == rowhead - ptr0 - (ptrdiff_t) (i1 - i0));
    qptr = nullptr;
    delete[] cols;
    cols = nullptr;

    /* The average dj is quite easy */
    H->dj_avg = ncols_t / (double) H->hdr->ncoeffs;

    verbose_printf(CADO_VERBOSE_PRINT_BWC_CACHE_BUILD,
            " w=%" PRIu64 ", avg dj=%.1f, max dj=%u, bucket block hit=1/%.1f\n",
            H->hdr->ncoeffs, H->dj_avg, H->dj_max,
            H->nlarge * H->dj_avg);

    if (0) {
        /* Don't keep -- well, here this never occurs */
        rowhead = ptr0;
        return 0;
    }

    split_huge_slice_in_vblocks(H, R, scratch2size);

    safe_assign(H->hdr->nchildren, H->vbl.size());

    return 1;
}

/* }}} */

/* {{{ vertical staircase slices */

static inline unsigned int when_flush(unsigned int k, unsigned int nv, unsigned int defer)
{
    k -= k % defer;
    k += defer - 1;
    k = min(k, nv-1);
    return k;
}

static inline int flush_here(unsigned int k, unsigned int nv, unsigned int defer) {
    return k == when_flush(k,nv,defer);
}

/* {{{ decide on our flushing periods (this step is data independent) */
static vector<unsigned int> flush_periods(unsigned int nvstrips)
{
    vector<unsigned int> fp;
    for(unsigned int nf = 1 ;  ; ) {
        unsigned int const quo = min(255u, iceildiv(nvstrips, nf));
        if (fp.empty() || quo < fp.back()) {
            fp.push_back(quo);
        }
        if (quo == 1)
            break;
        // next ?
        nf = max(nf + 1, (unsigned int) (nf * VSC_BLOCKS_FLUSHFREQ_RATIO));
    }
    std::reverse(fp.begin(), fp.end());
    if (verbose_enabled(CADO_VERBOSE_PRINT_BWC_CACHE_BUILD)) {
        fmt::print("Considering cells to be flushed every {} vstrips\n",
                join(fp, " "));
    }
    ASSERT_ALWAYS(fp.back() == min(255u, nvstrips));
    ASSERT_ALWAYS(fp.front() == 1);

    return fp;
}
/* }}} */

/* {{{ read the row weights */
template<typename Arith>
vector<unsigned long> 
builder<Arith>::rowblock_weights(struct vsc_slice * V)
{
    vector<unsigned long> blockweight;
    const uint32_t * ptr = rowhead;
    /* It's a bit unfortunate, but since here we do speedy parsing of the
     * row weights, we're taking a shortcut which is valid only if our
     * data set spans the entire column range */
    ASSERT_ALWAYS(V->hdr->j0 == 0 && V->hdr->j1 == ncols_t);
    for(uint32_t i = V->hdr->i0 ; i < V->hdr->i1 ; ) {
        unsigned long bw = 0;
        for(uint32_t di = 0 ; di < VSC_BLOCKS_ROW_BATCH && i < V->hdr->i1 ; i++, di++) {
            unsigned int const w = *ptr++;
            ptr += w;
            bw += w;
        }
        blockweight.push_back(bw);
    }
    return blockweight;
}
/* }}} */

/*{{{*/
template<typename Arith>
void builder<Arith>::compute_staircase(struct vsc_slice * V)
{
    unsigned int const nvstrips = V->dispatch.size();

    vector<unsigned long> blockweight = rowblock_weights(V);
    vector<unsigned int> flushperiod = flush_periods(nvstrips);

    /* {{{ first dispatching pass */
    uint32_t ii = V->hdr->i0;
    int s = 0;
    for(unsigned int k = 0 ; ; k++) {
        if (k == blockweight.size() ||
                (s < (int) flushperiod.size() - 1 &&
                 blockweight[k] < VSC_BLOCKS_ROW_BATCH * nvstrips / flushperiod[s]))
        {
            uint32_t const old_ii = ii;
            ii = min(V->hdr->i0 + k * VSC_BLOCKS_ROW_BATCH, V->hdr->i1);
            uint32_t const di = ii - old_ii;
            if (di) {
                vsc_slice::step step;
                step.defer = flushperiod[s];
                step.nrows = di;
                step.tbuf_space = 0;    // to be set later.
                step.density_upper_bound = (s == 0) ?
                    UINT_MAX : iceildiv(nvstrips, flushperiod[s-1]);
                V->steps.push_back(step);
            }
            s++;
        }
        if (k == blockweight.size())
            break;
    }
    /* }}} */
    /* {{{ second pass: do some merges for smallish areas */
    ii = V->hdr->i0;
    for(unsigned int k = 0 ; k < V->steps.size() ; ) {
        unsigned long const w = V->steps[k].nrows;
        uint32_t const old_ii = ii;
        int merge_with_next = 0;
        int merge_with_previous = 0;
        if ((k+1) < V->steps.size())
            merge_with_next = w <= VSC_BLOCKS_TOO_SMALL_CUTOFF; // || w <= V->steps[k+1].nrows / 10;
        if (k > 0)
            merge_with_previous = w <= VSC_BLOCKS_TOO_SMALL_CUTOFF; // || w <= V->steps[k-1].nrows / 10;
        if (!merge_with_previous && !merge_with_next) {
            ii += w;
            k++;
            continue;
        }
        /* Otherwise it's a bit ridiculous to have a separate block for
         * such a small number of rows. So by default, we choose to merge
         * it with the next block, or the previous one if we're reaching
         * the end. */
        V->steps.erase(V->steps.begin() + k);
        if (merge_with_next) {
            verbose_printf(CADO_VERBOSE_PRINT_BWC_CACHE_BUILD,
                    "strip %" PRIu32 "+%lu merged with next\n",
                    old_ii, w);
            ASSERT_ALWAYS(k < V->steps.size());
            V->steps[k].nrows += w;
            // don't increment k here.
            // don't increment ii either.
        } else {
            ASSERT_ALWAYS(merge_with_previous);
            verbose_printf(CADO_VERBOSE_PRINT_BWC_CACHE_BUILD,
                    "strip %" PRIu32 "+%lu merged with previous\n",
                    old_ii, w);
            V->steps[k-1].nrows += w;
            // don't increment k here.
            ii += w;
        }
    }
    /* }}} */
}
/*}}}*/

/*{{{*/
template<typename Arith>
void builder<Arith>::vsc_fill_buffers(vsc_slice * V)
{
    unsigned int const nvstrips = V->dispatch.size();
    uint32_t const width = iceildiv(V->hdr->j1 - V->hdr->j0, nvstrips);
    ASSERT_ALWAYS(width > 0);
    const uint32_t * ptr = rowhead;
    uint32_t i = V->hdr->i0;
    V->tbuf_space = 0;

    for(unsigned int s = 0 ; s < V->steps.size() ; s++) {
        /* This strip is scheduled to be flushed every V->steps[s].defer
         * vstrips. */
        unsigned int const defer = V->steps[s].defer;
        uint32_t const old_ii = i;
        uint32_t const ii = i + V->steps[s].nrows;
        ASSERT_ALWAYS(ii <= V->hdr->i1);
        for( ; i < ii ; i++) {
            unsigned int w = *ptr++;
            for( ; w-- ; ) {
                uint32_t const j = *ptr++;
                unsigned int const d = j / width;
                V->dispatch[d].sub[s].hdr->ncoeffs++;
                V->dispatch[d].sub[s].x.push_back(j % width);
                unsigned int const fidx = when_flush(d,nvstrips,defer);
                V->dispatch[fidx].sub[s].c.push_back(1 + (d % defer));
            }
            for(unsigned int d = 0 ; d < nvstrips ; d+= defer) {
                /* row change */
                unsigned int const fidx = when_flush(d,nvstrips,defer);
                V->dispatch[fidx].sub[s].c.push_back(0);
            }
        }
        /*
           if (verbose_enabled(CADO_VERBOSE_PRINT_BWC_CACHE_BUILD)) {
               fmt::print("Post counts\n");
                for(unsigned int d = 0 ; d < nvstrips ; d++) {
                    fmt::print(" (nrows={}) [{}].sub[{}] : {} coeffs ; xsize={}, csize={}\n", V->steps[s].nrows, d, s, V->dispatch[d].sub[s].hdr->ncoeffs, V->dispatch[d].sub[s].x.size(), V->dispatch[d].sub[s].c.size());
                }
            }
        */
#ifndef NDEBUG
        unsigned int acc = 0;
#endif
        for(unsigned int d = 0 ; d < nvstrips ; d++) {
#ifndef NDEBUG
            acc += V->dispatch[d].sub[s].hdr->ncoeffs;
#endif
            if (!flush_here(d,nvstrips,defer))
                continue;
#ifndef NDEBUG
            ASSERT(V->dispatch[d].sub[s].c.size() == acc + V->steps[s].nrows);
            acc = 0;
#endif
        }
        unsigned long m = 0;
        unsigned long cm = 0;
        for(unsigned int d = 0 ; d < nvstrips ; d++) {
            if (d % defer == 0) cm = 0;
            cm += V->dispatch[d].sub[s].hdr->ncoeffs;
            if (cm >= m) m = cm;
        }
        if (verbose_enabled(CADO_VERBOSE_PRINT_BWC_CACHE_BUILD)) {
            fmt::print("Rows {}+{}", old_ii, V->steps[s].nrows);
            if (V->steps[s].density_upper_bound != UINT_MAX) {
                fmt::print(": d < {}", V->steps[s].density_upper_bound);
            }
            fmt::print("; {} flushes (every {}), tbuf space {}\n", iceildiv(nvstrips, defer), defer, m);
        }
        m += 1; // add space for one dummy pointer.
        V->steps[s].tbuf_space = m;
        V->tbuf_space += m;
    }
    for(unsigned int d = 0 ; d < nvstrips ; d++) {
        for(unsigned int s = 0 ; s < V->steps.size() ; s++) {
            V->dispatch[d].hdr->ncoeffs+=V->dispatch[d].sub[s].hdr->ncoeffs;
        }
        V->hdr->ncoeffs+=V->dispatch[d].hdr->ncoeffs;
    }

    verbose_printf(CADO_VERBOSE_PRINT_BWC_CACHE_BUILD,
            "Total tbuf space %lu (%lu MB)\n",
            V->tbuf_space, (V->tbuf_space * sizeof(typename Arith::elt)) >> 20);
}
/*}}}*/

template<typename Arith>
int builder<Arith>::prepare_vsc_slices(struct vsc_slice * V, uint32_t i0, uint32_t imax)
{
    memset(V->hdr, 0, sizeof(slice_header_t));
    V->hdr->t = SLICE_TYPE_DEFER_ENVELOPE;
    V->hdr->i0 = i0;
    V->hdr->i1 = imax;
    V->hdr->j0 = 0;
    V->hdr->j1 = ncols_t;
    V->hdr->ncoeffs = 0;

    unsigned int const nvstrips = iceildiv(V->hdr->j1 - V->hdr->j0, 1 << 16);
    safe_assign(V->hdr->nchildren, nvstrips);

    uint32_t const width = iceildiv(V->hdr->j1 - V->hdr->j0, nvstrips);
    verbose_printf(CADO_VERBOSE_PRINT_BWC_CACHE_BUILD,
            "%u vstrips of width %" PRIu32 "\n", nvstrips, width);

    // prepare the dispatch slice headers
    for(uint32_t j = V->hdr->j0 ; j < V->hdr->j1 ; j += width) {
        vsc_slice::middle_slice v;
        memset(v.hdr, 0, sizeof(slice_header_t));
        v.hdr->t = SLICE_TYPE_DEFER_COLUMN;
        v.hdr->i0 = V->hdr->i0;
        v.hdr->i1 = V->hdr->i1;
        v.hdr->j0 = j;
        v.hdr->j1 = min(j + width, V->hdr->j1);
        /* ncoeffs is set later, as is the data array x */
        v.hdr->ncoeffs = 0;
        V->dispatch.push_back(v);
    }
    ASSERT(V->dispatch.size() == nvstrips);

    compute_staircase(V);
    for(unsigned int k = 0 ; k < V->dispatch.size() ; k++) {
        safe_assign(V->dispatch[k].hdr->nchildren, V->steps.size());
        uint32_t i0 = V->hdr->i0;
        for(unsigned int l = 0 ; l < V->steps.size() ; l++) {
            vsc_slice::middle_slice::sub_slice w;
            memset(w.hdr, 0, sizeof(slice_header_t));
            w.hdr->t = SLICE_TYPE_DEFER_DIS;
            w.hdr->i0 = i0;
            w.hdr->i1 = (i0 += V->steps[l].nrows);
            w.hdr->j0 = V->dispatch[k].hdr->j0;
            w.hdr->j1 = V->dispatch[k].hdr->j1;
            w.hdr->ncoeffs = 0;
            w.hdr->nchildren = 0;
            V->dispatch[k].sub.push_back(w);
        }
    }
    vsc_fill_buffers(V);

    return 0;
}

/* }}} */

/*************************************************************************/
/* Pushing slices to mm ; all these routines clear the given slice stack */

/* {{{ small slices */
template<typename Arith>
void builder<Arith>::push_small_slice(matmul_bucket<Arith> * mm, small_slice * S)
{
    unsigned int const ncols_t = mm->dim[!mm->store_transposed];
    if (S->is_small2) {
        slice_header_t hdr[1];
        memset(hdr, 0, sizeof(slice_header_t));
        hdr->t = SLICE_TYPE_SMALL2;
        hdr->i0 = S->i0;
        hdr->i1 = S->i1;
        hdr->j0 = 0;
        hdr->j1 = ncols_t;
        hdr->ncoeffs = S->ncoeffs;
        hdr->nchildren = 0;

        /* small2 slices are smaller, and faster -- but more
         * restrictive.  */

        /* Because a small2 slice has data in 16bits chunks, we'll
         * have to ensure proper alignment in the end. */
        ASSERT_ALWAYS((mm->t16.size() & (2-1)) == 0);

        unsigned int j = 0;
        uint32_t lu_offset = 0;
        for( ; ! S->Lu.empty() ; S->Lu.pop_front(), lu_offset+=1<<16) {
            small_slice::Lv_t lv;
            swap(lv, S->Lu.front());
            for(auto const & s : lv) {
                uint32_t const dj = (lu_offset + s.first) - j;
                ASSERT_ALWAYS(dj < (1 << SMALL_SLICES_DJ_BITS) );
                ASSERT_ALWAYS(s.second < (1 << SMALL_SLICES_I_BITS) );
                mm->t16.push_back((dj << SMALL_SLICES_I_BITS) | s.second);
                j = lu_offset + s.first;
            }
        }

        // align.
        for( ; (mm->t16.size() & (2-1)) ; )
            mm->t16.push_back(0);

        mm->headers.push_back(*hdr);
    } else {
        /* How many vertical parts -- this merely has to do with index
         * wraparound, since we're processing data column-major anyway.
         */
        slice_header_t ehdr[1];
        vector<slice_header_t> hdrs;
        memset(ehdr, 0, sizeof(slice_header_t));
        ehdr->t = SLICE_TYPE_SMALL1;
        ehdr->i0 = S->i0;
        ehdr->i1 = S->i1;
        ehdr->j0 = 0;
        ehdr->j1 = ncols_t;
        ehdr->ncoeffs = S->ncoeffs;

        unsigned int lu_offset = 0;
        for( ; ! S->Lu.empty() ; S->Lu.pop_front()) {
            small_slice::Lv_t lv;
            swap(lv, S->Lu.front());

            slice_header_t hdr[1];
            memset(hdr, 0, sizeof(slice_header_t));
            hdr->t = SLICE_TYPE_SMALL1_VBLOCK;
            hdr->i0 = S->i0;
            hdr->i1 = S->i1;
            hdr->j0 = lu_offset;
            hdr->nchildren = 0;
            hdr->ncoeffs = lv.size();
            for(auto const & s : lv) {
                mm->t16.push_back(s.first);
                mm->t16.push_back(s.second);
            }
            lu_offset += (1UL << 16);
            hdr->j1 = min(lu_offset, ncols_t);
            hdrs.push_back(*hdr);
        }
        safe_assign(ehdr->nchildren, hdrs.size());
        mm->headers.push_back(*ehdr);
        mm->headers.insert(mm->headers.end(), hdrs.begin(), hdrs.end());
    }
    mm->ncoeffs += S->ncoeffs;
}

/* }}} */

/* {{{ large slices */
template<typename Arith>
void builder<Arith>::push_large_slice(matmul_bucket<Arith> * mm, large_slice * L)
{
    mm->headers.push_back(*L->hdr);
    for( ; ! L->vbl.empty() ; L->vbl.pop_front()) {
        large_slice::vblock & V(L->vbl.front());
        mm->auxiliary.insert(mm->auxiliary.end(), V.auxc.begin(), V.auxc.end());
        unsigned const t8_size =  V.t8c.size();
        mm->t8.insert(mm->t8.end(), V.t8c.begin(), V.t8c.end());
        // large slices, but not huge slices, may put an odd number
        // of coefficients in t8
        if (t8_size & 1) { mm->t8.push_back(0); }
    }
    mm->ncoeffs += L->hdr->ncoeffs;
}

/* }}} */

/* {{{ huge slices */

template<typename Arith>
void builder<Arith>::push_huge_slice(matmul_bucket<Arith> * mm, huge_slice * H)
{
    mm->headers.push_back(*H->hdr);
    mm->auxiliary.push_back(H->nlarge);
    for( ; ! H->vbl.empty() ; H->vbl.pop_front()) {
        large_slice::vblock & V(H->vbl.front());
        mm->auxiliary.insert(mm->auxiliary.end(), V.auxc.begin(), V.auxc.end());
        unsigned const t8_size =  V.t8c.size();
        ASSERT_ALWAYS((t8_size & 1) == 0);
        mm->t8.insert(mm->t8.end(), V.t8c.begin(), V.t8c.end());
    }
    mm->ncoeffs += H->hdr->ncoeffs;
}
/* }}} */

/******************************************/
/* Iteratively call the building routines */

/* {{{ small slices */
template<typename Arith>
void builder<Arith>::do_all_small_slices(uint32_t * p_i0, uint32_t imax, unsigned int npack)
{
    /* npack is a guess for the expected size of small slices ; they are
     * arranged later to all have approximately equal size.
     */
    unsigned int s;
    uint32_t const i00 = *p_i0;
    unsigned int const nslices = iceildiv(imax - i00, npack);
    uint32_t done = 0;
    for(s = 0 ; s < nslices ; s++) {
        uint32_t const i0 = i00 +  s    * (uint64_t) (imax - i00) / nslices;
        uint32_t const i1 = i00 + (s+1) * (uint64_t) (imax - i00) / nslices;

        small_slice S[1];

        if (verbose_enabled(CADO_VERBOSE_PRINT_BWC_CACHE_BUILD)) {
            fmt::print("Ssl{}: {} {}+{}...", s, rowname, i0, i1-i0);
            fflush(stdout);
        }

        int const keep = do_small_slice(S, i0, i1);

        if (verbose_enabled(CADO_VERBOSE_PRINT_BWC_CACHE_BUILD)) {
            fmt::print(" w {} ; avg dj {:.1f} ; max dj {}{}\n",
                    S->ncoeffs, S->dj_avg, S->dj_max,
                    S->is_small2 ? " [packed]" : "");
            fflush(stdout);
        }

        if (!keep) {
            verbose_printf(CADO_VERBOSE_PRINT_BWC_CACHE_BUILD,
                    "Switching to large slices. Ssl%u to be redone\n", s);
            break;
        }
        push_small_slice(mm, S);
        // transfer(Sq, S);
        done += i1 - i0;
    }
    *p_i0 += done;
}
/* }}} */

/* {{{ large slices */
template<typename Arith>
void builder<Arith>::do_all_large_slices(uint32_t * p_i0, unsigned int imax, unsigned int scratch1size)
{
    unsigned int const rem_nrows = imax - *p_i0;
    unsigned int const nlarge_slices = iceildiv(rem_nrows, LSL_NBUCKETS_MAX * 256);
    uint32_t done = 0;
    for(unsigned int s = 0 ; s < nlarge_slices ; s++) {
        large_slice L[1];
        uint32_t const i0 = * p_i0 +  s      * (uint64_t) rem_nrows / nlarge_slices;
        uint32_t const i1 = * p_i0 + (s + 1) * (uint64_t) rem_nrows / nlarge_slices;

        if (verbose_enabled(CADO_VERBOSE_PRINT_BWC_CACHE_BUILD)) {
            fmt::print("Lsl{} {} {}+{}", s, rowname, i0, i1-i0);
            fflush(stdout);
        }

        int const keep = do_large_slice(L, i0, i1, imax, scratch1size);

        if (!keep) {
            verbose_printf(CADO_VERBOSE_PRINT_BWC_CACHE_BUILD,
                    "Switching to huge slices. Lsl%u to be redone\n", s);
            break;
        }

        // transfer(Lq, L);
        push_large_slice(mm, L);
        done += i1 - i0;
    }
    *p_i0 += done;
}
/* }}} */

/* {{{ huge slices */
template<typename Arith>
void builder<Arith>::do_all_huge_slices(uint32_t * p_i0, unsigned int imax, unsigned int scratch2size)
{
    unsigned int const rem_nrows = imax - *p_i0;
    unsigned int const nhuge_slices = iceildiv(rem_nrows, HUGE_MPLEX_MAX * LSL_NBUCKETS_MAX * 256);
    uint32_t done = 0;
    for(unsigned int s = 0 ; s < nhuge_slices ; s++) {
        huge_slice H[1];
        uint32_t const i0 = * p_i0 +  s      * (uint64_t) rem_nrows / nhuge_slices;
        uint32_t const i1 = * p_i0 + (s + 1) * (uint64_t) rem_nrows / nhuge_slices;

        if (verbose_enabled(CADO_VERBOSE_PRINT_BWC_CACHE_BUILD)) {
            fmt::print("Hsl{} {} {}+{}", s, rowname, i0, i1-i0);
            fflush(stdout);
        }
        do_huge_slice(H, i0, i1, scratch2size);
        push_huge_slice(mm, H);
        // transfer(Hq, H);
        done += i1 - i0;
    }
    *p_i0 += done;
}
/* }}} */

#define xxxCOMPRESS_COMBINERS_1
#define xxxCOMPRESS_COMBINERS_2
#define xxxCOMPRESS_COMBINERS_4

#if defined(COMPRESS_COMBINERS_1) || \
    defined(COMPRESS_COMBINERS_2) || \
    defined(COMPRESS_COMBINERS_4)
#ifdef CADO_MATMUL_SUB_VSC_COMBINE_H
#error "Please either fix or disable the assembly code for COMPRESS_COMBINERS"
/* Anyway this idea doesn't work so well and changes little to the
 * running time. */
#endif
#endif


/* {{{ staircase */
static unsigned long compressed_size(unsigned long s, unsigned int defer MAYBE_UNUSED)
{
    if (0) {
#ifdef COMPRESS_COMBINERS_1
    } else if (defer == 1) {
        return iceildiv(s, 8);
#endif
#ifdef COMPRESS_COMBINERS_2
    } else if (defer <= 3) {
        return iceildiv(s, 4);
#endif
#ifdef COMPRESS_COMBINERS_4
    } else if (defer <= 15) {
        return iceildiv(s, 2);
#endif
    }
    return s;
}

static void append_compressed(vector<uint8_t>& t8, vector<uint8_t> const& S, unsigned int defer MAYBE_UNUSED)
{
    if (0) {
#ifdef  COMPRESS_COMBINERS_1
    } else if (defer == 1) {
        const unsigned int nbits = 1;
        unsigned int c = S.size();
        t8.reserve(t8.size() + compressed_size(S.size(), defer));
        for(unsigned int i = 0 ; i < c ; i += 8 / nbits) {
            uint8_t x = 0;
            for(unsigned int j = 0 ; nbits * j < 8 && i + j < c ; j++) {
                x ^= S[i+j] << (nbits * j);
            }
            t8.push_back(x);
        }
#endif
#ifdef COMPRESS_COMBINERS_2
    } else if (defer <= 3) {
        const unsigned int nbits = 2;
        unsigned int c = S.size();
        t8.reserve(t8.size() + compressed_size(S.size(), defer));
        for(unsigned int i = 0 ; i < c ; i += 8 / nbits) {
            uint8_t x = 0;
            for(unsigned int j = 0 ; nbits * j < 8 && i + j < c ; j++) {
                x ^= S[i+j] << (nbits * j);
            }
            t8.push_back(x);
        }
#endif
#ifdef COMPRESS_COMBINERS_4
    } else if (defer <= 15) {
        const unsigned int nbits = 4;
        unsigned int c = S.size();
        t8.reserve(t8.size() + compressed_size(S.size(), defer));
        for(unsigned int i = 0 ; i < c ; i += 8 / nbits) {
            uint8_t x = 0;
            for(unsigned int j = 0 ; nbits * j < 8 && i + j < c ; j++) {
                x ^= S[i+j] << (nbits * j);
            }
            t8.push_back(x);
        }
#endif
    } else {
        t8.insert(t8.end(), S.begin(), S.end());
    }
}

template<typename Arith>
void builder<Arith>::push_vsc_slices(matmul_bucket<Arith> * mm, vsc_slice * V)
{
    verbose_printf(CADO_VERBOSE_PRINT_BWC_CACHE_BUILD,
            "Flushing staircase slices\n");

    mm->headers.push_back(*V->hdr);

    mm->auxiliary.push_back(V->dispatch.size());
    mm->auxiliary.push_back(V->steps.size());

    for(auto const & step : V->steps) {
        mm->auxiliary.push_back(step.defer);
        mm->auxiliary.push_back(step.tbuf_space);
    }

    /* all headers */
    for(auto const & D : V->dispatch) {
        mm->headers.push_back(*D.hdr);
    }
    for(auto const & D : V->dispatch) {
        ASSERT_ALWAYS(D.sub.size() == V->steps.size());
        for(auto const & S : D.sub) {
            mm->headers.push_back(*S.hdr);
        }
    }

    /* We also add header information for the combination operations,
     * even though the payload is stored alongside with the dispatching
     * stuff for convenience. We do need the combination headers for
     * properly keeping track of where the time is spent eventually.
     *
     * The combining headers are stored en route while flushing
     * combination points.
     */

    unsigned int const cidx = mm->headers.size();

    for(auto const & S : V->dispatch[0].sub) {
        slice_header_t chdr[1];
        memset(chdr, 0, sizeof(slice_header_t));
        chdr->t = SLICE_TYPE_DEFER_ROW;
        chdr->i0 = S.hdr->i0;
        chdr->i1 = S.hdr->i1;
        chdr->j0 = V->hdr->j0;
        chdr->j1 = V->hdr->j1;
        chdr->nchildren = 0;
        chdr->ncoeffs = 0;
        mm->headers.push_back(*chdr);
    }

    vector<pair<pair<unsigned int, unsigned int>, pair<uint64_t, uint64_t> > > csizes;

    // unsigned int o8 = mm->t8.size();

    /* dispatch data goes in dispatching order */
    for(unsigned int k = 0 ; k < V->dispatch.size() ; k++) {
        auto & D(V->dispatch[k]);

        /* The stream of u16 values corresponding to this strip */
        for(auto & S : D.sub) {
            // fmt::print("save dispatch({}), sum={}\n", S.x.size(), idiotic_sum((void*) ptrbegin(S.x), S.x.size() * sizeof(uint16_t)));
            mm->t16.insert(mm->t16.end(), S.x.begin(), S.x.end());
            S.x.clear();
        }
        /* Now the combining data */
        for(unsigned int l = 0 ; l < V->steps.size() ; l++) {
            if (!flush_here(k, V->dispatch.size(), V->steps[l].defer))
                continue;
            unsigned int k0 = k - (k % V->steps[l].defer);
            /* Think about the combining header ! */
            slice_header_t chdr[1];
            memset(chdr, 0, sizeof(slice_header_t));
            chdr->t = SLICE_TYPE_DEFER_CMB;
            chdr->i0 = V->dispatch[k0].sub[l].hdr->i0;
            chdr->i1 = V->dispatch[k0].sub[l].hdr->i1;
            chdr->j0 = V->dispatch[k0].sub[l].hdr->j0;
            chdr->j1 = V->dispatch[k].sub[l].hdr->j1;
            chdr->nchildren = 0;
            chdr->ncoeffs = 0;
            /* flush this one */
            vsc_slice::middle_slice::sub_slice & S(V->dispatch[k].sub[l]);
            for( ; k0 <= k ; k0++) {
                chdr->ncoeffs += V->dispatch[k0].sub[l].hdr->ncoeffs;
            }
            unsigned int const defer = V->steps[l].defer;

            ASSERT_ALWAYS(S.c.size() == chdr->ncoeffs + V->steps[l].nrows);

            csizes.emplace_back(
                        make_pair(k - (k % defer), l),
                        make_pair(0, S.c.size()));
            mm->headers[cidx + l].ncoeffs += chdr->ncoeffs;
            // fmt::print("save combine({}), sum={}\n", S.c.size(), idiotic_sum((void*) ptrbegin(S.c), S.c.size()));
            // fmt::print("m={}\n", *max_element(S.c.begin(), S.c.end()));
            append_compressed(mm->t8, S.c, defer);
            S.c.clear();
            mm->headers.push_back(*chdr);
        }
    }

    /* This fixes the reading order w.r.t. the positioning and size of
     * the combining slices */
    uint64_t p8 = 0;
    for(auto & cs : csizes) {
        // unsigned int k = cs.first.first;
        unsigned int const l = cs.first.second;
        unsigned long const s = cs.second.second;
        unsigned int const defer = V->steps[l].defer;
        unsigned long const count = compressed_size(s, defer);
        /*
           fmt::print("straight-order (defer {}, vstrip #{})"
                ": @{} combine({}), sum={:x}\n",
                V->steps[csizes[i].first.second].defer, csizes[i].first.first,
                csizes[i].second.first, csizes[i].second.second,
                idiotic_sum((void*) &*(mm->t8.begin() + o8 + csizes[i].second.first), csizes[i].second.second));
                */
        cs.second.first = p8;
        p8 += count;
    }
    sort(csizes.begin(), csizes.end());
    mm->auxiliary.push_back(2 * csizes.size() + 1);
    for(auto const & cs : csizes) {
        /*
           fmt::print("straight-order (defer {}, vstrip #{})"
                ": @{} combine({}), sum={:x}\n",
                V->steps[cs.first.second].defer, cs.first.first,
                cs.second.first, cs.second.second,
                idiotic_sum((void*) &*(mm->t8.begin() + o8 + cs.second.first), cs.second.second));
                */
        mm->auxiliary.push_back(cs.second.first);
        mm->auxiliary.push_back(cs.second.second);
    }
    mm->auxiliary.push_back(p8);

    mm->ncoeffs += V->hdr->ncoeffs;
}
/* }}} */

template<typename Arith>
void matmul_bucket<Arith>::build_cache(matrix_u32 && m)
{
    builder<Arith> mb(this, std::move(m));

    verbose_printf(CADO_VERBOSE_PRINT_BWC_CACHE_BUILD,
            "%u rows %u cols\n", dim[0], dim[1]);

    uint32_t main_i0 = 0;
    uint32_t const fence = mb.nrows_t;
    if (methods.small1 || methods.small2)
        mb.do_all_small_slices(&main_i0, fence, npack);

    if (methods.large)
        mb.do_all_large_slices(&main_i0, fence, scratch1size);

    /* Note that vsc and huge are exclusive ! */
    if (methods.vsc && main_i0 < fence) {
        vsc_slice V[1];
        mb.prepare_vsc_slices(V, main_i0, fence);
        scratch3size = MAX(scratch3size, V->tbuf_space);
        builder<Arith>::push_vsc_slices(this, V);
        main_i0 = fence;
    }
    if (methods.huge && main_i0 < fence) {
        mb.do_all_huge_slices(&main_i0, fence, scratch2size);
    }
    if (main_i0 < fence) {
        fmt::print(stderr, "ARGH ! only created a submatrix ({} < {}) !!\n", main_i0, fence);
        exit(1);
    }


    /* done, at last ! */

    finish_init();
}

template<typename Arith>
int matmul_bucket<Arith>::reload_cache_private()/* {{{ */
{
    auto f = matmul_common_reload_cache_fopen(sizeof(elt), *this, MM_MAGIC);
    if (!f) return 0;

    for( ;; ) {
        slice_header_t h {};
        MATMUL_COMMON_READ_ONE16(h.t, f.get());
        MATMUL_COMMON_READ_ONE16(h.nchildren, f.get());
        MATMUL_COMMON_READ_MANY8(h.pad, 4, f.get());
        MATMUL_COMMON_READ_ONE32(h.i0, f.get());
        MATMUL_COMMON_READ_ONE32(h.i1, f.get());
        MATMUL_COMMON_READ_ONE32(h.j0, f.get());
        MATMUL_COMMON_READ_ONE32(h.j1, f.get());
        MATMUL_COMMON_READ_ONE64(h.ncoeffs, f.get());
        if (h.t == SLICE_TYPE_NONE)
            break;
        headers.emplace_back(h);
    }

    MATMUL_COMMON_READ_ONE32(scratch1size, f.get());
    MATMUL_COMMON_READ_ONE32(scratch2size, f.get());
    MATMUL_COMMON_READ_ONE32(scratch3size, f.get());

    size_t n16, n8, naux;
    MATMUL_COMMON_READ_ONE32(n16, f.get());
    MATMUL_COMMON_READ_ONE32(n8, f.get());
    MATMUL_COMMON_READ_ONE32(naux, f.get());

    resize_and_check_meaningful(t16, n16, f.get());
    MATMUL_COMMON_READ_MANY16(ptrbegin(t16), n16, f.get());


    resize_and_check_meaningful(t8, n8, f.get());
    MATMUL_COMMON_READ_MANY8(ptrbegin(t8), n8, f.get());


    resize_and_check_meaningful(auxiliary, naux, f.get());
    MATMUL_COMMON_READ_MANY32(ptrbegin(auxiliary), naux, f.get());

    finish_init();

    return 1;
}/*}}}*/

template<typename Arith>
void matmul_bucket<Arith>::save_cache_private()/*{{{*/
{
    auto f = matmul_common_save_cache_fopen(sizeof(elt), *this, MM_MAGIC);
    if (!f) return;

    for(auto & h : headers) {
        MATMUL_COMMON_WRITE_ONE16(h.t, f.get());
        MATMUL_COMMON_WRITE_ONE16(h.nchildren, f.get());
        MATMUL_COMMON_WRITE_MANY8(h.pad, 4, f.get());
        MATMUL_COMMON_WRITE_ONE32(h.i0, f.get());
        MATMUL_COMMON_WRITE_ONE32(h.i1, f.get());
        MATMUL_COMMON_WRITE_ONE32(h.j0, f.get());
        MATMUL_COMMON_WRITE_ONE32(h.j1, f.get());
        MATMUL_COMMON_WRITE_ONE64(h.ncoeffs, f.get());
    }
    /* useful padding */
    MATMUL_COMMON_WRITE_ONE16(0, f.get());
    MATMUL_COMMON_WRITE_ONE16(0, f.get());
    MATMUL_COMMON_WRITE_ONE32(0, f.get());    // nchildren+padding
    MATMUL_COMMON_WRITE_ONE32(0, f.get());
    MATMUL_COMMON_WRITE_ONE32(0, f.get());
    MATMUL_COMMON_WRITE_ONE32(0, f.get());
    MATMUL_COMMON_WRITE_ONE32(0, f.get());
    MATMUL_COMMON_WRITE_ONE64(0, f.get());

    MATMUL_COMMON_WRITE_ONE32(scratch1size, f.get());
    MATMUL_COMMON_WRITE_ONE32(scratch2size, f.get());
    MATMUL_COMMON_WRITE_ONE32(scratch3size, f.get());
    size_t const n16 = t16.size();
    size_t const n8 = t8.size();
    size_t const naux = auxiliary.size();
    MATMUL_COMMON_WRITE_ONE32(n16, f.get());
    MATMUL_COMMON_WRITE_ONE32(n8, f.get());
    MATMUL_COMMON_WRITE_ONE32(naux, f.get());
    MATMUL_COMMON_WRITE_MANY16(ptrbegin(t16), n16, f.get());
    MATMUL_COMMON_WRITE_MANY8(ptrbegin(t8), n8, f.get());
    MATMUL_COMMON_WRITE_MANY32(ptrbegin(auxiliary), naux, f.get());
}/*}}}*/

// static inline uint32_t read32(uint16_t const * & q) /* {{{ */
// {
//     uint32_t res;
//     res = *q++;
//     res |= ((uint32_t) *q++) << 16;
//     return res;
// } /* }}} */

static inline uint64_t * cvt(arith_hard::elt * a) MAYBE_UNUSED;
static inline uint64_t * cvt(arith_hard::elt * a) {
    return reinterpret_cast<uint64_t *>(a);
}
static inline uint64_t const * cvt(arith_hard::elt const * a) MAYBE_UNUSED;
static inline uint64_t const * cvt(arith_hard::elt const * a) {
    return reinterpret_cast<uint64_t const *>(a);
}
static inline uint64_t ** cvt(arith_hard::elt ** a) MAYBE_UNUSED;
static inline uint64_t ** cvt(arith_hard::elt ** a) {
    return reinterpret_cast<uint64_t **>(a);
}
static inline uint64_t const ** cvt(arith_hard::elt const ** a) MAYBE_UNUSED;
static inline uint64_t const ** cvt(arith_hard::elt const ** a) {
    return reinterpret_cast<uint64_t const **>(a);
}

static const uint16_t * matmul_sub_small1(arith_hard * ab MAYBE_UNUSED, arith_hard::elt * where, arith_hard::elt const * from, const uint16_t * q, unsigned int count)
{
#ifdef CADO_MATMUL_SUB_SMALL1_H
    return matmul_sub_small1_asm(cvt(where), cvt(from), q, count);
#else
    for(uint32_t c = 0 ; c < count ; c++) {
        uint32_t const j = *q++;
        uint32_t const di = *q++;
        // ASSERT(j < (1UL << 16));
        // ASSERT(hdr->j0 + j < hdr->j1);
        // ASSERT(hdr->i0 + di < hdr->i1);
        ab->add(ab->vec_item(where, di), ab->vec_item(from, j));
    }
    return q;
#endif
}

static const uint16_t * matmul_sub_small1_tr(arith_hard * ab MAYBE_UNUSED, arith_hard::elt * where, arith_hard::elt const * from, const uint16_t * q, unsigned int count)
{
#ifdef CADO_MATMUL_SUB_SMALL1_TR_H
    return matmul_sub_small1_tr_asm(cvt(where), cvt(from), q, count);
#else
    for(uint32_t c = 0 ; c < count ; c++) {
        uint32_t const j = *q++;
        uint32_t const di = *q++;
        // ASSERT(j < (1UL << 16));
        // ASSERT(hdr->j0 + j < hdr->j1);
        // ASSERT(hdr->i0 + di < hdr->i1);
        ab->add(ab->vec_item(where, j), ab->vec_item(from, di));
    }
    return q;
#endif
}

static const uint16_t * matmul_sub_small2(arith_hard * ab MAYBE_UNUSED, arith_hard::elt * where, arith_hard::elt const * from, const uint16_t * q, unsigned int count)
{
#ifdef CADO_MATMUL_SUB_SMALL2_H
    return matmul_sub_small2_asm(cvt(where), cvt(from), q, count);
#else
    uint32_t j = 0;
    for(uint32_t c = 0 ; c < count ; c++) {
        uint16_t const h = *q++;
        j += h >> SMALL_SLICES_I_BITS;
        // ASSERT(j < pos->ncols_t);
        uint32_t const di = h & ((1UL << SMALL_SLICES_I_BITS) - 1);
        // ASSERT(hdr->j0 + j < hdr->j1);
        // ASSERT(hdr->i0 + di < hdr->i1);
        ab->add(ab->vec_item(where, di), ab->vec_item(from, j));
    }
    return q;
#endif
}

static const uint16_t * matmul_sub_small2_tr(arith_hard * ab MAYBE_UNUSED, arith_hard::elt * where, arith_hard::elt const * from, const uint16_t * q, unsigned int count)
{
#ifdef CADO_MATMUL_SUB_SMALL2_TR_H
    return matmul_sub_small2_tr_asm(cvt(where), cvt(from), q, count);
#else
    uint32_t j = 0;
    for(uint32_t c = 0 ; c < count ; c++) {
        uint16_t const h = *q++;
        j += h >> SMALL_SLICES_I_BITS;
        // ASSERT(j < pos->ncols_t);
        uint32_t const di = h & ((1UL << SMALL_SLICES_I_BITS) - 1);
        // ASSERT(hdr->j0 + j < hdr->j1);
        // ASSERT(hdr->i0 + di < hdr->i1);
        ab->add(ab->vec_item(where, j), ab->vec_item(from, di));
    }
    return q;
#endif
}

struct pos_desc {
    const uint16_t * q16;
    const uint8_t * q8;
    const unsigned int * ql;
    uint32_t i;
    uint32_t nrows_t;
    uint32_t ncols_t;
};

template<typename Arith>
static inline void matmul_bucket_mul_small1_vblock(matmul_bucket<Arith> * mm, slice_header_t * hdr, typename Arith::elt * dst, typename Arith::elt const * src, int d, struct pos_desc * pos)
{
    ASM_COMMENT("multiplication code -- small1 (dense) slices"); /* {{{ */
    int const usual = d == ! mm->store_transposed;
    arith_hard * ab = mm->xab;
    typename Arith::elt * where      = dst + (usual ? hdr->i0 : hdr->j0);
    typename Arith::elt const * from = src + (usual ? hdr->j0 : hdr->i0);
    
    if ((usual ? hdr->j0 : hdr->i0) == 0) { /* first to come, first to clear */
        ab->vec_set_zero(where, (usual ? (hdr->i1 - hdr->i0) : (hdr->j1 - hdr->j0)));
    }

    uint32_t const ncoeffs_slice = hdr->ncoeffs;
    ASSERT_ALWAYS(pos->i == hdr->i0);   // FIXME -- should disappear.

    if (usual) {
        pos->q16 = matmul_sub_small1(ab, where, from, pos->q16, ncoeffs_slice);
    } else {
        pos->q16 = matmul_sub_small1_tr(ab, where, from, pos->q16, ncoeffs_slice);
    }
    ASM_COMMENT("end of small1 (dense) slices"); /* }}} */
}

template<typename Arith>
static inline void matmul_bucket_mul_small2(matmul_bucket<Arith> * mm, slice_header_t * hdr, typename Arith::elt * dst, typename Arith::elt const * src, int d, struct pos_desc * pos)
{
    ASM_COMMENT("multiplication code -- small2 (dense) slices"); /* {{{ */
    int const usual = d == ! mm->store_transposed;
    arith_hard * ab = mm->xab;
    typename Arith::elt * where      = dst + (usual ? hdr->i0 : hdr->j0);
    typename Arith::elt const * from = src + (usual ? hdr->j0 : hdr->i0);

    if ((usual ? hdr->j0 : hdr->i0) == 0) { /* first to come, first to clear */
        ab->vec_set_zero(where, (usual ? (hdr->i1 - hdr->i0) : (hdr->j1 - hdr->j0)));
    }

    uint32_t const ncoeffs_slice = hdr->ncoeffs;
    ASSERT_ALWAYS(pos->i == hdr->i0);   // FIXME -- should disappear.

    if (usual) {
        pos->q16 = matmul_sub_small2(ab, where, from, pos->q16, ncoeffs_slice);
    } else {
        pos->q16 = matmul_sub_small2_tr(ab, where, from, pos->q16, ncoeffs_slice);
    }
    /* fix alignment in any case */
    pos->q16 += (2 - 1) & - ncoeffs_slice;
    ASM_COMMENT("end of small2 (dense) slices"); /* }}} */
}

// static inline void prepare_buckets_and_fences(arith_hard * x MAYBE_UNUSED, typename Arith::elt ** b, typename Arith::elt ** f, typename Arith::elt * z, const unsigned int * ql, unsigned int n)
// {
//     for(unsigned int k = 0 ; k < n ; k++) {
//         b[k] = z;
//         f[k] = b[k] + ql[k];
//         z += ql[k];
//     }
// }

static inline void prepare_buckets(arith_hard * ab MAYBE_UNUSED, arith_hard::elt ** b, arith_hard::elt * z, const unsigned int * ql, unsigned int n)
{
    for(unsigned int k = 0 ; k < n ; k++) {
        b[k] = z;
        z = ab->vec_subvec(z, ql[k]);
    }
}

static inline void matmul_sub_large_fbi(arith_hard * ab MAYBE_UNUSED, arith_hard::elt ** sb, const arith_hard::elt * z, const uint8_t * q, unsigned int n)
{
#ifdef CADO_MATMUL_SUB_LARGE_FBI_H
    matmul_sub_large_fbi_asm(cvt(sb), cvt(z), q, n);
#else
    /* Dispatch data found in z[0]...z[f(n-1)] such that z[f(i)] is in
     * array pointed to by sb[q[2*i+1]]. The function f(i) is given by
     * the sum q[0]+q[2]+...+q[2*(i-1)]. Exactly 2n coefficients are
     * expected in q[] All the sb[] pointers are increased */
    for(unsigned int c = 0 ; c < n ; c++) {
        z += *q;
        // we might receive zmax and do some checking (see caller)
        // ASSERT_ALWAYS(z < zmax);
        q++;
        ab->set(*sb[*q], *z);
        sb[*q] = ab->vec_subvec(sb[*q], 1);
        q++;
    }
#endif
}

static inline void
matmul_sub_large_fbi_tr(arith_hard * ab MAYBE_UNUSED, arith_hard::elt ** sb, arith_hard::elt * z, const uint8_t * q, unsigned int n)
{
#ifdef CADO_MATMUL_SUB_LARGE_FBI_TR_H
    matmul_sub_large_fbi_tr_asm(cvt(sb), cvt(z), q, n);
#else
    /* Does the converse of the above */
    for(unsigned int c = 0 ; c < n ; c++) {
        z = ab->vec_subvec(z, *q);
        q++;
        ab->add(*z, *sb[*q]);
        sb[*q] = ab->vec_subvec(sb[*q], 1);
        q++;
    }
#endif
}

// static void matmul_sub_large_asb(arith_hard::elt * dst, const arith_hard::elt * z, const uint8_t * q, const unsigned int * ql) __attribute__((__noinline__));
static void matmul_sub_large_asb(arith_hard * ab, arith_hard::elt * dst, const arith_hard::elt * z, const uint8_t * q, const unsigned int * ql)
{
#ifdef CADO_MATMUL_SUB_LARGE_ASB_H
    matmul_sub_large_asb(cvt(dst), cvt(z), q, ql);
#else
    /* This ``applies'' the LSL_NBUCKETS_MAX small buckets whose
     * respective lengths are given by ql[0] to ql[LSL_NBUCKETS_MAX-1].
     *
     * For ql[k] <= i < ql[k+1], z[i] is added to dst[k*256+qk[i]], with
     * qk = q + ql[0] + ... + ql[k-1].
     */
    for(int k = 0 ; k < LSL_NBUCKETS_MAX ; k++) {
        unsigned int l = ql[k];
        for( ; l-- ; )
            /* For padding coeffs, the assertion can fail if
             * we choose a row not equal to (0,0) -- first in
             * the first bucket.
             */
            ab->add(dst[*q++], *z++);
        dst += 256;
    }
#endif
}

static inline void matmul_sub_large_asb_tr(arith_hard * ab MAYBE_UNUSED, const arith_hard::elt * src, arith_hard::elt * z, const uint8_t * q, const unsigned int * ql)
{
#ifdef CADO_MATMUL_SUB_LARGE_ASB_TR_H
    matmul_sub_large_asb_tr(cvt(src), cvt(z), q, ql);
#else
    /* converse of the above */
    for(int k = 0 ; k < LSL_NBUCKETS_MAX ; k++) {
        unsigned int l = ql[k];
        for( ; l-- ; ) {
            /* For padding coeffs, the assertion can fail if
             * we choose a row not equal to (0,0) -- first in
             * the first bucket.
             */
            ab->set(*z, ab->vec_item(src, *q++));
            z = ab->vec_subvec(z, 1);
        }
        src += 256;
    }
#endif
}

static void
matmul_sub_large_fbd(arith_hard * ab MAYBE_UNUSED, arith_hard::elt ** sb, const arith_hard::elt * z, const uint8_t * q, unsigned int n)
{
#ifdef CADO_MATMUL_SUB_LARGE_FBD_H
    matmul_sub_large_fbd_asm(cvt(sb), cvt(z), q, n);
#else
    /* Dispatch data found in z[0]...z[n] such that z[i] is in array
     * pointed to by sb[q[i]]. Exactly n coefficients are expected. All
     * the sb[] pointers are increased */
    for(unsigned int c = 0 ; c < n ; c++) {
        ab->set(*sb[q[c]], *z);
        sb[q[c]] = ab->vec_subvec(sb[q[c]], 1);
        z = ab->vec_subvec(z, 1);
    }
#endif
}

static inline void
matmul_sub_large_fbd_tr(arith_hard * ab, arith_hard::elt ** sb, arith_hard::elt * z, const uint8_t * q, unsigned int n)
{
#ifdef CADO_MATMUL_SUB_LARGE_FBD_TR_H
    matmul_sub_large_fbd_tr_asm(cvt(sb), arith_hard::from_pointer(z), q, n);
#else
    /* Does the converse of the above */
    for(unsigned int c = 0 ; c < n ; c++) {
        ab->set(*z++, *sb[q[c]]++);
    }
#endif
}

template<typename Arith>
static inline void matmul_bucket_mul_large(matmul_bucket<Arith> * mm, slice_header_t * hdr, typename Arith::elt * dst, typename Arith::elt const * src, int d, struct pos_desc * pos)
{
    ASM_COMMENT("multiplication code -- large (sparse) slices"); /* {{{ */

    arith_hard * ab = mm->xab;

    int const usual = d == ! mm->store_transposed;

    typename Arith::elt * where      = dst + (usual ? hdr->i0 : hdr->j0);
    typename Arith::elt const * from = src + (usual ? hdr->j0 : hdr->i0);

    if ((usual ? hdr->j0 : hdr->i0) == 0) { /* first to come, first to clear */
        ab->vec_set_zero(where, (usual ? (hdr->i1 - hdr->i0) : (hdr->j1 - hdr->j0)));
    }

    typename Arith::elt * scratch = mm->scratch1;

    uint32_t j = 0;

    if (usual) {
        for( ; j < pos->ncols_t ; ) {
            uint32_t const j1 = j + *pos->ql++;
            uint32_t const n = *pos->ql++;
            typename Arith::elt * bucket[LSL_NBUCKETS_MAX];
            typename Arith::elt const * inp = from + j;
            prepare_buckets(ab, bucket,scratch,pos->ql,LSL_NBUCKETS_MAX);
            ASSERT_ALWAYS((((unsigned long)pos->q8)&1)==0);
            matmul_sub_large_fbi(ab, bucket, inp, pos->q8, n);
            // the (inp) variable in the call above never exceeds from +
            // j1, although we don't do the check in reality.
            // matmul_sub_large_fbi_boundschecked(x, bucket, inp, pos->q8, n, from + j1);
            matmul_sub_large_asb(ab, where, scratch, pos->q8+2*n, pos->ql);
            pos->q8 += 3*n;
            // fix alignment !
            pos->q8 += n & 1;
            /* FIXME. LSL_NBUCKETS_MAX is a compile-time constant, but it
             * might differ in the computed cache file !!! (occurs here
             * and in several other places as well) */
            pos->ql += LSL_NBUCKETS_MAX;
            j = j1;
        }
    } else {
        for( ; j < pos->ncols_t ; ) {
            uint32_t const j1 = j + *pos->ql++;
            uint32_t const n = *pos->ql++;
            typename Arith::elt * bucket[LSL_NBUCKETS_MAX];
            typename Arith::elt * outp = where + j;
            prepare_buckets(ab, bucket,scratch,pos->ql,LSL_NBUCKETS_MAX);
            matmul_sub_large_asb_tr(ab, from, scratch, pos->q8+2*n, pos->ql);
            ASSERT_ALWAYS((((unsigned long)pos->q8)&1)==0);
            matmul_sub_large_fbi_tr(ab, bucket, outp, pos->q8, n);
            pos->q8 += 3 * n;
            // fix alignment !
            pos->q8 += n & 1;
            pos->ql += LSL_NBUCKETS_MAX;
            j = j1;
        }
    }
    ASM_COMMENT("end of large (sparse) slices"); /* }}} */
}

template<typename Arith>
static inline void matmul_bucket_mul_huge(matmul_bucket<Arith> * mm, slice_header_t * hdr, typename Arith::elt * dst, typename Arith::elt const * src, int d, struct pos_desc * pos)
{
    ASM_COMMENT("multiplication code -- huge (very sparse) slices"); /* {{{ */

    arith_hard * ab = mm->xab;

    int const usual = d == ! mm->store_transposed;

    typename Arith::elt * where      = dst + (usual ? hdr->i0 : hdr->j0);
    typename Arith::elt const * from = src + (usual ? hdr->j0 : hdr->i0);
    if ((usual ? hdr->j0 : hdr->i0) == 0) { /* first to come, first to clear */
        ab->vec_set_zero(where, (usual ? (hdr->i1 - hdr->i0) : (hdr->j1 - hdr->j0)));
    }

    typename Arith::elt * scratch = mm->scratch1;

    uint32_t j = 0;

    uint32_t const di = hdr->i1 - hdr->i0;
    unsigned int const nlarge = *pos->ql++;
    uint32_t const di_sub = iceildiv(di, nlarge);

    if (usual) {
        /* ok -- hold your breath a second. */
        for( ; j < pos->ncols_t ; ) {
            uint32_t const j1 = j + *pos->ql++;
            unsigned int const n = *pos->ql++;
            ASSERT_ALWAYS(n <= mm->scratch2size);
            typename Arith::elt * scratch2 = mm->scratch2;
            typename Arith::elt const * inp = src + j;
            typename Arith::elt * bucket[HUGE_MPLEX_MAX];
            const unsigned int * Lsizes = pos->ql;
            prepare_buckets(ab, bucket,scratch2,pos->ql,nlarge);
            pos->ql += nlarge;
            ASSERT_ALWAYS((((unsigned long)pos->q8)&1)==0);
            matmul_sub_large_fbi(ab, bucket, inp, pos->q8, n);
            pos->q8 += 2 * n;
            for(unsigned int k = 0 ; k < nlarge ; k++) {
                typename Arith::elt * sbucket[LSL_NBUCKETS_MAX];
                prepare_buckets(ab, sbucket,scratch,pos->ql,LSL_NBUCKETS_MAX);
                bucket[k] -= Lsizes[k];
                matmul_sub_large_fbd(ab, sbucket, bucket[k], pos->q8, Lsizes[k]);
                pos->q8 += Lsizes[k];
                typename Arith::elt * outp = where + k * di_sub;
                matmul_sub_large_asb(ab, outp, scratch, pos->q8, pos->ql);
                pos->q8 += Lsizes[k];
                pos->ql += LSL_NBUCKETS_MAX;
            }
            j = j1;
        }
    } else {
        for( ; j < pos->ncols_t ; ) {
            uint32_t const j1 = j + *pos->ql++;
            unsigned int const n = *pos->ql++;
            typename Arith::elt * scratch2 = mm->scratch2;
            typename Arith::elt * outp = dst + j;
            typename Arith::elt * bucket[HUGE_MPLEX_MAX];
            const unsigned int * Lsizes = pos->ql;
            prepare_buckets(ab, bucket,scratch2,pos->ql,nlarge);
            pos->ql += nlarge;
            const uint8_t * q8_saved = pos->q8;
            pos->q8 += 2 * n;
            for(unsigned int k = 0 ; k < nlarge ; k++) {
                typename Arith::elt * sbucket[LSL_NBUCKETS_MAX];
                prepare_buckets(ab, sbucket,scratch,pos->ql,LSL_NBUCKETS_MAX);
                const uint8_t * fill = pos->q8;
                const uint8_t * apply = pos->q8 + Lsizes[k];
                const typename Arith::elt * inp = from + k * di_sub;
                matmul_sub_large_asb_tr(ab, inp, scratch, apply, pos->ql);

                matmul_sub_large_fbd_tr(ab, sbucket, bucket[k], fill, Lsizes[k]);

                pos->q8 += 2 * Lsizes[k];
                pos->ql += LSL_NBUCKETS_MAX;
            }
            swap(pos->q8, q8_saved);
            ASSERT_ALWAYS((((unsigned long)pos->q8)&1)==0);
            matmul_sub_large_fbi_tr(ab, bucket, outp, pos->q8, n);
            pos->q8 += 2 * n;
            swap(pos->q8, q8_saved);
            j = j1;
        }
    }
    ASM_COMMENT("end of huge (very sparse) slices"); /* }}} */
}

static inline void matmul_sub_vsc_dispatch(arith_hard * ab MAYBE_UNUSED, arith_hard::elt * dst, arith_hard::elt const * src, const uint16_t * q, unsigned int count)
{
#ifdef CADO_MATMUL_SUB_VSC_DISPATCH_H
    matmul_sub_vsc_dispatch_asm(cvt(dst), cvt(src), q, count);
#else
    // fmt::print("dispatch({}), sum={:x}\n", count, idiotic_sum((void*)q, count * sizeof(uint16_t)));
    for( ; count-- ; ) {
        ab->set(*dst, src[*q++]);
        dst = ab->vec_subvec(dst, 1);
    }
#endif
}

static inline void matmul_sub_vsc_combine(arith_hard * ab MAYBE_UNUSED, arith_hard::elt * dst, const arith_hard::elt * * mptrs, const uint8_t * q, unsigned int count, unsigned int defer MAYBE_UNUSED)
{
#ifdef CADO_MATMUL_SUB_VSC_COMBINE_H
    matmul_sub_vsc_combine_asm(cvt(dst), cvt(mptrs), q, count, defer);
#else
    // fmt::print("combine({}), defer {}\n", count, defer);
    // fmt::print("combine({}), sum={:x}\n", count, idiotic_sum((void*)q, compressed_size(count, defer)));
    if (0) {
#ifdef COMPRESS_COMBINERS_1
    } else if (defer == 1) {
        const unsigned int nbits = 1;
        for(unsigned int i = 0 ; i < count ; i += 8 / nbits) {
            uint8_t wx = *q++;
            for(unsigned int j = 0 ; nbits * j < 8 && i + j < count ; j++) {
                uint8_t c = wx & ((1 << nbits) - 1);
                ASSERT(c <= defer);
                wx >>= nbits;
                ab->add(*dst, *mptrs[c]);
                mptrs[c] = ab->vec_subvec(mptrs[c], c != 0);
                dst = ab->vec_subvec(dst, c == 0);
            }
        }
#endif
#ifdef COMPRESS_COMBINERS_2
    } else if (defer <= 3) {
        const unsigned int nbits = 2;
        for(unsigned int i = 0 ; i < count ; i += 8 / nbits) {
            uint8_t wx = *q++;
            for(unsigned int j = 0 ; nbits * j < 8 && i + j < count ; j++) {
                uint8_t c = wx & ((1 << nbits) - 1);
                ASSERT(c <= defer);
                wx >>= nbits;
                ab->add(*dst, *mptrs[c]);
                mptrs[c] = ab->vec_subvec(mptrs[c], c != 0);
                dst = ab->vec_subvec(dst, c == 0);
            }
        }
#endif
#ifdef COMPRESS_COMBINERS_4
    } else if (defer <= 15) {
        const unsigned int nbits = 4;
        for(unsigned int i = 0 ; i < count ; i += 8 / nbits) {
            uint8_t wx = *q++;
            for(unsigned int j = 0 ; nbits * j < 8 && i + j < count ; j++) {
                uint8_t c = wx & ((1 << nbits) - 1);
                ASSERT(c <= defer);
                wx >>= nbits;
                ab->add(*dst, *mptrs[c]);
                mptrs[c] = ab->vec_subvec(mptrs[c], c != 0);
                dst = ab->vec_subvec(dst, c == 0);
            }
        }
#endif
    } else {
        for( ; count-- ; ) {
            uint8_t const c = *q++;
            ASSERT(c <= defer);
            ab->add(*dst, *mptrs[c]);
            mptrs[c] = ab->vec_subvec(mptrs[c], c != 0);
            dst = ab->vec_subvec(dst, c == 0);
        }
    }
#endif
}

#ifndef CADO_MATMUL_SUB_VSC_COMBINE_TR_H
static inline void matmul_sub_vsc_combine_tr(arith_hard * ab MAYBE_UNUSED, arith_hard::elt ** mptrs, const arith_hard::elt * qw, const uint8_t * z, unsigned int count, unsigned int defer MAYBE_UNUSED)
{
    // fmt::print("uncombine({}), defer {}\n", count, defer);
    // fmt::print("uncombine({}), sum={:x}\n", count, idiotic_sum((void*)z, compressed_size(count, defer)));
    if (0) {
#ifdef COMPRESS_COMBINERS_1
    } else if (defer == 1) {
        const unsigned int nbits = 1;
        for(unsigned int i = 0 ; i < count ; i += 8 / nbits) {
            uint8_t wx = *z++;
            for(unsigned int j = 0 ; nbits * j < 8 && i + j < count ; j++) {
                uint8_t c = wx & ((1 << nbits) - 1);
                ASSERT(c <= defer);
                wx >>= nbits;
                ab->set(*mptrs[c], *qw);
                mptrs[c] = ab->vec_subvec(mptrs[c], c != 0);
                qw = ab->vec_subvec(qw, c == 0);
            }
        }
#endif
#ifdef COMPRESS_COMBINERS_2
    } else if (defer <= 3) {
        const unsigned int nbits = 2;
        for(unsigned int i = 0 ; i < count ; i += 8 / nbits) {
            uint8_t wx = *z++;
            for(unsigned int j = 0 ; nbits * j < 8 && i + j < count ; j++) {
                uint8_t c = wx & ((1 << nbits) - 1);
                ASSERT(c <= defer);
                wx >>= nbits;
                ab->set(*mptrs[c], *qw);
                mptrs[c] = ab->vec_subvec(mptrs[c], c != 0);
                qw = ab->vec_subvec(qw, c == 0);
            }
        }
#endif
#ifdef COMPRESS_COMBINERS_4
    } else if (defer <= 15) {
        const unsigned int nbits = 4;
        for(unsigned int i = 0 ; i < count ; i += 8 / nbits) {
            uint8_t wx = *z++;
            for(unsigned int j = 0 ; nbits * j < 8 && i + j < count ; j++) {
                uint8_t c = wx & ((1 << nbits) - 1);
                ASSERT(c <= defer);
                wx >>= nbits;
                ab->set(*mptrs[c], *qw);
                mptrs[c] = ab->vec_subvec(mptrs[c], c != 0);
                qw = ab->vec_subvec(qw, c == 0);
            }
        }
#endif
    } else {
        for( ; count-- ; ) {
            uint8_t const c = *z++;
            ASSERT(c <= defer);
            ab->set(*mptrs[c], *qw);
            mptrs[c] = ab->vec_subvec(mptrs[c], c != 0);
            qw = ab->vec_subvec(qw, c == 0);
        }
    }
}
#endif

#ifndef CADO_MATMUL_SUB_VSC_DISPATCH_TR_H
static inline void matmul_sub_vsc_dispatch_tr(arith_hard * ab MAYBE_UNUSED, arith_hard::elt * qr, const arith_hard::elt * q, const uint16_t * z, unsigned int count)
{
    // fmt::print("undispatch({}), sum={:x}\n", count, idiotic_sum((void*)z, count * sizeof(uint16_t)));
    for( ; count-- ; ) {
        ab->add(ab->vec_item(qr, *z++), *q);
        q = ab->vec_subvec(q, 1);
    }
}
#endif

/* in preparation for different steps, including the matrix
 * multiplication itself, we are rebuilding the vsc_slice object, or at
 * least its skeleton ; it makes it easy to run the control loops.  The
 * real data payload is of course not copied again ! */

template<typename Arith>
static inline void rebuild_vsc_slice_skeleton(
        matmul_bucket<Arith> * mm,
        vector<slice_header_t>::iterator & hdr, 
        vsc_slice * V, struct pos_desc * pos,
        unsigned int & nvstrips,
        unsigned int & Midx,
        unsigned int & Didx,
        unsigned int & Ridx,
        unsigned int & Cidx)
{
    *V->hdr = *hdr++;
    nvstrips = *pos->ql++;
    V->dispatch.resize(nvstrips);
    V->steps.resize(*pos->ql++);
    for(unsigned int l = 0 ; l < V->steps.size() ; l++) {
        V->steps[l].defer = *pos->ql++;
        V->steps[l].tbuf_space = *pos->ql++;
    }

    Midx = hdr - mm->headers.begin();
    for(unsigned int k = 0 ; k < V->dispatch.size() ; k++) {
        vsc_slice::middle_slice &D (V->dispatch[k]);
        *D.hdr = *hdr++;
    }
    Didx = hdr - mm->headers.begin();
    for(unsigned int k = 0 ; k < V->dispatch.size() ; k++) {
        vsc_slice::middle_slice &D (V->dispatch[k]);
        D.sub.resize(V->steps.size());
        for(unsigned int l = 0 ; l < V->steps.size() ; l++) {
            vsc_slice::middle_slice::sub_slice & S(D.sub[l]);
            *S.hdr = *hdr++;
            V->steps[l].nrows = S.hdr->i1 - S.hdr->i0;
        }
    }
    /* Skip over the combining headers as well */
    Ridx = hdr - mm->headers.begin();
    for(unsigned int l = 0 ; l < V->steps.size() ; l++) {
        hdr++;
    }
    Cidx = hdr - mm->headers.begin();
    for(unsigned int k = 0 ; k < V->dispatch.size() ; k++) {
        for(unsigned int l = 0 ; l < V->steps.size() ; l++) {
            if (!flush_here(k, V->dispatch.size(), V->steps[l].defer))
                continue;
            hdr++;
        }
    }
}
template<typename Arith>
static inline void matmul_bucket_mul_vsc(matmul_bucket<Arith> * mm, vector<slice_header_t>::iterator & hdr, typename Arith::elt * dst, typename Arith::elt const * src, int d, struct pos_desc * pos)
{
    arith_hard * ab = mm->xab;

    int const usual = d == ! mm->store_transposed;

    typename Arith::elt * where      = dst + (usual ? hdr->i0 : hdr->j0);
    // typename Arith::elt const * from = src + (usual ? hdr->j0 : hdr->i0);
    if ((usual ? hdr->j0 : hdr->i0) == 0) { /* first to come, first to clear */
        ab->vec_set_zero(where, (usual ? (hdr->i1 - hdr->i0) : (hdr->j1 - hdr->j0)));
    }

    typename Arith::elt * scratch = mm->scratch3;

    /* {{{ */

    vsc_slice V[1];
    unsigned int nvstrips;
    unsigned int Midx;
    unsigned int Didx;
    unsigned int Ridx;
    unsigned int Cidx;
    rebuild_vsc_slice_skeleton(mm, hdr, V, pos, nvstrips, Midx, Didx, Ridx, Cidx);

    /*}}}*/

    ASM_COMMENT("multiplication code -- vertical staircase");/*{{{*/

    /* now prepare pointers */
    /* NOTE: we've been getting gcc warnings here for quite a while, and
     * that is explained in the bug report below.
     * https://gcc.gnu.org/bugzilla/show_bug.cgi?id=97222
     */
    vector<typename Arith::elt *> base_ptrs;
    vector<typename Arith::elt *> cptrs;
    typename Arith::elt * q0 = scratch;
    typename Arith::elt * dummy = q0;
    ab->set_zero(*dummy);
    q0++;
    for(unsigned int l = 0 ; l < V->steps.size() ; l++) {
        base_ptrs.push_back(q0);
        cptrs.push_back(q0);
        q0 += V->steps[l].tbuf_space;
    }
    base_ptrs.push_back(q0);

    if (usual) {
        uint32_t const skipover = *pos->ql++;
        pos->ql += skipover;
        for(unsigned int k = 0 ; k < V->dispatch.size() ; k++) {
            vsc_slice::middle_slice const & D(V->dispatch[k]);
            const typename Arith::elt * qr = src + D.hdr->j0;
            mm->slice_timings[Midx].t -= wct_seconds();
            for(unsigned int l = 0 ; l < V->steps.size() ; l++) {
                vsc_slice::middle_slice::sub_slice const & S(D.sub[l]);
                typename Arith::elt * q = cptrs[l];
                unsigned int const count = S.hdr->ncoeffs;
                ASSERT(base_ptrs[l] <= q);
                ASSERT(q <= base_ptrs[l+1]);
                mm->slice_timings[Didx].t -= wct_seconds();
                matmul_sub_vsc_dispatch(ab, q, qr, pos->q16, count);
                q += count;
                pos->q16 += count;
                mm->slice_timings[Didx].t += wct_seconds();
                Didx++;
                ASSERT(q <= base_ptrs[l+1]);
                cptrs[l] = q;
            }
            mm->slice_timings[Midx].t += wct_seconds();
            Midx++;

            for(unsigned int l = 0 ; l < V->steps.size() ; l++) {
                vsc_slice::middle_slice::sub_slice const & S(D.sub[l]);
                unsigned int const defer = V->steps[l].defer;
                if (!flush_here(k,nvstrips,defer))
                    continue;

                /* our different read pointers */
                vector<typename Arith::elt const *> mptrs;
                mptrs.reserve(defer + 1);
                typename Arith::elt const * q = base_ptrs[l];
                mptrs.push_back(dummy);
                for(unsigned int k0 = k - k % defer ; k0 <= k ; k0++) {
                    mptrs.push_back(q);
                    q += V->dispatch[k0].sub[l].hdr->ncoeffs;
                }

                typename Arith::elt * qw = dst + S.hdr->i0;
                unsigned int const count = mm->headers[Cidx].ncoeffs + V->steps[l].nrows;
                ASSERT(V->steps[l].nrows == S.hdr->i1 - S.hdr->i0);
                ASSERT(q - base_ptrs[l] == (ptrdiff_t) mm->headers[Cidx].ncoeffs);

                double t = wct_seconds();
                mm->slice_timings[Cidx].t -= t;
                mm->slice_timings[Ridx+l].t -= t;
                matmul_sub_vsc_combine(ab, qw, ptrbegin(mptrs), pos->q8, count, defer);
                pos->q8 += compressed_size(count, defer);
                t = wct_seconds();
                mm->slice_timings[Cidx].t += t;
                mm->slice_timings[Ridx+l].t += t;
                Cidx++;
                cptrs[l]=base_ptrs[l];
            }
        }
        ASSERT((ptrdiff_t) Cidx == hdr - mm->headers.begin());
    } else {
        /* There's quite an annoying difficulty here. Combining buffers
         * are not read in the same order here, because ``combining'' (whose
         * meaning turns into actually dispatching) occurs of course
         * _earlier_ in the operation.
         */
        pos->ql++;      // only the count, for fast skipover.
        for(unsigned int k = 0 ; k < V->dispatch.size() ; k++) {
            vsc_slice::middle_slice const & D(V->dispatch[k]);
            typename Arith::elt * qr = dst + D.hdr->j0;
            for(unsigned int l = 0 ; l < V->steps.size() ; l++) {
                vsc_slice::middle_slice::sub_slice const & S(D.sub[l]);
                unsigned int const defer = V->steps[l].defer;

                /* Here the test is not the same. We have to fill the
                 * scratch buffers ahead of time. */
                if (k % defer)
                    continue;

                /* our different _write_ pointers */
                vector<typename Arith::elt *> mptrs;
                mptrs.reserve(defer + 1);
                typename Arith::elt * q = base_ptrs[l];
                cptrs[l] = q;
                mptrs.push_back(dummy);
                for(unsigned int k0 = k ; k0 <= when_flush(k,nvstrips,defer) ; k0++) {
                    mptrs.push_back(q);
                    q += V->dispatch[k0].sub[l].hdr->ncoeffs;
                }

                const typename Arith::elt * qw = src + S.hdr->i0;
                // unsigned int count = mm->headers[Cidx].ncoeffs + V->steps[l].nrows;
                const uint8_t * z = pos->q8 + *pos->ql++;
                unsigned int const count = *pos->ql++;

                double t = wct_seconds();
                // mm->slice_timings[Cidx].t -= t;
                mm->slice_timings[Ridx+l].t -= t;
                matmul_sub_vsc_combine_tr(ab, ptrbegin(mptrs), qw, z, count, defer);
                z += compressed_size(count, defer);
                t = wct_seconds();
                // mm->slice_timings[Cidx].t += t;
                mm->slice_timings[Ridx+l].t += t;
                // Cidx++;
            }

            mm->slice_timings[Midx].t -= wct_seconds();
            for(unsigned int l = 0 ; l < V->steps.size() ; l++) {
                vsc_slice::middle_slice::sub_slice const & S(D.sub[l]);
                typename Arith::elt * q = cptrs[l];
                unsigned int const count = S.hdr->ncoeffs;
                ASSERT(base_ptrs[l] <= q);
                ASSERT(q <= base_ptrs[l+1]);
                mm->slice_timings[Didx].t -= wct_seconds();
                matmul_sub_vsc_dispatch_tr(ab, qr, q, pos->q16, count);
                cptrs[l] += count;
                pos->q16 += count;
                mm->slice_timings[Didx].t += wct_seconds();
                Didx++;
                ASSERT(q <= base_ptrs[l+1]);
            }
            mm->slice_timings[Midx].t += wct_seconds();
            Midx++;
        }
        // ASSERT(Cidx == hdr - mm->headers.begin());
        pos->q8 += *pos->ql++;
    }

    ASM_COMMENT("end of vertical staircase");/*}}}*/

    /* It's an envelope type, so make sure the iterator is properly
     * placed eventually */
    hdr--;
}

///////////////////////////////////////////////////////////////////////
// just count how many times an iteration schedules a coefficient in
// the fbi/fbd/asb routines.

template<typename Arith>
void matmul_bucket<Arith>::finish_init()
{
    struct pos_desc pos[1];

    pos->q16 = ptrbegin(t16);
    pos->q8 = ptrbegin(t8);
    pos->ql = ptrbegin(auxiliary);
    pos->i = 0;
    pos->nrows_t = dim[ store_transposed];
    pos->ncols_t = dim[!store_transposed];

    /*
    for(uint16_t h = 0 ; h < headers.size() ; h++) {
        slice_header_t * hdr = & (headers[h]);
        fmt::print("block {} {} [{}..{}[ x [{}..{}[, {} cld, {} coeffs\n",
                h, slice_name(hdr->t),
                hdr->i0, hdr->i1,
                hdr->j0, hdr->j1,
                hdr->nchildren, hdr->ncoeffs);
    }
    */
    for(auto hdr = headers.begin() ; hdr != headers.end() ; ++hdr) {
        switch(hdr->t) {
            case SLICE_TYPE_SMALL1:
                break;
            case SLICE_TYPE_SMALL1_VBLOCK:
                // mms1_ncoeffs += hdr->ncoeffs;
                pos->q16 += 2 * hdr->ncoeffs;
                break;
            case SLICE_TYPE_SMALL2:
                pos->q16 += hdr->ncoeffs;
                pos->q16 += (2 - 1) & - hdr->ncoeffs;
                // mms2_ncoeffs += hdr->ncoeffs;
                break;
            case SLICE_TYPE_LARGE_ENVELOPE:
                {
                    uint32_t j = 0;
                    for( ; j < pos->ncols_t ; ) {
                        uint32_t const j1 = j + *pos->ql++;
                        uint32_t const n = *pos->ql++;
                        // fbi_ncoeffs += n;
                        // asb_ncoeffs += n;
                        pos->q8 += 3*n;
                        pos->q8 += n & 1;
                        pos->ql += LSL_NBUCKETS_MAX;
                        j = j1;
                    }
                }
                break;
            case SLICE_TYPE_HUGE_ENVELOPE:
                {
                    uint32_t j = 0;
                    unsigned int const nlarge = *pos->ql++;
                    for( ; j < pos->ncols_t ; ) {
                        uint32_t const j1 = j + *pos->ql++;
                        unsigned int const n = *pos->ql++;
                        // fbi_ncoeffs += n;
                        // fbd_ncoeffs += n;
                        // asb_ncoeffs += n;
                        pos->ql += nlarge;
                        pos->ql += nlarge * LSL_NBUCKETS_MAX;
                        pos->q8 += 4*n;
                        j = j1;
                    }
                }
                break;
            case SLICE_TYPE_DEFER_ENVELOPE:
                {
                    vsc_slice V[1];
                    unsigned int nvstrips;
                    unsigned int Midx;
                    unsigned int Didx;
                    unsigned int Ridx;
                    unsigned int Cidx;
                    rebuild_vsc_slice_skeleton(this, hdr, V, pos, nvstrips, Midx, Didx, Ridx, Cidx);
                    hdr--;
                    /*
                    for(unsigned int s = 0 ; s < V->steps.size() ; s++) {
                        unsigned int defer = V->steps[s].defer;
                        fmt::print("Rows {}+{}", headers[Ridx++].i0, V->steps[s].nrows);
                        if (V->steps[s].density_upper_bound != UINT_MAX) {
                            fmt::print(": d < {}", V->steps[s].density_upper_bound);
                        }
                        fmt::print("; {} flushes (every {}), tbuf space {}\n", iceildiv(nvstrips, defer), defer, V->steps[s].tbuf_space);
                    }
                    */
                }
                break;
            default:
                fmt::print(stderr, "Unexpected block {} encountered\n",
                        slice_name(hdr->t));
                break;
        }
        if (hdr->j1 == pos->ncols_t) {
            pos->i = hdr->i1;
        }
    }

#if 0
    for(uint16_t h = 0 ; h < headers.size() ; h++) {
        slice_header_t * hdr = & (headers[h]);
        ASSERT_ALWAYS(pos->i == hdr->i0);
        switch(hdr->t) {
            case SLICE_TYPE_SMALL1_VBLOCK:
                mms1_ncoeffs += hdr->ncoeffs;
                pos->q16 += 2 * hdr->ncoeffs;
                break;
            case SLICE_TYPE_SMALL2:
                pos->q16 += hdr->ncoeffs;
                pos->q16 += (2 - 1) & - hdr->ncoeffs;
                mms2_ncoeffs += hdr->ncoeffs;
                break;
            case SLICE_TYPE_LARGE_ENVELOPE:
                {
                    uint32_t j = 0;
                    for( ; j < pos->ncols_t ; ) {
                        uint32_t j1 = j + *pos->ql++;
                        uint32_t n = *pos->ql++;
                        fbi_ncoeffs += n;
                        asb_ncoeffs += n;
                        // pos->q8 += 3*n;
                        pos->ql += LSL_NBUCKETS_MAX;
                        j = j1;
                    }
                }
                break;
            case SLICE_TYPE_HUGE_ENVELOPE:
                {
                    uint32_t j = 0;
                    unsigned int nlarge = *pos->ql++;
                    for( ; j < pos->ncols_t ; ) {
                        uint32_t j1 = j + *pos->ql++;
                        unsigned int n = *pos->ql++;
                        fbi_ncoeffs += n;
                        fbd_ncoeffs += n;
                        asb_ncoeffs += n;
                        pos->ql += nlarge;
                        pos->ql += nlarge * LSL_NBUCKETS_MAX;
                        // pos->q8 += 4*n;
                        j = j1;
                    }
                }
            case SLICE_TYPE_DEFER_ENVELOPE:
            default:
                break;
        }
        if (hdr->j1 == pos->ncols_t) {
            pos->i = hdr->i1;
        }
    }
#endif

    arith_hard * ab = xab;
    // NOLINTBEGIN(readability-static-accessed-through-instance)
    scratch1 = ab->alloc(scratch1size);
    scratch2 = ab->alloc(scratch2size);
    scratch3 = ab->alloc(scratch3size);
    // NOLINTEND(readability-static-accessed-through-instance)

    slice_timings.resize(headers.size());
}

template<typename Arith>
static void matmul_bucket_zero_stats(matmul_bucket<Arith> * mm)
{
    vector<slice_header_t>::iterator hdr;
    for(hdr = mm->headers.begin() ; hdr != mm->headers.end() ; hdr++) {
        unsigned int const hidx = hdr - mm->headers.begin();
        mm->slice_timings[hidx].t = 0;
    }
}

template<typename Arith>
static inline void matmul_bucket_mul_small1(matmul_bucket<Arith> * mm, vector<slice_header_t>::iterator & hdr, typename Arith::elt * dst, typename Arith::elt const * src, int d, struct pos_desc * pos)
{
    slice_header_t  const& h(*hdr++);
    for(unsigned int i = 0 ; i < h.nchildren ; i++, hdr++) {
        ASSERT(hdr != mm->headers.end());
        ASSERT(hdr->t == SLICE_TYPE_SMALL1_VBLOCK);
        mm->slice_timings[hdr - mm->headers.begin()].t -= wct_seconds();
        matmul_bucket_mul_small1_vblock(mm, &*hdr, dst, src, d, pos);
        mm->slice_timings[hdr - mm->headers.begin()].t += wct_seconds();
    }
    hdr--;
}

template<typename Arith>
static inline void matmul_bucket_mul_loop(matmul_bucket<Arith> * mm, typename Arith::elt * dst, typename Arith::elt const * src, int d, struct pos_desc * pos)
{
    vector<slice_header_t>::iterator hdr;

    for(hdr = mm->headers.begin() ; hdr != mm->headers.end() ; hdr++) {
        unsigned int const hidx = hdr - mm->headers.begin();
        mm->slice_timings[hidx].t -= wct_seconds();
        switch(hdr->t) {
            case SLICE_TYPE_SMALL2:
                matmul_bucket_mul_small2(mm, &*hdr, dst, src, d, pos);
                break;
            case SLICE_TYPE_SMALL1:
                matmul_bucket_mul_small1(mm, hdr, dst, src, d, pos);
                break;
            case SLICE_TYPE_LARGE_ENVELOPE:
                matmul_bucket_mul_large(mm, &*hdr, dst, src, d, pos);
                break;
            case SLICE_TYPE_HUGE_ENVELOPE:
                matmul_bucket_mul_huge(mm, &*hdr, dst, src, d, pos);
                break;
            case SLICE_TYPE_DEFER_ENVELOPE:
                matmul_bucket_mul_vsc(mm, hdr, dst, src, d, pos);
                break;
            /* some slice types are ``contained'', and we should never
             * see them. NOTE that this implies in particular that the
             * timings for type "defer" (DEFER_ENVELOPE) is counted here,
             * thus including all overhead, while the timing for the
             * inner objects is counted from within the subroutines. */
            case SLICE_TYPE_SMALL1_VBLOCK:
            case SLICE_TYPE_DEFER_COLUMN:
            case SLICE_TYPE_DEFER_ROW:
            case SLICE_TYPE_DEFER_CMB:
            case SLICE_TYPE_DEFER_DIS:
                break;
            /* The current implementation no longer accepts data not
             * obeying the header structure, so the default branch also
             * aborts. */
            default:
                fmt::print(stderr, "Bogus slice type seen: {}\n", hdr->t);
                ASSERT_ALWAYS(0);
        }
        if (hdr->j1 == pos->ncols_t) { pos->i = hdr->i1; }
        mm->slice_timings[hidx].t += wct_seconds();
    }
}

template<typename Arith>
void matmul_bucket<Arith>::mul(void * xdst, void const * xsrc, int d)
{
    struct pos_desc pos[1];

    auto * dst = (elt *) xdst;
    auto const * src = (elt const *) xsrc;
    pos->q16 = ptrbegin(t16);
    pos->q8 = ptrbegin(t8);
    pos->ql = ptrbegin(auxiliary);
    pos->i = 0;
    pos->nrows_t = dim[ store_transposed];
    pos->ncols_t = dim[!store_transposed];

    arith_hard * ab = xab;

    main_timing.t -= wct_seconds();

    if (d == !store_transposed) {
        /* This is the ``normal'' case (matrix times vector). */
    } else {
        /* d == store_transposed */
        /* BEWARE, it's a priori sub-optimal ! In practice, the
         * difference isn't so striking though. */
        if (iteration[d] == 10) {
            fmt::print(stderr, "Warning: Doing many iterations with bad code\n");
        }
        /* We zero out the dst area beforehand */
        ab->vec_set_zero(dst, pos->ncols_t);
    }
    matmul_bucket_mul_loop(this, dst, src, d, pos);

    main_timing.t += wct_seconds();

    iteration[d]++;
}

template<typename Arith>
static std::ostream& matmul_bucket_report_vsc(std::ostream& os, matmul_bucket<Arith> * mm, double scale, vector<slice_header_t>::iterator & hdr, double * p_t_total)
{
    uint64_t scale0;
    scale0 = (mm->iteration[0] + mm->iteration[1]);
    unsigned int const nvstrips = hdr->nchildren;
    hdr++;
    unsigned int const nsteps = hdr->nchildren;
    vector<pair<uint64_t, double> > dtime(nsteps);
    vector<pair<uint64_t, double> > ctime(nsteps);
    for(unsigned int k = 0 ; k < nvstrips ; k++) {
        ASSERT_ALWAYS(hdr->t == SLICE_TYPE_DEFER_COLUMN);
        hdr++;
    }
    // hdr+=nvstrips;
    for(unsigned int k = 0 ; k < nvstrips ; k++) {
        for(unsigned int l = 0 ; l < nsteps ; l++) {
            ASSERT_ALWAYS(hdr->t == SLICE_TYPE_DEFER_DIS);
            double const t = mm->slice_timings[hdr - mm->headers.begin()].t;
            uint64_t const nc = hdr->ncoeffs;
            dtime[l].first += nc;
            dtime[l].second += t;
            hdr++;
        }
    }
    /*
    double total_from_defer_rows = 0;
    double total_from_defer_cmbs = 0;
    */
    for(unsigned int l = 0 ; l < nsteps ; l++) {
        ASSERT_ALWAYS(hdr->t == SLICE_TYPE_DEFER_ROW);
        double const t = mm->slice_timings[hdr - mm->headers.begin()].t;
        uint64_t const nc = hdr->ncoeffs;
        ctime[l].first += nc;
        ctime[l].second += t;
        // total_from_defer_rows+=t;
        hdr++;
    }
    /* Skip the combining blocks, because they're accounted for already
     * by the row blocks */
    /*
    for( ; hdr != mm->headers.end() && hdr->t == SLICE_TYPE_DEFER_CMB ; hdr++) {
        double t = mm->slice_timings[hdr - mm->headers.begin()].t;
        total_from_defer_cmbs+=t;
    }
    */
    /* Some jitter will appear if transposed mults are performed, because
     * for the moment transposed mults don't properly store timing info
     */
    /*
    fmt::print("jitter {:.2f} - {:.2f} = {:.2f}\n",
            total_from_defer_rows / scale0, total_from_defer_cmbs / scale0,
            (total_from_defer_rows - total_from_defer_cmbs) / scale0);
            */
    for(unsigned int l = 0 ; l < nsteps ; l++) {
        double t;
        uint64_t nc;
        double a;
        nc = dtime[l].first;
        t = dtime[l].second / scale0;
        a = 1.0e9 * t / nc;
        *p_t_total += t;
        os << fmt::format("defer\t{:.2f}s         ; n={:>9d} ; {:5.2f} ns/c ;"
            " scaled*{:.2f} : {:5.2f}/c\n",
            t, nc, a, scale, a * scale);
        nc = ctime[l].first;
        t = ctime[l].second / scale0;
        a = 1.0e9 * t / nc;
        *p_t_total += t;
        os << fmt::format("      + {:.2f}s [{:.2f}s] ; n={:>9d} ; {:5.2f} ns/c ;"
            " scaled*{:.2f} : {:5.2f}/c\n",
            t, *p_t_total, nc, a, scale, a * scale);
    }
    hdr--;
    return os;
}


template<typename Arith>
void matmul_bucket<Arith>::report(double scale)
{
    uint64_t scale0;
    scale0 = (iteration[0] + iteration[1]);

    std::ostringstream os;

    vector<slice_header_t>::iterator hdr;

    os << fmt::format("n {} {:.3f}s/iter (wct of cpu-bound loop)\n",
            ncoeffs,
            main_timing.t / scale0
            );

    double t_total = 0;
    for(hdr = headers.begin() ; hdr != headers.end() ; hdr++) {
        if (hdr->t == SLICE_TYPE_SMALL1_VBLOCK) continue;
        if (hdr->t == SLICE_TYPE_DEFER_ENVELOPE) {
            matmul_bucket_report_vsc(os, this, scale, hdr, &t_total);
            continue;
        }
        double t = slice_timings[hdr - headers.begin()].t;
        uint64_t const nc = hdr->ncoeffs;
        t /= scale0;
        t_total += t;
        double const a = 1.0e9 * t / nc;
        os << fmt::format("{}\t{:.2f}s [{:.2f}s] ; n={:>9d} ; {:5.2f} ns/c ;"
            " scaled*{:.2f} : {:5.2f}/c\n",
            slice_name(hdr->t), t, t_total,
            nc, a, scale, a * scale);
    }
    for(int i = 0 ; i < 40 ; i++) os << '-';
    os << '\n';
    for(int i = 0 ; i < SLICE_TYPE_MAX ; i++) {
        if (i == SLICE_TYPE_SMALL1_VBLOCK) continue;
        double t = 0;
        uint64_t nc = 0;
        for(hdr = headers.begin() ; hdr != headers.end() ; hdr++) {
            if (hdr->t != i) continue;
            nc += hdr->ncoeffs;
            t += slice_timings[hdr - headers.begin()].t;
        }
        if (nc == 0) continue;
        t /= scale0;
        double const a = 1.0e9 * t / nc;
        os << fmt::format("{}\t{:.2f}s ; n={:>9d} ; {:5.2f} ns/c ;"
            " scaled*{:.2f} : {:5.2f}/c\n",
            slice_name(i), t,
            nc, a, scale, a * scale);
    }
    report_string = os.str();
}

template<typename Arith>
void matmul_bucket<Arith>::aux(int op, ...)
{
    va_list ap;
    va_start(ap, op);
    if (op == MATMUL_AUX_ZERO_STATS) {
        matmul_bucket_zero_stats(this);
        iteration[0] = 0;
        iteration[1] = 0;
    }
#if 0
    if (op == MATMUL_AUX_GET_READAHEAD) {
        unsigned int * res = va_arg(ap, unsigned int *);
        ++*res;
    }
#endif
    va_end(ap);
}

// NOLINTNEXTLINE(misc-use-internal-linkage)
matmul_interface * CADO_CONCATENATE4(new_matmul_, ARITH_LAYER, _, MM_IMPL)(
        matmul_public && P,
        arith_generic * arith,
        cxx_param_list & pl,
        int optimized_direction)
{
    return new matmul_bucket<arith_hard>(std::move(P), arith->concrete(), pl, optimized_direction);
}

/* vim: set sw=4: */
