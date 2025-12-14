#include "cado.h" // IWYU pragma: keep

#include <cstdarg>         // for va_list, va_end, va_start
#include <cstddef>         // for ptrdiff_t
#include <cstdint>          // for uint16_t, uint32_t, int32_t, uint64_t
#include <string>           // for basic_string
#include <cstdio>
#include <cstdlib>
#include <cinttypes>
#include <cstring>

#include <vector>
#include <utility>
#include <sstream>
#include <algorithm>
#include <map>
#include <type_traits>  // for std::is_same (C++11)

#include <pthread.h>
#include <gmp.h>

#include "macros.h"
#include "matmul-common.hpp"
#include "matmul.hpp"
#include "matmul_facade.hpp"
#include "memory.h"         // for malloc_aligned
#include "arith-hard.hpp"
#include "params.h"

template<typename Arith> struct fast_default {
    // uncomment the following line (and comment out the next one) to
    // enable the SSE2/AVX2 carry-save code.
    // TODO: revive, maybe?
    // using type = arith_modp::fast_type<Arith>;
    using type = Arith;
};


/* This extension is used to distinguish between several possible
 * implementations of the product */
#define MM_EXTENSION   "-zone"
#define MM_MAGIC_FAMILY        0xb002UL
#define MM_MAGIC_VERSION       0x1001UL
#define MM_MAGIC (MM_MAGIC_FAMILY << 16 | MM_MAGIC_VERSION)

/* This selects the default behaviour as to which is our best code
 * for multiplying. If this flag is 1, then a multiplication matrix times
 * vector (direction==1) performs best if the in-memory structure
 * reflects the non-transposed matrix. Similarly, a vector times matrix
 * multiplication (direction==0) performs best if the in-memory structure
 * reflects the transposed matrix. When the flag is 1, the converse
 * happens.
 * This flag depends on the implementation, and possibly even on the cpu
 * type under certain circumstances.
 */
#define MM_DIR0_PREFERS_TRANSP_MULT   1

using namespace std;

/* Config block */

/* The implementation of the "dispatchers and combiners" here is rotten.
 * I mean to delete it at some point. Maybe replace it with something
 * better. The point is that the current code does not do what it should
 * do, namely write in several buckets simultaneously. So in essence,
 * it's just wasted energy.
 */
#define xxxDISPATCHERS_AND_COMBINERS

/* Size of the row blocks. This impacts both the immediate blocks, and
 * the dispatch/combine blocks */
static size_t rowbatch = 2048;

/* Immediate blocks have this many columns. */
static size_t colbatch0 = 65536;

/* Coefficients whose absolute value is below this bound are stored as
 * repetition of the column index several times.
 */
static int coeff_repeat_bound = 4;

static int debug_print = 0;

#ifdef DISPATCHERS_AND_COMBINERS
/* column indices below col_col_dispatcher_cutoff are treated as immediate
 * blocks.  Setting this cutoff to 1 is
 * sufficient to force the first vertical strip to be processed as
 * immediate blocks only, which is generally something we want
 * because there are so many coefficients there. */
static size_t col_dispatcher_cutoff = 262144;

/* For the dispatcher/combiner split, we treat this number of columns at
 * a time.  */

/* 8k 0.82 cpu @100
 * 16k 0.74 cpu @100
 * 32k 0.71 cpu @100
 * 64k 0.70 cpu @100
 */
static size_t colbatch1 = 65536;
#endif

/* {{{ some boilerplate related to saving cache files */

struct cachefile {
    FILE * f;
    cachefile(FILE * f) : f(f) {}
    template<typename T> struct basic_seq {
        static constexpr bool rigid_element_size = false;
        cachefile& out(cachefile& c, T const * t, size_t n) const {
            for(size_t i = 0 ; i < n ; i++) c << t[i];
            return c;
        }
        cachefile& in(cachefile& c, T * t, size_t n) {
            for(size_t i = 0 ; i < n ; i++) c >> t[i];
            return c;
        }
    };
    template<typename T> struct seq : public basic_seq<T> {
    };
};

template<> struct cachefile::seq<uint16_t> : public cachefile::basic_seq<uint16_t> {
    using T = uint16_t;
    static constexpr bool rigid_element_size = true;
    cachefile& in(cachefile& c, T * t, size_t n) {
        static_assert(sizeof(T) == sizeof(uint16_t), "please fix struct padding");
        MATMUL_COMMON_READ_MANY16(t, n, c.f);
        return c;
    }
    cachefile& out(cachefile& c, T const * t, size_t n) const {
        static_assert(sizeof(T) == sizeof(uint16_t), "please fix struct padding");
        MATMUL_COMMON_WRITE_MANY16(t, n, c.f);
        return c;
    }
};

template<> struct cachefile::seq<uint32_t> : public cachefile::basic_seq<uint32_t> {
    using T = uint32_t;
    static constexpr bool rigid_element_size = true;
    cachefile& in(cachefile& c, T * t, size_t n) {
        static_assert(sizeof(T) == sizeof(uint32_t), "please fix struct padding");
        MATMUL_COMMON_READ_MANY32(t, n, c.f);
        return c;
    }
    cachefile& out(cachefile& c, T const * t, size_t n) const {
        static_assert(sizeof(T) == sizeof(uint32_t), "please fix struct padding");
        MATMUL_COMMON_WRITE_MANY32(t, n, c.f);
        return c;
    }
};

template<> struct cachefile::seq<pair<int, uint32_t>> : public cachefile::basic_seq<pair<int, uint32_t>> {
    using T = pair<int, uint32_t>;
    static constexpr bool rigid_element_size = true;
    cachefile& in(cachefile& c, T * t, size_t n) {
        static_assert(sizeof(T) == 2 * sizeof(uint32_t), "please fix struct padding");
        MATMUL_COMMON_READ_MANY32(t, 2*n, c.f);
        return c;
    }
    cachefile& out(cachefile& c, T const * t, size_t n) const {
        static_assert(sizeof(T) == 2 * sizeof(uint32_t), "please fix struct padding");
        MATMUL_COMMON_WRITE_MANY32(t, 2*n, c.f);
        return c;
    }
};

template<> struct cachefile::seq<pair<uint16_t, int32_t>> : public cachefile::basic_seq<pair<uint16_t, int32_t>> {
#if 0
    using T = pair<uint16_t, int32_t>;
    static constexpr bool rigid_element_size = true;
    cachefile& in(cachefile& c, T * t, size_t n) {
        static_assert(sizeof(T) == 3 * sizeof(uint16_t), "please fix struct padding");
        MATMUL_COMMON_READ_MANY16(t, 3*n, c.f);
        return c;
    }
    cachefile& out(cachefile& c, T const * t, size_t n) const {
        static_assert(sizeof(T) == 3 * sizeof(uint16_t), "please fix struct padding");
        MATMUL_COMMON_WRITE_MANY16(t, 3*n, c.f);
        return c;
    }
#endif
};

template<> struct cachefile::seq<pair<uint16_t, uint16_t>> : public cachefile::basic_seq<pair<uint16_t, uint16_t>> {
    using T = pair<uint16_t, uint16_t>;
    static constexpr bool rigid_element_size = true;
    cachefile& in(cachefile& c, T * t, size_t n) {
        static_assert(sizeof(T) == 2 * sizeof(uint16_t), "please fix struct padding");
        MATMUL_COMMON_READ_MANY16(t, 2*n, c.f);
        return c;
    }
    cachefile& out(cachefile& c, T const * t, size_t n) const {
        static_assert(sizeof(T) == 2 * sizeof(uint16_t), "please fix struct padding");
        MATMUL_COMMON_WRITE_MANY16(t, 2*n, c.f);
        return c;
    }
};

template<> struct cachefile::seq<pair<int, pair<uint32_t, int32_t>>> : public cachefile::basic_seq<pair<int, pair<uint32_t, int32_t>>> {
    using T = pair<int, pair<uint32_t, int32_t>>;
    static constexpr bool rigid_element_size = true;
    cachefile& in(cachefile& c, T * t, size_t n) {
        static_assert(sizeof(T) == 3 * sizeof(uint32_t), "please fix struct padding");
        MATMUL_COMMON_READ_MANY32(t, 3*n, c.f);
        return c;
    }
    cachefile& out(cachefile& c, T const * t, size_t n) const {
        static_assert(sizeof(T) == 3 * sizeof(uint32_t), "please fix struct padding");
        MATMUL_COMMON_WRITE_MANY32(t, 3*n, c.f);
        return c;
    }
};

struct triple_161632 {
    uint16_t first;
    uint16_t second;
    uint32_t third;
    triple_161632() { first = second = third = 0; }
    triple_161632(uint16_t a, uint16_t b, uint32_t c) : first(a), second(b), third(c) {}
};

template<> struct cachefile::seq<triple_161632> : public cachefile::basic_seq<triple_161632> {
    using T = triple_161632;
    static constexpr bool rigid_element_size = true;
    cachefile& in(cachefile& c, T * t, size_t n) {
        static_assert(sizeof(T) == 4 * sizeof(uint16_t), "please fix struct padding");
        MATMUL_COMMON_READ_MANY16(t, 4*n, c.f);
        return c;
    }
    cachefile& out(cachefile& c, T const * t, size_t n) const {
        static_assert(sizeof(T) == 4 * sizeof(uint16_t), "please fix struct padding");
        MATMUL_COMMON_WRITE_MANY16(t, 4*n, c.f);
        return c;
    }
};


static cachefile& operator>>(cachefile & c, uint32_t & x) {
    MATMUL_COMMON_READ_ONE32(x, c.f);
    return c;
}

static cachefile& operator<<(cachefile & c, uint32_t const & x) {
    MATMUL_COMMON_WRITE_ONE32(x, c.f);
    return c;
}

#if 0
static cachefile& operator>>(cachefile & c, int32_t & x) {
    MATMUL_COMMON_READ_ONE32(x, c.f);
    return c;
}

static cachefile& operator<<(cachefile & c, int32_t const & x) {
    MATMUL_COMMON_WRITE_ONE32(x, c.f);
    return c;
}
#endif

static cachefile& operator>>(cachefile & c, uint64_t & x) {
    MATMUL_COMMON_READ_ONE64(x, c.f);
    return c;
}

static cachefile& operator<<(cachefile & c, uint64_t const & x) {
    MATMUL_COMMON_WRITE_ONE64(x, c.f);
    return c;
}

#if 0
static cachefile& operator>>(cachefile & c, int64_t & x) {
    MATMUL_COMMON_READ_ONE64(x, c.f);
    return c;
}

static cachefile& operator<<(cachefile & c, int64_t const & x) {
    MATMUL_COMMON_WRITE_ONE64(x, c.f);
    return c;
}
#endif

template<typename T>
static cachefile& operator>>(cachefile & c, cachefile::seq<T> & s);
template<typename T>
static cachefile& operator<<(cachefile & c, cachefile::seq<T> const & s);

template<typename T>
static cachefile& operator>>(cachefile & c, vector<T>&v)
{
    uint64_t size;
    c >> size;
    ASSERT_ALWAYS(v.empty());
    if (cachefile::seq<T>::rigid_element_size)
        resize_and_check_meaningful(v, size, c.f);
    else
        v.resize(size);
    return cachefile::seq<T>().in(c, v.data(), v.size());
}
template<typename T>
static cachefile& operator<<(cachefile & c, vector<T> const &v)
{
    c << (uint64_t) v.size();
    return cachefile::seq<T>().out(c, v.data(), v.size());
}
template<typename T>
static inline cachefile& operator>>(cachefile& c, T& z) { return z.cachefile_load(c); }
template<typename T>
static inline cachefile& operator<<(cachefile& c, T const & z) { return z.cachefile_save(c); }
/* }}} */

struct placed_block {/*{{{*/
    unsigned int i0, j0;
    placed_block() { i0 = j0 = 0; }
    placed_block(unsigned int i0, unsigned int j0) : i0(i0), j0(j0) {}
    struct rowmajor_sorter {/*{{{*/
        bool operator()(placed_block const& a, placed_block const& b) const {
            return a.i0 < b.i0 || (a.i0 == b.i0 && a.j0 < b.j0);
        }
    };/*}}}*/
    struct colmajor_sorter {/*{{{*/
        bool operator()(placed_block const& a, placed_block const& b) const {
            return a.j0 < b.j0 || (a.j0 == b.j0 && a.i0 < b.i0);
        }
    };/*}}}*/
    /* NOTE that the cachefile_load and _save functions below *must* be
     * overloaded in children classes, otherwise we'll end up saving only
     * *OUR* data, that's it.
     */
    cachefile& cachefile_load(cachefile& c) {/*{{{*/
        return c >> i0 >> j0;
    }/*}}}*/
    cachefile& cachefile_save(cachefile& c) const {/*{{{*/
        return c << i0 << j0;
    }/*}}}*/
};
/*}}}*/

template<typename Arith, typename fast_gfp = typename fast_default<Arith>::type >
class zone : public placed_block { /* {{{ (immediate zones) */
    using elt = typename Arith::elt;
    using fast_elt = typename fast_gfp::elt;
    using fast_elt_ur_for_add = typename fast_gfp::elt_ur_for_add;
public:
    using qpm_t = vector<pair<uint16_t, uint16_t>>;
    using qg_t = vector<triple_161632>;
    qpm_t qp, qm;
    qg_t qg;
    zone() {}
    zone(unsigned int i0, unsigned int j0) : placed_block(i0, j0) {}
    bool empty() const { return qp.empty() && qm.empty() && qg.empty(); }
    size_t size() const { return qp.size() + qm.size() + qg.size(); }
    void mul(Arith const *, fast_elt_ur_for_add *, const fast_elt *) const;
    void tmul(Arith const *, fast_elt_ur_for_add *, const fast_elt *) const;

    struct sort_qpm {
        bool operator()(qpm_t::value_type const& a, qpm_t::value_type const& b) const {
            return a.second < b.second;
        }
    };

    struct sort_qg {
        bool operator()(qg_t::value_type const& a, qg_t::value_type const& b) const {
            return a.second < b.second;
        }
    };

    void sort() {
        std::ranges::sort(qp, sort_qpm());
        std::ranges::sort(qm, sort_qpm());
        std::ranges::sort(qg, sort_qg());
    }
    cachefile& cachefile_load(cachefile& c) {/*{{{*/
        return c >> (placed_block&) *this >> qp >> qm >> qg;
    }/*}}}*/
    cachefile& cachefile_save(cachefile& c) const {/*{{{*/
        return c << (placed_block const&) *this << qp << qm << qg;
    }/*}}}*/
};
/*}}}*/
#ifdef DISPATCHERS_AND_COMBINERS
struct dispatcher : public vector<uint16_t>, public placed_block {/*{{{*/
    using super = vector<uint16_t>;
    /* a dispatcher contains: (col id)* */
    cachefile& cachefile_load(cachefile& c) {/*{{{*/
        return c >> (placed_block&)*this >> (super&)*this;
    }/*}}}*/
    cachefile& cachefile_save(cachefile& c) const {/*{{{*/
        return c << (placed_block const&)*this << (super const&)*this;
    }/*}}}*/
};
/*}}}*/
struct combiner : public placed_block {/*{{{*/
    /* which buffer will we read from, and how many values, from
     * where ?  */
    size_t index, offset, count;
    /* a combiner contains: (dest row id)* (index range is duplicated, in
     * order to account for the fact that we offset the index by a full
     * row batch to account for negative coefficients) */
    using main_value_type = uint16_t;
    vector<main_value_type> main;
    vector<pair<uint16_t, int32_t>> aux;
    cachefile& cachefile_load(cachefile& c) {/*{{{*/
        return c >> (placed_block&)*this
            >> index >> offset >> count
            >> main >> aux;
    }/*}}}*/
    cachefile& cachefile_save(cachefile& c) const {/*{{{*/
        return c << (placed_block const&)*this
            << index << offset << count
            << main << aux;
    }/*}}}*/
};
/*}}}*/
/* {{{ Temporary buffers which are used when reading the dispatchers.  */
template<typename Arith, typename fast_gfp = typename fast_default<Arith>::type >
class temp_buffer {
    using fast_elt = typename fast_gfp::elt;
public:
    unsigned int j0;
    using vec_type = vector<fast_elt, aligned_allocator<fast_elt, fast_elt::alignment> >;  
    vec_type v;
    vec_type::iterator ptr;
    temp_buffer(unsigned int j0, size_t n) : j0(j0), v(n), ptr(v.begin()) {}
    void rewind() { ptr = v.begin(); }
};
/* }}} */
#endif

template<typename Arith, typename fast_gfp = typename fast_default<Arith>::type >
struct block_of_rows : public placed_block {/*{{{*/
    /* A block of rows contains:
     *
     *  - immediate zones, presumably for leftmost blocks.
     *  - dispatcher blocks
     *  - combiner blocks
     *
     * It is processed in the order found in matmul_zone::mul
     */
    vector<zone<Arith, fast_gfp> > Z;
#ifdef DISPATCHERS_AND_COMBINERS
    vector<dispatcher> D;       /* can be empty */
    vector<combiner> C;
#endif
    block_of_rows() {}
    block_of_rows(unsigned int i0) : placed_block(i0, 0) {}
    cachefile& cachefile_load(cachefile& c) {/*{{{*/
        c >> (placed_block&)*this >> Z;
#ifdef DISPATCHERS_AND_COMBINERS
        c >> D >> C;
#endif
        return c;
    }/*}}}*/
    cachefile& cachefile_save(cachefile& c) const {/*{{{*/
        c << (placed_block const&)*this << Z;
#ifdef DISPATCHERS_AND_COMBINERS
        c << D << C;
#endif
        return c;
    }/*}}}*/
};/*}}}*/

template<typename Arith, typename fast_gfp = typename fast_default<Arith>::type >
class matmul_zone : public matmul_interface {/*{{{*/
    using elt = typename Arith::elt;
    using fast_elt = typename fast_gfp::elt;
    using fast_elt_ur_for_add = typename fast_gfp::elt_ur_for_add;
public:
    /* now our private fields */
    Arith * xab;

    vector<block_of_rows<Arith, fast_gfp> > blocks;

    static_assert(std::is_same<elt, fast_elt>::value, "if alternate representation of the source data is needed, please amend this code");
    // vector<fast_elt, aligned_allocator<fast_elt, fast_elt::alignment> > alternate[2];

#ifdef DISPATCHERS_AND_COMBINERS
    size_t maxmaxw = 0;
#endif

    /* {{{ timing data */
    struct twn {
        uint64_t tt;
        size_t w;
        unsigned int n;
        // twn() : tt(0), w(0), n(0) {}
    };

    using tmap_t = map<unsigned int, twn>;
    tmap_t tmap;
    /* }}} */

    void build_cache(matrix_u32 &&) override;
    int reload_cache_private() override;
    void save_cache_private() override;
    void mul(void *, const void *, int) override;
    void report(double) override;
    ~matmul_zone() override = default;

    matmul_zone(matmul_public &&, arith_concrete_base *, cxx_param_list &, int);

    matmul_zone(matmul_zone const &) = delete;
    matmul_zone& operator=(matmul_zone const &) = delete;
    matmul_zone(matmul_zone &&) noexcept = default;
    matmul_zone& operator=(matmul_zone &&) noexcept = default;
};
/*}}}*/
/**************************************************************************/


template<typename Arith, typename fast_gfp>
matmul_zone<Arith, fast_gfp>::matmul_zone(matmul_public && P, arith_concrete_base * pxx, cxx_param_list & pl, int optimized_direction) /*{{{*/
    : matmul_interface(std::move(P))
    , xab((Arith *) pxx) // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
{
    int const suggest = optimized_direction ^ MM_DIR0_PREFERS_TRANSP_MULT;
    store_transposed = suggest;
    param_list_parse(pl, "mm_store_transposed", store_transposed);
    if (store_transposed != suggest) {
        fprintf(stderr, "Warning, mm_store_transposed"
                " overrides suggested matrix storage ordering\n");
    }
}
/*}}}*/

struct coeff_stats {/*{{{*/
    int bound;
    vector<uint64_t> ccount;
    coeff_stats(int c) : bound(c), ccount(2*c+1, 0) {}
    void operator()(int c) {
        if (c < 0 && c >= -bound) {
            ccount[bound + c]++;
        } else if (c > 0 && c <= bound) {
            ccount[bound + c]++;
        } else {
            ccount[0]++;
        }
    }
    void report(std::ostream& o, size_t nrows) {
        o << "coeff stats per row";
        int maxnz = bound;
        for( ; maxnz > 1 ; maxnz--)
            if (ccount[bound + maxnz] + ccount[bound - maxnz]) break;
        for(int i = 1 ; i <= maxnz ; i++) {
            o << " ±" << i << ":" << (double) (ccount[bound + i] + ccount[bound - i]) / nrows;
        }
        if (ccount[bound]) o << " other:" << (double) ccount[bound] / nrows;   
        o << "\n";
    }
};
/*}}}*/
struct sort_jc {/*{{{*/
    bool operator()(pair<uint32_t, int32_t> const& a, pair<uint32_t,int32_t> const& b) const {
        return a.first < b.first;
    }
};
/*}}}*/

#ifdef DISPATCHERS_AND_COMBINERS
void merge_dispatchers(vector<dispatcher>& all, size_t maxmaxw)/*{{{*/
{
    using dispatcher = dispatcher;
    vector<dispatcher> merged;
    for(size_t k0 = 0, k1; k0 < all.size() ; ) {
        /* recompute maxw, since we haven't kept track */
        size_t maxw = 0;
        unsigned int j0 = all[k0].j0;
        for(k1 = k0 ; k1 < all.size() && all[k1].j0 == all[k0].j0 ; k1++) {
            maxw = max(maxw, all[k1].size());
        }
        /* the range [k0..k1[ has dispatchers for the same set of
         * columns. */
        unsigned int group = maxmaxw / maxw;
        unsigned int nm = 0;
        size_t old_k0 = k0;
        for( ; k0 < k1 ; ) {
            dispatcher D;
            D.i0 = all[k0].i0;
            D.j0 = all[k0].j0;
            /* concatenate. */
            for(unsigned int k2 = 0 ; k2 < group && k0 < k1 ; k2++, k0++) {
                D.insert(D.end(), all[k0].begin(), all[k0].end());
                /* free this memory as early as we can */
                all[k0].clear();
                /* XXX at this point, combiner must get the info that it
                 * should read in the buffer from index D.size()
                 */
            }
            merged.push_back(std::move(D));
            nm++;
        }
        if (debug_print)
            printf("j0=%u: maxw = %zu ; group %u together. Split %zd dispatchers [%zu..%zu[ into %u dispatchers\n", j0, maxw, group, k1 - old_k0, old_k0, k1, nm);
    }
    if (debug_print)
        printf("We now have %zu dispatchers, combined from %zu original\n",
                merged.size(), all.size());
    all.swap(merged);
}
/*}}}*/

pair<dispatcher, combiner> create_dispatcher_and_combiner(zone const& q)/*{{{*/
{
    dispatcher D;
    combiner C;
    D.i0 = C.i0 = q.i0;
    D.j0 = C.j0 = q.j0;
    for(auto const& ij : q.qp) {
        uint16_t destrow_id = ij.first;
        uint16_t col_id = ij.second;
        D.push_back(col_id);
        C.main.push_back(destrow_id);
    }
    for(auto const& ij : q.qm) {
        uint16_t destrow_id = ij.first;
        uint16_t col_id = ij.second;
        D.push_back(col_id);
        /* specific offset for negative coefficients */
        C.main.push_back(destrow_id + rowbatch);
    }
    for(auto const& ijc : q.qg) {
        uint16_t col_id = ijc.second;
        uint16_t destrow_id = ijc.first;
        int32_t coeff = ijc.third;
        D.push_back(col_id);
        C.aux.emplace_back(destrow_id, coeff);
    }
    return make_pair(D, C);
}/*}}}*/
#endif

template<typename Arith, typename fast_gfp>
void matmul_zone<Arith, fast_gfp>::build_cache(matrix_u32 && m)/*{{{*/
{
    ASSERT_ALWAYS(!m.p.empty());

    unsigned int const nrows_t = dim[ store_transposed];
    unsigned int const ncols_t = dim[!store_transposed];

    uint32_t * ptr = m.p.data();
    size_t size = m.p.size();

    /* count coefficients */
    ncoeffs = 0;

    double zavg = 0;

#ifdef DISPATCHERS_AND_COMBINERS
    /* immediate zones as well as combiners are stored right inside the
     * block_of_rows structures. As for dispatchers, we'll put them there
     * in a second pass.
     */
    vector<dispatcher> D;
#endif

    size_t n_immediate = 0;

    coeff_stats cstats(coeff_repeat_bound);

    unsigned int maxrow = 0;

    for(unsigned int i0 = 0 ; i0 < nrows_t ; i0 += rowbatch) {
        /* Create the data structures for the horizontal strip starting
         * at row i0, column 0.
         */
        block_of_rows<Arith, fast_gfp> B(i0);

        /* Because this horizontal strip will be split in many blocks, we
         * need to have a batch of pointers for reading each row. */
        uint32_t * pp[rowbatch + 1];
        uint32_t * cc[rowbatch + 1];
        pp[0] = ptr;
        for(unsigned int k = 0 ; k < rowbatch ; k++) {
            cc[k] = pp[k] + 1;

            if (pp[k] == m.p.data() + size) {
                /* reached the end of our data stream. We have empty
                 * padding rows, we must treat them accordingly */
                pp[k+1] = pp[k];
            } else if (i0 + k >= nrows_t) {
                pp[k+1] = pp[k];
            } else {
                maxrow = i0 + k + 1;
                ASSERT_ALWAYS((pp[k] - m.p.data()) < (ptrdiff_t) size);
                uint32_t const weight = *pp[k];
                pp[k+1] = pp[k] + 1 + 2*weight;
                ncoeffs += weight;
                /* This is very important. We must sort rows before
                 * processing. */
                auto * cb = (pair<uint32_t, int32_t> *) cc[k];
                pair<uint32_t, int32_t> * ce = cb + weight;
                sort(cb, ce, sort_jc());
            }
        }
        ptr = pp[rowbatch];
        ASSERT_ALWAYS((ptr - m.p.data()) <= (ptrdiff_t) size);
        for(unsigned int j0 = 0, colbatch ; j0 < ncols_t ; j0 += colbatch) {
#ifdef DISPATCHERS_AND_COMBINERS
            colbatch = (j0 < col_dispatcher_cutoff) ? colbatch0 : colbatch1;
#else
            colbatch = colbatch0;
#endif
            zone<Arith, fast_gfp> z(i0, j0);
            for(unsigned int k = 0 ; k < rowbatch ; k++) {
                for( ; cc[k] < pp[k+1] ; cc[k] += 2) {
                    uint32_t const j = cc[k][0] - j0;
                    auto c = int32_t(cc[k][1]);
                    if (j >= colbatch) break;
                    cstats(c);
                    if (c < 0 && c >= -coeff_repeat_bound) {
                        for( ; c++ ; ) {
                            z.qm.emplace_back(k, j);
                        }
                    } else if (c > 0 && c <= coeff_repeat_bound) {
                        for( ; c-- ; ) {
                            z.qp.emplace_back(k, j);
                        }
                    } else {
                        z.qg.push_back(triple_161632(k, j, c));
                    }
                }
            }
            z.sort();
            // printf("Zone %zu @[%u,%u]: %zu+%zu+%zu\n", B.Z.size(), z.i0, z.j0, z.qp.size(), z.qm.size(), z.qg.size());
            if (z.empty()) continue;

#ifdef DISPATCHERS_AND_COMBINERS
            if (j0 >= col_dispatcher_cutoff) {
                /* blocks which are not in the first strip go for
                 * the dispatcher/combiner split. */
                if (z.size() > maxmaxw)
                    maxmaxw = z.size();
                auto DC = create_dispatcher_and_combiner(z);
                D.push_back(std::move(DC.first));
                B.C.push_back(std::move(DC.second));
            } else
#endif
            {
                zavg += z.size();
                B.Z.push_back(std::move(z));
                n_immediate++;
            }
        }
        blocks.push_back(std::move(B));
    }
    ASSERT_ALWAYS(maxrow <= nrows_t);
    ASSERT_ALWAYS(ptr - m.p.data() == (ptrdiff_t) (maxrow + 2 * ncoeffs));
#ifdef DISPATCHERS_AND_COMBINERS
    /* We now merge dispatchers */
    sort(D.begin(), D.end(), dispatcher::colmajor_sorter());
    merge_dispatchers(D, maxmaxw);
    sort(D.begin(), D.end(), dispatcher::rowmajor_sorter());
    /* And put the dispatchers to the proper place */
    for(size_t i = 0, j = 0; i < D.size() ; i++) {
        for ( ; D[i].i0 > blocks[j].i0 ; j++) ;
        blocks[j].D.push_back(std::move(D[i]));
    }
#endif

    ostringstream os;
    if (debug_print) {
        cstats.report(os, nrows_t);
        printf("Stats: [%" PRIu64 " coeffs]:"
                " %zu immediate zones of average weight %.1f\n%s",
                ncoeffs,
                n_immediate, zavg / n_immediate,
                os.str().c_str());
    }
}
/*}}}*/
/* cache load and save {{{ */
template<typename Arith, typename fast_gfp>
int matmul_zone<Arith, fast_gfp>::reload_cache_private()
{
    auto f = matmul_common_reload_cache_fopen(sizeof(typename Arith::elt), *this, MM_MAGIC);
    if (!f) return 0;
    cachefile c(f.get());
    c >> blocks;
#ifdef DISPATCHERS_AND_COMBINERS
    maxmaxw = 0;
    for(auto const& B : blocks) {
        for(auto const& D : B.D) {
            maxmaxw = std::max(D.size(), maxmaxw);
        }
    }
#endif

    return 1;
}

template<typename Arith, typename fast_gfp>
void matmul_zone<Arith, fast_gfp>::save_cache_private()
{
    auto f = matmul_common_save_cache_fopen(sizeof(typename Arith::elt), *this, MM_MAGIC);
    if (!f) return;
    cachefile c(f.get());
    c << blocks;
}
/* }}} */

#if 0
extern "C" {
extern void gfp3_dispatch_add(void * tdst, const void * tsrc, const void * p, size_t size);
extern void gfp3_dispatch_sub(void * tdst, const void * tsrc, const void * p, size_t size);
}

void __attribute__((noinline)) gfp3_dispatch_add (void * tdst, const void * tsrc, const void * p, size_t size)
{
    Arith::elt_ur_for_add * tdst0 = (Arith::elt_ur_for_add *) tdst; 
    const Arith::elt * tsrc0 = (const Arith::elt *) tsrc; 
    const zone::qpm_t::value_type * q = (const zone::qpm_t::value_type *) p;
    for( ; size-- ; q++) {
        x->add(tdst0[q->first], tsrc0[q->second]);
    }
}

void __attribute__((noinline)) gfp3_dispatch_sub (void * tdst, const void * tsrc, const void * p, size_t size)
{
    Arith::elt_ur_for_add * tdst0 = (Arith::elt_ur_for_add *) tdst; 
    const Arith::elt * tsrc0 = (const Arith::elt *) tsrc; 
    const zone::qpm_t::value_type * q = (const zone::qpm_t::value_type *) p;
    for( ; size-- ; q++) {
        x->sub(tdst0[q->first], tsrc0[q->second]);
    }
}
#endif
#if 0
extern "C" {
extern void gfp3_dispatch_add(void * tdst, const void * tsrc, const void * p, size_t size);
extern void gfp3_dispatch_sub(void * tdst, const void * tsrc, const void * p, size_t size);
}

void __attribute__((noinline)) gfp3_dispatch_add (void * tdst, const void * tsrc, const void * p, size_t size)
{
    fast_elt_ur_for_add * tdst0 = (fast_elt_ur_for_add *) tdst; 
    const fast_elt * tsrc0 = (const fast_elt *) tsrc; 
    const zone::qpm_t::value_type * q = (const zone::qpm_t::value_type *) p;
    for( ; size-- ; q++) {
        x->add(tdst0[q->first], tsrc0[q->second]);
    }
}

void __attribute__((noinline)) gfp3_dispatch_sub (void * tdst, const void * tsrc, const void * p, size_t size)
{
    fast_elt_ur_for_add * tdst0 = (fast_elt_ur_for_add *) tdst; 
    const fast_elt * tsrc0 = (const fast_elt *) tsrc; 
    const zone::qpm_t::value_type * q = (const zone::qpm_t::value_type *) p;
    for( ; size-- ; q++) {
        x->sub(tdst0[q->first], tsrc0[q->second]);
    }
}
#endif

template<typename Arith, typename fast_gfp>
void zone<Arith, fast_gfp>::mul(Arith const * x, fast_elt_ur_for_add * tdst, const fast_elt * tsrc)
        const
{
#if 1
    for(auto const& ij : qp) 
        x->add(x->vec_item(tdst, ij.first), x->vec_item(tsrc, ij.second));
    for(auto const& ij : qm) 
        x->sub(x->vec_item(tdst, ij.first), x->vec_item(tsrc, ij.second));
#else
    gfp3_dispatch_add((void*)tdst, (const void*)tsrc, (const void*)(&*qp.begin()), qp.size());
    gfp3_dispatch_sub((void*)tdst, (const void*)tsrc, (const void*)(&*qm.begin()), qm.size());
#endif
    for(auto const& ijc : qg) {
        uint16_t const i = ijc.first;
        uint16_t const j = ijc.second;
        int32_t const c = int32_t(ijc.third);
        if (c>0) {
            x->addmul_ui(x->vec_item(tdst, i), x->vec_item(tsrc, j), c);
        } else {
            x->submul_ui(x->vec_item(tdst, i), x->vec_item(tsrc, j), -c);
        }
    }
}

template<typename Arith, typename fast_gfp>
void zone<Arith, fast_gfp>::tmul(Arith const * x, fast_elt_ur_for_add * tdst, const fast_elt * tsrc) const
{
#if 1
    for(auto const& ij : qp) 
        x->add(x->vec_item(tdst, ij.second), x->vec_item(tsrc, ij.first));
    for(auto const& ij : qm) 
        x->sub(x->vec_item(tdst, ij.second), x->vec_item(tsrc, ij.first));
#else
    gfp3_dispatch_add((void*)tdst, (const void*)tsrc, (const void*)(&*qp.begin()), qp.size());
    gfp3_dispatch_sub((void*)tdst, (const void*)tsrc, (const void*)(&*qm.begin()), qm.size());
#endif
    for(auto const& ijc : qg) {
        uint16_t const i = ijc.first;
        uint16_t const j = ijc.second;
        auto const c = int32_t(ijc.third);
        if (c>0) {
            x->addmul_ui(x->vec_item(tdst, j), x->vec_item(tsrc, i), c);
        } else {
            x->submul_ui(x->vec_item(tdst, j), x->vec_item(tsrc, i), -c);
        }
    }
}

template<typename Arith, typename fast_gfp>
void matmul_zone<Arith, fast_gfp>::mul(void * xdst, void const * xsrc, int d)
{
    ASM_COMMENT("multiplication code");
    Arith * x = xab;

    const fast_elt * src;
    fast_elt * dst;

    // size_t nsrc = dim[d];
    size_t const ndst = dim[!d];

    /* ef7268528 had some remains of support for changing the source
     * vector to alternate representation. However, this does
     * not play well with adding the pz layer, so we take this feature
     * out for the moment. (updating would be doable, but since we have
     * zero coverage for it, who cares...)
     */
    static_assert(std::is_same<elt, fast_elt>::value, "if alternate representation of the source data is needed, please amend this code");
    src = (const fast_elt *) xsrc;
    dst = (fast_elt *) xdst;
#if 0
    if (std::is_same<elt, fast_elt>::value) {
        src = (const fast_elt *) xsrc;
        dst = (fast_elt *) xdst;
    } else {
        for(int j = 0 ; j < 2 ; j++) {
            if (alternate[j].empty())
                alternate[j].assign(dim[j], fast_elt());
            ASSERT_ALWAYS(alternate[j].size() == dim[j]);
        }
        /* we read items in xsrc exactly as they are, which is Arith::elt's.
         * And because those are convertible to fast_gfp::elt's, we'll
         * get our vector.
         */
        const elt * begin = (const elt *) xsrc;
        const elt * end = begin + nsrc;
        alternate[d].assign(begin, end);
        src = alternate[d].data();
        dst = alternate[!d].data();
    }
#endif


    /* d == 1: matrix times vector product */
    /* d == 0: vector times matrix product */

    /* However the matrix may be stored either row-major
     * (store_transposed == 0) or column-major (store_transposed == 1)
     */

    x->vec_set_zero(dst, ndst);

    /* Processing order:
     *
     *  - vector data corresponding to the combiner blocks is
     *  read from the buffers, and stored to a buffer which stores a
     *  write window to the destination
     *  - immediate blocks are read an applied
     *  - the write window to the destination undergoes reduction, and
     *  final store.
     *  - dispatcher blocks (if any) are processed, thereby refilling
     *  buffers which were emptied when processing their last combiner.
     */

#ifdef DISPATCHERS_AND_COMBINERS
    vector<temp_buffer<Arith, fast_gfp> > buffers;
    unsigned int ncols_t = dim[!store_transposed];

    /* replay the sequence of starting column indices */
    for(unsigned int j0 = 0, colbatch ; j0 < ncols_t ; j0 += colbatch) {
        if (j0 < col_dispatcher_cutoff) {
            colbatch = colbatch0;
        } else {
            colbatch = colbatch1;
            buffers.push_back(std::move(temp_buffer(j0, maxmaxw)));
        }
    }
#endif

    /* TODO: missing in mpfq elt_ur_{add,sub}_elt */
    if (d == !store_transposed) {
#ifdef DISPATCHERS_AND_COMBINERS
        /* Doubling rowbatch is because we do a nasty trick with the
         * indices for the negative coefficients. Oddly enough, we don't
         * seem to do so currently for the immediate zones...
         */
        auto * tdst = x->template alloc<fast_elt_ur_for_add>(2*rowbatch);
#else
        auto * tdst = x->template alloc<fast_elt_ur_for_add>(rowbatch);
#endif

        ASM_COMMENT("critical loop");

        for(auto const& B : blocks) {
#ifdef DISPATCHERS_AND_COMBINERS
            x->vec_set_zero(tdst, 2 * rowbatch);
#else
            x->vec_set_zero(tdst, rowbatch);
#endif

#ifdef DISPATCHERS_AND_COMBINERS
            /* loop through all dispatcher strips, and pre-fill buffers if
             * we happen to need it. */
            for(size_t j = 0, k = 0; j < B.D.size() ; j++) {
                auto const& d(B.D[j]);
                ASSERT_ALWAYS(d.i0 == B.i0);
                for( ; buffers[k].j0 < d.j0 ; k++) ;

                /* Fill that buffer now ! */
                auto ptr = buffers[k].v.begin();
                for(auto const& x : d) {
                    /* It should use a vec_subvec instead of
                     * pointer/iterator arithmetic. And anyway we must
                     * not use vectors for the temp buffers!!
                     */
                    static_assert(0,"FIXME");
                    fast_gfp::stream_store(&*ptr++, x->vec_item(src, d.j0 + x));
                    // *ptr++ = x->vec_item(src, d.j0 + x);
                }
                buffers[k].rewind();
            }
#endif

            /* process immediate zones */
            for(auto const & z : B.Z) {
                z.mul(x, tdst, x->vec_subvec(src, z.j0));
            }

#ifdef DISPATCHERS_AND_COMBINERS
            /* Read all combiner data, and apply it */
            for(size_t j = 0, k = 0 ; j < B.C.size() ; j++, k++) {
                for( ; buffers[k].j0 < B.C[j].j0 ; k++) ;
                for(auto const& i : B.C[j].main) {
                    ASSERT(buffers[k].ptr < buffers[k].v.end());
                    x->add(x->vec_item(tdst, i), *(Arith::elt *)&*buffers[k].ptr++);
                }
                /* Almost surely this loop will never run */
                for(auto const& xc : B.C[j].aux) {
                    auto i = xc.first;
                    auto c = xc.second;
                    if (c > 0) {
                        x->addmul_ui(x->vec_item(tdst, i), *buffers[k].ptr++, c);
                    } else {
                        x->submul_ui(x->vec_item(tdst, i), *buffers[k].ptr++, -c);
                    }
                }
            }
#endif

            /* reduce last batch. It could possibly be incomplete */
            size_t const active = std::min(rowbatch, (size_t) (ndst - B.i0));
            for(size_t i = 0 ; i < active ; i++) {
#ifdef DISPATCHERS_AND_COMBINERS
                fast_gfp::sub_ur(x->vec_item(tdst, i), x->vec_item(tdst, i + rowbatch));
#endif
                x->reduce(x->vec_item(dst, B.i0 + i), x->vec_item(tdst, i));
            }
        }

        ASM_COMMENT("end of critical loop");

#ifdef DISPATCHERS_AND_COMBINERS
        // x->template free<fast_elt_ur_for_add>(tdst, 2*rowbatch);
#else
        // x->template free<fast_elt_ur_for_add>(tdst, rowbatch);
#endif
        x->template free<fast_elt_ur_for_add>(tdst);
    } else {
#ifdef DISPATCHERS_AND_COMBINERS
        fprintf(stderr, "transposed product not yet implemented for the dispatch-combine trick");
        ASSERT_ALWAYS(0);
#endif
        auto * tdst = x->template alloc<fast_elt_ur_for_add>(ndst);
        x->vec_set_zero(tdst, ndst);

        ASM_COMMENT("critical loop");

        for(auto const& B : blocks) {
            /* process immediate zones */
            for(auto const & z : B.Z) {
                z.tmul(x, x->vec_subvec(tdst, z.j0), x->vec_subvec(src, B.i0));
            }
        }
        for(size_t j = 0 ; j < ndst ; j++) {
            x->reduce(x->vec_item(dst, j), x->vec_item(tdst, j));
        }

        ASM_COMMENT("end of critical loop");
        x->template free<fast_elt_ur_for_add>(tdst);
    }

    static_assert(std::is_same_v<elt, fast_elt>, "if alternate representation of the source data is needed, please amend this code");
    /*
    if (!std::is_same<elt, fast_elt>::value) {
        std::ranges::copy(alternate[!d], (elt*) xdst);
    }
    */

    ASM_COMMENT("end of multiplication code");

    iteration[d]++;
}

template<typename Arith, typename fast_gfp>
void matmul_zone<Arith, fast_gfp>::report(double scale MAYBE_UNUSED)
{
    static pthread_mutex_t lk = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&lk);
    unsigned int const niter = iteration[0] + iteration[1];
    for(auto const & it : tmap) {
        printf("j0=%u [%u zones]: avg %.1f cycles/c [%.1f coeffs avg] - %.1f Mcycles/iter\n",
                it.first, it.second.n / niter, it.second.tt / scale / it.second.w, (double) it.second.w / it.second.n, (double) it.second.tt / niter * 1.0e-6);
    }
    pthread_mutex_unlock(&lk);
}

// NOLINTNEXTLINE(misc-use-internal-linkage)
matmul_interface * CADO_CONCATENATE4(new_matmul_, ARITH_LAYER, _, MM_IMPL)(
        matmul_public && P,
        arith_generic * arith,
        cxx_param_list & pl,
        int optimized_direction)
{
    return new matmul_zone<arith_hard>(std::move(P), arith->concrete(), pl, optimized_direction);
}

/* vim: set sw=4: */
