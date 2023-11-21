#ifndef LINGEN_MATPOLY_BINARY_HPP_
#define LINGEN_MATPOLY_BINARY_HPP_

/* The outer interface of the matpoly_binary type is exactly the same as
 * for the matpoly type. This is enforced by the test program. In
 * particular, we copy even the "arith_hard" argument that is passed
 * everywhere, with the slight catch that the code *here* does not use it
 * at all.
 *
 *
 * All the "arith_hard" stuff in here is replaced by (nearly) empty proxies.
 */

#include <cstdlib>
#include <gmp.h>
#include <algorithm>
#include "cado_config.h"
#include "macros.h"
#include "lingen_memory_pool.hpp"
#include "submatrix_range.hpp"
// #include "arith-hard.hpp"
#include "tree_stats.hpp"
#include "lingen_call_companion.hpp"
#include "cxx_mpz.hpp"

class matpoly {
    friend class bigmatpoly;

    // static_assert(arith_hard::elt::compatible_with<unsigned long>::value, "This file requires that the arithmetic layer be based on unsigned long");

public:
    /*
    typedef arith_hard::elt elt;
    typedef arith_hard::elt * ptr;
    typedef arith_hard::elt const * srcptr;
    */
    typedef unsigned long elt;
    typedef unsigned long * ptr;
    typedef unsigned long const * srcptr;

    struct arith_hard {
        static size_t elt_stride() { return sizeof(elt); }
        static mpz_srcptr characteristic() {
            static mp_limb_t limbs[1] = {2};
            static __mpz_struct a = { 1, 1, limbs };
            return &a;
        }
        arith_hard(cxx_mpz const &, unsigned int) {}
        static bool is_zero(elt x) { return x == 0; }
        public:
        /* Note that here, k is in numbers of bits, and it must always be
         * a multiple of ULONG_BITS
         */
        static elt * vec_subvec(ptr p, size_t k) {
            return p + k/ULONG_BITS;
        }
        static elt const * vec_subvec(srcptr p, size_t k) {
            return p + k/ULONG_BITS;
        }
        static void vec_add(elt * u, elt const * x, size_t n) {
            for(size_t i = 0 ; i < n/ULONG_BITS ; i++) {
                u[i] ^= x[i];
            }
        }
        static void vec_set(elt * u, elt const * x, size_t n) {
            std::copy_n(x, n/ULONG_BITS, u);
        }
    };

    /* it's only needed when we have code, such as in lingen_matpoly_ft,
     * which wants to pass matpoly init data from earlier stuff */
    static constexpr arith_hard * ab = nullptr;

private:
    typedef memory_pool_wrapper<ptr, true> memory_pool_type;
    static memory_pool_type memory;

public:
    struct memory_guard : private memory_pool_type::guard_base {
        memory_guard(size_t s) : memory_pool_type::guard_base(memory, s) {}
        ~memory_guard() { memory_pool_type::guard_base::pre_dtor(memory); }
    };

    static constexpr bool over_gf2 = true;
    // static void add_to_main_memory_pool(size_t s);
    // arith_hard * ab = NULL;
    /* Note that in a sense, we're perverting the b64 layer anyway, since
     * what we're really doing is _NOT_ simd at all.
     */
    // static_assert(arith_hard::simd_groupsize == 64, "This code only works with the b64 layer");
    unsigned int m = 0;
    unsigned int n = 0;
    /* alloc_words is the number of unsigned longs used to store each
     * coefficient */
private:
    size_t size = 0;    /* in bits */
    size_t alloc_words = 0;
    ptr x = NULL;
#define BITS_TO_WORDS(B,W)      iceildiv((B),(W))
    static inline size_t b2w_x(size_t n) {
        /* We always use an even number of words. It seems stupid, but
         * some of the routines that play an important role in
         * block-based lingen basecase really want 64*64 matrices of
         * 64-bit word-based polynomials. In a way it's a shortcoming,
         * yes.
         */
        static_assert(64 % ULONG_BITS == 0, "ULONG_BITS must divide 64");
        return BITS_TO_WORDS(n, 64) * (64 / ULONG_BITS);
    }
    static inline size_t b2w(size_t n) {
        return BITS_TO_WORDS(n, ULONG_BITS);
    }
    /*{{{*/
    // inline size_t colstride() const { return nrows() * stride(); }/*}}}*/
public:
    inline size_t capacity() const { return alloc_words * ULONG_BITS; }
    const void * data_area() const { return x; }
    bool is_tight() const { return alloc_words == b2w(size); }
    size_t data_entry_size_in_bytes() const {
        return b2w(size) * sizeof(unsigned long);
    }
    size_t data_entry_size_in_words() const {
        return b2w(size);
    }
    size_t data_size_in_bytes() const {
        return m * n * data_entry_size_in_bytes();
    }
    size_t data_entry_alloc_size_in_words() const {
        return alloc_words;
    }
    size_t data_alloc_size_in_words() const {
        return m * n * data_entry_alloc_size_in_words();
    }
    size_t data_entry_alloc_size_in_bytes() const {
        return alloc_words * sizeof(unsigned long);
    }
    size_t data_alloc_size_in_bytes() const {
        return m * n * data_entry_alloc_size_in_bytes();
    }
    inline unsigned int nrows() const { return m; }
    inline unsigned int ncols() const { return n; }
    inline size_t get_size() const { return size; }
    void set_size(size_t s) { size = s; }
    size_t get_true_nonzero_size() const;
    unsigned int valuation() const;

    matpoly() { m=n=0; size=0; alloc_words=0; x=NULL; }
    matpoly(arith_hard * ab, unsigned int m, unsigned int n, int len);
    matpoly(matpoly const&) = delete;
    matpoly& operator=(matpoly const&) = delete;
    matpoly& set(matpoly const&);
    matpoly(matpoly &&);
    matpoly& operator=(matpoly &&);
    ~matpoly();
    matpoly similar_shell() const { return matpoly(nullptr, m, n, 0); }
    bool check_pre_init() const ATTRIBUTE_WARN_UNUSED_RESULT {
        return x == NULL;
    }
    void realloc(size_t new_number_of_coeffs);
    inline void shrink_to_fit() { realloc(size); }
    void zero();
    void clear() { *this = matpoly(); }

    /* {{{ access interface for matpoly */
    /* part_head does not _promise_ to point to coefficient k exactly. In
     * the simd case (lingen_matpoly_binary.hpp) it points to the word
     * where coefficient k can be found */
    inline ptr part(unsigned int i, unsigned int j) {
        return x + (i*n+j)*alloc_words;
    }
    // inline elt& coeff(unsigned int i, unsigned int j) { return *part(i, j); }
    inline srcptr part(unsigned int i, unsigned int j) const {
        return x + (i*n+j)*alloc_words;
    }
    public:
    inline ptr part_head(unsigned int i, unsigned int j, unsigned int k) {
        /* k is a number of bits. Only k/ULONG_BITS matter (see
         * ab->vec_subvec)
         */
        return ab->vec_subvec(part(i, j), k);
    }
    inline srcptr part_head(unsigned int i, unsigned int j, unsigned int k) const {
        /* k is a number of bits. Only k/ULONG_BITS matter (see
         * ab->vec_subvec)
         */
        return ab->vec_subvec(part(i, j), k);
    }
    inline elt coeff(unsigned int i, unsigned int j, unsigned int k) const {
        /* k is a number of bits */
        return (*part_head(i, j, k) >> (k % ULONG_BITS)) != 0;
    }
    struct coeff_accessor_proxy {
        unsigned long & p;
        unsigned int kr;
        coeff_accessor_proxy(matpoly& F, unsigned int i,
                unsigned int j, unsigned int k)
            : p(*F.part_head(i, j, k))
            , kr(k % ULONG_BITS)
        {
        }
        operator elt() const { return (p>>kr) & 1; }
        coeff_accessor_proxy& operator+=(unsigned long x) {
            p ^= (x&1) << kr;
            return *this;
        }
    };
    inline coeff_accessor_proxy coeff(unsigned int i, unsigned int j, unsigned int k) {
        return coeff_accessor_proxy(*this, i, j, k);
    }
    public:
    /* }}} */

    /* The interfaces below used to exist for the old binary "polmat"
     * type, and we wish to do away with them.
     */
    void addpoly(unsigned int i, unsigned int j, matpoly const& y, unsigned int iy, unsigned int jy) __attribute__((deprecated));
    void xmul_poly(unsigned int i, unsigned int j, unsigned long s) __attribute__((deprecated));
    ptr poly(unsigned int i, unsigned int j) __attribute__((deprecated)) { return part(i, j); }
    srcptr poly(unsigned int i, unsigned int j) const __attribute__((deprecated)) { return part(i, j); }

    void set_constant_ui(unsigned long e);
    void set_constant(unsigned long e) { set_constant_ui(e); }

    /* Note that this does not affect the size field */
    void fill_random(unsigned int k0, unsigned int k1, gmp_randstate_t rstate);

    void clear_and_set_random(unsigned int len, gmp_randstate_t rstate)
    {
        if (len > capacity())
            zero_pad(len);
        zero_pad(capacity());
        fill_random(0, len, rstate);
        set_size(len);
    }
    int cmp(matpoly const & b) const;
    void multiply_column_by_x(unsigned int j, unsigned int size);
    void divide_column_by_x(unsigned int j, unsigned int size);
    void truncate(matpoly const & src, unsigned int size);
    void truncate(unsigned int size) { truncate(*this, size); }
    /* This checks that coefficients of degree k to size-1 are zero.
     */
    int tail_is_zero(unsigned int k) const;
public:
    /* not to be confused with the former. a priori this is an
     * implementation detail. At times, we want to assert that.
     */
    bool high_word_is_clear() const;
private:
    void clear_high_word_common(unsigned int length);
public:
    inline void clear_high_word() { clear_high_word_common(size); }
public:
    /* This changes size to nsize, and fills [size..nsize[ with zeroes */
    void zero_pad(unsigned int nsize);

    void zero_with_size(size_t size) { set_size(0); zero_pad(size); }

    void extract_column(
        unsigned int jdst, unsigned int kdst,
        matpoly const & src, unsigned int jsrc, unsigned int ksrc);
    void zero_column(unsigned int jdst, unsigned int kdst);
    void rshift(matpoly const &, unsigned int k);
    void rshift(unsigned int k);

    void add(matpoly const & a, matpoly const & b);
    void sub(matpoly const & a, matpoly const & b);
    void add(matpoly const & a) { add(*this, a); }
    void sub(matpoly const & a) { sub(*this, a); }

    static matpoly mul(matpoly const & a, matpoly const & b);
    static matpoly mp(matpoly const & a, matpoly const & c);
    void addmul(matpoly const & a, matpoly const & b);
    void addmp(matpoly const & a, matpoly const & c);

    /* We have nothing terrific to identify as substeps here, and there's
     * no schedule information that we intend to exploit */
    static matpoly mp(tree_stats & stats, matpoly const & a, matpoly const & b, lingen_call_companion::mul_or_mp_times * M)
    {
        if (!M) return mp(a, b);
        tree_stats::smallstep_sentinel dummy(stats, M->step_name());
        return mp(a, b);
    }

    static matpoly mul(tree_stats & stats, matpoly const & a, matpoly const & b, lingen_call_companion::mul_or_mp_times * M) {
        if (!M) return mul(a, b);
        tree_stats::smallstep_sentinel dummy(stats, M->step_name());
        return mul(a, b);
    }
 

    // void set_polymat(polymat const & src);
    int coeff_is_zero(unsigned int k) const;
    void coeff_set_zero(unsigned int k);

    struct view_t : public submatrix_range {
        matpoly & M;
        view_t(matpoly & M, submatrix_range S) : submatrix_range(S), M(M) {}
        view_t(matpoly & M) : submatrix_range(M), M(M) {}
        inline ptr part(unsigned int i, unsigned int j) {
            return M.part(i0+i, j0+j);
        }
        inline srcptr part(unsigned int i, unsigned int j) const {
            return M.part(i0+i, j0+j);
        }
        void zero();
    };

    struct const_view_t : public submatrix_range {
        matpoly const & M;
        const_view_t(matpoly const & M, submatrix_range S) : submatrix_range(S), M(M) {}
        const_view_t(matpoly const & M) : submatrix_range(M), M(M) {}
        const_view_t(view_t const & V) : submatrix_range(V), M(V.M) {}
        inline srcptr part(unsigned int i, unsigned int j) const {
            return M.part(i0+i, j0+j);
        }
    };
    view_t view(submatrix_range S) { ASSERT_ALWAYS(S.valid(*this)); return view_t(*this, S); }
    const_view_t view(submatrix_range S) const { ASSERT_ALWAYS(S.valid(*this)); return const_view_t(*this, S); }
    view_t view() { return view_t(*this); }
    const_view_t view() const { return const_view_t(*this); }

    static void copy(view_t t, const_view_t a);
    static void addmul(view_t t, const_view_t t0, const_view_t t1);
    static void addmp(view_t t, const_view_t t0, const_view_t t1);

    matpoly truncate_and_rshift(unsigned int truncated_size, unsigned int rshift);
};

#endif	/* LINGEN_MATPOLY_BINARY_HPP_ */
