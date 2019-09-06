#ifndef LINGEN_MATPOLY_BINARY_HPP_
#define LINGEN_MATPOLY_BINARY_HPP_

/* The outer interface of the matpoly_binary type is exactly the same as
 * for the matpoly type. This is enforced by the test program. In
 * particular, we copy even the "abdst_field" argument that is passed
 * everywhere, with the slight catch that the code *here* does not use it
 * at all.
 *
 * We have two options:
 *
 * For convenience, consider this binary matpoly type as
 * linked to the u64k1 type
 *
 * Replace all the "ab" stuff by empty proxies.
 *
 * I haven't made my mind yet as to which is best.
 */

#include <cstdlib>
#include <gmp.h>

#include <mutex>

#include "cado_config.h"
#include "macros.h"
#include "lingen_memory_pool.hpp"
#include "lingen_submatrix.hpp"
#include "mpfq_fake.hpp"

class matpoly {
    static memory_pool_loose memory;
    friend decltype(memory)::guard<matpoly>;
public:
    typedef decltype(memory)::guard<matpoly> memory_guard;
    static constexpr bool over_gf2 = true;
    // static void add_to_main_memory_pool(size_t s);
    abdst_field ab = NULL;
    unsigned int m = 0;
    unsigned int n = 0;
    size_t size = 0;
    /* alloc_words is the number of unsigned longs used to store each
     * coefficient */
private:
    size_t alloc_words = 0;
    unsigned long * x = NULL;
#define BITS_TO_WORDS(B,W)      iceildiv((B),(W))
    static inline size_t b2w(size_t n) { return BITS_TO_WORDS(n, ULONG_BITS); }/*{{{*/
    // inline size_t colstride() const { return nrows() * stride(); }/*}}}*/
    size_t alloc_size_words() const { return nrows() * ncols() * alloc_words; }
public:
    inline size_t capacity() const { return alloc_words * ULONG_BITS; }
    inline unsigned int nrows() const { return m; }
    inline unsigned int ncols() const { return n; }
    matpoly() { m=n=0; size=0; alloc_words=0; ab=NULL; x=NULL; }
    matpoly(abdst_field ab, unsigned int m, unsigned int n, int len);
    matpoly(matpoly const&) = delete;
    matpoly& operator=(matpoly const&) = delete;
    matpoly& set(matpoly const&);
    matpoly(matpoly &&);
    matpoly& operator=(matpoly &&);
    ~matpoly();
    bool check_pre_init() const ATTRIBUTE_WARN_UNUSED_RESULT { return x == NULL; }
    void realloc(size_t new_number_of_coeffs);
    inline void shrink_to_fit() { realloc(size); }
    void zero();
    /* {{{ access interface for matpoly */
    inline abdst_vec part(unsigned int i, unsigned int j, unsigned int k=0) {
        ASSERT_ALWAYS((k % ULONG_BITS) == 0);
        return x + (i*n+j)*alloc_words+k / ULONG_BITS;
    }
    inline abdst_elt coeff(unsigned int i, unsigned int j, unsigned int k=0) {
        return part(i,j,k);
    }
    inline absrc_vec part(unsigned int i, unsigned int j, unsigned int k=0) const {
        ASSERT_ALWAYS((k % ULONG_BITS) == 0);
        return x + (i*n+j)*alloc_words+k / ULONG_BITS;
    }
    inline absrc_elt coeff(unsigned int i, unsigned int j, unsigned int k=0) const {
        return part(i,j,k);
    }
    /* }}} */

    /* The interfaces below used to exist for the old binary "polmat"
     * type, and we wish to do away with them.
     */
    void addpoly(unsigned int i, unsigned int j, matpoly const& y, unsigned int iy, unsigned int jy) __attribute__((deprecated));
    void xmul_poly(unsigned int i, unsigned int j, unsigned long s) __attribute__((deprecated));
    unsigned long * poly(unsigned int i, unsigned int j) __attribute__((deprecated)) { return part(i, j); }
    const unsigned long * poly(unsigned int i, unsigned int j) const __attribute__((deprecated)) { return part(i, j); }

    void set_constant_ui(unsigned long e);
    void set_constant(absrc_elt e) { set_constant_ui(*e); }
    void fill_random(unsigned int size, gmp_randstate_t rstate);
    int cmp(matpoly const & b) const;
    void multiply_column_by_x(unsigned int j, unsigned int size);
    void divide_column_by_x(unsigned int j, unsigned int size);
    void truncate(matpoly const & src, unsigned int size);
    void truncate(unsigned int size) { truncate(*this, size); }
    /* This checks that coefficients of degree k to size-1 are zero.
     */
    int tail_is_zero(unsigned int k) const;
private:
    /* not to be confused with the former. a priori this is an
     * implementation detail. At times, we want to assert that.
     */
    bool high_word_is_clear() const;
    void clear_high_word();
public:
    void zero_pad(unsigned int nsize); /* changes size to nsize */
    void extract_column(
        unsigned int jdst, unsigned int kdst,
        matpoly const & src, unsigned int jsrc, unsigned int ksrc);
    void zero_column(unsigned int jdst, unsigned int kdst);
    void rshift(matpoly const &, unsigned int k);

    void add(matpoly const & a, matpoly const & b);
    void sub(matpoly const & a, matpoly const & b);
    void add(matpoly const & a) { add(*this, a); }
    void sub(matpoly const & a) { sub(*this, a); }

    /* It is probably wise to avoid the mul and mp functions below. The
     * first-class citizens are the caching alternatives.
     */
    void mul(matpoly const & a, matpoly const & b);
    void mp(matpoly const & a, matpoly const & c);
    void addmul(matpoly const & a, matpoly const & b);
    void addmp(matpoly const & a, matpoly const & c);

    // void set_polymat(polymat const & src);
    int coeff_is_zero(unsigned int k) const;
    void coeff_set_zero(unsigned int k);
    struct view_t;
    struct const_view_t;

    struct view_t : public submatrix_range {
        matpoly & M;
        view_t(matpoly & M, submatrix_range S) : submatrix_range(S), M(M) {}
        view_t(matpoly & M) : submatrix_range(M), M(M) {}
        inline abdst_vec part(unsigned int i, unsigned int j) {
            return M.part(i0+i, j0+j);
        }
        inline absrc_vec part(unsigned int i, unsigned int j) const {
            return M.part(i0+i, j0+j);
        }
    };

    struct const_view_t : public submatrix_range {
        matpoly const & M;
        const_view_t(matpoly const & M, submatrix_range S) : submatrix_range(S), M(M) {}
        const_view_t(matpoly const & M) : submatrix_range(M), M(M) {}
        const_view_t(view_t const & V) : submatrix_range(V), M(V.M) {}
        inline absrc_vec part(unsigned int i, unsigned int j) const {
            return M.part(i0+i, j0+j);
        }
    };
    view_t view(submatrix_range S) { ASSERT_ALWAYS(S.valid(*this)); return view_t(*this, S); }
    const_view_t view(submatrix_range S) const { ASSERT_ALWAYS(S.valid(*this)); return const_view_t(*this, S); }
    view_t view() { return view_t(*this); }
    const_view_t view() const { return const_view_t(*this); }
    matpoly truncate_and_rshift(unsigned int truncated_size, unsigned int rshift);
};

#endif	/* LINGEN_MATPOLY_BINARY_HPP_ */
