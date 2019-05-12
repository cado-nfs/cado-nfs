#ifndef LINGEN_MATPOLY_H_
#define LINGEN_MATPOLY_H_

#include "mpfq_layer.h"

struct matpoly;
struct polymat;

#include "lingen-polymat.hpp"

/* This is used only for plingen. */

/* We use abvec because this offers the possibility of having flat data
 *
 * Note that this ends up being exactly the same data type as polymat.
 * The difference here is that the stride is not the same.
 */
struct matpoly {
    abdst_field ab = NULL;
    unsigned int m = 0;
    unsigned int n = 0;
    size_t size = 0;
    size_t alloc = 0;
    abvec x = NULL;

    matpoly() { m=n=0; size=alloc=0; ab=NULL; x=NULL; }
    matpoly(abdst_field ab, unsigned int m, unsigned int n, int len);
    matpoly(matpoly const&) = delete;
    matpoly& operator=(matpoly const&) = delete;
    matpoly(matpoly &&);
    matpoly& operator=(matpoly &&);
    ~matpoly();
    int check_pre_init() const;
    void realloc(size_t);
    void zero();
    /* {{{ access interface for matpoly */
    inline abdst_vec part(unsigned int i, unsigned int j, unsigned int k) {
        ASSERT_ALWAYS(size);
        return abvec_subvec(ab, x, (i*n+j)*alloc+k);
    }
    inline abdst_elt coeff(unsigned int i, unsigned int j, unsigned int k) {
        return abvec_coeff_ptr(ab, part(i,j,k), 0);
    }
    inline absrc_vec part(unsigned int i, unsigned int j, unsigned int k) const {
        ASSERT_ALWAYS(size);
        return abvec_subvec_const(ab, x, (i*n+j)*alloc+k);
    }
    inline absrc_elt coeff(unsigned int i, unsigned int j, unsigned int k) const {
        return abvec_coeff_ptr_const(ab, part(i,j,k), 0);
    }
    /* }}} */
    void set_constant_ui(unsigned long e);
    void set_constant(absrc_elt e);
#if 0
    void swap(matpoly&);
#endif
    void fill_random(unsigned int size, gmp_randstate_t rstate);
    int cmp(matpoly const & b) const;
    void multiply_column_by_x(unsigned int j, unsigned int size);
    void divide_column_by_x(unsigned int j, unsigned int size);
    void truncate(matpoly const & src, unsigned int size);
    void extract_column(
        unsigned int jdst, unsigned int kdst,
        matpoly const & src, unsigned int jsrc, unsigned int ksrc);
    void transpose_dumb(matpoly const & src);
    void zero_column(unsigned int jdst, unsigned int kdst);
    void extract_row_fragment(unsigned int i1, unsigned int j1,
        matpoly const & src, unsigned int i0, unsigned int j0,
        unsigned int n);
    void rshift(matpoly const &, unsigned int k);
    void addmul(matpoly const & a, matpoly const & b);
    void mul(matpoly const & a, matpoly const & b);
    void addmp(matpoly const & a, matpoly const & c);
    void mp(matpoly const & a, matpoly const & c);
    void set_polymat(polymat const & src);
    int coeff_is_zero(unsigned int k) const;
    void coeff_set_zero(unsigned int k);
};

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __cplusplus
}
#endif


#endif	/* LINGEN_MATPOLY_H_ */
