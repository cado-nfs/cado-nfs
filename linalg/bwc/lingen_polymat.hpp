#ifndef LINGEN_POLYMAT_HPP_
#define LINGEN_POLYMAT_HPP_

#ifdef LINGEN_BINARY
#error "lingen_polymat does not work with binary (maybe use bblas instead -- could end up being the same interface, IDK)"
#endif

#include <cstddef>      // for size_t, NULL
#include <gmp.h>         // for gmp_randstate_t
#include "arith-hard.hpp"
#include "macros.h"

class matpoly;

/* This is used only for lingen. */
struct polymat {
    arith_hard * ab = NULL;
    unsigned int m = 0;
    unsigned int n = 0;
    private:
    size_t size = 0;
    size_t alloc = 0;
    public:
    inline size_t capacity() const { return alloc; }
    inline size_t get_size() const { return size; }
    void set_size(size_t s) { size = s; }
    typedef arith_hard::elt elt;
    typedef elt * ptr;
    typedef elt const * srcptr;
    ptr x = NULL;

    polymat() : polymat(NULL, 0,0,0) {}
    polymat(polymat const&) = delete;
    polymat& operator=(polymat const&) = delete;
    polymat(polymat &&);
    polymat& operator=(polymat &&);
    ~polymat();
    polymat(arith_hard * ab, unsigned int m, unsigned int n, int len);
    int check_pre_init() const;
    void realloc(size_t newalloc);
    void zero();
#if 0
    void swap(polymat & b);
#endif
    void fill_random(unsigned int nsize, gmp_randstate_t rstate);
    int cmp(polymat const & b);

    /* {{{ access interface for polymat */
    inline ptr part(unsigned int i, unsigned int j, unsigned int k) {
        /* Assume row-major in all circumstances. Old code used to support
         * various orderings, here we don't */
        ASSERT_ALWAYS(size);
        return ab->vec_subvec(x, (k*m+i)*n+j);
    }
    inline elt & coeff(unsigned int i, unsigned int j, unsigned int k) {
        return ab->vec_item(part(i,j,k),0);
    }
    inline srcptr part(unsigned int i, unsigned int j, unsigned int k) const {
        /* Assume row-major in all circumstances. Old code used to support
         * various orderings, here we don't */
        ASSERT_ALWAYS(size);
        return ab->vec_subvec(x,(k*m+i)*n+j);
    }
    inline elt const & coeff(unsigned int i, unsigned int j, unsigned int k) const {
        return ab->vec_item(part(i,j,k),0);
    }
    /* }}} */

    void addmat(
        unsigned int kc,
        polymat const & a, unsigned int ka,
        polymat const & b, unsigned int kb);
    void submat(
        unsigned int kc,
        polymat const & a, unsigned int ka,
        polymat const & b, unsigned int kb);
    void mulmat(
        unsigned int kc,
        polymat const & a, unsigned int ka,
        polymat const & b, unsigned int kb);
    void addmulmat(
        unsigned int kc,
        polymat const & a, unsigned int ka,
        polymat const & b, unsigned int kb);
    void multiply_column_by_x(unsigned int j, unsigned int size);
    void truncate(polymat const & src, unsigned int nsize);
    void extract_column(
        unsigned int jdst, unsigned int kdst,
        polymat const & src, unsigned int jsrc, unsigned int ksrc);
    void extract_row_fragment(
        unsigned int i1, unsigned int j1,
        polymat const & src, unsigned int i0, unsigned int j0,
        unsigned int n);
    void mul(polymat const & a, polymat const & b);
    void addmul(polymat const & a, polymat const & b);
    void mp(polymat const & a, polymat const & b);
    void addmp(polymat const & a, polymat const & b);

    void set_matpoly(matpoly const &);
    void rshift(polymat const & src, unsigned int k);

};

/* {{{ cut-off structures */
/* This structure is used to decide which algorithm to use for a given
 * input length. This essentially is a function from N*N to a finite set of
 * choices (so far, {0,1} only). The value returned for a balanced input
 * length x*x is:
 * if x >= cut: 1
 * if x < cut:
 *      if table == NULL:
 *              return 0
 *      find last (s,a) in table such that s <= x
 *      return a
 * For unbalanced input length x*y, we compare MIN(x,y) with subdivide.
 * At the moment there is no additional table designed to improve the
 * choice.
 *
 * Cutoffs are used only by lingen-polymat.c and lingen-bigpolymat.c -- however the
 * tuning program also uses it, so we expose it here too.
 */
/* Such a cutoff structure is used for finding which algorithm to used
 * for:
 *  - a balanced x*x product. -> basecase(0) or karatsuba(1)
 *  - an unbalanced x*y product. -> basecase(0) or split into balanced
 *  sizes(1)
 */

struct polymat_cutoff_info {
    unsigned int cut;
    unsigned int subdivide;
    unsigned int (*table)[2];
    unsigned int table_size;
};
void polymat_cutoff_info_init(struct polymat_cutoff_info * c);
void polymat_cutoff_info_clear(struct polymat_cutoff_info * c);
void polymat_cutoff_add_step(struct polymat_cutoff_info * c, unsigned int size, unsigned int alg);
void polymat_set_mul_kara_cutoff(const struct polymat_cutoff_info * new_cutoff, struct polymat_cutoff_info * old_cutoff);
void polymat_set_mp_kara_cutoff(const struct polymat_cutoff_info * new_cutoff, struct polymat_cutoff_info * old_cutoff);
/* }}} */

#endif	/* LINGEN_POLYMAT_HPP_ */
