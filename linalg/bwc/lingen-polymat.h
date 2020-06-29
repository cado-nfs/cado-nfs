#ifndef POLYMAT_H_
#define POLYMAT_H_

#include <stddef.h>      // for size_t
#include <gmp.h>         // for gmp_randstate_t
#include "macros.h"      // for ASSERT_ALWAYS
#include "mpfq_layer.h"  // for abdst_field, mpfq_pz_dst_field, MAYBE_UNUSED

/* This is used only for plingen. */

struct polymat_s {
    unsigned int m;
    unsigned int n;
    size_t size;
    size_t alloc;
    abvec x;
};
typedef struct polymat_s polymat[1];
#ifndef MATPOLY_H_
typedef struct polymat_s * polymat_ptr;
typedef const struct polymat_s * polymat_srcptr;
typedef struct matpoly_s * matpoly_ptr;
typedef const struct matpoly_s * matpoly_srcptr;
#endif

#include "lingen-matpoly.h"

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
#ifdef __cplusplus
extern "C" {
#endif
void polymat_cutoff_info_init(struct polymat_cutoff_info * c);
void polymat_cutoff_info_clear(struct polymat_cutoff_info * c);
void polymat_cutoff_add_step(struct polymat_cutoff_info * c, unsigned int size, unsigned int alg);
void polymat_set_mul_kara_cutoff(const struct polymat_cutoff_info * new_cutoff, struct polymat_cutoff_info * old_cutoff);
void polymat_set_mp_kara_cutoff(const struct polymat_cutoff_info * new_cutoff, struct polymat_cutoff_info * old_cutoff);
#ifdef __cplusplus
}
#endif
/* }}} */

#ifdef __cplusplus
extern "C" {
#endif
void polymat_init(abdst_field ab, polymat_ptr p, unsigned int m, unsigned int n, int len);
int polymat_check_pre_init(polymat_srcptr p);
void polymat_realloc(abdst_field ab, polymat_ptr p, size_t newalloc);
void polymat_zero(abdst_field ab, polymat_ptr p);
void polymat_clear(abdst_field ab, polymat_ptr p);
void polymat_fill_random(abdst_field ab MAYBE_UNUSED, polymat_ptr a, unsigned int size, gmp_randstate_t rstate);
void polymat_swap(polymat_ptr a, polymat_ptr b);
int polymat_cmp(abdst_field ab MAYBE_UNUSED, polymat_srcptr a, polymat_srcptr b);
static inline abdst_vec polymat_part(abdst_field ab, polymat_ptr p, unsigned int i, unsigned int j, unsigned int k);
static inline abdst_elt polymat_coeff(abdst_field ab, polymat_ptr p, unsigned int i, unsigned int j, unsigned int k);

void polymat_set_matpoly(abdst_field ab MAYBE_UNUSED, polymat_ptr a, matpoly_srcptr b);



void polymat_addmat(abdst_field ab,
        polymat c, unsigned int kc,
        polymat a, unsigned int ka,
        polymat b, unsigned int kb);

void polymat_submat(abdst_field ab,
        polymat c, unsigned int kc,
        polymat a, unsigned int ka,
        polymat b, unsigned int kb);

void polymat_truncate(abdst_field ab, polymat_ptr dst, polymat_ptr src, unsigned int size);
void polymat_multiply_column_by_x(abdst_field ab, polymat_ptr pi, unsigned int j, unsigned int size);
void polymat_extract_column(abdst_field ab,
        polymat_ptr dst, unsigned int jdst, unsigned int kdst,
        polymat_ptr src, unsigned int jsrc, unsigned int ksrc);
void polymat_extract_row_fragment(abdst_field ab,
        polymat_ptr dst, unsigned int i1, unsigned int j1,
        polymat_ptr src, unsigned int i0, unsigned int j0,
        unsigned int n);
void polymat_rshift(abdst_field ab, polymat_ptr dst, polymat_ptr src, unsigned int k);


void polymat_addmul(abdst_field ab, polymat c, polymat a, polymat b);
void polymat_addmp(abdst_field ab, polymat c, polymat a, polymat b);
void polymat_mul(abdst_field ab, polymat c, polymat a, polymat b);
void polymat_mp(abdst_field ab, polymat c, polymat a, polymat b);


#ifdef __cplusplus
}
#endif

/* {{{ access interface for polymat */
static inline abdst_vec polymat_part(abdst_field ab, polymat_ptr p, unsigned int i, unsigned int j, unsigned int k) {
    /* Assume row-major in all circumstances. Old code used to support
     * various orderings, here we don't */
    ASSERT_ALWAYS(p->size);
    return abvec_subvec(ab, p->x, (k*p->m+i)*p->n+j);
}
static inline abdst_elt polymat_coeff(abdst_field ab, polymat_ptr p, unsigned int i, unsigned int j, unsigned int k) {
    return abvec_coeff_ptr(ab, polymat_part(ab, p,i,j,k),0);
}
static inline absrc_vec polymat_part_const(abdst_field ab, polymat_srcptr p, unsigned int i, unsigned int j, unsigned int k) {
    /* Assume row-major in all circumstances. Old code used to support
     * various orderings, here we don't */
    ASSERT_ALWAYS(p->size);
    return abvec_subvec_const(ab, (absrc_vec)p->x,(k*p->m+i)*p->n+j);
}
static inline absrc_elt polymat_coeff_const(abdst_field ab, polymat_srcptr p, unsigned int i, unsigned int j, unsigned int k) {
    return abvec_coeff_ptr_const(ab, polymat_part_const(ab, p,i,j,k),0);
}
/* }}} */



#endif	/* POLYMAT_H_ */
