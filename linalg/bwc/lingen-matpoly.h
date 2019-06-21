#ifndef LINGEN_MATPOLY_H_
#define LINGEN_MATPOLY_H_

#include "mpfq_layer.h"

/* This is used only for plingen. */

/* We use abvec because this offers the possibility of having flat data
 *
 * Note that this ends up being exactly the same data type as polymat.
 * The difference here is that the stride is not the same.
 */
struct matpoly_s {
    unsigned int m;
    unsigned int n;
    size_t size;
    size_t alloc;
    abvec x;
};
typedef struct matpoly_s matpoly[1];
typedef struct matpoly_s * matpoly_ptr;
typedef const struct matpoly_s * matpoly_srcptr;

#include "lingen-polymat.h"

#ifdef __cplusplus
extern "C" {
#endif

void matpoly_init(abdst_field ab, matpoly_ptr p, unsigned int m, unsigned int n, int len);
int matpoly_check_pre_init(matpoly_srcptr p);
void matpoly_realloc(abdst_field ab, matpoly_ptr p, size_t newalloc);
void matpoly_zero(abdst_field ab, matpoly_ptr p);
void matpoly_clear(abdst_field ab, matpoly_ptr p);
void matpoly_fill_random(abdst_field ab MAYBE_UNUSED, matpoly_ptr a, unsigned int size, gmp_randstate_t rstate);
void matpoly_swap(matpoly_ptr a, matpoly_ptr b);
int matpoly_cmp(abdst_field ab MAYBE_UNUSED, matpoly_srcptr a, matpoly_srcptr b);
static inline abdst_vec matpoly_part(abdst_field ab, matpoly_ptr p, unsigned int i, unsigned int j, unsigned int k);
static inline abdst_elt matpoly_coeff(abdst_field ab, matpoly_ptr p, unsigned int i, unsigned int j, unsigned int k);
void matpoly_set_polymat(abdst_field ab MAYBE_UNUSED, matpoly_ptr dst, polymat_srcptr src);
void matpoly_set_constant_ui(abdst_field ab, matpoly_ptr p, unsigned long e);
void matpoly_set_constant(abdst_field ab, matpoly_ptr p, absrc_elt e);

#if 0   /* nonexistent ? */
void matpoly_addmat(abdst_field ab,
        matpoly c, unsigned int kc,
        matpoly a, unsigned int ka,
        matpoly b, unsigned int kb);

void matpoly_submat(abdst_field ab,
        matpoly c, unsigned int kc,
        matpoly a, unsigned int ka,
        matpoly b, unsigned int kb);
#endif

void matpoly_truncate(abdst_field ab, matpoly_ptr dst, matpoly_srcptr src, unsigned int size);
void matpoly_multiply_column_by_x(abdst_field ab, matpoly_ptr pi, unsigned int j, unsigned int size);
void matpoly_divide_column_by_x(abdst_field ab, matpoly_ptr pi, unsigned int j, unsigned int size);
void matpoly_zero_column(abdst_field ab,
        matpoly_ptr dst, unsigned int jdst, unsigned int kdst);
void matpoly_extract_column(abdst_field ab,
        matpoly_ptr dst, unsigned int jdst, unsigned int kdst,
        matpoly_srcptr src, unsigned int jsrc, unsigned int ksrc);
void matpoly_extract_row_fragment(abdst_field ab,
        matpoly_ptr dst, unsigned int i1, unsigned int j1,
        matpoly_srcptr src, unsigned int i0, unsigned int j0,
        unsigned int n);
void matpoly_rshift(abdst_field ab, matpoly_ptr dst, matpoly_srcptr src, unsigned int k);

void matpoly_transpose_dumb(abdst_field ab, matpoly_ptr dst, matpoly_srcptr src);

double matpoly_addmul(abdst_field ab, matpoly c, matpoly a, matpoly b, int draft);
double matpoly_addmp(abdst_field ab, matpoly c, matpoly a, matpoly b, int draft);
double matpoly_mul(abdst_field ab, matpoly c, matpoly a, matpoly b, int draft);
double matpoly_mp(abdst_field ab, matpoly c, matpoly a, matpoly b, int draft);

int matpoly_coeff_is_zero(abdst_field ab, matpoly_srcptr pi, unsigned int k);
void matpoly_coeff_set_zero(abdst_field ab, matpoly_ptr pi, unsigned int k);

/* {{{ access interface for matpoly */
static inline abdst_vec matpoly_part(abdst_field ab, matpoly_ptr p, unsigned int i, unsigned int j, unsigned int k) {
    ASSERT_ALWAYS(p->size);
    return abvec_subvec(ab, p->x, (i*p->n+j)*p->alloc+k);
}
static inline abdst_elt matpoly_coeff(abdst_field ab, matpoly_ptr p, unsigned int i, unsigned int j, unsigned int k) {
    return abvec_coeff_ptr(ab, matpoly_part(ab, p,i,j,k), 0);
}
static inline absrc_vec matpoly_part_const(abdst_field ab, matpoly_srcptr p, unsigned int i, unsigned int j, unsigned int k) {
    ASSERT_ALWAYS(p->size);
    return abvec_subvec_const(ab, p->x, (i*p->n+j)*p->alloc+k);
}
static inline absrc_elt matpoly_coeff_const(abdst_field ab, matpoly_srcptr p, unsigned int i, unsigned int j, unsigned int k) {
    return abvec_coeff_ptr_const(ab, matpoly_part_const(ab, p,i,j,k), 0);
}
/* }}} */

#ifdef __cplusplus
}
#endif


#endif	/* LINGEN_MATPOLY_H_ */
