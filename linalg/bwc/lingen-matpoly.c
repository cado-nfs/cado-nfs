#include "cado.h" // IWYU pragma: keep
#include "mpfq_layer.h"
#include <stdlib.h>
#include <gmp.h>
#include "flint-fft/fft.h"
#include "lingen-matpoly.h"
#include "macros.h"
#include "timing.h"     // seconds

int matpoly_check_pre_init(matpoly_srcptr p)
{
    if (p->m && p->n && p->alloc)
        return 0;
    if (!p->m && !p->n && !p->alloc)
        return 1;
    abort();
    return 0;
}

/* with the exception of matpoly_realloc, all functions here are exactly
 * identical to those in lingen-polymat.c */
/* {{{ init/zero/clear interface for matpoly */
void matpoly_init(abdst_field ab, matpoly_ptr p, unsigned int m, unsigned int n, int len) {
    memset(p, 0, sizeof(matpoly));
    /* As a special case, we allow a pre-init state with m==n==len==0 */
    ASSERT(!m == !n);
    ASSERT(!m == !len);
    if (!m) return;
    p->m = m;
    p->n = n;
    p->alloc = len;
    p->size = 0;
    abvec_init(ab, &(p->x), m*n*p->alloc);
    abvec_set_zero(ab, p->x, m*n*p->alloc);
}
void matpoly_realloc(abdst_field ab, matpoly_ptr p, size_t newalloc) {
    matpoly_check_pre_init(p);
    /* zero out the newly added data */
    if (newalloc > p->alloc) {
        /* allocate new space, then deflate */
        abvec_reinit(ab, &(p->x), p->m*p->n*p->alloc, p->m*p->n*newalloc);
        abvec rhead = abvec_subvec(ab, p->x, p->m*p->n*p->alloc);
        abvec whead = abvec_subvec(ab, p->x, p->m*p->n*newalloc);
        for(unsigned int i = p->m ; i-- ; ) {
            for(unsigned int j = p->n ; j-- ; ) {
                whead = abvec_subvec(ab, whead, -newalloc);
                rhead = abvec_subvec(ab, rhead, -p->alloc);
                abvec_set(ab, whead, rhead, p->alloc);
                abvec_set_zero(ab,
                        abvec_subvec(ab, whead, p->alloc),
                        newalloc - p->alloc);
            }
        }
    } else {
        /* inflate, then free space */
        ASSERT_ALWAYS(p->size <= newalloc);
        abvec rhead = p->x;
        abvec whead = p->x;
        for(unsigned int i = 0 ; i < p->m ; i++) {
            for(unsigned int j = 0 ; j < p->n ; j++) {
                abvec_set(ab, whead, rhead, newalloc);
                whead = abvec_subvec(ab, whead, newalloc);
                rhead = abvec_subvec(ab, rhead, p->alloc);
            }
        }
        abvec_reinit(ab, &(p->x), p->m*p->n*p->alloc, p->m*p->n*newalloc);
    }
    p->alloc = newalloc;
}
void matpoly_zero(abdst_field ab, matpoly_ptr p) {
    p->size = 0;
    abvec_set_zero(ab, p->x, p->m*p->n*p->alloc);
}
void matpoly_set_constant_ui(abdst_field ab, matpoly_ptr p, unsigned long e) {
    ASSERT_ALWAYS(p->m == p->n);
    p->size = 0;
    if (p->alloc == 0 && e)
        matpoly_realloc(ab, p, 1);
    matpoly_zero(ab, p);
    if (!e) return;
    p->size = 1;
    for(unsigned int i = 0 ; i < p->m ; ++i)
        abset_ui(ab, matpoly_coeff(ab, p, i, i, 0), e);
}
void matpoly_set_constant(abdst_field ab, matpoly_ptr p, absrc_elt e) {
    ASSERT_ALWAYS(p->m == p->n);
    p->size = 0;
    if (p->alloc == 0 && abcmp_ui(ab, e, 0) != 0)
        matpoly_realloc(ab, p, 1);
    matpoly_zero(ab, p);
    if (abcmp_ui(ab, e, 0) == 0) return;
    p->size = 1;
    for(unsigned int i = 0 ; i < p->m ; ++i)
        abset(ab, matpoly_coeff(ab, p, i, i, 0), e);
}
void matpoly_clear(abdst_field ab, matpoly_ptr p) {
    abvec_clear(ab, &(p->x), p->m*p->n*p->alloc);
    memset(p, 0, sizeof(matpoly));
}
void matpoly_swap(matpoly_ptr a, matpoly_ptr b)
{
    matpoly x;
    memcpy(x, a, sizeof(matpoly));
    memcpy(a, b, sizeof(matpoly));
    memcpy(b, x, sizeof(matpoly));
}
/* }}} */

void matpoly_fill_random(abdst_field ab MAYBE_UNUSED, matpoly_ptr a, unsigned int size, gmp_randstate_t rstate)
{
    ASSERT_ALWAYS(size <= a->alloc);
    a->size = size;
    abvec_random(ab, a->x, a->m*a->n*size, rstate);
}

int matpoly_cmp(abdst_field ab MAYBE_UNUSED, matpoly_srcptr a, matpoly_srcptr b)
{
    ASSERT_ALWAYS(a->n == b->n);
    ASSERT_ALWAYS(a->m == b->m);
    if (a->size != b->size) return (a->size > b->size) - (b->size > a->size);
    return abvec_cmp(ab, a->x, b->x, a->m*a->n*a->size);
}


/* polymat didn't have a proper add(), which is somewhat curious */
void matpoly_add(abdst_field ab,
        matpoly c,
        matpoly a,
        matpoly b)
{
    ASSERT_ALWAYS(c->size == a->size);
    ASSERT_ALWAYS(c->size == b->size);
    for(unsigned int i = 0 ; i < a->m ; i++) {
        for(unsigned int j = 0 ; j < a->n ; j++) {
            abvec_add(ab,
                    matpoly_part(ab, c, i, j, 0),
                    matpoly_part(ab, a, i, j, 0),
                    matpoly_part(ab, b, i, j, 0), c->size);
        }
    }
}

/* shift by a multiplication by x all coefficients of degree k and above
 * in column j. (that is, whole polynomial entries if k==0)
 */
void matpoly_multiply_column_by_x(abdst_field ab, matpoly_ptr pi, unsigned int j, unsigned int size)/*{{{*/
{
    ASSERT_ALWAYS((size + 1) <= pi->alloc);
    for(unsigned int i = 0 ; i < pi->m ; i++) {
        memmove(matpoly_part(ab, pi, i, j, 1), matpoly_part(ab, pi, i, j, 0), 
                size * abvec_elt_stride(ab, 1));
        abset_ui(ab, matpoly_coeff(ab, pi, i, j, 0), 0);
    }
}/*}}}*/
void matpoly_divide_column_by_x(abdst_field ab, matpoly_ptr pi, unsigned int j, unsigned int size)/*{{{*/
{
    ASSERT_ALWAYS(size <= pi->alloc);
    for(unsigned int i = 0 ; i < pi->m ; i++) {
        memmove(matpoly_part(ab, pi, i, j, 0), matpoly_part(ab, pi, i, j, 1), 
                (size-1) * abvec_elt_stride(ab, 1));
        abset_ui(ab, matpoly_coeff(ab, pi, i, j, size-1), 0);
    }
}/*}}}*/

void matpoly_truncate(abdst_field ab, matpoly_ptr dst, matpoly_srcptr src, unsigned int size)/*{{{*/
{
    ASSERT_ALWAYS(size <= src->alloc);
    if (dst == src) {
        /* Never used by the code, so far. We're leaving garbage coeffs
         * on top, could this be a problem ? */
        dst->size = size;
        return;
    }
    if (matpoly_check_pre_init(dst)) {
        matpoly_init(ab, dst, src->m, src->n, size);
    }
    ASSERT_ALWAYS(dst->m == src->m);
    ASSERT_ALWAYS(dst->n == src->n);
    ASSERT_ALWAYS(size <= dst->alloc);
    dst->size = size;
    /* XXX Much more cumbersome here than for polymat, of course */
    for(unsigned int i = 0 ; i < src->m ; i++) {
        for(unsigned int j = 0 ; j < src->n ; j++) {
            abvec_set(ab,
                    matpoly_part(ab, dst, i, j, 0),
                    matpoly_part_const(ab, src, i, j, 0),
                    size);
        }
    }
}/*}}}*/

/* XXX compared to polymat, our diffferent stride has a consequence,
 * clearly ! */
void matpoly_extract_column(abdst_field ab,/*{{{*/
        matpoly_ptr dst, unsigned int jdst, unsigned int kdst,
        matpoly_srcptr src, unsigned int jsrc, unsigned int ksrc)
{
    ASSERT_ALWAYS(dst->m == src->m);
    for(unsigned int i = 0 ; i < src->m ; i++)
        abset(ab,
            matpoly_coeff(ab, dst, i, jdst, kdst),
            matpoly_coeff_const(ab, src, i, jsrc, ksrc));
}/*}}}*/

void matpoly_transpose_dumb(abdst_field ab, matpoly_ptr dst, matpoly_srcptr src) /*{{{*/
{
    if (dst == src) {
        matpoly tmp;
        matpoly_init(ab, tmp, 0, 0, 0);  /* pre-init for now */
        matpoly_transpose_dumb(ab, tmp, src);
        matpoly_swap(dst, tmp);
        matpoly_clear(ab, tmp);
        return;
    }
    if (!matpoly_check_pre_init(dst))
        matpoly_clear(ab, dst);
    matpoly_init(ab, dst, src->n, src->m, src->size);
    dst->size = src->size;
    for(unsigned int i = 0 ; i < src->m ; i++) {
        for(unsigned int j = 0 ; j < src->n ; j++) {
            for(unsigned int k = 0 ; k < src->size ; k++) {
                abset(ab,
                        matpoly_coeff(ab, dst, j, i, k),
                        matpoly_coeff_const(ab, src, i, j, k));
            }
        }
    }
}/*}}}*/

void matpoly_zero_column(abdst_field ab,/*{{{*/
        matpoly_ptr dst, unsigned int jdst, unsigned int kdst)
{
    for(unsigned int i = 0 ; i < dst->m ; i++)
        abset_zero(ab,
            matpoly_coeff(ab, dst, i, jdst, kdst));
}/*}}}*/

void matpoly_extract_row_fragment(abdst_field ab,/*{{{*/
        matpoly_ptr dst, unsigned int i1, unsigned int j1,
        matpoly_srcptr src, unsigned int i0, unsigned int j0,
        unsigned int n)
{
    ASSERT_ALWAYS(src->size <= dst->alloc);
    ASSERT_ALWAYS(dst->size == src->size);
    for(unsigned int k = 0 ; k < n ; k++)
        abvec_set(ab,
                matpoly_part(ab, dst, i1, j1 + k, 0),
                matpoly_part_const(ab, src, i0, j0 + k, 0), dst->size);
}/*}}}*/

void matpoly_rshift(abdst_field ab, matpoly_ptr dst, matpoly_srcptr src, unsigned int k)/*{{{*/
{
    ASSERT_ALWAYS(k <= src->size);
    unsigned int newsize = src->size - k;
    if (matpoly_check_pre_init(dst)) {
        matpoly_init(ab, dst, src->m, src->n, newsize);
    }
    ASSERT_ALWAYS(dst->m == src->m);
    ASSERT_ALWAYS(dst->n == src->n);
    ASSERT_ALWAYS(newsize <= dst->alloc);
    dst->size = newsize;
    for(unsigned int i = 0 ; i < src->m ; i++) {
        for(unsigned int j = 0 ; j < src->n ; j++) {
            abvec_set(ab,
                    matpoly_part(ab, dst, i, j, 0),
                    matpoly_part_const(ab, src, i, j, k),
                    newsize);
        }
    }
    if (dst == src)
        matpoly_realloc(ab, dst, newsize);
}/*}}}*/

double matpoly_addmul(abdst_field ab, matpoly c, matpoly a, matpoly b, int draft)/*{{{*/
{
    size_t csize = a->size + b->size; csize -= (csize > 0);
    if (c == a || c == b) {
        matpoly tc;
        matpoly_init(ab, tc, a->m, b->n, csize);
        double t = matpoly_addmul(ab, tc, a, b, draft);
        matpoly_swap(tc, c);
        matpoly_clear(ab, tc);
        return t;
    }
    ASSERT_ALWAYS(a->n == b->m);
    if (matpoly_check_pre_init(c)) {
        matpoly_init(ab, c, a->m, b->n, csize);
    }
    ASSERT_ALWAYS(c->m == a->m);
    ASSERT_ALWAYS(c->n == b->n);
    ASSERT_ALWAYS(c->alloc >= csize);

    for(unsigned int i = 0 ; i < c->m ; ++i) {
        for(unsigned int j = 0 ; j < c->n ; ++j) {
            abvec_set_zero(ab,
                    matpoly_part(ab, c, i, j, c->size),
                    csize - c->size);
        }
    }

    c->size = csize;
    abvec_ur tmp[2];
    abvec_ur_init(ab, &tmp[0], c->size);
    abvec_ur_init(ab, &tmp[1], c->size);

    double tt0=0, tt1=0, tt=0;  /* placate gcc */
    unsigned long c0 = 0, c1 = 0;
    int set = 0;

    for(unsigned int i = 0 ; i < a->m && !set; i++) {
        for(unsigned int j = 0 ; j < b->n&& !set ; j++) {
            if (draft) { tt=wct_seconds(); tt0 -= tt; }
            abvec_ur_set_vec(ab, tmp[1], matpoly_part(ab, c, i, j, 0), c->size);
            if (draft) { tt=wct_seconds(); tt0 += tt; tt1 -= tt; }
            for(unsigned int k = 0 ; k < a->n && !set; k++) {
                abvec_conv_ur(ab, tmp[0],
                        matpoly_part(ab, a, i, k, 0), a->size,
                        matpoly_part(ab, b, k, j, 0), b->size);
                abvec_ur_add(ab, tmp[1], tmp[1], tmp[0], c->size);
                if (draft) { ++c1; set = (tt1 + (tt = wct_seconds())) >= draft; }
            }
            if (draft) { tt1 += tt; tt0 -= tt; }
            abvec_reduce(ab, matpoly_part(ab, c, i, j, 0), tmp[1], c->size);
            if (draft) { ++c0; tt0 += wct_seconds(); set = set || tt0 >= draft; }
        }
    }
    abvec_ur_clear(ab, &tmp[0], c->size);
    abvec_ur_clear(ab, &tmp[1], c->size);

    if (!set) return 0;

    /* At this point we've done c0 load & reduce, and c1 products.
     * How much extra time would we have spent here, had we done the
     * computation in full ?
     */
    unsigned long n0 = (unsigned long) a->m * (unsigned long) c->n;
    unsigned long n1 = n0 * (unsigned long) a->n;
    double x0 = tt0 / c0 * (n0 - c0);
    double x1 = tt1 / c1 * (n1 - c1);
    return x0 + x1;
}/*}}}*/

double matpoly_mul(abdst_field ab, matpoly c, matpoly a, matpoly b, int draft)/*{{{*/
{
    size_t csize = a->size + b->size; csize -= (csize > 0);

    if (c == a || c == b) {
        matpoly tc;
        matpoly_init(ab, tc, a->m, b->n, csize);
        double t = matpoly_mul(ab, tc, a, b, draft);
        matpoly_swap(tc, c);
        matpoly_clear(ab, tc);
        return t;
    }
        
    ASSERT_ALWAYS(a->n == b->m);
    if (matpoly_check_pre_init(c)) {
        matpoly_init(ab, c, a->m, b->n, csize);
    }
    ASSERT_ALWAYS(c->m == a->m);
    ASSERT_ALWAYS(c->n == b->n);
    ASSERT_ALWAYS(c->alloc >= csize);
    c->size = csize;
    abvec_set_zero(ab, c->x, c->m * c->n * c->size);
    return matpoly_addmul(ab, c, a, b, draft);
}/*}}}*/

double matpoly_addmp(abdst_field ab, matpoly b, matpoly a, matpoly c, int draft)/*{{{*/
{
    size_t fullsize = a->size + c->size; fullsize -= (fullsize > 0);
    unsigned int nb = MAX(a->size, c->size) - MIN(a->size, c->size) + 1;
    ASSERT_ALWAYS(a->n == c->m);
    if (matpoly_check_pre_init(b)) {
        matpoly_init(ab, b, a->m, c->n, nb);
    }
    ASSERT_ALWAYS(b->m == a->m);
    ASSERT_ALWAYS(b->n == c->n);
    ASSERT_ALWAYS(b->alloc >= nb);
    b->size = nb;

    /* We are going to make it completely stupid for a beginning. */
    /* XXX XXX XXX FIXME !!! */
    abvec_ur tmp[2];
    abvec_ur_init(ab, &tmp[0], fullsize);
    abvec_ur_init(ab, &tmp[1], b->size);

    double tt0=0, tt1=0, tt=0;  /* placate gcc */
    unsigned long c0 = 0, c1 = 0;
    int set = 0;

    for(unsigned int i = 0 ; i < a->m && !set ; i++) {
        for(unsigned int j = 0 ; j < c->n && !set ; j++) {
            if (draft) { tt=wct_seconds(); tt0 -= tt; }
            abvec_ur_set_vec(ab, tmp[1], matpoly_part(ab, b, i, j, 0), b->size);
            if (draft) { tt=wct_seconds(); tt0 += tt; tt1 -= tt; }
            for(unsigned int k = 0 ; k < a->n && !set ; k++) {
                abvec_conv_ur(ab, tmp[0],
                        matpoly_part(ab, a, i, k, 0), a->size,
                        matpoly_part(ab, c, k, j, 0), c->size);
                abvec_ur_add(ab, tmp[1], tmp[1],
                        abvec_ur_subvec(ab, tmp[0], MIN(a->size, c->size) - 1),
                        nb);
                if (draft) { ++c1; set = (tt1 + (tt = wct_seconds())) >= 1; }
            }
            if (draft) { tt1 += tt; tt0 -= tt; }
            abvec_reduce(ab, matpoly_part(ab, b, i, j, 0), tmp[1], b->size);
            if (draft) { ++c0; tt0 += wct_seconds(); set = set || tt0 >= 1; }
        }
    }
    abvec_ur_clear(ab, &tmp[0], fullsize);
    abvec_ur_clear(ab, &tmp[1], b->size);

    if (!set) return 0;

    /* At this point we've done c0 load & reduce, and c1 products.
     * How much extra time would we have spent here, had we done the
     * computation in full ?
     */
    unsigned long n0 = (unsigned long) a->m * (unsigned long) c->n;
    unsigned long n1 = n0 * (unsigned long) a->n;
    double x0 = tt0 / c0 * (n0 - c0);
    double x1 = tt1 / c1 * (n1 - c1);
    return x0 + x1;
}/*}}}*/

double matpoly_mp(abdst_field ab, matpoly b, matpoly a, matpoly c, int draft)/*{{{*/
{
    unsigned int nb = MAX(a->size, c->size) - MIN(a->size, c->size) + 1;
    ASSERT_ALWAYS(a->n == c->m);
    if (matpoly_check_pre_init(b)) {
        matpoly_init(ab, b, a->m, c->n, nb);
    }
    ASSERT_ALWAYS(b->m == a->m);
    ASSERT_ALWAYS(b->n == c->n);
    ASSERT_ALWAYS(b->alloc >= nb);
    b->size = nb;
    abvec_set_zero(ab, b->x, b->m * b->n * b->size);
    return matpoly_addmp(ab, b, a, c, draft);
}/*}}}*/


void matpoly_set_polymat(abdst_field ab MAYBE_UNUSED, matpoly_ptr dst, polymat_srcptr src)
{
    if (!matpoly_check_pre_init(dst))
        matpoly_clear(ab, dst);
    matpoly_init(ab, dst, src->m, src->n, src->size);
    dst->size = src->size;

    for(unsigned int i = 0 ; i < src->m ; i++) {
        for(unsigned int j = 0 ; j < src->n ; j++) {
            for(unsigned int k = 0 ; k < src->size ; k++) {
                abset(ab,
                        matpoly_coeff(ab, dst, i, j, k),
                        polymat_coeff_const(ab, src, i, j, k));
            }
        }
    }
}

int matpoly_coeff_is_zero(abdst_field ab, matpoly_srcptr pi, unsigned int k)
{
    for(unsigned int j = 0; j < pi->m; j++)
        for(unsigned int i = 0 ; i < pi->n ; i++)
            if (!abis_zero(ab, matpoly_coeff_const(ab, pi, i, j, k)))
                return 0;
    return 1;
}
void matpoly_coeff_set_zero(abdst_field ab, matpoly_ptr pi, unsigned int k)
{
    for(unsigned int j = 0; j < pi->m; j++)
        for(unsigned int i = 0 ; i < pi->n ; i++)
            abset_zero(ab, matpoly_coeff(ab, pi, i, j, k));
}
