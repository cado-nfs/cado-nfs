#include "cado.h" // IWYU pragma: keep
#include "mpfq_layer.h"
#include <stdlib.h>
#include <limits.h>
#include <stdio.h>
#include <gmp.h>
#include "macros.h"
#include "lingen-polymat.h"

#define POLYMAT_MUL_KARA_CUTOFF_DEFAULT { .cut = 10, .subdivide = 10, .table = NULL, .table_size = 0}
#define POLYMAT_MP_KARA_CUTOFF_DEFAULT { .cut = 10, .subdivide = 10, .table = NULL, .table_size = 0}

/*{{{ cutoffs and such */
struct polymat_cutoff_info polymat_mul_kara_cutoff = POLYMAT_MUL_KARA_CUTOFF_DEFAULT;

/* For the middle product, this is similar, but applies to the small
 * operand length and to the result length
 * (A balanced middle product is n times 2*n - 1, and gives n result
 * coefficients. Here the important data is n).
 * */
/* FIXME: polymat_mul_kara_threshold and polymat_mp_kara_threshold should
 * not differ, really. The only thing is that dimensions for mp and mul
 * are not the same, hence the difference...
 */
struct polymat_cutoff_info polymat_mp_kara_cutoff = POLYMAT_MP_KARA_CUTOFF_DEFAULT;


/* This sets the cutoffs. The table, if present in new_cutoff, is copied,
 * and hence must be freed by the caller who allocated it in the first
 * place.
 * The old_cutoff (if not NULL) value should be initialized already, and
 * receives initialized data.
 */

void polymat_cutoff_info_init(struct polymat_cutoff_info * c)
{
    memset(c, 0, sizeof(*c));
    c->cut = UINT_MAX;
}

void polymat_cutoff_info_clear(struct polymat_cutoff_info * c)
{
    if (c->table) free(c->table);
    c->cut = UINT_MAX;
}

void polymat_cutoff_add_step(struct polymat_cutoff_info * c, unsigned int size, unsigned int alg)
{
    c->table = realloc(c->table, (c->table_size + 1) * sizeof(unsigned int [2]));
    if (c->table) ASSERT_ALWAYS(size > c->table[c->table_size-1][0]);
    c->table[c->table_size][0] = size;
    c->table[c->table_size][1] = alg;
    c->table_size++;
}

static void polymat_set_generic_cutoff(struct polymat_cutoff_info * slot, const struct polymat_cutoff_info * new_cutoff, struct polymat_cutoff_info * old_cutoff)
{
    if (old_cutoff) {
        if (old_cutoff->table) free(old_cutoff->table);
        memcpy(old_cutoff, slot, sizeof(struct polymat_cutoff_info));
    } else if (slot->table) {
        free(slot->table);
    }
    memcpy(slot, new_cutoff, sizeof(struct polymat_cutoff_info));
    if (slot->table) {
        slot->table = malloc(slot->table_size * sizeof(unsigned int [2]));
        memcpy(slot->table, new_cutoff->table, slot->table_size * sizeof(unsigned int [2]));
    }
}

void polymat_set_mul_kara_cutoff(const struct polymat_cutoff_info * new_cutoff, struct polymat_cutoff_info * old_cutoff)
{
    polymat_set_generic_cutoff(&polymat_mul_kara_cutoff, new_cutoff, old_cutoff);
}

void polymat_set_mp_kara_cutoff(const struct polymat_cutoff_info * new_cutoff, struct polymat_cutoff_info * old_cutoff)
{
    polymat_set_generic_cutoff(&polymat_mp_kara_cutoff, new_cutoff, old_cutoff);
}

int polymat_cutoff_get_alg_b(const struct polymat_cutoff_info * cutoff, unsigned int s)
{
    if (s >= cutoff->cut) return 1;
    if (!cutoff->table) return 0;
    unsigned int a = 0;
    unsigned int b = cutoff->table_size;
    for( ; b-a > 1 ; ) {
        unsigned int c = (a+b)/2;
        if (s >= cutoff->table[c][0]) {
            a = c;
        } else {
            b = c;
        }
    }
    return cutoff->table[a][1];
}

int polymat_cutoff_get_subdivide_ub(const struct polymat_cutoff_info * cutoff, unsigned int s0, unsigned int s1)
{
    return cutoff->table != NULL && MIN(s0, s1) >= cutoff->subdivide;
}

/*}}}*/

/* {{{ unreduced interface, for the implementation only */
struct polymat_ur_s {
    unsigned int m;
    unsigned int n;
    size_t size;
    size_t alloc;
    abvec_ur x;
};
typedef struct polymat_ur_s polymat_ur[1];
typedef struct polymat_ur_s * polymat_ur_ptr;
typedef const struct polymat_ur_s * polymat_ur_srcptr;

static void polymat_ur_init(abdst_field ab, polymat_ur_ptr p, unsigned int m, unsigned int n, int len);
static void polymat_ur_zero(abdst_field ab, polymat_ur_ptr p);
static void polymat_ur_clear(abdst_field ab, polymat_ur_ptr p);
static inline abvec_ur polymat_ur_part(abdst_field ab, polymat_ur_ptr p, unsigned int i, unsigned int j, unsigned int k);
static inline abdst_elt_ur polymat_ur_coeff(abdst_field ab, polymat_ur_ptr p, unsigned int i, unsigned int j, unsigned int k);

/* {{{ access interface for polymat_ur */
static inline abdst_vec_ur polymat_ur_part(abdst_field ab, polymat_ur_ptr p, unsigned int i, unsigned int j, unsigned int k) {
    /* Assume row-major in all circumstances. Old code used to support
     * various orderings, here we don't */
    return abvec_ur_subvec(ab, p->x, (k*p->m+i)*p->n+j);
}
static inline abdst_elt_ur polymat_ur_coeff(abdst_field ab, polymat_ur_ptr p, unsigned int i, unsigned int j, unsigned int k) {
    return abvec_ur_coeff_ptr(ab, polymat_ur_part(ab, p,i,j,k), 0);
}
void polymat_reducemat(abdst_field ab,
        polymat c, unsigned int kc,
        polymat_ur a, unsigned int ka);

/* }}} */
/* }}} */

int polymat_check_pre_init(polymat_srcptr p)
{
    if (p->m && p->n && p->alloc)
        return 0;
    if (!p->m && !p->n && !p->alloc)
        return 1;
    abort();
    return 0;
}

/* {{{ init/zero/clear interface for polymat */
void polymat_init(abdst_field ab, polymat_ptr p, unsigned int m, unsigned int n, int len) {
    memset(p, 0, sizeof(polymat));
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
void polymat_realloc(abdst_field ab, polymat_ptr p, size_t newalloc) {
    polymat_check_pre_init(p);
    abvec_reinit(ab, &(p->x), p->m*p->n*p->alloc, p->m*p->n*newalloc);
    /* zero out the newly added data */
    if (newalloc > p->alloc) {
        abvec_set_zero(ab, p->x + p->m*p->n*p->alloc, p->m*p->n*(newalloc - p->alloc));
    } else {
        ASSERT_ALWAYS(p->size <= newalloc);
    }
    p->alloc = newalloc;
}
void polymat_zero(abdst_field ab, polymat_ptr p) {
    p->size = 0;
    abvec_set_zero(ab, p->x, p->m*p->n*p->alloc);
}
void polymat_clear(abdst_field ab, polymat_ptr p) {
    abvec_clear(ab, &(p->x), p->m*p->n*p->alloc);
    memset(p, 0, sizeof(polymat));
}
void polymat_swap(polymat_ptr a, polymat_ptr b)
{
    polymat x;
    memcpy(x, a, sizeof(polymat));
    memcpy(a, b, sizeof(polymat));
    memcpy(b, x, sizeof(polymat));
}

void polymat_fill_random(abdst_field ab MAYBE_UNUSED, polymat_ptr a, unsigned int size, gmp_randstate_t rstate)
{
    ASSERT_ALWAYS(size <= a->alloc);
    a->size = size;
    abvec_random(ab, a->x, a->m*a->n*size, rstate);
}

int polymat_cmp(abdst_field ab MAYBE_UNUSED, polymat_srcptr a, polymat_srcptr b)
{
    ASSERT_ALWAYS(a->n == b->n);
    ASSERT_ALWAYS(a->m == b->m);
    if (a->size != b->size) return (a->size > b->size) - (b->size > a->size);
    return abvec_cmp(ab, a->x, b->x, a->m*a->n*a->size);
}


/* }}} */

/* {{{ init/zero/clear interface for polymat_ur */
static void polymat_ur_init(abdst_field ab, polymat_ur_ptr p, unsigned int m, unsigned int n, int len) {
    memset(p, 0, sizeof(polymat_ur));
    p->m = m;
    p->n = n;
    p->alloc = len;
    p->size = 0;
    abvec_ur_init(ab, &(p->x), m*n*p->alloc);
    abvec_ur_set_zero(ab, p->x, m*n*p->alloc);
}
static void polymat_ur_zero(abdst_field ab, polymat_ur_ptr p) {
    p->size = 0;
    abvec_ur_set_zero(ab, p->x, p->m*p->n*p->alloc);
}
static void polymat_ur_clear(abdst_field ab, polymat_ur_ptr p) {
    abvec_ur_clear(ab, &(p->x), p->m*p->n*p->alloc);
    memset(p, 0, sizeof(polymat_ur));
}
/* }}} */

/* It's used from lingen-bigpolymat.c as well */
void bwmat_copy_coeffs(abdst_field ab MAYBE_UNUSED, abvec x0, int stride0, abvec x1, int stride1, unsigned int n)
{
    for(unsigned int i = 0 ; i < n ; i++) {
        abset(ab,
                    abvec_coeff_ptr(ab, x0, i * stride0),
                    abvec_coeff_ptr(ab, x1, i * stride1));
    }
}
void bwmat_zero_coeffs(abdst_field ab MAYBE_UNUSED, abvec x0, int stride0, unsigned int n)
{
    for(unsigned int i = 0 ; i < n ; i++) {
        abset_zero(ab, 
                abvec_coeff_ptr(ab, x0, i * stride0));
    }
}
void bwmat_move_coeffs(abdst_field ab MAYBE_UNUSED, abvec x0, int stride0, abvec x1, int stride1, unsigned int n)
{
    ASSERT_ALWAYS(stride0 == stride1); /* Otherwise there's probably no point */
    if (x0 < x1) {
        for(unsigned int i = 0 ; i < n ; i++) {
            abset(ab,
                    abvec_coeff_ptr(ab, x0, i * stride0),
                    abvec_coeff_ptr(ab, x1, i * stride1));
        }
    } else {
        for(unsigned int i = n ; i-- ; ) {
            abset(ab,
                    abvec_coeff_ptr(ab, x0, i * stride0),
                    abvec_coeff_ptr(ab, x1, i * stride1));
        }
    }
}

/* {{{ Elementary operations on the matrices, which are the coefficients
 * of our polynomials */
void polymat_addmat(abdst_field ab,
        polymat c, unsigned int kc,
        polymat a, unsigned int ka,
        polymat b, unsigned int kb)
{
    for(unsigned int i = 0 ; i < a->m ; i++) {
        for(unsigned int j = 0 ; j < a->n ; j++) {
            abadd(ab,
                    polymat_coeff(ab, c, i, j, kc),
                    polymat_coeff(ab, a, i, j, ka),
                    polymat_coeff(ab, b, i, j, kb));
        }
    }
}

void polymat_submat(abdst_field ab,
        polymat c, unsigned int kc,
        polymat a, unsigned int ka,
        polymat b, unsigned int kb)
{
    for(unsigned int i = 0 ; i < a->m ; i++) {
        for(unsigned int j = 0 ; j < a->n ; j++) {
            absub(ab,
                    polymat_coeff(ab, c, i, j, kc),
                    polymat_coeff(ab, a, i, j, ka),
                    polymat_coeff(ab, b, i, j, kb));
        }
    }
}

void polymat_addmulmat_ur(abdst_field ab,
        polymat_ur c, unsigned int kc,
        polymat a, unsigned int ka,
        polymat b, unsigned int kb)
{
    abelt_ur tmp_ur;
    abelt_ur_init(ab, &tmp_ur);
    ASSERT_ALWAYS(a->n == b->m);
    for(unsigned int i = 0 ; i < a->m ; i++) {
        for(unsigned int j = 0 ; j < b->n ; j++) {
            for(unsigned int k = 0 ; k < a->n ; k++) {
                absrc_elt ea = polymat_coeff(ab, a, i, k, ka);
                absrc_elt eb = polymat_coeff(ab, b, k, j, kb);
                abmul_ur(ab, tmp_ur, ea, eb);
                abelt_ur_add(ab,
                        polymat_ur_coeff(ab, c, i, j, kc),
                        polymat_ur_coeff(ab, c, i, j, kc),
                        tmp_ur);
            }
        }
    }
    abelt_ur_clear(ab, &tmp_ur);
}

void polymat_reducemat(abdst_field ab,
        polymat c, unsigned int kc,
        polymat_ur a, unsigned int ka)
{
    for(unsigned int i = 0 ; i < a->m ; i++) {
        for(unsigned int j = 0 ; j < a->n ; j++) {
            abreduce(ab, 
                    polymat_coeff(ab, c, i, j, kc),
                    polymat_ur_coeff(ab, a, i, j, ka));
        }
    }
}

void polymat_mulmat(abdst_field ab,
        polymat c, unsigned int kc,
        polymat a, unsigned int ka,
        polymat b, unsigned int kb)
{
    polymat_ur cx;
    polymat_ur_init(ab, cx, c->m, c->n, 1);
    polymat_addmulmat_ur(ab, cx, 0, a, ka, b, kb);
    polymat_reducemat(ab, c, kc, cx, 0);
    polymat_ur_clear(ab, cx);
}

void polymat_addmulmat(abdst_field ab,
        polymat c, unsigned int kc,
        polymat a, unsigned int ka,
        polymat b, unsigned int kb)
{
    polymat cc;
    polymat_init(ab, cc, c->m, c->n, 1);
    polymat_mulmat(ab,cc,0,a,ka,b,kb);
    polymat_addmat(ab,c,kc,c,kc,cc,0);
    polymat_clear(ab, cc);
}

/* }}} */

void polymat_multiply_column_by_x(abdst_field ab, polymat_ptr pi, unsigned int j, unsigned int size)/*{{{*/
{
    ASSERT_ALWAYS((size + 1) <= pi->alloc);
    bwmat_move_coeffs(ab,
            polymat_part(ab, pi, 0, j, 1), pi->n,
            polymat_part(ab, pi, 0, j, 0), pi->n,
            pi->m * size);
    for(unsigned int i = 0 ; i < pi->m ; i++) {
        abset_ui(ab, polymat_coeff(ab, pi, i, j, 0), 0);
    }
}/*}}}*/

void polymat_truncate(abdst_field ab, polymat_ptr dst, polymat_ptr src, unsigned int size)/*{{{*/
{
    ASSERT_ALWAYS(size <= src->alloc);
    if (dst == src) {
        /* Never used by the code, so far. We're leaving garbage coeffs
         * on top, could this be a problem ? */
        abort(); /* so that programmer becomes aware of the gotcha */
        dst->size = size;
        return;
    }
    if (polymat_check_pre_init(dst)) {
        polymat_init(ab, dst, src->m, src->n, size);
    }
    ASSERT_ALWAYS(dst->m == src->m);
    ASSERT_ALWAYS(dst->n == src->n);
    ASSERT_ALWAYS(size <= dst->alloc);
    dst->size = size;
    bwmat_copy_coeffs(ab,
            polymat_part(ab, dst,0,0,0),1,
            polymat_part(ab, src,0,0,0),1,
            src->m*src->n*size);
}/*}}}*/

void polymat_extract_column(abdst_field ab,/*{{{*/
        polymat_ptr dst, unsigned int jdst, unsigned int kdst,
        polymat_ptr src, unsigned int jsrc, unsigned int ksrc)
{
    ASSERT_ALWAYS(dst->m == src->m);
    bwmat_copy_coeffs(ab,
            polymat_part(ab, dst, 0, jdst, kdst), dst->n,
            polymat_part(ab, src, 0, jsrc, ksrc), src->n, src->m);
}/*}}}*/

void polymat_extract_row_fragment(abdst_field ab,/*{{{*/
        polymat_ptr dst, unsigned int i1, unsigned int j1,
        polymat_ptr src, unsigned int i0, unsigned int j0,
        unsigned int n)
{
    ASSERT_ALWAYS(src->size <= dst->alloc);
    for(unsigned int k = 0 ; k < src->size ; k++) {
        bwmat_copy_coeffs(ab,
                polymat_part(ab, dst, i1, j1, k), 1,
                polymat_part(ab, src, i0, j0, k), 1,
                n);
    }
}/*}}}*/


void polymat_rshift(abdst_field ab, polymat_ptr dst, polymat_ptr src, unsigned int k)/*{{{*/
{
    ASSERT_ALWAYS(k <= src->size);
    unsigned int newsize = src->size - k;
    if (dst != src) {
        if (polymat_check_pre_init(dst)) {
            polymat_init(ab, dst, src->m, src->n, newsize);
        }
        ASSERT_ALWAYS(dst->m == src->m);
        ASSERT_ALWAYS(dst->n == src->n);
        ASSERT_ALWAYS(newsize <= dst->alloc);
        dst->size = newsize;
        bwmat_copy_coeffs(ab,
                polymat_part(ab, dst,0,0,0),1,
                polymat_part(ab, src,0,0,k),1,
                src->m*src->n*newsize);
    } else {
        dst->size = newsize;
        bwmat_move_coeffs(ab,
                polymat_part(ab, dst,0,0,0),1,
                polymat_part(ab, src,0,0,k),1,
                src->m*src->n*newsize);
    }
    polymat_realloc(ab, dst, newsize);
}/*}}}*/


void polymat_mul_raw_basecase(abdst_field ab,/*{{{*/
        polymat c, unsigned int xc,
        polymat a, unsigned int xa, unsigned int na,
        polymat b, unsigned int xb, unsigned int nb,
        int transpose, int add)
{
    unsigned int nc = na + nb - 1;
    ASSERT_ALWAYS(c->alloc >= xc + nc);

    polymat_ur tmat_ur;
    polymat_ur_init(ab, tmat_ur, c->m, c->n, 1);

    for(unsigned int k = 0 ; k < nc ; k++) {
        unsigned int i0 = k >= nb ? k + 1 - nb : 0;
        unsigned int i1 = k + 1 < na ? k + 1 : na;
        polymat_ur_zero(ab, tmat_ur);
        for(unsigned int i = i0 ; i < i1 ; i++) {
            unsigned int j = k - i;
            if (!transpose) {
                polymat_addmulmat_ur(ab, tmat_ur, 0, a, xa + i, b, xb + j);
            } else {
                polymat_addmulmat_ur(ab, tmat_ur, 0, b, xb + j, a, xa + i);
            }
        }
        if (!add) {
            polymat_reducemat(ab, c, xc + k, tmat_ur, 0);
        } else {
            polymat tmat;
            polymat_init(ab, tmat, c->m, c->n, 1);
            tmat->size = 1;
            polymat_reducemat(ab, tmat, 0, tmat_ur, 0);
            polymat_addmat(ab, c, xc + k, c, xc + k, tmat, 0);
            polymat_clear(ab, tmat);
        }
    }
    polymat_ur_clear(ab, tmat_ur);
}/*}}}*/

void polymat_mul_raw_kara(abdst_field ab,/*{{{*/
        polymat c, unsigned int xc,
        polymat a, unsigned int xa, unsigned int na,
        polymat b, unsigned int xb, unsigned int nb,
        int transpose, int add)
{
    /* This code works __only__ for kara-sized products */
    ASSERT_ALWAYS(nb == na);

    if (polymat_cutoff_get_alg_b(&polymat_mul_kara_cutoff, MIN(na, nb)) == 0) {
        polymat_mul_raw_basecase(ab, c, xc, a, xa, na, b, xb, nb, transpose, add);
        return;
    }

    unsigned int nc = na + nb - 1;
    ASSERT_ALWAYS(c->alloc >= xc + nc);
    if (!add) {
        bwmat_zero_coeffs(ab, polymat_part(ab, c, 0, 0, xc), 1, c->m*c->n*nc);
    }
    unsigned int m0 = na / 2;
    unsigned int m1 = na - m0;

    polymat q0, q2, s, t;
    polymat_init(ab, s, a->m, a->n, m1);
    polymat_init(ab, t, b->m, b->n, m1);
    s->size = m1;
    t->size = m1;
    for(unsigned int k = 0 ; k < m0 ; k++) {
        polymat_addmat(ab, s, k, a, xa + k, a, xa + m0 + k);
        polymat_addmat(ab, t, k, b, xb + k, b, xb + m0 + k);
    }
    if (m1 > m0) {
        polymat_addmat(ab, s, m0, s, m0, a, xa + m0 + m0);
        polymat_addmat(ab, t, m0, t, m0, b, xb + m0 + m0);
    }
    polymat_mul_raw_kara(ab, c, xc+m0, s, 0, m1, t, 0, m1, transpose, 1);
    polymat_clear(ab, s);
    polymat_clear(ab, t);

    polymat_init(ab, q0, c->m, c->n, 2*m0-1);
    q0->size = 2*m0-1;
    polymat_init(ab, q2, c->m, c->n, 2*m1-1);
    q2->size = 2*m1-1;
    polymat_mul_raw_kara(ab, q0, 0, a, xa, m0, b, xb, m0, transpose, 0);
    polymat_mul_raw_kara(ab, q2, 0, a, xa+m0, m1, b, xb+m0, m1, transpose, 0);
    for(unsigned int k = 0 ; k < 2*m0 - 1 ; k++) {
        polymat_addmat(ab, c, xc + k, c, xc + k, q0, k);
        polymat_submat(ab, c, xc + m0 + k, c, xc + m0 + k, q0, k);
    }
    for(unsigned int k = 0 ; k < 2*m1 - 1 ; k++) {
        polymat_addmat(ab, c, xc + 2 * m0 + k, c, xc + 2 * m0 + k, q2, k);
        polymat_submat(ab, c, xc + m0 + k, c, xc + m0 + k, q2, k);
    }
    polymat_clear(ab, q0);
    polymat_clear(ab, q2);
} /* }}} */

void polymat_mul_raw_subdivide(abdst_field ab,/*{{{*/
        polymat c, unsigned int xc,
        polymat a, unsigned int xa, unsigned int na,
        polymat b, unsigned int xb, unsigned int nb,
        int transpose, int add)
{
    if (!na || !nb) {
        return;
    }
    /* we consider using karatsuba only when the smallest of the two
     * operands is large enough.
     */
    if (polymat_cutoff_get_subdivide_ub(&polymat_mul_kara_cutoff, na, nb) == 0) {
        polymat_mul_raw_basecase(ab, c, xc, a, xa, na, b, xb, nb, transpose, add);
        return;
    }
    unsigned int nc = na + nb - 1;
    if (!add) {
        bwmat_zero_coeffs(ab, polymat_part(ab, c, 0, 0, xc), 1, c->m*c->n*nc);
    }

    for( ; nb >= na ; ) {
        polymat_mul_raw_kara(ab, c, xc, a, xa, na, b, xb, na, transpose, 1);
        xb += na;
        xc += na;
        nb -= na;
    }
    /* see below comment in #if0'ed block. Seems to me that the code
     * theere is bogus */
    ASSERT_ALWAYS(nb < na);
    polymat_mul_raw_subdivide(ab, c, xc, a, xa, na, b, xb, nb, transpose, 1);
#if 0
    if (!nb) {
        /* This condition also catches the Karatsuba case n, n */
        return;
    }
    /* FIXME: nb might be very small here. Why do we do
     * polymat_mul_raw_kara ? This defeats the semantics of the subdivide
     * check, apparently.
     */
    ASSERT_ALWAYS(nb < na);
    /* Fixup needed, now. Treat the largest possible subblock. We want
     * the intervals [0,na-k[ and [0, nb[ to match a kara scheme, i.e.
     * have k=nb-na.
     */
    unsigned int chop = na - nb;
    polymat_mul_raw_kara(ab, c, xc,
            a, xa, nb, b, xb, nb, transpose, 1);
    polymat_mul_raw_subdivide(ab, c, xc + na - chop,
            a, xa + nb, chop, b, xb, nb, transpose, 1);
#endif
}/*}}}*/

void polymat_mul_raw(abdst_field ab,/*{{{*/
        polymat c, unsigned int xc,
        polymat a, unsigned int xa, unsigned int na,
        polymat b, unsigned int xb, unsigned int nb,
        int transpose, int add)
{
    polymat_mul_raw_subdivide(ab, c, xc, a, xa, na, b, xb, nb, transpose, add);
}/*}}}*/

void polymat_mul(abdst_field ab, polymat c, polymat a, polymat b)/*{{{*/
{
    size_t csize = a->size + b->size; csize -= (csize > 0);
    ASSERT_ALWAYS(a->n == b->m);
    if (polymat_check_pre_init(c)) {
        polymat_init(ab, c, a->m, b->n, csize);
    }
    ASSERT_ALWAYS(c->m == a->m);
    ASSERT_ALWAYS(c->n == b->n);
    ASSERT_ALWAYS(c->alloc >= csize);
    c->size = csize;
    polymat_mul_raw(ab, c, 0, a, 0, a->size, b, 0, b->size, 0, 0);
}/*}}}*/

void polymat_addmul(abdst_field ab, polymat c, polymat a, polymat b)/*{{{*/
{
    size_t csize = a->size + b->size; csize -= (csize > 0);
    ASSERT_ALWAYS(a->n == b->m);
    if (polymat_check_pre_init(c)) {
        polymat_init(ab, c, a->m, b->n, csize);
    }
    ASSERT_ALWAYS(c->m == a->m);
    ASSERT_ALWAYS(c->n == b->n);
    ASSERT_ALWAYS(c->alloc >= csize);
    c->size = csize;
    polymat_mul_raw(ab, c, 0, a, 0, a->size, b, 0, b->size, 0, 1);
}/*}}}*/

/* Middle product b of a and c. This is the part of the product where
 * there are the largest number of summands used to accumulate the terms
 * of the result.
 *
 * Note that the first coefficient computed is stored as the coefficient
 * of degree 0 in b, even though it corresponds to something of larger
 * degree in a*c.
 *
 * We require a to be the small operand, and c the larger one (we support
 * the other case, though, by swapping arguments).
 *
 * Exactly nc-na+1 coefficients are stored in b (so that nc = na + nb -
 * 1). Those are the coefficients of degrees [na-1..nc[ in the product
 * a*c.
 *
 * Some examples:
 *
 *      na==nc==n+1 --> nb==1: [n]
 *      na=n+1, nc=2n --> nb==n: [n,2n[
 */
void polymat_mp_raw_basecase(abdst_field ab,/*{{{*/
        polymat b, unsigned int xb,
        polymat a, unsigned int xa, unsigned int na,
        polymat c, unsigned int xc, unsigned int nc,
        int transpose, int add)
{
    if (nc < na) {
        polymat_mp_raw_basecase(ab, b, xb, c, xc, nc, a, xa, na, !transpose, add);
        return;
    }
    unsigned int nb = nc - na + 1;
    ASSERT_ALWAYS(b->alloc >= xb + nb);

    polymat_ur tmat_ur;
    polymat_ur_init(ab, tmat_ur, b->m, b->n, 1);
    tmat_ur->size = 1;

    for(unsigned int j = 0 ; j < nb ; j++) {
        polymat_ur_zero(ab, tmat_ur);
        for(unsigned int i = 0 ; i < na ; i++) {
            unsigned int k = j + na - 1 - i;
            if (!transpose) {
                polymat_addmulmat_ur(ab, tmat_ur, 0, a, xa + i, c, xc + k);
            } else {
                polymat_addmulmat_ur(ab, tmat_ur, 0, c, xc + k, a, xa + i);
            }
        }
        if (!add) {
            polymat_reducemat(ab, b, xb + j, tmat_ur, 0);
        } else {
            polymat tmat;
            polymat_init(ab, tmat, b->m, b->n, 1);
            tmat->size = 1;
            polymat_reducemat(ab, tmat, 0, tmat_ur, 0);
            polymat_addmat(ab, b, xb + j, b, xb + j, tmat, 0);
            polymat_clear(ab, tmat);
        }
    }
    polymat_ur_clear(ab, tmat_ur);
}/*}}}*/

void polymat_mp_raw_kara(abdst_field ab,/*{{{*/
        polymat b, unsigned int xb,
        polymat a, unsigned int xa, unsigned int na,
        polymat c, unsigned int xc, unsigned int nc,
        int transpose, int add)
{
    /* This code works __only__ for kara-sized middle products */
    ASSERT_ALWAYS(nc == 2*na - 1 || na == 2*nc - 1);
    // printf("mp_kara(%d,%d)%b\n",na,nc,transpose?'T':' ');
    if (polymat_cutoff_get_alg_b(&polymat_mp_kara_cutoff, MIN(na, nc)) == 0) {
        polymat_mp_raw_basecase(ab, b, xb, a, xa, na, c, xc, nc, transpose, add);
        return;
    }
    if (nc < na) {
        polymat_mp_raw_kara(ab, b, xb, c, xc, nc, a, xa, na, !transpose, add);
        return;
    }
    unsigned int m0 = na / 2;
    unsigned int m1 = na - m0;
    /* Spans of different chunks, in (offset, length) format */
    unsigned int span_a0[2] =  {xa,        m0};
    unsigned int span_a1[2] =  {xa+m0,     m1};
    unsigned int span_c0[2] =  {xc,      2*m1-1};
    unsigned int span_c10[2] = {xc+m1,   2*m0-1};
    unsigned int span_c11[2] = {xc+m1,   2*m1-1};
    unsigned int span_c2[2] =  {xc+2*m1, 2*m0-1};
    /* Dammit, we need some temp storage, of course */
    polymat q1, s, t;
    /* The polymats s and t are used to build temporaries, respectively
     * from a and c, of respective maximum lengths m1 and 2*m1-1 */
    /* q0 is MP(a1, b0+b11), i.e. MP(m1, 2*m1-1), which
     * produces m1 coefficients. This goes to offset 0 (well, xb) in b. */
    polymat_init(ab, t, c->m, c->n, 2*m1-1);
    t->size = 2*m1-1;
    for(unsigned int k = 0 ; k < 2*m1 - 1 ; k++) {
        polymat_addmat(ab, t, k, c, span_c0[0] + k, c, span_c11[0] + k);
    }
    polymat_mp_raw_kara(ab, b, xb,
            a, span_a1[0], span_a1[1],
            t, 0, 2*m1-1,
            transpose, add);
    polymat_clear(ab, t);
    /* q1 is MP(a1-a0, b11), i.e. MP(m1, 2*m1-1) again */
    polymat_init(ab, q1, b->m, b->n, m1);
    polymat_init(ab, s, a->m, a->n, m1);
    q1->size = m1;
    s->size = m1;
    bwmat_copy_coeffs(ab,
            polymat_part(ab,s,0,0,0),1,
            polymat_part(ab,a,0,0,span_a1[0]),1,
            m1*a->m*a->n);
    for(unsigned int k = 0 ; k < m0 ; k++) {
        /* We've already copied the a1 coefficients above */
        polymat_submat(ab,
                s, k + m1 - m0,
                s, k + m1 - m0,
                a, span_a0[0] + k);
    }
    polymat_mp_raw_kara(ab, q1, 0,
            s, 0, m1,
            c, span_c11[0], span_c11[1],
            transpose, 0);
    polymat_clear(ab, s);
    /* q2 is MP(a0, b10+b2), i.e. MP(m0, 2*m0-1) */
    /* This goes to offset m1 in b */
    polymat_init(ab, t, c->m, c->n, 2*m0-1);
    t->size = 2*m0-1;
    for(unsigned int k = 0 ; k < 2*m0 - 1 ; k++) {
        polymat_addmat(ab, t, k, c, span_c10[0] + k, c, span_c2[0] + k);
    }
    polymat_mp_raw_kara(ab, b, xb + m1,
            a, span_a0[0], span_a0[1],
            t, 0, 2*m0-1,
            transpose, add);
    polymat_clear(ab, t);

    /* We now have to append the coefficients to form the result. Not all
     * coefficients in q1 are read for the high part of the middle
     * product. */
    {
        for(unsigned int k = 0 ; k < m1 ; k++) {
            polymat_submat(ab, b, xb + k, b, xb + k, q1, k);
        }
        for(unsigned int k = 0 ; k < m0 ; k++) {
            polymat_addmat(ab, b, xb + m1 + k, b, xb + m1 + k, q1, k);
        }
    }
    polymat_clear(ab, q1);
}/*}}}*/

void polymat_mp_raw_subdivide(abdst_field ab,/*{{{*/
        polymat b, unsigned int xb,
        polymat a, unsigned int xa, unsigned int na,
        polymat c, unsigned int xc, unsigned int nc,
        int transpose, int add)
{
    // printf("mp_subdivide(%d,%d)%b\n",na,nc,transpose?'T':' ');
    /* 
     * In block Wiedemann, nc/na is typically 1+2n/m, which may well be 3
     * for instance if m=n=1. In such a case, the present code will
     * reduce our mp to 2*3 = 6 balanced middle products
     * (Karatsuba-friendly).
     *
     * Other options could be to rely on e.g. a polymat_mp_raw_toom24,
     * which would use the following operations to fall back on only 5
     * balanced middle products:
     *
     * Formulae for MP(2k, 6k)-> 3k
     *  t0 = 1/2*(2*x0 - x1 - 2*x2 + x3)*y1
     *  t1 = 1/6*(y0 - y1)*(2*x1 - 3*x2 + x3)
     *  t2 = (2*x1 - x2 - 2*x3 + x4)*y0
     *  t3 = 1/2*(y0 + y1)*(2*x1 + x2 - x3)
     *  t4 = -1/6*(x1 - x3)*(2*y0 + y1)
     *  z0 = t0 + t1 + t3 + t4
     *  z1 = -t1 + t3 + 2*t4
     *  z2 = t1 + t3 + 4*t4
     *  z3 = -t1 + t2 + t3 + 8*t4
     *
     * The different options, for reaching balanced products of size
     * k/2,2*(k/2) from something of size k,3k are:
     *  - as we do here: Karatsuba twice: 2*3=6M
     *  - using toom24: 5M (but note the divisions by constants !)
     */
    if (!na || !nc) {
        return;
    }
    if (nc < na) {
        polymat_mp_raw_subdivide(ab, b, xb, c, xc, nc, a, xa, na, !transpose, add);
        return;
    }
    if (polymat_cutoff_get_subdivide_ub(&polymat_mp_kara_cutoff, na, nc) == 0) {
        polymat_mp_raw_basecase(ab, b, xb, a, xa, na, c, xc, nc, transpose, add);
        return;
    }
    ASSERT(nc >= na);
    for( ; nc >= 2 * na - 1 ; ) {
        polymat_mp_raw_kara(ab, b, xb,
                a, xa, na, c, xc, 2*na-1, transpose, add);
        xc += na;
        xb += na;
        nc -= na;
    }
    if (nc < na) {
        /* This condition also catches the Karatsuba case n, 2n-1 */
        return;
    }
    // printf("mp_subdivide_tail(%d,%d)%b\n",na,nc,transpose?'T':' ');
    ASSERT_ALWAYS(nc >= na);
    /* Fixup needed, now. Treat the largest possible subblock. We want
     * the intervals [0,na-k[ and [k, nc[ to match a kara scheme, i.e.
     * have the constraints 0<=k<=na, and nc-k = 2*(na-k)-1, i.e. k =
     * 2*na-1-nc.  By assumption, the expression above, with nc>=na,
     * guarantees that the first condition 0<=k<=na is satisfied.
     */
    unsigned int chop = 2*na-1-nc;
    polymat_mp_raw_kara(ab, b, xb,
            a, xa, na - chop, c, xc + chop, nc - chop, transpose, add);
    polymat_mp_raw_subdivide(ab, b, xb,
            a, xa + na - chop, chop, c, xc, na - 1, transpose, 1);
}/*}}}*/

void polymat_mp_raw(abdst_field ab,/*{{{*/
        polymat b, unsigned int xb,
        polymat a, unsigned int xa, unsigned int na,
        polymat c, unsigned int xc, unsigned int nc,
        int transpose, int add)
{
    polymat_mp_raw_subdivide(ab, b, xb, a, xa, na, c, xc, nc, transpose, add);
}/*}}}*/

void polymat_mp(abdst_field ab, polymat b, polymat a, polymat c)/*{{{*/
{
    unsigned int nb = MAX(a->size, c->size) - MIN(a->size, c->size) + 1;
    ASSERT_ALWAYS(a->n == c->m);
    if (polymat_check_pre_init(b)) {
        polymat_init(ab, b, a->m, c->n, nb);
    }
    ASSERT_ALWAYS(b->m == a->m);
    ASSERT_ALWAYS(b->n == c->n);
    ASSERT_ALWAYS(b->alloc >= nb);
    b->size = nb;

    polymat_mp_raw(ab, b, 0, a, 0, a->size, c, 0, c->size, 0, 0);
}/*}}}*/

void polymat_addmp(abdst_field ab, polymat b, polymat a, polymat c)/*{{{*/
{
    unsigned int nb = MAX(a->size, c->size) - MIN(a->size, c->size) + 1;
    ASSERT_ALWAYS(a->n == c->m);
    if (polymat_check_pre_init(b)) {
        polymat_init(ab, b, a->m, c->n, nb);
    }
    ASSERT_ALWAYS(b->m == a->m);
    ASSERT_ALWAYS(b->n == c->n);
    ASSERT_ALWAYS(b->alloc >= nb);
    b->size = nb;

    polymat_mp_raw(ab, b, 0, a, 0, a->size, c, 0, c->size, 0, 1);
}/*}}}*/

void polymat_set_matpoly(abdst_field ab MAYBE_UNUSED, polymat_ptr dst, matpoly_srcptr src)
{
    if (!polymat_check_pre_init(dst))
        polymat_clear(ab, dst);
    polymat_init(ab, dst, src->m, src->n, src->size);
    dst->size = src->size;

    for(unsigned int i = 0 ; i < src->m ; i++) {
        for(unsigned int j = 0 ; j < src->n ; j++) {
            for(unsigned int k = 0 ; k < src->size ; k++) {
                abset(ab,
                        polymat_coeff(ab, dst, i, j, k),
                        matpoly_coeff_const(ab, src, i, j, k));
            }
        }
    }
}
