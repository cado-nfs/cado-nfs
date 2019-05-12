#include "cado.h"
#include "mpfq_layer.h"
#include <stdlib.h>
#include <gmp.h>
#include "portability.h"
#include "macros.h"
#include "utils.h"
#include "lingen-matpoly.hpp"
#include "flint-fft/fft.h"

int matpoly::check_pre_init() const
{
    if (m && n && alloc)
        return 0;
    if (!m && !n && !alloc)
        return 1;
    abort();
    return 0;
}

/* with the exception of matpoly_realloc, all functions here are exactly
 * identical to those in lingen-polymat.c */
/* {{{ init/zero/clear interface for matpoly */
matpoly::matpoly(abdst_field ab, unsigned int m, unsigned int n, int len) : ab(ab), m(m), n(n), alloc(len) {
    /* As a special case, we allow a pre-init state with m==n==len==0 */
    ASSERT(!m == !n);
    ASSERT(!m == !len);
    if (!m) return;
    size = 0;
    abvec_init(ab, &(x), m*n*alloc);
    abvec_set_zero(ab, x, m*n*alloc);
}
matpoly::~matpoly() {
    abvec_clear(ab, &(x), m*n*alloc);
}
matpoly::matpoly(matpoly && a)
    : ab(a.ab), m(a.m), n(a.n), alloc(a.alloc)
{
    size=a.size;
    x=a.x;
    a.x=NULL;
    a.m=a.n=a.size=a.alloc=0;
    a.ab=NULL;
}
matpoly& matpoly::operator=(matpoly&& a)
{
    if (m)
        abvec_clear(ab, &(x), m*n*alloc);
    ab = a.ab;
    m = a.m;
    n = a.n;
    alloc = a.alloc;
    size = a.size;
    x=a.x;
    a.x=NULL;
    a.m=a.n=a.size=a.alloc=0;
    a.ab=NULL;
    return *this;
}

void matpoly::realloc(size_t newalloc) {
    check_pre_init();
    /* zero out the newly added data */
    if (newalloc > alloc) {
        /* allocate new space, then deflate */
        abvec_reinit(ab, &(x), m*n*alloc, m*n*newalloc);
        abvec rhead = abvec_subvec(ab, x, m*n*alloc);
        abvec whead = abvec_subvec(ab, x, m*n*newalloc);
        for(unsigned int i = m ; i-- ; ) {
            for(unsigned int j = n ; j-- ; ) {
                whead = abvec_subvec(ab, whead, -newalloc);
                rhead = abvec_subvec(ab, rhead, -alloc);
                abvec_set(ab, whead, rhead, alloc);
                abvec_set_zero(ab,
                        abvec_subvec(ab, whead, alloc),
                        newalloc - alloc);
            }
        }
    } else {
        /* inflate, then free space */
        ASSERT_ALWAYS(size <= newalloc);
        abvec rhead = x;
        abvec whead = x;
        for(unsigned int i = 0 ; i < m ; i++) {
            for(unsigned int j = 0 ; j < n ; j++) {
                abvec_set(ab, whead, rhead, newalloc);
                whead = abvec_subvec(ab, whead, newalloc);
                rhead = abvec_subvec(ab, rhead, alloc);
            }
        }
        abvec_reinit(ab, &(x), m*n*alloc, m*n*newalloc);
    }
    alloc = newalloc;
}
void matpoly::zero() {
    size = 0;
    abvec_set_zero(ab, x, m*n*alloc);
}

void matpoly::set_constant_ui(unsigned long e) {
    ASSERT_ALWAYS(m == n);
    size = 0;
    if (alloc == 0 && e)
        realloc(1);
    zero();
    if (!e) return;
    size = 1;
    for(unsigned int i = 0 ; i < m ; ++i)
        abset_ui(ab, coeff(i, i, 0), e);
}
void matpoly::set_constant(absrc_elt e) {
    ASSERT_ALWAYS(m == n);
    size = 0;
    if (alloc == 0 && abcmp_ui(ab, e, 0) != 0)
        realloc(1);
    zero();
    if (abcmp_ui(ab, e, 0) == 0) return;
    size = 1;
    for(unsigned int i = 0 ; i < m ; ++i)
        abset(ab, coeff(i, i, 0), e);
}
#if 0
void matpoly::swap(matpoly & a)
{
    matpoly x = std::move(*this);
    *this = std::move(a);
    a = std::move(x);
}
#endif
/* }}} */

void matpoly::fill_random(unsigned int nsize, gmp_randstate_t rstate)
{
    ASSERT_ALWAYS(nsize <= alloc);
    size = nsize;
    abvec_random(ab, x, m*n*size, rstate);
}

int matpoly::cmp(matpoly const& b) const
{
    ASSERT_ALWAYS(n == b.n);
    ASSERT_ALWAYS(m == b.m);
    if (size != b.size) return (size > b.size) - (b.size > size);
    return abvec_cmp(ab, x, b.x, m*n*size);
}

#if 0   /* unused ? */
/* polymat didn't have a proper add(), which is somewhat curious */
void matpoly::add(matpoly const & a, matpoly const & b)
{
    ASSERT_ALWAYS(size == a.size);
    ASSERT_ALWAYS(size == b.size);
    for(unsigned int i = 0 ; i < a.m ; i++) {
        for(unsigned int j = 0 ; j < a.n ; j++) {
            abvec_add(ab,
                      part(i, j, 0),
                    a.part(i, j, 0),
                    b.part(i, j, 0), size);
        }
    }
}
#endif

/* shift by a multiplication by x all coefficients of degree k and above
 * in column j. (that is, whole polynomial entries if k==0)
 */
void matpoly::multiply_column_by_x(unsigned int j, unsigned int nsize)/*{{{*/
{
    ASSERT_ALWAYS((nsize + 1) <= alloc);
    for(unsigned int i = 0 ; i < m ; i++) {
        memmove(part(i, j, 1), part(i, j, 0), nsize * abvec_elt_stride(ab, 1));
        abset_ui(ab, coeff(i, j, 0), 0);
    }
}/*}}}*/
void matpoly::divide_column_by_x(unsigned int j, unsigned int nsize)/*{{{*/
{
    ASSERT_ALWAYS(nsize <= alloc);
    for(unsigned int i = 0 ; i < m ; i++) {
        memmove(part(i, j, 0), part(i, j, 1), 
                (nsize-1) * abvec_elt_stride(ab, 1));
        abset_ui(ab, coeff(i, j, nsize-1), 0);
    }
}/*}}}*/

void matpoly::truncate(matpoly const & src, unsigned int nsize)/*{{{*/
{
    ASSERT_ALWAYS(nsize <= src.alloc);
    if (this == &src) {
        /* Never used by the code, so far. We're leaving garbage coeffs
         * on top, could this be a problem ? */
        size = nsize;
        return;
    }
    if (check_pre_init()) {
        /* Need to call ctor... */
        *this = matpoly(src.ab, src.m, src.n, nsize);
    }
    ASSERT_ALWAYS(m == src.m);
    ASSERT_ALWAYS(n == src.n);
    ASSERT_ALWAYS(nsize <= alloc);
    size = nsize;
    /* XXX Much more cumbersome here than for polymat, of course */
    for(unsigned int i = 0 ; i < src.m ; i++) {
        for(unsigned int j = 0 ; j < src.n ; j++) {
            abvec_set(ab, part(i, j, 0), src.part(i, j, 0), nsize);
        }
    }
}/*}}}*/

/* XXX compared to polymat, our diffferent stride has a consequence,
 * clearly ! */
void matpoly::extract_column( /*{{{*/
        unsigned int jdst, unsigned int kdst,
        matpoly const & src, unsigned int jsrc, unsigned int ksrc)
{
    ASSERT_ALWAYS(m == src.m);
    for(unsigned int i = 0 ; i < src.m ; i++)
        abset(ab, coeff(i, jdst, kdst), src.coeff(i, jsrc, ksrc));
}/*}}}*/

void matpoly::transpose_dumb(matpoly const & src) /*{{{*/
{
    if (this == &src) {
        matpoly tmp;
        tmp.transpose_dumb(src);
        *this = std::move(tmp);
        return;
    }
    if (!check_pre_init()) {
        *this = matpoly();
    }
    matpoly tmp(src.ab, src.n, src.m, src.size);
    tmp.size = src.size;
    for(unsigned int i = 0 ; i < src.m ; i++) {
        for(unsigned int j = 0 ; j < src.n ; j++) {
            for(unsigned int k = 0 ; k < src.size ; k++) {
                abset(ab, tmp.coeff(j, i, k), src.coeff(i, j, k));
            }
        }
    }
    *this = std::move(tmp);
}/*}}}*/

void matpoly::zero_column(unsigned int jdst, unsigned int kdst) /*{{{*/
{
    for(unsigned int i = 0 ; i < m ; i++)
        abset_zero(ab, coeff(i, jdst, kdst));
}/*}}}*/

void matpoly::extract_row_fragment(/*{{{*/
        unsigned int i1, unsigned int j1,
        matpoly const & src, unsigned int i0, unsigned int j0,
        unsigned int n)
{
    ASSERT_ALWAYS(src.size <= alloc);
    ASSERT_ALWAYS(size == src.size);
    for(unsigned int k = 0 ; k < n ; k++)
        abvec_set(ab, part(i1, j1 + k, 0), src.part(i0, j0 + k, 0), size);
}/*}}}*/

void matpoly::rshift(matpoly const & src, unsigned int k)/*{{{*/
{
    ASSERT_ALWAYS(k <= src.size);
    unsigned int newsize = src.size - k;
    if (check_pre_init()) {
        *this = matpoly(src.ab, src.m, src.n, newsize);
    }
    ASSERT_ALWAYS(m == src.m);
    ASSERT_ALWAYS(n == src.n);
    ASSERT_ALWAYS(newsize <= alloc);
    size = newsize;
    for(unsigned int i = 0 ; i < src.m ; i++) {
        for(unsigned int j = 0 ; j < src.n ; j++) {
            abvec_set(ab, part(i, j, 0), src.part(i, j, k), newsize);
        }
    }
    if (this == &src) realloc(newsize);
}/*}}}*/

void matpoly::addmul(matpoly const & a, matpoly const & b)/*{{{*/
{
    size_t csize = a.size + b.size; csize -= (csize > 0);
    if (this == &a || this == &b) {
        matpoly tc(a.ab, a.m, b.n, csize);
        tc.addmul(a, b);
        *this = std::move(tc);
        return;
    }
    ASSERT_ALWAYS(a.n == b.m);
    if (check_pre_init()) {
        *this = matpoly(a.ab, a.m, b.n, csize);
    }
    ASSERT_ALWAYS(m == a.m);
    ASSERT_ALWAYS(n == b.n);
    ASSERT_ALWAYS(alloc >= csize);

    for(unsigned int i = 0 ; i < m ; ++i) {
        for(unsigned int j = 0 ; j < n ; ++j) {
            abvec_set_zero(ab, part(i, j, size), csize - size);
        }
    }

    size = csize;
    abvec_ur tmp[2];
    abvec_ur_init(ab, &tmp[0], size);
    abvec_ur_init(ab, &tmp[1], size);

    for(unsigned int i = 0 ; i < a.m ; i++) {
        for(unsigned int j = 0 ; j < b.n ; j++) {
            abvec_ur_set_vec(ab, tmp[1], part(i, j, 0), size);
            for(unsigned int k = 0 ; k < a.n; k++) {
                abvec_conv_ur(ab, tmp[0],
                        a.part(i, k, 0), a.size,
                        b.part(k, j, 0), b.size);
                abvec_ur_add(ab, tmp[1], tmp[1], tmp[0], size);
            }
            abvec_reduce(ab, part(i, j, 0), tmp[1], size);
        }
    }
    abvec_ur_clear(ab, &tmp[0], size);
    abvec_ur_clear(ab, &tmp[1], size);
}/*}}}*/

void matpoly::mul(matpoly const & a, matpoly const & b)/*{{{*/
{
    size_t csize = a.size + b.size; csize -= (csize > 0);

    if (this == &a || this == &b) {
        matpoly tc(a.ab, a.m, b.n, csize);
        tc.mul(a, b);
        *this = std::move(tc);
        return;
    }
    ASSERT_ALWAYS(a.n == b.m);
    if (check_pre_init()) {
        *this = matpoly(a.ab, a.m, b.n, csize);
    }
    ASSERT_ALWAYS(m == a.m);
    ASSERT_ALWAYS(n == b.n);
    ASSERT_ALWAYS(alloc >= csize);
    size = csize;
    abvec_set_zero(ab, x, m * n * size);
    addmul(a, b);
}/*}}}*/

void matpoly::addmp(matpoly const & a, matpoly const & c)/*{{{*/
{
    size_t fullsize = a.size + c.size; fullsize -= (fullsize > 0);
    unsigned int nb = MAX(a.size, c.size) - MIN(a.size, c.size) + 1;
    ASSERT_ALWAYS(a.n == c.m);
    if (check_pre_init()) {
        *this = matpoly(a.ab, a.m, c.n, nb);
    }
    ASSERT_ALWAYS(m == a.m);
    ASSERT_ALWAYS(n == c.n);
    ASSERT_ALWAYS(alloc >= nb);
    size = nb;

    /* We are going to make it completely stupid for a beginning. */
    /* XXX XXX XXX FIXME !!! */
    abvec_ur tmp[2];
    abvec_ur_init(ab, &tmp[0], fullsize);
    abvec_ur_init(ab, &tmp[1], size);

    for(unsigned int i = 0 ; i < a.m ; i++) {
        for(unsigned int j = 0 ; j < c.n ; j++) {
            abvec_ur_set_vec(ab, tmp[1], part(i, j, 0), size);
            for(unsigned int k = 0 ; k < a.n ; k++) {
                abvec_conv_ur(ab, tmp[0],
                        a.part(i, k, 0), a.size,
                        c.part(k, j, 0), c.size);
                abvec_ur_add(ab, tmp[1], tmp[1],
                        abvec_ur_subvec(ab, tmp[0], MIN(a.size, c.size) - 1),
                        nb);
            }
            abvec_reduce(ab, part(i, j, 0), tmp[1], size);
        }
    }
    abvec_ur_clear(ab, &tmp[0], fullsize);
    abvec_ur_clear(ab, &tmp[1], size);
}/*}}}*/

void matpoly::mp(matpoly const & a, matpoly const & c)/*{{{*/
{
    unsigned int nb = MAX(a.size, c.size) - MIN(a.size, c.size) + 1;
    ASSERT_ALWAYS(a.n == c.m);
    if (check_pre_init()) {
        *this = matpoly(a.ab, a.m, c.n, nb);
    }
    ASSERT_ALWAYS(m == a.m);
    ASSERT_ALWAYS(n == c.n);
    ASSERT_ALWAYS(alloc >= nb);
    size = nb;
    abvec_set_zero(ab, x, m * n * size);
    addmp(a, c);
}/*}}}*/


void matpoly::set_polymat(polymat const & src)
{
    *this = matpoly(src.ab, src.m, src.n, src.size);
    size = src.size;

    for(unsigned int i = 0 ; i < src.m ; i++) {
        for(unsigned int j = 0 ; j < src.n ; j++) {
            for(unsigned int k = 0 ; k < src.size ; k++) {
                abset(ab, coeff(i, j, k), src.coeff(i, j, k));
            }
        }
    }
}

int matpoly::coeff_is_zero(unsigned int k) const
{
    for(unsigned int j = 0; j < m; j++)
        for(unsigned int i = 0 ; i < n ; i++)
            if (!abis_zero(ab, coeff(i, j, k)))
                return 0;
    return 1;
}
void matpoly::coeff_set_zero(unsigned int k)
{
    for(unsigned int j = 0; j < m; j++)
        for(unsigned int i = 0 ; i < n ; i++)
            abset_zero(ab, coeff(i, j, k));
}
