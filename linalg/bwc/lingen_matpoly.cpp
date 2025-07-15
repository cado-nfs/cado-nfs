#include "cado.h" // IWYU pragma: keep

#include <algorithm>
#include <utility>

#include <gmp.h>

#include "gmp_aux.h"
#include "lingen_matpoly_select.hpp"
#include "arith-hard.hpp"
#ifdef HAVE_OPENMP
#include "omp_proxy.h"
#endif
#include "macros.h"
#include "lingen_matpoly.hpp"
#include "lingen_polymat.hpp"

matpoly<false>::memory_pool_type matpoly<false>::memory;

/* with the exception of matpoly_realloc, all functions here are exactly
 * identical to those in lingen-polymat.c */
/* {{{ init/zero/clear interface for matpoly */
matpoly<false>::matpoly(arith_hard * ab, unsigned int m, unsigned int n, size_t len) : ab(ab), m(m), n(n), alloc(len) {
    /* As a special case, we allow a pre-init state with m==n==len==0 */
    /* Note that because we want to handle homogenous and non-homogenous
     * cases the same way, we support matrices of size 0*n, so that is
     * not technically a pre-init state */
    if (!m && !n) {
        ASSERT_ALWAYS(!len);
        return;
    }
    if (alloc) {
        if (data_alloc_size_in_bytes()) {
            x = (arith_hard::elt *) memory.alloc(data_alloc_size_in_bytes());
            ab->vec_set_zero(x, m*n*alloc);
        } else {
            x = nullptr;
        }
    }
}

matpoly<false>::~matpoly() {
    if (x) {
        memory.free(x, data_alloc_size_in_bytes());
        // ab->vec_clear(&(x), m*n*alloc);
    }
}
matpoly<false>::matpoly(matpoly && a) noexcept
    : ab(a.ab)
    , m(a.m)
    , n(a.n)
    , size(a.size)
    , alloc(a.alloc)
    , x(a.x)
{
    a.x = nullptr;
    a.m = a.n = a.size = a.alloc = 0;
    a.ab = nullptr;
}

matpoly<false>& matpoly<false>::operator=(matpoly&& a) noexcept
{
    if (x) {
        memory.free(x, data_alloc_size_in_bytes());
        // ab->vec_clear(&(x), m*n*alloc);
    }
    ab = a.ab;
    m = a.m;
    n = a.n;
    alloc = a.alloc;
    size = a.size;
    x = a.x;
    a.x = nullptr;
    a.m = a.n = a.size = a.alloc = 0;
    a.ab = nullptr;
    return *this;
}

matpoly<false>& matpoly<false>::set(matpoly const& a)
{
    if (x) {
        memory.free(x, data_alloc_size_in_bytes());
        // ab->vec_clear(&(x), m*n*alloc);
    }
    ab = a.ab;
    m = a.m;
    n = a.n;
    alloc = a.alloc;
    size = a.size;
    // ab->vec_init(&(x), m*n*alloc);
    x = (arith_hard::elt *) memory.alloc(data_alloc_size_in_bytes());
    ab->vec_set(x, a.x, m*n*alloc);
    return *this;
}

/* For structures in pre_init state, this is equivalent to an initial
 * allocation.
 *
 * If data is grown, old data is kept, and the 'size' field is unchanged.
 *
 * If data is shrunk below the previous value of 'size', then 'size' is
 * set to zero.
 *
 * The contents of the data area above 'size' on return is unspecified.
 */

void matpoly<false>::realloc(size_t newalloc) {
    ASSERT_ALWAYS(size <= alloc);
    size_t const oldmem = data_alloc_size_in_bytes();
    size_t const newmem = data_alloc_size_in_bytes(newalloc);

    /* zero out the newly added data */
    if (newalloc > alloc) {
        /* allocate new space, then inflate */
        // ab->vec_reinit(&(x), m*n*alloc, m*n*newalloc);
        x = (arith_hard::elt *) memory.realloc(x, oldmem, newmem);
        arith_hard::elt const * rhead = ab->vec_subvec(x, m*n*alloc);
        arith_hard::elt * whead = ab->vec_subvec(x, m*n*newalloc);
        if (size)
            for(unsigned int i = m ; i-- ; ) {
                for(unsigned int j = n ; j-- ; ) {
                    whead = ab->vec_subvec(whead, -newalloc);
                    rhead = ab->vec_subvec(rhead, -alloc);
                    ab->vec_set(whead, rhead, size);
                    ab->vec_set_zero(
                            ab->vec_subvec(whead, alloc),
                            newalloc - alloc);
                }
            }
    } else {
        if (size > newalloc)
            size = 0;
        /* deflate, then free space */
        ASSERT_ALWAYS(size <= newalloc);
        arith_hard::elt const * rhead = x;
        arith_hard::elt * whead = x;
        if (size)
            for(unsigned int i = 0 ; i < m ; i++) {
                for(unsigned int j = 0 ; j < n ; j++) {
                    ab->vec_set(whead, rhead, size);
                    whead = ab->vec_subvec(whead, newalloc);
                    rhead = ab->vec_subvec(rhead, alloc);
                }
            }
        // ab->vec_reinit(&(x), m*n*alloc, m*n*newalloc);
        x =(arith_hard::elt *)  memory.realloc(x, oldmem, newmem);
    }
    // if (!size) ab->vec_set_zero(x, m*n*alloc);
    alloc = newalloc;
}

size_t matpoly<false>::get_true_nonzero_size() const
{
    size_t lb = 0;
    size_t const ub = get_size();
    for(unsigned int ij = 0 ; ij < m*n && lb < ub ; ij++) {
        unsigned int const i = ij / n;
        unsigned int const j = ij % n;
        /* Find the last nonzero in the range [lb, ub[ */
        for(size_t k = ub ; k > lb ; k--) {
            if (!ab->is_zero(coeff(i, j, k-1))) {
                lb = k;
                break;
            }
        }
    }
    return lb;
}

void matpoly<false>::zero() {
    size = 0;
    ab->vec_set_zero(x, m*n*alloc);
}

void matpoly<false>::set_constant_ui(unsigned long e) {
    ASSERT_ALWAYS(m == n);
    size = 0;
    if (alloc == 0 && e)
        realloc(1);
    zero();
    if (!e) return;
    size = 1;
    for(unsigned int i = 0 ; i < m ; ++i)
        ab->set(coeff(i, i, 0), e);
}

void matpoly<false>::set_constant(arith_hard::elt const & e) {
    ASSERT_ALWAYS(m == n);
    size = 0;
    if (alloc == 0 && !ab->is_zero(e))
        realloc(1);
    zero();
    if (ab->is_zero(e)) return;
    size = 1;
    for(unsigned int i = 0 ; i < m ; ++i)
        ab->set(coeff(i, i, 0), e);
}
/* }}} */

void matpoly<false>::fill_random(size_t k0, size_t k1, cxx_gmp_randstate & rstate)
{
    ASSERT_ALWAYS(k1 <= alloc);
    if (k0 == 0 && k1 == alloc) {
        ab->vec_set_random(x, m*n*k1, rstate);
    } else if (k0 < k1) {
        for(unsigned int i = 0 ; i < m ; i++) {
            for(unsigned int j = 0 ; j < n ; j++) {
                ab->vec_set_random(part_head(i, j, k0), k1 - k0, rstate);
            }
        }
    }
}

int matpoly<false>::cmp(matpoly const& b) const
{
    ASSERT_ALWAYS(n == b.n);
    ASSERT_ALWAYS(m == b.m);
    if (size != b.size) {
        return (size > b.size) - (b.size > size);
    } else if (size == 0 && b.size == 0) {
        /* This is for the "intermediary pre-init" state. size is zero,
         * but a priori the dimension fields are ok. x might be NULL or
         * not, it depends on the previous life of the object */
        return 0;
    } else if (size == alloc && b.size == b.alloc) {
        return ab->vec_cmp(x, b.x, m*n*size);
    } else {
        for(unsigned int i = 0 ; i < m ; i++) {
            for(unsigned int j = 0 ; j < n ; j++) {
                int const r = ab->vec_cmp(part(i, j), b.part(i, j), size);
                if (r) return r;
            }
        }
        return 0;
    }
}

/* shift by a multiplication by x all coefficients of degree less than
 * colsize in column j. This results in a new column of length colsize+1.
 * Allocation must be sufficient for this length to fit.  What happens to
 * coefficients of degree larger than colsize in the result is unspecified.
 *
 * It is often relevant to "colsize++" right after this call, since the
 * coefficient of degree colsize is well-defined on output
 */
void matpoly<false>::multiply_column_by_x(unsigned int j, size_t colsize)/*{{{*/
{
    ASSERT_ALWAYS((colsize + 1) <= alloc);
    for(unsigned int i = 0 ; i < m ; i++) {
        ab->vec_set(part_head(i, j, 1), part(i, j), colsize);
        ab->set_zero(coeff(i, j, 0));
    }
}/*}}}*/

/* shift right (by a division by x) all coefficients of degree less than
 * colsize in column j. This results in a new column of length colsize-1,
 * with a zero coefficient shifted in at degree colsize-1.
 *
 * It is often relevant to "colsize--" right after this call.
 */
void matpoly<false>::divide_column_by_x(unsigned int j, size_t colsize)/*{{{*/
{
    if (!colsize) return;
    ASSERT_ALWAYS(colsize <= alloc);
    for(unsigned int i = 0 ; i < m ; i++) {
        ab->vec_set(part(i, j), part_head(i, j, 1), colsize-1);
        ab->set_zero(coeff(i, j, colsize-1));
    }
}/*}}}*/

void matpoly<false>::truncate(matpoly const & src, size_t nsize)/*{{{*/
{
    ASSERT_ALWAYS(nsize <= src.alloc);
    if (check_pre_init()) {
        /* Need to call ctor... */
        *this = matpoly(src.ab, src.m, src.n, nsize);
    }
    ASSERT_ALWAYS(m == src.m);
    ASSERT_ALWAYS(n == src.n);
    ASSERT_ALWAYS(nsize <= alloc);
    ASSERT_ALWAYS(nsize <= src.size);
    size = nsize;
    if (this == &src) return;
    /* XXX Much more cumbersome here than for polymat, of course */
    for(unsigned int i = 0 ; i < src.m ; i++) {
        for(unsigned int j = 0 ; j < src.n ; j++) {
            ab->vec_set(part(i, j), src.part(i, j), nsize);
        }
    }
}/*}}}*/

int matpoly<false>::tail_is_zero(size_t size0)/*{{{*/
{
    ASSERT_ALWAYS(size0 <= size);
    for(unsigned int i = 0 ; i < m ; i++) {
        for(unsigned int j = 0 ; j < n ; j++) {
            if (!ab->vec_is_zero(part_head(i, j, size0), size - size0))
                return 0;
        }
    }
    return 1;
}/*}}}*/

void matpoly<false>::zero_pad(size_t nsize)/*{{{*/
{
    ASSERT_ALWAYS(nsize >= size);
    if (check_pre_init() || nsize > alloc)
        realloc(nsize);
    for(unsigned int i = 0 ; i < m ; i++) {
        for(unsigned int j = 0 ; j < n ; j++) {
            ab->vec_set_zero(part_head(i, j, size), nsize - size);
        }
    }
    size = nsize;
}/*}}}*/


/* This takes coefficient ksrc of column jsrc, and copies it to
 * coefficient kdst of column jdst
 */
/* XXX compared to polymat, our diffferent stride has a consequence,
 * clearly ! */

void matpoly<false>::extract_column( /*{{{*/
        unsigned int jdst, size_t kdst,
        matpoly const & src, unsigned int jsrc, size_t ksrc)
{
    ASSERT_ALWAYS(m == src.m);
    for(unsigned int i = 0 ; i < src.m ; i++)
        ab->set(coeff(i, jdst, kdst), src.coeff(i, jsrc, ksrc));
}/*}}}*/

#if 0
void matpoly<false>::transpose_dumb(matpoly const & src) /*{{{*/
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
            for(size_t k = 0 ; k < src.size ; k++) {
                ab->set(tmp.coeff(j, i, k), src.coeff(i, j, k));
            }
        }
    }
    *this = std::move(tmp);
}/*}}}*/
#endif

void matpoly<false>::zero_column(unsigned int jdst, size_t kdst) /*{{{*/
{
    for(unsigned int i = 0 ; i < m ; i++)
        ab->set_zero(coeff(i, jdst, kdst));
}/*}}}*/

#if 0
void matpoly<false>::extract_row_fragment(/*{{{*/
        unsigned int i1, unsigned int j1,
        matpoly const & src, unsigned int i0, unsigned int j0,
        size_t n)
{
    ASSERT_ALWAYS(src.size <= alloc);
    ASSERT_ALWAYS(size == src.size);
    for(size_t k = 0 ; k < n ; k++)
        ab->vec_set(part(i1, j1 + k, 0), src.part(i0, j0 + k, 0), size);
}/*}}}*/
#endif

void matpoly<false>::view_t::zero() { /*{{{*/
    unsigned int const nrows = this->nrows();
    unsigned int const ncols = this->ncols();
#ifdef HAVE_OPENMP
#pragma omp parallel for collapse(2)
#endif
    for(unsigned int i = 0 ; i < nrows ; i++) {
        for(unsigned int j = 0 ; j < ncols ; j++) {
            M.ab->vec_set_zero(part(i,j), M.get_size());
        }
    }
}/*}}}*/

void matpoly<false>::rshift(matpoly const & src, size_t k)/*{{{*/
{
    ASSERT_ALWAYS(k <= src.size);
    size_t const newsize = src.size - k;
    if (check_pre_init()) {
        *this = matpoly(src.ab, src.m, src.n, newsize);
    }
    ASSERT_ALWAYS(m == src.m);
    ASSERT_ALWAYS(n == src.n);
    ASSERT_ALWAYS(newsize <= alloc);
    size = newsize;
    for(unsigned int i = 0 ; i < src.m ; i++) {
        for(unsigned int j = 0 ; j < src.n ; j++) {
            ab->vec_set(part(i, j), src.part_head(i, j, k), newsize);
        }
    }
}/*}}}*/

void matpoly<false>::rshift(size_t k)/*{{{*/
{
    ASSERT_ALWAYS(k <= size);
    size_t const newsize = size - k;
    ASSERT_ALWAYS(newsize <= alloc);
    size = newsize;
    for(unsigned int i = 0 ; i < m ; i++) {
        for(unsigned int j = 0 ; j < n ; j++) {
            /* can't use abvec_set because memcpy does not accept overlap */
            for(unsigned s = 0 ; s < newsize ; s++) {
                ab->set(coeff(i, j, s), coeff(i, j, s + k));
            }
        }
    }
}/*}}}*/

void matpoly<false>::add(matpoly const & a, matpoly const & b)/*{{{*/
{
    size_t const csize = std::max(a.size, b.size);
    ASSERT_ALWAYS(a.m == b.m);
    ASSERT_ALWAYS(a.n == b.n);
    if (check_pre_init()) {
        *this = matpoly(a.ab, a.m, a.n, csize);
    }
    if (alloc < csize)
        realloc(csize);
    ASSERT_ALWAYS(alloc >= csize);

    for(unsigned int i = 0 ; i < m ; ++i) {
        for(unsigned int j = 0 ; j < n ; ++j) {
            size_t const s0 = std::min(a.size, b.size);
            ab->vec_add_and_reduce(part(i, j), a.part(i, j), b.part(i, j), s0);
            if (a.size > s0)
                ab->vec_set(part_head(i, j, s0), a.part_head(i, j, s0), a.size - s0);
            if (b.size > s0)
                ab->vec_set(part_head(i, j, s0), b.part_head(i, j, s0), b.size - s0);
        }
    }
    size = csize;
}/*}}}*/

void matpoly<false>::sub(matpoly const & a, matpoly const & b)/*{{{*/
{
    size_t const csize = std::max(a.size, b.size);
    ASSERT_ALWAYS(a.m == b.m);
    ASSERT_ALWAYS(a.n == b.n);
    if (check_pre_init()) {
        *this = matpoly(a.ab, a.m, a.n, csize);
    }
    if (alloc < csize)
        realloc(csize);
    ASSERT_ALWAYS(alloc >= csize);

    for(unsigned int i = 0 ; i < m ; ++i) {
        for(unsigned int j = 0 ; j < n ; ++j) {
            size_t const s0 = std::min(a.size, b.size);
            ab->vec_sub_and_reduce(part(i, j), a.part(i, j), b.part(i, j), s0);
            if (a.size > s0)
                ab->vec_set(part_head(i, j, s0), a.part_head(i, j, s0), a.size - s0);
            if (b.size > s0)
                ab->vec_neg(part_head(i, j, s0), b.part_head(i, j, s0), b.size - s0);
        }
    }
    size = csize;
}/*}}}*/

void matpoly<false>::addmul(matpoly const & a, matpoly const & b)/*{{{*/
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

    if (a.size == 0 || b.size == 0)
        return;

    for(unsigned int i = 0 ; i < m ; ++i) {
        for(unsigned int j = 0 ; j < n ; ++j) {
            ab->vec_set_zero(part_head(i, j, size), csize - size);
        }
    }

    if (csize >= size)
        zero_pad(csize);

    matpoly<false>::addmul(*this, a, b);

    size = csize;
}/*}}}*/

void matpoly<false>::copy(matpoly<false>::view_t t, matpoly<false>::const_view_t a)/*{{{*/
{
    unsigned int const nrows = a.nrows();
    unsigned int const ncols = a.ncols();
    ASSERT_ALWAYS(t.nrows() == nrows);
    ASSERT_ALWAYS(t.ncols() == ncols);
#ifdef HAVE_OPENMP
    unsigned int T = std::min((unsigned int) omp_get_max_threads(), nrows*ncols);
#pragma omp parallel num_threads(T)
#endif
    {
#ifdef HAVE_OPENMP
#pragma omp for collapse(2)
#endif
        for(unsigned int i = 0 ; i < nrows ; i++) {
            for(unsigned int j = 0 ; j < ncols ; j++) {
                ptr tij = t.part(i, j);
                srcptr aij = a.part(i, j);
                a.M.ab->vec_set(tij, aij, a.M.get_size());
            }
        }
    }
}/*}}}*/

void matpoly<false>::addmul(matpoly<false>::view_t t, matpoly<false>::const_view_t t0, matpoly<false>::const_view_t t1)/*{{{*/
{
    unsigned int const nrows = t.nrows();
    unsigned int const ncols = t.ncols();
    unsigned int const nadd = t0.ncols();
    ASSERT_ALWAYS(&t.M != &t0.M);
    ASSERT_ALWAYS(&t.M != &t1.M);
    ASSERT_ALWAYS(&t0.M != &t1.M);
    ASSERT_ALWAYS(t0.nrows() == nrows);
    ASSERT_ALWAYS(t1.ncols() == ncols);
    ASSERT_ALWAYS(t0.ncols() == nadd);
    ASSERT_ALWAYS(t1.nrows() == nadd);
    // ASSERT(t0.check());
    // ASSERT(t1.check());
    // ASSERT(t.check());
    size_t csize = t0.M.size + t1.M.size; csize -= (csize > 0);
    arith_hard * ab = t0.M.ab;
    if (t0.M.size == 0 || t1.M.size == 0)
        return;
    /* Attention: we don't want to count on t.M.size being set yet. */
#ifdef HAVE_OPENMP
    unsigned int T = std::min((unsigned int) omp_get_max_threads(), nrows*ncols);
#pragma omp parallel num_threads(T)
#endif
    {
        auto * tmp0 = ab->alloc<arith_hard::elt_ur_for_addmul>(csize);
        auto * tmp1 = ab->alloc<arith_hard::elt_ur_for_addmul>(csize);

#ifdef HAVE_OPENMP
#pragma omp for collapse(2)
#endif
        for(unsigned int i = 0 ; i < nrows ; i++) {
            for(unsigned int j = 0 ; j < ncols ; j++) {
                ab->vec_set(tmp1, t.part(i, j), csize);
                for(size_t k = 0 ; k < nadd ; k++) {
                    ab->vec_conv_ur(tmp0,
                            t0.part(i, k), t0.M.size,
                            t1.part(k, j), t1.M.size);
                    ab->vec_add(tmp1, tmp0, csize);
                }
                ab->vec_reduce(t.part(i, j), tmp1, csize);
            }
        }
        ab->free(tmp0);
        ab->free(tmp1);
    }
}/*}}}*/

void matpoly<false>::addmp(matpoly<false>::view_t t, matpoly<false>::const_view_t t0, matpoly<false>::const_view_t t1)/*{{{*/
{
    unsigned int const nrows = t.nrows();
    unsigned int const ncols = t.ncols();
    unsigned int const nadd = t0.ncols();
    ASSERT_ALWAYS(&t.M != &t0.M);
    ASSERT_ALWAYS(&t.M != &t1.M);
    ASSERT_ALWAYS(&t0.M != &t1.M);
    ASSERT_ALWAYS(t0.nrows() == nrows);
    ASSERT_ALWAYS(t1.ncols() == ncols);
    ASSERT_ALWAYS(t0.ncols() == nadd);
    ASSERT_ALWAYS(t1.nrows() == nadd);
    // ASSERT(t0.check());
    // ASSERT(t1.check());
    // ASSERT(t.check());
    size_t fullsize = t0.M.size + t1.M.size; fullsize -= (fullsize > 0);
    size_t const nb = std::max(t0.M.size, t1.M.size) - std::min(t0.M.size, t1.M.size) + 1;
    arith_hard * ab = t0.M.ab;
    if (t0.M.size == 0 || t1.M.size == 0)
        return;
    /* Attention: we don't want to count on t.M.size being set yet. */

    /* We are going to make it completely stupid for a beginning. */
    /* Note that the caching code that is used for larger sizes does a
     * _real_ middle product.
     */
#ifdef HAVE_OPENMP
    unsigned int T = std::min((unsigned int) omp_get_max_threads(), nrows*ncols);
#pragma omp parallel num_threads(T)
#endif
    {
        auto * tmp0 = ab->alloc<arith_hard::elt_ur_for_addmul>(fullsize);
        auto * tmp1 = ab->alloc<arith_hard::elt_ur_for_addmul>(nb);

#ifdef HAVE_OPENMP
#pragma omp for collapse(2)
#endif
        for(unsigned int i = 0 ; i < nrows ; i++) {
            for(unsigned int j = 0 ; j < ncols ; j++) {
            ab->vec_set(tmp1, t.part(i, j), nb);
                for(unsigned int k = 0 ; k < nadd ; k++) {
                    ab->vec_conv_ur(tmp0,
                            t0.part(i, k), t0.M.size,
                            t1.part(k, j), t1.M.size);
                    ab->vec_add(tmp1,
                            ab->vec_subvec(tmp0, MIN(t0.M.size, t1.M.size) - 1), nb);
                }
                    ab->vec_reduce(t.part(i, j), tmp1, nb);
            }
        }
        ab->free(tmp0);
        ab->free(tmp1);
    }
}/*}}}*/

matpoly<false> matpoly<false>::mul(matpoly const & a, matpoly const & b)/*{{{*/
{
    size_t csize = a.size + b.size; csize -= (csize > 0);

    matpoly tc(a.ab, a.m, b.n, csize);

    ASSERT_ALWAYS(a.n == b.m);
    tc.set_size(csize);
    tc.zero();
    tc.addmul(a, b);
    return tc;
}/*}}}*/

void matpoly<false>::addmp(matpoly const & a, matpoly const & c)/*{{{*/
{
    // size_t fullsize = a.size + c.size; fullsize -= (fullsize > 0);
    size_t const nb = std::max(a.size, c.size) - std::min(a.size, c.size) + 1;
    ASSERT_ALWAYS(a.n == c.m);
    if (check_pre_init()) {
        *this = matpoly(a.ab, a.m, c.n, nb);
    }
    ASSERT_ALWAYS(m == a.m);
    ASSERT_ALWAYS(n == c.n);
    ASSERT_ALWAYS(alloc >= nb);

    if (a.size == 0 || c.size == 0)
        return;

    if (nb >= size)
        zero_pad(nb);

    matpoly<false>::addmp(*this, a, c);
}/*}}}*/

matpoly<false> matpoly<false>::mp(matpoly const & a, matpoly const & c)/*{{{*/
{
    size_t const nb = std::max(a.size, c.size) - std::min(a.size, c.size) + 1;
    ASSERT_ALWAYS(a.n == c.m);
    matpoly b(a.ab, a.m, c.n, nb);
    b.zero();
    b.set_size(nb);
    b.addmp(a, c);
    return b;
}/*}}}*/


void matpoly<false>::set_polymat(polymat const & src)
{
    *this = matpoly(src.ab, src.m, src.n, src.get_size());
    set_size(src.get_size());

    for(unsigned int i = 0 ; i < src.m ; i++) {
        for(unsigned int j = 0 ; j < src.n ; j++) {
            for(size_t k = 0 ; k < size ; k++) {
                ab->set(coeff(i, j, k), src.coeff(i, j, k));
            }
        }
    }
}

int matpoly<false>::coeff_is_zero(size_t k) const
{
    for(unsigned int j = 0; j < ncols(); j++)
        for(unsigned int i = 0 ; i < nrows() ; i++)
            if (!ab->is_zero(coeff(i, j, k)))
                return 0;
    return 1;
}

void matpoly<false>::coeff_set_zero(size_t k)
{
    for(unsigned int j = 0; j < ncols(); j++)
        for(unsigned int i = 0 ; i < nrows() ; i++)
            ab->set_zero(coeff(i, j, k));
}

matpoly<false> matpoly<false>::truncate_and_rshift(size_t truncated_size, size_t shiftcount)
{
    matpoly other(ab, m, n, size - shiftcount);
    other.rshift(*this, shiftcount);
    truncate(*this, truncated_size);
    shrink_to_fit();
    std::swap(*this, other);
    return other;
}
