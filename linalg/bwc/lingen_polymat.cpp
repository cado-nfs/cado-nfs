#include "cado.h" // IWYU pragma: keep

#include <climits>
#include <cstdlib>
#include <cstring>

#include <algorithm>

#include <gmp.h>

#include "arith-hard.hpp"
#include "gmp_aux.h"
#include "lingen_matpoly_select.hpp"
#include "lingen_polymat.hpp"
#include "macros.h"

#define POLYMAT_MUL_KARA_CUTOFF_DEFAULT { 10, 10 }
#define POLYMAT_MP_KARA_CUTOFF_DEFAULT  { 10, 10 }

/*{{{ cutoffs and such */
static struct polymat_cutoff_info polymat_mul_kara_cutoff =
    POLYMAT_MUL_KARA_CUTOFF_DEFAULT;

/* For the middle product, this is similar, but applies to the small
 * operand length and to the result length
 * (A balanced middle product is n times 2*n - 1, and gives n result
 * coefficients. Here the important data is n).
 * */
/* FIXME: polymat_mul_kara_threshold and polymat_mp_kara_threshold should
 * not differ, really. The only thing is that dimensions for mp and mul
 * are not the same, hence the difference...
 */
static struct polymat_cutoff_info polymat_mp_kara_cutoff =
    POLYMAT_MP_KARA_CUTOFF_DEFAULT;

/* This sets the cutoffs. The table, if present in new_cutoff, is copied,
 * and hence must be freed by the caller who allocated it in the first
 * place.
 * The old_cutoff (if not NULL) value should be initialized already, and
 * receives initialized data.
 */

void polymat_cutoff_add_step(struct polymat_cutoff_info * c, size_t size,
                             int alg)
{
    c->table.emplace_back(size, alg);
}

static void
polymat_set_generic_cutoff(struct polymat_cutoff_info * slot,
                           const struct polymat_cutoff_info * new_cutoff,
                           struct polymat_cutoff_info * old_cutoff)
{
    if (old_cutoff)
        *old_cutoff = *slot;
    *slot = *new_cutoff;
}

void polymat_set_mul_kara_cutoff(const struct polymat_cutoff_info * new_cutoff,
                                 struct polymat_cutoff_info * old_cutoff)
{
    polymat_set_generic_cutoff(&polymat_mul_kara_cutoff, new_cutoff,
                               old_cutoff);
}

void polymat_set_mp_kara_cutoff(const struct polymat_cutoff_info * new_cutoff,
                                struct polymat_cutoff_info * old_cutoff)
{
    polymat_set_generic_cutoff(&polymat_mp_kara_cutoff, new_cutoff, old_cutoff);
}

static int
polymat_cutoff_get_alg_b(const struct polymat_cutoff_info * cutoff,
                         size_t s)
{
    if (s >= cutoff->cut)
        return 1;
    if (cutoff->table.empty())
        return 0;
    size_t a = 0;
    size_t b = cutoff->table.size();
    for (; b - a > 1;) {
        size_t const c = (a + b) / 2;
        if (s >= cutoff->table[c].first) {
            a = c;
        } else {
            b = c;
        }
    }
    return cutoff->table[a].second;
}

static int
polymat_cutoff_get_subdivide_ub(const struct polymat_cutoff_info * cutoff,
                                size_t s0, size_t s1)
{
    return !cutoff->table.empty() && std::min(s0, s1) >= cutoff->subdivide;
}

/*}}}*/

/* {{{ unreduced interface, for the implementation only */

template <typename X> struct polymat_ur {
    arith_hard * ab = nullptr;
    unsigned int m = 0;
    unsigned int n = 0;

  private:
    size_t size = 0;
    size_t alloc = 0;

  public:
    polymat_ur(polymat_ur const &) = delete;
    polymat_ur(polymat_ur &&) = delete;
    polymat_ur& operator=(polymat_ur const &) = delete;
    polymat_ur& operator=(polymat_ur &&) = delete;
    size_t capacity() const { return alloc; }
    size_t get_size() const { return size; }
    void set_size(size_t s) { size = s; }

    using ptr = X *;
    using srcptr = X const *;

    ptr x = nullptr;

    polymat_ur() = default;
    ~polymat_ur();
    polymat_ur(arith_hard * ab, unsigned int m, unsigned int n, size_t len);
    void zero();
    /* {{{ access interface for polymat_ur */
    ptr part(unsigned int i, unsigned int j, size_t k)
    {
        /* Assume row-major in all circumstances. Old code used to support
         * various orderings, here we don't */
        return ab->vec_subvec(x, (k * m + i) * n + j);
    }
    X & coeff(unsigned int i, unsigned int j, size_t k)
    {
        return ab->vec_item(part(i, j, k), 0);
    }
#if 0
    inline srcptr part(unsigned int i, unsigned int j, size_t k) const {
        /* Assume row-major in all circumstances. Old code used to support
         * various orderings, here we don't */
        return ab->vec_subvec(x, (k*m+i)*n+j);
    }
    inline X & coeff(unsigned int i, unsigned int j, size_t k) const {
        return ab->vec_item(part(i,j,k), 0);
    }
#endif
    /* }}} */
    void reducemat(polymat & c, size_t kc, size_t ka);
    void addmulmat(size_t kc, polymat const & a, size_t ka,
                   polymat const & b, size_t kb);
};
/* }}} */

int polymat::check_pre_init() const
{
    if (m && n && alloc)
        return 0;
    if (!m && !n && !alloc)
        return 1;
    abort();
    return 0;
}

/* {{{ init/zero/clear interface for polymat */
polymat::polymat(arith_hard * ab, unsigned int m, unsigned int n, size_t len)
    : ab(ab)
    , m(m)
    , n(n)
    , alloc(len)
{
    if (!m)
        return;
    x = ab->alloc(m * n * alloc);
    ab->vec_set_zero(x, m * n * alloc);
}

void polymat::realloc(size_t newalloc)
{
    ASSERT_ALWAYS(!check_pre_init());
    x = ab->realloc(x, m * n * alloc, m * n * newalloc);
    /* zero out the newly added data */
    if (newalloc > alloc) {
        ab->vec_set_zero(x + m * n * alloc, m * n * (newalloc - alloc));
    } else {
        ASSERT_ALWAYS(size <= newalloc);
    }
    alloc = newalloc;
}
void polymat::zero()
{
    size = 0;
    ab->vec_set_zero(x, m * n * alloc);
}
polymat::~polymat()
{
    if (ab)
        ab->free(x);
}
polymat::polymat(polymat && a) noexcept
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
polymat & polymat::operator=(polymat && a) noexcept
{
    if (m)
        ab->free(x);
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

#if 0
void polymat::swap(polymat & b)
{
    polymat x;
    memcpy(&x, this, sizeof(polymat));
    memcpy(this, &b, sizeof(polymat));
    memcpy(&b,   &x, sizeof(polymat));
}
#endif

void polymat::fill_random(unsigned int nsize, cxx_gmp_randstate & rstate)
{
    ASSERT_ALWAYS(nsize <= alloc);
    size = nsize;
    ab->vec_set_random(x, m * n * size, rstate);
}

int polymat::cmp(polymat const & b) const
{
    ASSERT_ALWAYS(n == b.n);
    ASSERT_ALWAYS(m == b.m);
    if (size != b.size)
        return (size > b.size) - (b.size > size);
    return ab->vec_cmp(x, b.x, m * n * size);
}

/* }}} */

/* {{{ init/zero/clear interface for polymat_ur */
template <typename X>
polymat_ur<X>::polymat_ur(arith_hard * ab, unsigned int m, unsigned int n,
                          size_t len)
    : ab(ab)
    , m(m)
    , n(n)
    , alloc(len)
    , x(nullptr)
{
    if (!m)
        return;
    x = ab->alloc<X>(m * n * alloc);
    ab->vec_set_zero(x, m * n * alloc);
}
template <typename X> void polymat_ur<X>::zero()
{
    size = 0;
    ab->vec_set_zero(x, m * n * alloc);
}
template <typename X> polymat_ur<X>::~polymat_ur()
{
    ab->free(x);
}
/* }}} */

/* It's used from lingen-bigpolymat.c as well */
static void bwmat_copy_coeffs(arith_hard * ab MAYBE_UNUSED,
                              arith_hard::elt * x0, unsigned int stride0,
                              arith_hard::elt const * x1, unsigned int stride1,
                              size_t n)
{
    for (size_t i = 0; i < n; i++) {
        ab->set(ab->vec_item(x0, i * stride0), ab->vec_item(x1, i * stride1));
    }
}
static void bwmat_zero_coeffs(arith_hard * ab MAYBE_UNUSED,
                              arith_hard::elt * x0, unsigned int stride0, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        ab->set_zero(ab->vec_item(x0, i * stride0));
    }
}
static void bwmat_move_coeffs(arith_hard * ab MAYBE_UNUSED,
                              arith_hard::elt * x0, unsigned int stride0,
                              arith_hard::elt const * x1, unsigned int stride1,
                              size_t n)
{
    ASSERT_ALWAYS(stride0 == stride1); /* Otherwise there's probably no point */
    if (x0 < x1) {
        for (size_t i = 0; i < n; i++) {
            ab->set(ab->vec_item(x0, i * stride0),
                    ab->vec_item(x1, i * stride1));
        }
    } else {
        for (size_t i = n; i--;) {
            ab->set(ab->vec_item(x0, i * stride0),
                    ab->vec_item(x1, i * stride1));
        }
    }
}

/* {{{ Elementary operations on the matrices, which are the coefficients
 * of our polynomials */
void polymat::addmat(size_t kc, polymat const & a, size_t ka,
                     polymat const & b, size_t kb)
{
    for (unsigned int i = 0; i < a.m; i++) {
        for (unsigned int j = 0; j < a.n; j++) {
            ab->add_and_reduce(coeff(i, j, kc), a.coeff(i, j, ka),
                               b.coeff(i, j, kb));
        }
    }
}

void polymat::submat(size_t kc, polymat const & a, size_t ka,
                     polymat const & b, size_t kb)
{
    for (unsigned int i = 0; i < a.m; i++) {
        for (unsigned int j = 0; j < a.n; j++) {
            ab->sub_and_reduce(coeff(i, j, kc), a.coeff(i, j, ka),
                               b.coeff(i, j, kb));
        }
    }
}

template <typename X>
void polymat_ur<X>::addmulmat(size_t kc, polymat const & a,
                              size_t ka, polymat const & b,
                              size_t kb)
{
    ASSERT_ALWAYS(a.n == b.m);
    for (unsigned int i = 0; i < a.m; i++) {
        for (unsigned int j = 0; j < b.n; j++) {
            for (unsigned int k = 0; k < a.n; k++) {
                polymat::elt const & ea = a.coeff(i, k, ka);
                polymat::elt const & eb = b.coeff(k, j, kb);
                ab->addmul_ur(coeff(i, j, kc), ea, eb);
            }
        }
    }
}

template <typename X>
void polymat_ur<X>::reducemat(polymat & c, size_t kc, size_t ka)
{
    for (unsigned int i = 0; i < m; i++) {
        for (unsigned int j = 0; j < n; j++) {
            ab->reduce(c.coeff(i, j, kc), coeff(i, j, ka));
        }
    }
}

void polymat::mulmat(size_t kc, polymat const & a, size_t ka,
                     polymat const & b, size_t kb)
{
    polymat_ur<arith_hard::elt_ur_for_addmul> cx(ab, m, n, 1);
    cx.addmulmat(0, a, ka, b, kb);
    cx.reducemat(*this, kc, 0);
}

void polymat::addmulmat(size_t kc, polymat const & a, size_t ka,
                        polymat const & b, size_t kb)
{
    polymat cc(ab, m, n, 1);
    cc.mulmat(0, a, ka, b, kb);
    addmat(kc, *this, kc, cc, 0);
}

/* }}} */

void polymat::multiply_column_by_x(unsigned int j, size_t nsize) /*{{{*/
{
    ASSERT_ALWAYS((nsize + 1) <= alloc);
    bwmat_move_coeffs(ab, part(0, j, 1), n, part(0, j, 0), n, m * nsize);
    for (unsigned int i = 0; i < m; i++)
        ab->set_zero(coeff(i, j, 0));
} /*}}}*/

void polymat::truncate(polymat const & src, size_t nsize) /*{{{*/
{
    ASSERT_ALWAYS(nsize <= src.alloc);
    if (this == &src) {
        /* Never used by the code, so far. We're leaving garbage coeffs
         * on top, could this be a problem ? */
        abort(); /* so that programmer becomes aware of the gotcha */
        size = nsize;
        return;
    }
    if (check_pre_init())
        *this = polymat(ab, src.m, src.n, nsize);
    ASSERT_ALWAYS(m == src.m);
    ASSERT_ALWAYS(n == src.n);
    ASSERT_ALWAYS(nsize <= alloc);
    size = nsize;
    bwmat_copy_coeffs(ab, part(0, 0, 0), 1, src.part(0, 0, 0), 1, m * n * size);
} /*}}}*/

void polymat::extract_column(/*{{{*/
                             unsigned int jdst, size_t kdst,
                             polymat const & src, unsigned int jsrc,
                             size_t ksrc)
{
    ASSERT_ALWAYS(m == src.m);
    bwmat_copy_coeffs(ab, part(0, jdst, kdst), n, src.part(0, jsrc, ksrc),
                      src.n, src.m);
} /*}}}*/

void polymat::extract_row_fragment(/*{{{*/
                                   unsigned int i1, unsigned int j1,
                                   polymat const & src, unsigned int i0,
                                   unsigned int j0, size_t nn)
{
    ASSERT_ALWAYS(src.size <= alloc);
    for (size_t k = 0; k < src.size; k++) {
        bwmat_copy_coeffs(ab, part(i1, j1, k), 1, src.part(i0, j0, k), 1, nn);
    }
} /*}}}*/

void polymat::rshift(polymat const & src, size_t k) /*{{{*/
{
    ASSERT_ALWAYS(k <= src.size);
    size_t const newsize = src.size - k;
    ASSERT_ALWAYS(m == src.m);
    ASSERT_ALWAYS(n == src.n);
    if (this != &src) {
        if (check_pre_init()) {
            *this = polymat(ab, src.m, src.n, newsize);
        }
        ASSERT_ALWAYS(newsize <= alloc);
        size = newsize;
        bwmat_copy_coeffs(ab, part(0, 0, 0), 1, src.part(0, 0, k), 1,
                          m * n * size);
    } else {
        size = newsize;
        bwmat_move_coeffs(ab, part(0, 0, 0), 1, src.part(0, 0, k), 1,
                          m * n * size);
    }
    realloc(size);
} /*}}}*/

static void polymat_mul_raw_basecase(/*{{{*/
                                     polymat & c, size_t xc,
                                     polymat const & a, size_t xa, size_t na,
                                     polymat const & b, size_t xb, size_t nb,
                                     int transpose, int add)
{
    size_t const nc = na + nb - 1;
    ASSERT_ALWAYS(c.capacity() >= xc + nc);

    polymat_ur<arith_hard::elt_ur_for_addmul> tmat_ur(c.ab, c.m, c.n, 1);

    for (size_t k = 0; k < nc; k++) {
        size_t const i0 = k >= nb ? k + 1 - nb : 0;
        size_t const i1 = k + 1 < na ? k + 1 : na;
        tmat_ur.zero();
        for (size_t i = i0; i < i1; i++) {
            size_t const j = k - i;
            if (!transpose) {
                tmat_ur.addmulmat(0, a, xa + i, b, xb + j);
            } else {
                tmat_ur.addmulmat(0, b, xb + j, a, xa + i);
            }
        }
        if (!add) {
            tmat_ur.reducemat(c, xc + k, 0);
        } else {
            polymat tmat(c.ab, c.m, c.n, 1);
            tmat.set_size(1);
            tmat_ur.reducemat(tmat, 0, 0);
            c.addmat(xc + k, c, xc + k, tmat, 0);
        }
    }
} /*}}}*/

static void polymat_mul_raw_kara(/*{{{*/
                                 polymat & c, size_t xc,
                                 polymat const & a, size_t xa, size_t na,
                                 polymat const & b, size_t xb, size_t nb,
                                 int transpose, int add)
{
    /* This code works __only__ for kara-sized products */
    ASSERT_ALWAYS(nb == na);

    if (polymat_cutoff_get_alg_b(&polymat_mul_kara_cutoff, std::min(na, nb)) == 0) {
        polymat_mul_raw_basecase(c, xc, a, xa, na, b, xb, nb, transpose, add);
        return;
    }

    size_t const nc = na + nb - 1;
    ASSERT_ALWAYS(c.capacity() >= xc + nc);
    if (!add) {
        bwmat_zero_coeffs(c.ab, c.part(0, 0, xc), 1, c.m * c.n * nc);
    }
    size_t const m0 = na / 2;
    size_t const m1 = na - m0;

    {
        polymat s(c.ab, a.m, a.n, m1);
        polymat t(c.ab, b.m, b.n, m1);
        s.set_size(m1);
        t.set_size(m1);
        for (size_t k = 0; k < m0; k++) {
            s.addmat(k, a, xa + k, a, xa + m0 + k);
            t.addmat(k, b, xb + k, b, xb + m0 + k);
        }
        if (m1 > m0) {
            s.addmat(m0, s, m0, a, xa + m0 + m0);
            t.addmat(m0, t, m0, b, xb + m0 + m0);
        }
        polymat_mul_raw_kara(c, xc + m0, s, 0, m1, t, 0, m1, transpose, 1);
    }

    {
        polymat q0(c.ab, c.m, c.n, 2 * m0 - 1);
        polymat q2(c.ab, c.m, c.n, 2 * m1 - 1);
        q0.set_size(2 * m0 - 1);
        q2.set_size(2 * m1 - 1);
        polymat_mul_raw_kara(q0, 0, a, xa, m0, b, xb, m0, transpose, 0);
        polymat_mul_raw_kara(q2, 0, a, xa + m0, m1, b, xb + m0, m1, transpose,
                             0);
        for (size_t k = 0; k < 2 * m0 - 1; k++) {
            c.addmat(xc + k, c, xc + k, q0, k);
            c.submat(xc + m0 + k, c, xc + m0 + k, q0, k);
        }
        for (size_t k = 0; k < 2 * m1 - 1; k++) {
            c.addmat(xc + 2 * m0 + k, c, xc + 2 * m0 + k, q2, k);
            c.submat(xc + m0 + k, c, xc + m0 + k, q2, k);
        }
    }
} /* }}} */

static void polymat_mul_raw_subdivide(/*{{{*/
                                      polymat & c, size_t xc,
                                      polymat const & a, size_t xa, size_t na,
                                      polymat const & b, size_t xb, size_t nb,
                                      int transpose, int add)
{
    if (!na || !nb) {
        return;
    }
    /* we consider using karatsuba only when the smallest of the two
     * operands is large enough.
     */
    if (polymat_cutoff_get_subdivide_ub(&polymat_mul_kara_cutoff, na, nb) ==
        0) {
        polymat_mul_raw_basecase(c, xc, a, xa, na, b, xb, nb, transpose, add);
        return;
    }
    size_t const nc = na + nb - 1;
    if (!add) {
        bwmat_zero_coeffs(c.ab, c.part(0, 0, xc), 1, c.m * c.n * nc);
    }

    for (; nb >= na;) {
        polymat_mul_raw_kara(c, xc, a, xa, na, b, xb, na, transpose, 1);
        xb += na;
        xc += na;
        nb -= na;
    }
    /* see below comment in #if0'ed block. Seems to me that the code
     * theere is bogus */
    ASSERT_ALWAYS(nb < na);
    polymat_mul_raw_subdivide(c, xc, a, xa, na, b, xb, nb, transpose, 1);
#if 0
    if (!nb) {
        /* This condition also catches the Karatsuba case n, n */
        return;
    }
    /* FIXME: nb might be very small here. Why do we do
     * polymat::mul_raw_kara ? This defeats the semantics of the subdivide
     * check, apparently.
     */
    ASSERT_ALWAYS(nb < na);
    /* Fixup needed, now. Treat the largest possible subblock. We want
     * the intervals [0,na-k[ and [0, nb[ to match a kara scheme, i.e.
     * have k=nb-na.
     */
    size_t chop = na - nb;
    polymat::mul_raw_kara(ab, c, xc,
            a, xa, nb, b, xb, nb, transpose, 1);
    polymat::mul_raw_subdivide(ab, c, xc + na - chop,
            a, xa + nb, chop, b, xb, nb, transpose, 1);
#endif
} /*}}}*/

static void polymat_mul_raw(/*{{{*/
                            polymat & c, size_t xc,
                            polymat const & a, size_t xa, size_t na,
                            polymat const & b, size_t xb, size_t nb,
                            int transpose, int add)
{
    polymat_mul_raw_subdivide(c, xc, a, xa, na, b, xb, nb, transpose, add);
} /*}}}*/

void polymat::mul(polymat const & a, polymat const & b) /*{{{*/
{
    size_t csize = a.size + b.size;
    csize -= (csize > 0);
    ASSERT_ALWAYS(a.n == b.m);
    if (check_pre_init()) {
        *this = polymat(a.ab, a.m, b.n, csize);
    }
    ASSERT_ALWAYS(m == a.m);
    ASSERT_ALWAYS(n == b.n);
    ASSERT_ALWAYS(alloc >= csize);
    size = csize;
    polymat_mul_raw(*this, 0, a, 0, a.size, b, 0, b.size, 0, 0);
} /*}}}*/

void polymat::addmul(polymat const & a, polymat const & b) /*{{{*/
{
    size_t csize = a.size + b.size;
    csize -= (csize > 0);
    ASSERT_ALWAYS(a.n == b.m);
    if (check_pre_init()) {
        *this = polymat(a.ab, a.m, b.n, csize);
    }
    ASSERT_ALWAYS(m == a.m);
    ASSERT_ALWAYS(n == b.n);
    ASSERT_ALWAYS(alloc >= csize);
    size = csize;
    polymat_mul_raw(*this, 0, a, 0, a.size, b, 0, b.size, 0, 1);
} /*}}}*/

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
 *      na==nc==n+1 -. nb==1: [n]
 *      na=n+1, nc=2n -. nb==n: [n,2n[
 */
static void polymat_mp_raw_basecase(/*{{{*/
                                    polymat & b, size_t xb,
                                    polymat const & a, size_t xa, size_t na,
                                    polymat const & c, size_t xc, size_t nc,
                                    int transpose, int add)
{
    if (nc < na) {
        polymat_mp_raw_basecase(b, xb, c, xc, nc, a, xa, na, !transpose, add);
        return;
    }
    size_t const nb = nc - na + 1;
    ASSERT_ALWAYS(b.capacity() >= xb + nb);

    polymat_ur<arith_hard::elt_ur_for_addmul> tmat_ur(b.ab, b.m, b.n, 1);
    tmat_ur.set_size(1);

    for (size_t j = 0; j < nb; j++) {
        tmat_ur.zero();
        for (size_t i = 0; i < na; i++) {
            size_t const k = j + na - 1 - i;
            if (!transpose) {
                tmat_ur.addmulmat(0, a, xa + i, c, xc + k);
            } else {
                tmat_ur.addmulmat(0, c, xc + k, a, xa + i);
            }
        }
        if (!add) {
            tmat_ur.reducemat(b, xb + j, 0);
        } else {
            polymat tmat(b.ab, b.m, b.n, 1);
            tmat.set_size(1);
            tmat_ur.reducemat(tmat, 0, 0);
            b.addmat(xb + j, b, xb + j, tmat, 0);
        }
    }
} /*}}}*/

static void polymat_mp_raw_kara(/*{{{*/
                                polymat & b, size_t xb,
                                polymat const & a, size_t xa, size_t na,
                                polymat const & c, size_t xc, size_t nc,
                                int transpose, int add)
{
    /* This code works __only__ for kara-sized middle products */
    ASSERT_ALWAYS(nc == 2 * na - 1 || na == 2 * nc - 1);
    // printf("mp_kara(%d,%d)%b\n",na,nc,transpose?'T':' ');
    if (polymat_cutoff_get_alg_b(&polymat_mp_kara_cutoff, std::min(na, nc)) == 0) {
        polymat_mp_raw_basecase(b, xb, a, xa, na, c, xc, nc, transpose, add);
        return;
    }
    if (nc < na) {
        polymat_mp_raw_kara(b, xb, c, xc, nc, a, xa, na, !transpose, add);
        return;
    }
    size_t const m0 = na / 2;
    size_t const m1 = na - m0;
    /* Spans of different chunks, in (offset, length) format */
    size_t const span_a0[2] = {xa, m0};
    size_t const span_a1[2] = {xa + m0, m1};
    size_t const span_c0[2] = {xc, 2 * m1 - 1};
    size_t const span_c10[2] = {xc + m1, 2 * m0 - 1};
    size_t const span_c11[2] = {xc + m1, 2 * m1 - 1};
    size_t const span_c2[2] = {xc + 2 * m1, 2 * m0 - 1};
    /* Dammit, we need some temp storage, of course */
    /* The polymats s and t are used to build temporaries, respectively
     * from a and c, of respective maximum lengths m1 and 2*m1-1 */
    /* q0 is MP(a1, b0+b11), i.e. MP(m1, 2*m1-1), which
     * produces m1 coefficients. This goes to offset 0 (well, xb) in b. */
    {
        polymat t(c.ab, c.m, c.n, 2 * m1 - 1);
        t.set_size(2 * m1 - 1);
        for (size_t k = 0; k < 2 * m1 - 1; k++) {
            t.addmat(k, c, span_c0[0] + k, c, span_c11[0] + k);
        }
        polymat_mp_raw_kara(b, xb, a, span_a1[0], span_a1[1], t, 0, 2 * m1 - 1,
                            transpose, add);
    }
    {
        /* q1 is MP(a1-a0, b11), i.e. MP(m1, 2*m1-1) again */
        polymat q1(b.ab, b.m, b.n, m1);
        q1.set_size(m1);
        {
            polymat s(b.ab, a.m, a.n, m1);
            s.set_size(m1);
            bwmat_copy_coeffs(b.ab, s.part(0, 0, 0), 1,
                              a.part(0, 0, span_a1[0]), 1, m1 * a.m * a.n);
            for (size_t k = 0; k < m0; k++) {
                /* We've already copied the a1 coefficients above */
                s.submat(k + m1 - m0, s, k + m1 - m0, a, span_a0[0] + k);
            }
            polymat_mp_raw_kara(q1, 0, s, 0, m1, c, span_c11[0], span_c11[1],
                                transpose, 0);
        }
        /* q2 is MP(a0, b10+b2), i.e. MP(m0, 2*m0-1) */
        {
            /* This goes to offset m1 in b */
            polymat t(c.ab, c.m, c.n, 2 * m0 - 1);
            t.set_size(2 * m0 - 1);
            for (size_t k = 0; k < 2 * m0 - 1; k++) {
                t.addmat(k, c, span_c10[0] + k, c, span_c2[0] + k);
            }
            polymat_mp_raw_kara(b, xb + m1, a, span_a0[0], span_a0[1], t, 0,
                                2 * m0 - 1, transpose, add);
        }

        /* We now have to append the coefficients to form the result. Not all
         * coefficients in q1 are read for the high part of the middle
         * product. */
        {
            for (size_t k = 0; k < m1; k++) {
                b.submat(xb + k, b, xb + k, q1, k);
            }
            for (size_t k = 0; k < m0; k++) {
                b.addmat(xb + m1 + k, b, xb + m1 + k, q1, k);
            }
        }
    }
} /*}}}*/

static void polymat_mp_raw_subdivide(/*{{{*/
                                     polymat & b, size_t xb,
                                     polymat const & a, size_t xa, size_t na,
                                     polymat const & c, size_t xc, size_t nc,
                                     int transpose, int add)
{
    // printf("mp_subdivide(%d,%d)%b\n",na,nc,transpose?'T':' ');
    /*
     * In block Wiedemann, nc/na is typically 1+2n/m, which may well be 3
     * for instance if m=n=1. In such a case, the present code will
     * reduce our mp to 2*3 = 6 balanced middle products
     * (Karatsuba-friendly).
     *
     * Other options could be to rely on e.g. a polymat::mp_raw_toom24,
     * which would use the following operations to fall back on only 5
     * balanced middle products:
     *
     * Formulae for MP(2k, 6k). 3k
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
        polymat_mp_raw_subdivide(b, xb, c, xc, nc, a, xa, na, !transpose, add);
        return;
    }
    if (polymat_cutoff_get_subdivide_ub(&polymat_mp_kara_cutoff, na, nc) == 0) {
        polymat_mp_raw_basecase(b, xb, a, xa, na, c, xc, nc, transpose, add);
        return;
    }
    ASSERT(nc >= na);
    for (; nc >= 2 * na - 1;) {
        polymat_mp_raw_kara(b, xb, a, xa, na, c, xc, 2 * na - 1, transpose,
                            add);
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
    size_t const chop = 2 * na - 1 - nc;
    polymat_mp_raw_kara(b, xb, a, xa, na - chop, c, xc + chop, nc - chop,
                        transpose, add);
    polymat_mp_raw_subdivide(b, xb, a, xa + na - chop, chop, c, xc, na - 1,
                             transpose, 1);
} /*}}}*/

static void polymat_mp_raw(/*{{{*/
                           polymat & b, size_t xb,
                           polymat const & a, size_t xa, size_t na,
                           polymat const & c, size_t xc, size_t nc,
                           int transpose,
                           int add)
{
    polymat_mp_raw_subdivide(b, xb, a, xa, na, c, xc, nc, transpose, add);
} /*}}}*/

void polymat::mp(polymat const & a, polymat const & c) /*{{{*/
{
    size_t const nb = std::max(a.size, c.size) - std::min(a.size, c.size) + 1;
    ASSERT_ALWAYS(a.n == c.m);
    if (check_pre_init()) {
        *this = polymat(a.ab, a.m, c.n, nb);
    }
    ASSERT_ALWAYS(m == a.m);
    ASSERT_ALWAYS(n == c.n);
    ASSERT_ALWAYS(alloc >= nb);
    size = nb;

    polymat_mp_raw(*this, 0, a, 0, a.size, c, 0, c.size, 0, 0);
} /*}}}*/

void polymat::addmp(polymat const & a, polymat const & c) /*{{{*/
{
    size_t const nb = std::max(a.size, c.size) - std::min(a.size, c.size) + 1;
    ASSERT_ALWAYS(a.n == c.m);
    if (check_pre_init()) {
        *this = polymat(a.ab, a.m, c.n, nb);
    }
    ASSERT_ALWAYS(m == a.m);
    ASSERT_ALWAYS(n == c.n);
    ASSERT_ALWAYS(alloc >= nb);
    size = nb;

    polymat_mp_raw(*this, 0, a, 0, a.size, c, 0, c.size, 0, 1);
} /*}}}*/

void polymat::set_matpoly(matpoly<false> const & src)
{
    *this = polymat(src.ab, src.m, src.n, src.get_size());

    size = src.get_size();

    for (unsigned int i = 0; i < src.m; i++) {
        for (unsigned int j = 0; j < src.n; j++) {
            for (size_t k = 0; k < src.get_size(); k++) {
                ab->set(coeff(i, j, k), src.coeff(i, j, k));
            }
        }
    }
}
