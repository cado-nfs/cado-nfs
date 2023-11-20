#include "cado.h" // IWYU pragma: keep
#include <climits>                   // for UINT_MAX
#include <cstdlib>
#include <cstring>
#include <algorithm>                  // for min, max
#include <type_traits>                // for integral_constant<>::value, is_...
#include <utility>                    // for move, swap
#include <gmp.h>
#include "gmp_aux.h"                  // for mpn_randomb
#include "lingen_matpoly_binary.hpp"  // for matpoly, matpoly::const_view_t
#include "omp_proxy.h"
#include "macros.h"

#include "arith-hard.hpp"
#include "misc.h"       // cado_ctzl
#include "gf2x.h"

matpoly::memory_pool_type matpoly::memory;

static_assert(std::is_same<unsigned long, mp_limb_t>::value, "need fix for mp_limb_t != unsigned long");

/* {{{ init/zero/clear interface for matpoly */
matpoly::matpoly(matpoly::arith_hard *, unsigned int m, unsigned int n, int len) : m(m), n(n), alloc_words(b2w_x(len)) {
    /* As a special case, we allow a pre-init state with m==n==len==0 */
    /* Note that because we want to handle homogenous and non-homogenous
     * cases the same way, we support matrices of size 0*n, so that is
     * not technically a pre-init state */
    if (!m && !n) {
        ASSERT_ALWAYS(!len);
        return;
    }
    if (alloc_words) {
        if (data_alloc_size_in_bytes()) {
            x = (unsigned long *) memory.alloc(data_alloc_size_in_bytes());
            memset(x, 0, data_alloc_size_in_bytes());
        } else {
            x = NULL;
        }
    }
}
matpoly::~matpoly() {
    if (x)
        memory.free(x, data_alloc_size_in_bytes());
}
matpoly::matpoly(matpoly && a)
    : m(a.m), n(a.n), alloc_words(a.alloc_words)
{
    size=a.size;
    x=a.x;
    a.x=NULL;
    a.m=a.n=a.size=a.alloc_words=0;
    // a.ab=NULL;
}
matpoly& matpoly::operator=(matpoly&& a)
{
    if (x)
        memory.free(x, data_alloc_size_in_bytes());
    // ab = a.ab;
    m = a.m;
    n = a.n;
    alloc_words = a.alloc_words;
    size = a.size;
    x=a.x;
    a.x=NULL;
    a.m=a.n=a.size=a.alloc_words=0;
    // a.ab=NULL;
    return *this;
}
matpoly& matpoly::set(matpoly const& a)
{
    if (x)
        memory.free(x, data_alloc_size_in_bytes());
    // ab = a.ab;
    m = a.m;
    n = a.n;
    alloc_words = a.alloc_words;
    size = a.size;
    // abvec_init(ab, &(x), m*n*alloc);
    x = (unsigned long *) memory.alloc(data_alloc_size_in_bytes());
    memcpy(x, a.x, data_alloc_size_in_bytes());
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
void matpoly::realloc(size_t new_ncoeffs) {
    size_t newalloc_words = b2w_x(new_ncoeffs);
    ASSERT_ALWAYS(b2w(size) <= alloc_words);
    size_t oldmem = m * n * alloc_words * sizeof(unsigned long);
    size_t newmem = m * n * newalloc_words * sizeof(unsigned long);

    /* zero out the newly added data */
    if (newalloc_words > alloc_words) {
        newalloc_words = MAX(newalloc_words, alloc_words + alloc_words / 8);
#if ULONG_BITS < 64
        if (newalloc_words % (64 / ULONG_BITS)) {
            newalloc_words = b2w_x(newalloc_words * ULONG_BITS);
        }
#endif
        /* allocate new space, then inflate */
        x = (unsigned long *) memory.realloc(x, oldmem, newmem);
        const unsigned long * rhead = x + m * n * alloc_words;
        unsigned long * whead = x + m * n * newalloc_words;
        if (size)
            for(unsigned int i = m ; i-- ; ) {
                for(unsigned int j = n ; j-- ; ) {
                    whead -= newalloc_words;
                    rhead -= alloc_words;
                    mpn_copyd(whead, rhead, alloc_words);
                    mpn_zero(whead + alloc_words, newalloc_words - alloc_words);
                }
            }
    } else {
        if (size > new_ncoeffs)
            size = 0;
        /* deflate, then free space */
        ASSERT_ALWAYS(b2w(size) <= newalloc_words);
        const unsigned long * rhead = x;
        unsigned long * whead = x;
        if (size)
            for(unsigned int i = 0 ; i < m ; i++) {
                for(unsigned int j = 0 ; j < n ; j++) {
                    mpn_copyi(whead, rhead, newalloc_words);
                    whead += newalloc_words;
                    rhead += alloc_words;
                }
            }
        x = (unsigned long *) memory.realloc(x, oldmem, newmem);
    }
    alloc_words = newalloc_words;
}
void matpoly::zero() {
    size = 0;
    memset(x, 0, data_alloc_size_in_bytes());
}
size_t matpoly::get_true_nonzero_size() const
{
    size_t lb = 0;
    size_t ub = get_size();
    for(unsigned int ij = 0 ; ij < m*n && lb < ub ; ij++) {
        unsigned int i = ij / n;
        unsigned int j = ij % n;
        /* Find the last nonzero in the range [lb, ub[ */
        for(unsigned int k = ub ; k > lb ; k--) {
            unsigned long x = part_head(i, j, (k-1))[0] >> ((k-1) % ULONG_BITS);
            if (x&1) {
                lb = k;
                break;
            }
        }
    }
    return lb;
}

void matpoly::set_constant_ui(unsigned long e) {
    ASSERT_ALWAYS(m == n);
    size = 0;
    if (alloc_words == 0 && e)
        realloc(1);
    zero();
    if (!e) return;
    size = 1;
    for(unsigned int i = 0 ; i < m ; ++i)
        *part(i, i) = e;
}
/* }}} */

void matpoly::fill_random(unsigned int k0, unsigned int k1, gmp_randstate_t rstate)
{
    ASSERT_ALWAYS(b2w(k1) <= alloc_words);
    size_t nw0 = b2w(k0);
    size_t nw1 = b2w(k1);
    unsigned long nmask1 = ((1UL << (k1 % ULONG_BITS)) - 1);
    unsigned long nmask0 = ((1UL << (k0 % ULONG_BITS)) - 1);
    if (k0 == 0 && k1 == size && nw1 == alloc_words) {
        mpn_randomb(x, rstate, m*n*alloc_words);
        clear_high_word();
    } else if (k0 >= k1) {
        return;
    } else {
        if (k0 % ULONG_BITS) {
            /* Then if nw1==nw0, then k0 can't be a multiple of ULONG_BITS.
             * (but k1 could be)
             */
            unsigned long mask = ~nmask0;
            if (nw1 == nw0 && (k1 % ULONG_BITS))
                mask = mask & nmask1;
            /* put random bits at the right place in the low words, not
             * touching what might be there presently */
            for(unsigned int i = 0 ; i < m ; i++) {
                for(unsigned int j = 0 ; j < n ; j++) {
                    unsigned long * pa = part(i, j);
                    pa[nw0-1] &= ~mask;
                    pa[nw0-1] ^= gmp_urandomb_ui(rstate, ULONG_BITS) & mask;
                }
            }
            k0 = nw0 * ULONG_BITS;
        }
        if (nw1 == nw0)
            return;
        if (nw0 < nw1-1) {
            for(unsigned int i = 0 ; i < m ; i++) {
                for(unsigned int j = 0 ; j < n ; j++) {
                    unsigned long * pa = part(i, j);
                    mpn_randomb(pa + nw0, rstate, nw1 - 1 - nw0);
                }
            }
        }
        for(unsigned int i = 0 ; i < m ; i++) {
            for(unsigned int j = 0 ; j < n ; j++) {
                unsigned long * pa = part(i, j);
                pa[nw1-1] &= ~nmask1;
                pa[nw1-1] ^= gmp_urandomb_ui(rstate, ULONG_BITS) & nmask1;
            }
        }
    }
}

int matpoly::cmp(matpoly const& b) const
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
    } else {
        size_t nw = b2w(size);
        unsigned long mask = (1UL << (size % ULONG_BITS)) - 1;
        ASSERT_ALWAYS(nw > 0);
        for(unsigned int i = 0 ; i < m ; i++) {
            for(unsigned int j = 0 ; j < n ; j++) {
                const unsigned long * pa = part(i,j);
                const unsigned long * pb = b.part(i,j);
                int r = 0;
                if (nw > 1)
                    r = memcmp(pa, pb, (nw-1) * sizeof(unsigned long));
                if (r) return r;
                if (!mask) continue;
                unsigned long aa = pa[nw-1] & mask;
                unsigned long bb = pb[nw-1] & mask;
                r = (aa > bb) - (bb > aa);
                if (r) return r;
            }
        }
        return 0;
    }
}

void matpoly::addpoly(unsigned int i, unsigned int j, matpoly const& y, unsigned int iy, unsigned int jy) {/*{{{*/
    ASSERT(i < nrows());
    ASSERT(j < ncols());
    ASSERT(iy < y.nrows());
    ASSERT(jy < y.ncols());
    unsigned long * dst = part(i, j);
    const unsigned long * src = y.part(iy, jy);
    ASSERT_ALWAYS(y.high_word_is_clear());
    mpn_xor_n(dst, dst, src, std::min(alloc_words, y.alloc_words));
    clear_high_word();
}/*}}}*/
void matpoly::xmul_poly(unsigned int i, unsigned int j, unsigned long s)/*{{{*/
{
    unsigned long * dst = part(i, j);
    ASSERT(1 <= s && s <= GMP_LIMB_BITS-1);
    mpn_lshift(dst, dst, alloc_words, s);
}/*}}}*/

/* shift by a multiplication by x all coefficients of degree less than
 * colsize in column j. This results in a new column of length colsize+1.
 * Allocation must be sufficient for this length to fit.  What happens to
 * coefficients of degree larger than colsize in the result is unspecified.
 *
 * It is often relevant to "colsize++" right after this call, since the
 * coefficient of degree colsize is well-defined on output
 *
 * Other columns are unaffected.
 */
void matpoly::multiply_column_by_x(unsigned int j, unsigned int colsize)/*{{{*/
{
    unsigned int k = colsize + 1;
    size_t kw = b2w(k);
    ASSERT_ALWAYS(kw <= alloc_words);
    for(unsigned int i = 0 ; i < m ; i++) {
        unsigned long * pa = part(i, j);
        mpn_lshift(pa, pa, kw, 1);
    }
}/*}}}*/

/* shift right (by a division by x) all coefficients of degree less than
 * colsize in column j. This results in a new column of length colsize-1,
 * with a zero coefficient shifted in at degree colsize-1.
 *
 * It is often relevant to "colsize--" right after this call.
 */
void matpoly::divide_column_by_x(unsigned int j, unsigned int colsize)/*{{{*/
{
    size_t nw = b2w(colsize);
    ASSERT_ALWAYS(nw <= alloc_words);
    ASSERT_ALWAYS(colsize);
    unsigned int k = colsize - 1;
    size_t kw = b2w(k);
    unsigned long kmask = ((1UL << (k % ULONG_BITS)) - 1);
    for(unsigned int i = 0 ; i < m ; i++) {
        unsigned long * pa = part(i, j);
        mpn_rshift(pa, pa, nw, 1);
        if (kmask)
            pa[kw - 1] &= kmask;
    }
}/*}}}*/

void matpoly::truncate(matpoly const & src, unsigned int new_ncoeffs)/*{{{*/
{
    ASSERT_ALWAYS(b2w(new_ncoeffs) <= src.alloc_words);
    if (check_pre_init()) {
        /* Need to call ctor... */
        *this = matpoly(src.ab, src.m, src.n, new_ncoeffs);
    }
    ASSERT_ALWAYS(m == src.m);
    ASSERT_ALWAYS(n == src.n);
    ASSERT_ALWAYS(b2w(new_ncoeffs) <= alloc_words);
    ASSERT_ALWAYS(new_ncoeffs <= src.size);
    size = new_ncoeffs;
    if (this != &src) {
        for(unsigned int i = 0 ; i < src.m ; i++) {
            for(unsigned int j = 0 ; j < src.n ; j++) {
                mpn_copyi(part(i, j), src.part(i, j), b2w(new_ncoeffs));
            }
        }
    }
    clear_high_word();
}/*}}}*/
int matpoly::tail_is_zero(unsigned int k) const /*{{{*/
{
    if (k == size) return 1;
    ASSERT_ALWAYS(k < size);
    size_t kw = b2w(k);
    unsigned long kmask = ((1UL << (k % ULONG_BITS)) - 1);
    size_t nw = b2w(size);
    unsigned long nmask = ((1UL << (size % ULONG_BITS)) - 1);
    for(unsigned int i = 0 ; i < m ; i++) {
        for(unsigned int j = 0 ; j < n ; j++) {
            const unsigned long * pa = part(i, j);
            if (kw < nw) {
                if (kmask && (pa[kw - 1] & ~kmask)) return 0;
                for(unsigned int s =  kw ; s < nw - 1 ; s++)
                    if (pa[s]) return 0;
                if (nmask && (pa[nw - 1] & nmask)) return 0;
            } else {
                /* kw == nw, which implies that if kmask can't be zero
                 * (since we would have size > k = kw*64). We may, on
                 * the other hand, have nmask==0. */
                ASSERT_ALWAYS(kw == nw);
                ASSERT_ALWAYS(kmask);
                unsigned long z = pa[kw - 1] & ~kmask;
                if (nmask) z &= nmask;
                if (z) return 0;
            }
        }
    }
    return 1;
}/*}}}*/
void matpoly::clear_high_word_common(unsigned int length)/*{{{*/
{
    ASSERT_ALWAYS(length <= capacity());
    size_t nw = b2w(length);
    unsigned long nmask = ((1UL << (length % ULONG_BITS)) - 1);
    if (nmask == 0) return;
    for(unsigned int i = 0 ; i < m ; i++) {
        for(unsigned int j = 0 ; j < n ; j++) {
            unsigned long * pa = part(i, j);
            pa[nw - 1] &= nmask;
        }
    }
}/*}}}*/
bool matpoly::high_word_is_clear() const/*{{{*/
{
    size_t nw = b2w(size);
    unsigned long nmask = ((1UL << (size % ULONG_BITS)) - 1);
    if (nmask == 0) return true;
    for(unsigned int i = 0 ; i < m ; i++) {
        for(unsigned int j = 0 ; j < n ; j++) {
            const unsigned long * pa = part(i, j);
            if (pa[nw - 1] & ~nmask) return false;
        }
    }
    return true;
}/*}}}*/
void matpoly::zero_pad(unsigned int k)/*{{{*/
{
    ASSERT_ALWAYS(k >= size);
    if (check_pre_init() || b2w(k) > alloc_words)
        realloc(k);
    size_t kw = b2w(k);
    size_t nw = b2w(size);
    unsigned long nmask = ((1UL << (size % ULONG_BITS)) - 1);
    for(unsigned int i = 0 ; i < m ; i++) {
        for(unsigned int j = 0 ; j < n ; j++) {
            unsigned long * pa = part(i, j);
            if (nmask) pa[nw - 1] &= nmask;
            mpn_zero(pa + nw, kw - nw);
        }
    }
    size = k;
}/*}}}*/

/* This takes coefficient ksrc of column jsrc, and copies it to
 * coefficient kdst of column jdst
 */
void matpoly::extract_column( /*{{{*/
        unsigned int jdst, unsigned int kdst,
        matpoly const & src, unsigned int jsrc, unsigned int ksrc)
{
    ASSERT_ALWAYS(m == src.m);
    size_t sw = ksrc / ULONG_BITS;
    unsigned long sbit = 1UL << (ksrc % ULONG_BITS);
    size_t dw = kdst / ULONG_BITS;
    unsigned long dbit = 1UL << (kdst % ULONG_BITS);
    ASSERT_ALWAYS(sw <= src.alloc_words);
    ASSERT_ALWAYS(dw <= alloc_words);

    for(unsigned int i = 0 ; i < m ; i++) {
        const unsigned long * ps = src.part(i, jsrc);
        unsigned long * pd = part(i, jdst);
        pd[dw] = (pd[dw] & ~dbit) | (dbit & -((ps[sw] & sbit) != 0));
    }
}/*}}}*/

void matpoly::zero_column(unsigned int jdst, unsigned int kdst) /*{{{*/
{
    size_t dw = kdst / ULONG_BITS;
    unsigned long dbit = 1UL << (kdst % ULONG_BITS);
    ASSERT_ALWAYS(dw <= alloc_words);

    for(unsigned int i = 0 ; i < m ; i++) {
        unsigned long * pd = part(i, jdst);
        pd[dw] = pd[dw] & ~dbit;
    }
}/*}}}*/

/* Copy nb bits from c1 to c, starting from bit (shift) */
static inline void CopyBitsRsh(unsigned long * c, const unsigned long * c1, size_t bits_c, size_t shift)
{
    size_t t = BITS_TO_WORDS(bits_c, ULONG_BITS);
    size_t words_full = BITS_TO_WORDS(bits_c + shift, ULONG_BITS);
    size_t pick = shift / ULONG_BITS;
    size_t cnt = shift % ULONG_BITS;
    size_t tnc = ULONG_BITS - cnt;
    if (cnt) {
        mpn_rshift(c, c1 + pick, t, cnt);
    } else {
        mpn_copyi(c, c1 + pick, t);
    }
    /* words_full - pick - t is either 0 or 1 */
    if (words_full - pick == t + 1)
        c[t - 1] |= c1[pick + t] << tnc;
    if (bits_c % ULONG_BITS)
        c[bits_c / ULONG_BITS] &= (1UL << (bits_c % ULONG_BITS)) - 1;
}

void matpoly::rshift(matpoly const & src, unsigned int k)/*{{{*/
{
    ASSERT_ALWAYS(k <= src.size);
    unsigned int newsize = src.size - k;
    if (check_pre_init()) {
        *this = matpoly(src.ab, src.m, src.n, newsize);
    }
    ASSERT_ALWAYS(m == src.m);
    ASSERT_ALWAYS(n == src.n);
    for(unsigned int i = 0 ; i < src.m ; i++) {
        for(unsigned int j = 0 ; j < src.n ; j++) {
            const unsigned long * ps = src.part(i, j);
            unsigned long * pd = part(i, j);
            CopyBitsRsh(pd, ps, src.size - k, k);
        }
    }
    size = newsize;
}/*}}}*/
void matpoly::rshift(unsigned int k)/*{{{*/
{
    ASSERT_ALWAYS(k <= size);
    unsigned int newsize = size - k;
    if (newsize) {
        for(unsigned int i = 0 ; i < m ; i++) {
            for(unsigned int j = 0 ; j < n ; j++) {
                const unsigned long * ps = part(i, j);
                unsigned long * pd = part(i, j);
                CopyBitsRsh(pd, ps, size - k, k);
            }
        }
    }
    size = newsize;
}/*}}}*/
unsigned int matpoly::valuation() const /*{{{*/
{
    auto isnz = [this](unsigned int k) {
        for(unsigned int i = 0 ; i < m ; ++i) {
            for(unsigned int j = 0 ; j < n ; ++j) {
                const mp_limb_t * z = part(i, j);
                if (z[k]) return 1;
            }
        }
        return 0;
    };
    unsigned int k = 0;
    for( ; k < b2w(size) ; k++) {
        if (isnz(k)) break;
    }
    if (k >= b2w(size)) return UINT_MAX;
    mp_limb_t x = 0;
    for(unsigned int i = 0 ; i < m ; ++i) {
        for(unsigned int j = 0 ; j < n ; ++j) {
            const mp_limb_t * z = part(i, j);
            x |= z[k];
        }
    }
    k = k * ULONG_BITS + cado_ctzl(x);
    if (k >= size) return UINT_MAX;
    return k;
}/*}}}*/

void matpoly::view_t::zero() { /*{{{*/
    unsigned int nrows = this->nrows();
    unsigned int ncols = this->ncols();
#ifdef HAVE_OPENMP
#pragma omp parallel for collapse(2)
#endif
    for(unsigned int i = 0 ; i < nrows ; i++) {
        for(unsigned int j = 0 ; j < ncols ; j++) {
            memset(part(i,j), 0, M.data_entry_alloc_size_in_bytes());
        }
    }
}/*}}}*/

void matpoly::add(matpoly const & a, matpoly const & b)/*{{{*/
{
    size_t new_ncoeffs = std::max(a.size, b.size);
    ASSERT_ALWAYS(a.m == b.m);
    ASSERT_ALWAYS(a.n == b.n);
    if (check_pre_init()) {
        *this = matpoly(a.ab, a.m, a.n, new_ncoeffs);
    }
    if (alloc_words < b2w(new_ncoeffs))
        realloc(new_ncoeffs);
    ASSERT_ALWAYS(b2w(new_ncoeffs) <= alloc_words);

    for(unsigned int i = 0 ; i < m ; ++i) {
        for(unsigned int j = 0 ; j < n ; ++j) {
            size_t s0 = std::min(a.size, b.size);
            size_t s0q = s0 / ULONG_BITS;
            size_t s0w = b2w(s0);
            size_t aw = b2w(a.size);
            size_t bw = b2w(b.size);
            unsigned long s0mask = ((1UL << (s0 % ULONG_BITS)) - 1);
            mp_limb_t * z = part(i, j);
            const mp_limb_t * az = a.part(i, j);
            const mp_limb_t * bz = b.part(i, j);
            if (s0q) mpn_xor_n(z, az, bz, s0q);
            /* Note that we haven't xored the high bits yet ! */
            if (a.size > s0) {
                if (s0mask) z[s0q] = az[s0q] ^ (bz[s0q] & s0mask);
                if (aw > s0w) mpn_copyi(z + s0w, az + s0w, aw - s0w);
            } else if (b.size > s0) {
                if (s0mask) z[s0q] = (az[s0q] & s0mask) ^ bz[s0q];
                if (bw > s0w) mpn_copyi(z + s0w, bz + s0w, bw - s0w);
            } else {
                if (s0mask) z[s0q] = (az[s0q] ^ bz[s0q]) & s0mask;
            }
        }
    }
    size = new_ncoeffs;
}/*}}}*/
void matpoly::sub(matpoly const & a, matpoly const & b)/*{{{*/
{
    add(a, b);
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
    ASSERT_ALWAYS(b2w(csize) <= alloc_words);

    if (a.size == 0 || b.size == 0)
        return;

    if (csize >= size) zero_pad(csize);

    ASSERT_ALWAYS(a.high_word_is_clear());
    ASSERT_ALWAYS(b.high_word_is_clear());

    matpoly::addmul(*this, a, b);

    size = std::max(size, csize);
}/*}}}*/

void matpoly::copy(matpoly::view_t t, matpoly::const_view_t a)
{
    unsigned int nrows = a.nrows();
    unsigned int ncols = a.ncols();
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
                memcpy(tij, aij, a.M.data_entry_alloc_size_in_bytes());
            }
        }
    }
}

void matpoly::addmul(matpoly::view_t t, matpoly::const_view_t t0, matpoly::const_view_t t1)/*{{{*/
{
    unsigned int nrows = t0.nrows();
    unsigned int ncols = t1.ncols();
    size_t csize = t0.M.size + t1.M.size; csize -= (csize > 0);
    ASSERT_ALWAYS(t0.ncols() == t1.nrows());
    ASSERT_ALWAYS(t.nrows() == nrows);
    ASSERT_ALWAYS(t.ncols() == ncols);
    ASSERT_ALWAYS(b2w(csize) <= t.M.alloc_words);

    if (t0.M.size == 0 || t1.M.size == 0)
        return;

#ifdef HAVE_OPENMP
    unsigned int T = std::min((unsigned int) omp_get_max_threads(), t0.nrows() * t1.ncols());
#pragma omp parallel num_threads(T)
#endif
    {
        unsigned long * tmp = (unsigned long *) malloc((b2w(t0.M.size) + b2w(t1.M.size)) * sizeof(unsigned long));
#ifdef HAVE_OPENMP
#pragma omp for collapse(2)
#endif
        for(unsigned int i = 0 ; i < nrows ; i++) {
            for(unsigned int j = 0 ; j < ncols ; j++) {
                for(unsigned int k = 0 ; k < t0.ncols() ; k++) {
                    gf2x_mul(tmp,
                            t0.part(i, k), b2w(t0.M.size),
                            t1.part(k, j), b2w(t1.M.size));
                    mpn_xor_n(t.part(i, j), t.part(i, j), tmp, b2w(csize));
                }
            }
        }
        free(tmp);
    }
}/*}}}*/

matpoly matpoly::mul(matpoly const & a, matpoly const & b)/*{{{*/
{
    size_t csize = a.size + b.size; csize -= (csize > 0);
    matpoly c(a.ab, a.m, b.n, csize);
    ASSERT_ALWAYS(a.n == b.m);
    c.zero();
    c.addmul(a, b);
    return c;
}/*}}}*/

void matpoly::addmp(matpoly const & a, matpoly const & c)/*{{{*/
{
    size_t fullsize = a.size + c.size; fullsize -= (fullsize > 0);
    size_t nb = MAX(a.size, c.size) - MIN(a.size, c.size) + 1;
    ASSERT_ALWAYS(a.n == c.m);
    if (check_pre_init()) {
        *this = matpoly(a.ab, a.m, c.n, nb);
    }
    ASSERT_ALWAYS(m == a.m);
    ASSERT_ALWAYS(n == c.n);
    ASSERT_ALWAYS(b2w(nb) <= alloc_words);

    if (a.size == 0 || c.size == 0)
        return;

    if (nb >= size) zero_pad(nb);

    matpoly::addmp(*this, a, c);

    size = std::max(size, nb);
}/*}}}*/

void matpoly::addmp(matpoly::view_t t, matpoly::const_view_t t0, matpoly::const_view_t t1)/*{{{*/
{
    unsigned int nrows = t0.nrows();
    unsigned int ncols = t1.ncols();
    size_t fullsize = t0.M.size + t1.M.size; fullsize -= (fullsize > 0);
    size_t shift = MIN(t0.M.size, t1.M.size) - 1;
    size_t nb = MAX(t0.M.size, t1.M.size) - MIN(t0.M.size, t1.M.size) + 1;
    ASSERT_ALWAYS(t0.ncols() == t1.nrows());
    ASSERT_ALWAYS(t.nrows() == nrows);
    ASSERT_ALWAYS(t.ncols() == ncols);
    ASSERT_ALWAYS(b2w(nb) <= t.M.alloc_words);

    if (t0.M.size == 0 || t1.M.size == 0)
        return;

#ifdef HAVE_OPENMP
    unsigned int T = std::min((unsigned int) omp_get_max_threads(), t0.nrows() * t1.ncols());
#pragma omp parallel num_threads(T)
#endif
    {
        unsigned long * tmp = (unsigned long *) malloc((b2w(t0.M.size) + b2w(t1.M.size)) * sizeof(unsigned long));
#ifdef HAVE_OPENMP
#pragma omp for collapse(2)
#endif
        for(unsigned int i = 0 ; i < nrows ; i++) {
            for(unsigned int j = 0 ; j < ncols ; j++) {
                for(unsigned int k = 0 ; k < t0.ncols() ; k++) {
                    gf2x_mul(tmp,
                            t0.part(i, k), b2w(t0.M.size),
                            t1.part(k, j), b2w(t1.M.size));
                    CopyBitsRsh(tmp, tmp, nb, shift);
                    mpn_xor_n(t.part(i, j), t.part(i, j), tmp, b2w(nb));
                }
            }
        }
        free(tmp);
    }
}/*}}}*/

matpoly matpoly::mp(matpoly const & a, matpoly const & c)/*{{{*/
{
    unsigned int nb = MAX(a.size, c.size) - MIN(a.size, c.size) + 1;
    ASSERT_ALWAYS(a.n == c.m);
    matpoly b(a.ab, a.m, c.n, nb);
    b.zero();
    b.addmp(a, c);
    return b;
}/*}}}*/

int matpoly::coeff_is_zero(unsigned int k) const
{
    unsigned int kq = k / ULONG_BITS;
    unsigned long kbit = 1UL << (k % ULONG_BITS);
    for(unsigned int i = 0 ; i < m ; i++)
        for(unsigned int j = 0; j < n; j++)
            if (part(i, j)[kq] & kbit)
                return 0;
    return 1;
}
void matpoly::coeff_set_zero(unsigned int k)
{
    unsigned int kq = k / ULONG_BITS;
    unsigned long kbit = 1UL << (k % ULONG_BITS);
    for(unsigned int i = 0 ; i < m ; i++)
        for(unsigned int j = 0; j < n; j++)
            part(i, j)[kq] &= ~kbit;
}

matpoly matpoly::truncate_and_rshift(unsigned int truncated_size, unsigned int shiftcount)
{
    matpoly other(ab, m, n, size - shiftcount);
    other.rshift(*this, shiftcount);
    truncate(*this, truncated_size);
    shrink_to_fit();
    std::swap(*this, other);
    return other;
}
