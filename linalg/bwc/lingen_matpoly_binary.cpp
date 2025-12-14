#include "cado.h" // IWYU pragma: keep

#include <climits>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <type_traits>
#include <memory>
#include <utility>

#include <gmp.h>

#include "gmp_aux.h"
#include "lingen_matpoly_binary.hpp"
#include "omp_proxy.h"
#include "macros.h"
#include "arith-hard.hpp"
#include "misc.h"
#include "gf2x.h"
#include "runtime_numeric_cast.hpp"

// NOLINTBEGIN(readability-static-accessed-through-instance)

matpoly<true>::memory_pool_type matpoly<true>::memory;

static_assert(std::is_same_v<unsigned long, mp_limb_t>, "need fix for mp_limb_t != unsigned long");

/* {{{ init/zero/clear interface for matpoly */
matpoly<true>::matpoly(matpoly<true>::arith_hard *, unsigned int m, unsigned int n, size_t len) : m(m), n(n), alloc_words(b2w_x(len)) {
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
            x = memory.alloc(data_alloc_size_in_bytes());
            memset(x, 0, data_alloc_size_in_bytes());
        } else {
            x = nullptr;
        }
    }
}
matpoly<true>::~matpoly() {
    if (x)
        memory.free(x, data_alloc_size_in_bytes());
}
matpoly<true>::matpoly(matpoly && a) noexcept
    : m(a.m)
    , n(a.n)
    , size(a.size)
    , alloc_words(a.alloc_words)
    , x(a.x)
{
    a.x=nullptr;
    a.m=a.n=a.size=a.alloc_words=0;
    // a.ab=nullptr;
}
matpoly<true>& matpoly<true>::operator=(matpoly&& a) noexcept
{
    if (x)
        memory.free(x, data_alloc_size_in_bytes());
    // ab = a.ab;
    m = a.m;
    n = a.n;
    alloc_words = a.alloc_words;
    size = a.size;
    x=a.x;
    a.x=nullptr;
    a.m=a.n=a.size=a.alloc_words=0;
    // a.ab=nullptr;
    return *this;
}
matpoly<true>& matpoly<true>::set(matpoly const& a)
{
    if (x)
        memory.free(x, data_alloc_size_in_bytes());
    // ab = a.ab;
    m = a.m;
    n = a.n;
    alloc_words = a.alloc_words;
    size = a.size;
    // abvec_init(ab, &(x), m*n*alloc);
    x = memory.alloc(data_alloc_size_in_bytes());
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
void matpoly<true>::realloc(size_t new_ncoeffs) {
    size_t newalloc_words = b2w_x(new_ncoeffs);
    ASSERT_ALWAYS(b2w(size) <= alloc_words);
    size_t const oldmem = m * n * alloc_words * sizeof(unsigned long);
    size_t const newmem = m * n * newalloc_words * sizeof(unsigned long);

    /* zero out the newly added data */
    if (newalloc_words > alloc_words) {
        newalloc_words = MAX(newalloc_words, alloc_words + alloc_words / 8);
#if ULONG_BITS < 64
        if (newalloc_words % (64 / ULONG_BITS)) {
            newalloc_words = b2w_x(newalloc_words * ULONG_BITS);
        }
#endif
        /* allocate new space, then inflate */
        x = memory.realloc(x, oldmem, newmem);
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
        x = memory.realloc(x, oldmem, newmem);
    }
    alloc_words = newalloc_words;
}
void matpoly<true>::zero() {
    size = 0;
    memset(x, 0, data_alloc_size_in_bytes());
}
size_t matpoly<true>::get_true_nonzero_size() const
{
    size_t lb = 0;
    size_t const ub = get_size();
    for(unsigned int ij = 0 ; ij < m*n && lb < ub ; ij++) {
        unsigned int const i = ij / n;
        unsigned int const j = ij % n;
        /* Find the last nonzero in the range [lb, ub[ */
        for(size_t k = ub ; k > lb ; k--) {
            unsigned long const x = part_head(i, j, (k-1))[0] >> ((k-1) % ULONG_BITS);
            if (x&1) {
                lb = k;
                break;
            }
        }
    }
    return lb;
}

void matpoly<true>::set_constant_ui(unsigned long e) {
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

void matpoly<true>::fill_random(size_t k0, size_t k1, cxx_gmp_randstate & rstate)
{
    ASSERT_ALWAYS(b2w(k1) <= alloc_words);
    size_t const nw0 = b2w(k0);
    size_t const nw1 = b2w(k1);
    unsigned long const nmask1 = ((1UL << (k1 % ULONG_BITS)) - 1);
    unsigned long const nmask0 = ((1UL << (k0 % ULONG_BITS)) - 1);
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

int matpoly<true>::cmp(matpoly const& b) const
{
    ASSERT_ALWAYS(n == b.n);
    ASSERT_ALWAYS(m == b.m);
    if (size != b.size) {
        return (size > b.size) - (b.size > size);
    } else if (size == 0 && b.size == 0) {
        /* This is for the "intermediary pre-init" state. size is zero,
         * but a priori the dimension fields are ok. x might be nullptr or
         * not, it depends on the previous life of the object */
        return 0;
    } else {
        size_t const nw = b2w(size);
        unsigned long const mask = (1UL << (size % ULONG_BITS)) - 1;
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
                unsigned long const aa = pa[nw-1] & mask;
                unsigned long const bb = pb[nw-1] & mask;
                r = (aa > bb) - (bb > aa);
                if (r) return r;
            }
        }
        return 0;
    }
}

void matpoly<true>::addpoly(unsigned int i, unsigned int j, matpoly const& y, unsigned int iy, unsigned int jy) {/*{{{*/
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
void matpoly<true>::xmul_poly(unsigned int i, unsigned int j, unsigned long s)/*{{{*/
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
void matpoly<true>::multiply_column_by_x(unsigned int j, size_t colsize)/*{{{*/
{
    size_t const k = colsize + 1;
    size_t const kw = b2w(k);
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
void matpoly<true>::divide_column_by_x(unsigned int j, size_t colsize)/*{{{*/
{
    size_t const nw = b2w(colsize);
    ASSERT_ALWAYS(nw <= alloc_words);
    ASSERT_ALWAYS(colsize);
    size_t const k = colsize - 1;
    size_t const kw = b2w(k);
    unsigned long const kmask = ((1UL << (k % ULONG_BITS)) - 1);
    for(unsigned int i = 0 ; i < m ; i++) {
        unsigned long * pa = part(i, j);
        mpn_rshift(pa, pa, nw, 1);
        if (kmask)
            pa[kw - 1] &= kmask;
    }
}/*}}}*/

void matpoly<true>::truncate(matpoly const & src, size_t new_ncoeffs)/*{{{*/
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
int matpoly<true>::tail_is_zero(size_t k) const /*{{{*/
{
    if (k == size) return 1;
    ASSERT_ALWAYS(k < size);
    size_t const kw = b2w(k);
    unsigned long const kmask = ((1UL << (k % ULONG_BITS)) - 1);
    size_t const nw = b2w(size);
    unsigned long const nmask = ((1UL << (size % ULONG_BITS)) - 1);
    for(unsigned int i = 0 ; i < m ; i++) {
        for(unsigned int j = 0 ; j < n ; j++) {
            const unsigned long * pa = part(i, j);
            if (kw < nw) {
                if (kmask && (pa[kw - 1] & ~kmask)) return 0;
                for(size_t s =  kw ; s < nw - 1 ; s++)
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
void matpoly<true>::clear_high_word_common(size_t length)/*{{{*/
{
    ASSERT_ALWAYS(length <= capacity());
    size_t const nw = b2w(length);
    unsigned long const nmask = ((1UL << (length % ULONG_BITS)) - 1);
    if (nmask == 0) return;
    for(unsigned int i = 0 ; i < m ; i++) {
        for(unsigned int j = 0 ; j < n ; j++) {
            unsigned long * pa = part(i, j);
            pa[nw - 1] &= nmask;
        }
    }
}/*}}}*/
bool matpoly<true>::high_word_is_clear() const/*{{{*/
{
    size_t const nw = b2w(size);
    unsigned long const nmask = ((1UL << (size % ULONG_BITS)) - 1);
    if (nmask == 0) return true;
    for(unsigned int i = 0 ; i < m ; i++) {
        for(unsigned int j = 0 ; j < n ; j++) {
            const unsigned long * pa = part(i, j);
            if (pa[nw - 1] & ~nmask) return false;
        }
    }
    return true;
}/*}}}*/
void matpoly<true>::zero_pad(size_t k)/*{{{*/
{
    ASSERT_ALWAYS(k >= size);
    if (check_pre_init() || b2w(k) > alloc_words)
        realloc(k);
    size_t const kw = b2w(k);
    size_t const nw = b2w(size);
    unsigned long const nmask = ((1UL << (size % ULONG_BITS)) - 1);
    for(unsigned int i = 0 ; i < m ; i++) {
        for(unsigned int j = 0 ; j < n ; j++) {
            unsigned long * pa = part(i, j);
            if (nmask) pa[nw - 1] &= nmask;
            mpn_zero(pa + nw, runtime_numeric_cast<mp_size_t>(kw - nw));
        }
    }
    size = k;
}/*}}}*/

/* This takes coefficient ksrc of column jsrc, and copies it to
 * coefficient kdst of column jdst
 */
void matpoly<true>::extract_column( /*{{{*/
        unsigned int jdst, size_t kdst,
        matpoly const & src, unsigned int jsrc, size_t ksrc)
{
    ASSERT_ALWAYS(m == src.m);
    size_t const sw = ksrc / ULONG_BITS;
    unsigned long const sbit = 1UL << (ksrc % ULONG_BITS);
    size_t const dw = kdst / ULONG_BITS;
    unsigned long const dbit = 1UL << (kdst % ULONG_BITS);
    ASSERT_ALWAYS(sw <= src.alloc_words);
    ASSERT_ALWAYS(dw <= alloc_words);

    for(unsigned int i = 0 ; i < m ; i++) {
        const unsigned long * ps = src.part(i, jsrc);
        unsigned long * pd = part(i, jdst);
        pd[dw] = (pd[dw] & ~dbit) | (dbit & -((ps[sw] & sbit) != 0));
    }
}/*}}}*/

void matpoly<true>::zero_column(unsigned int jdst, size_t kdst) /*{{{*/
{
    size_t const dw = kdst / ULONG_BITS;
    unsigned long const dbit = 1UL << (kdst % ULONG_BITS);
    ASSERT_ALWAYS(dw <= alloc_words);

    for(unsigned int i = 0 ; i < m ; i++) {
        unsigned long * pd = part(i, jdst);
        pd[dw] = pd[dw] & ~dbit;
    }
}/*}}}*/

/* Copy nb bits from c1 to c, starting from bit (shift) */
static inline void CopyBitsRsh(unsigned long * c, const unsigned long * c1, size_t bits_c, size_t shift)
{
    size_t const t = BITS_TO_WORDS(bits_c, ULONG_BITS);
    size_t const words_full = BITS_TO_WORDS(bits_c + shift, ULONG_BITS);
    size_t const pick = shift / ULONG_BITS;
    size_t const cnt = shift % ULONG_BITS;
    size_t const tnc = ULONG_BITS - cnt;
    if (cnt) {
        mpn_rshift(c, c1 + pick, 
                runtime_numeric_cast<mp_size_t>(t), cnt);
        /* words_full - pick - t is either 0 or 1 */
        if (words_full - pick == t + 1)
            c[t - 1] |= c1[pick + t] << tnc;
    } else {
        mpn_copyi(c, c1 + pick, 
                runtime_numeric_cast<mp_size_t>(t));
        ASSERT(words_full == pick + t);
    }
    if (bits_c % ULONG_BITS)
        c[bits_c / ULONG_BITS] &= (1UL << (bits_c % ULONG_BITS)) - 1;
}

void matpoly<true>::rshift(matpoly const & src, size_t k)/*{{{*/
{
    ASSERT_ALWAYS(k <= src.size);
    size_t const newsize = src.size - k;
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
void matpoly<true>::rshift(size_t k)/*{{{*/
{
    ASSERT_ALWAYS(k <= size);
    size_t const newsize = size - k;
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
size_t matpoly<true>::valuation() const /*{{{*/
{
    auto isnz = [this](size_t k) {
        for(unsigned int i = 0 ; i < m ; ++i) {
            for(unsigned int j = 0 ; j < n ; ++j) {
                const mp_limb_t * z = part(i, j);
                if (z[k]) return 1;
            }
        }
        return 0;
    };
    size_t k = 0;
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

void matpoly<true>::view_t::zero() { /*{{{*/
    unsigned int const nrows = this->nrows();
    unsigned int const ncols = this->ncols();
#ifdef HAVE_OPENMP
#pragma omp parallel for collapse(2)
#endif
    for(unsigned int i = 0 ; i < nrows ; i++) {
        for(unsigned int j = 0 ; j < ncols ; j++) {
            memset(part(i,j), 0, M.data_entry_alloc_size_in_bytes());
        }
    }
}/*}}}*/

void matpoly<true>::add(matpoly const & a, matpoly const & b)/*{{{*/
{
    size_t const new_ncoeffs = std::max(a.size, b.size);
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
            size_t const s0 = std::min(a.size, b.size);
            size_t const s0q = s0 / ULONG_BITS;
            size_t const s0w = b2w(s0);
            size_t const aw = b2w(a.size);
            size_t const bw = b2w(b.size);
            unsigned long const s0mask = ((1UL << (s0 % ULONG_BITS)) - 1);
            mp_limb_t * z = part(i, j);
            const mp_limb_t * az = a.part(i, j);
            const mp_limb_t * bz = b.part(i, j);
            if (s0q) mpn_xor_n(z, az, bz,
                    runtime_numeric_cast<mp_size_t>(s0q));
            /* Note that we haven't xored the high bits yet ! */
            if (a.size > s0) {
                if (s0mask) z[s0q] = az[s0q] ^ (bz[s0q] & s0mask);
                if (aw > s0w) mpn_copyi(z + s0w, az + s0w,
                        runtime_numeric_cast<mp_size_t>(aw - s0w));
            } else if (b.size > s0) {
                if (s0mask) z[s0q] = (az[s0q] & s0mask) ^ bz[s0q];
                if (bw > s0w) mpn_copyi(z + s0w, bz + s0w,
                        runtime_numeric_cast<mp_size_t>(bw - s0w));
            } else {
                if (s0mask) z[s0q] = (az[s0q] ^ bz[s0q]) & s0mask;
            }
        }
    }
    size = new_ncoeffs;
}/*}}}*/
void matpoly<true>::sub(matpoly const & a, matpoly const & b)/*{{{*/
{
    add(a, b);
}/*}}}*/

void matpoly<true>::addmul(matpoly const & a, matpoly const & b)/*{{{*/
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

    matpoly<true>::addmul(*this, a, b);

    size = std::max(size, csize);
}/*}}}*/

void matpoly<true>::copy(matpoly<true>::view_t t, matpoly<true>::const_view_t a)
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
                memcpy(tij, aij, a.M.data_entry_alloc_size_in_bytes());
            }
        }
    }
}

void matpoly<true>::addmul(matpoly<true>::view_t t, matpoly<true>::const_view_t t0, matpoly<true>::const_view_t t1)/*{{{*/
{
    unsigned int const nrows = t0.nrows();
    unsigned int const ncols = t1.ncols();
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
        const std::unique_ptr<unsigned long[]> tmp(
                new unsigned long[b2w(t0.M.size) + b2w(t1.M.size)]);
#ifdef HAVE_OPENMP
#pragma omp for collapse(2)
#endif
        for(unsigned int i = 0 ; i < nrows ; i++) {
            for(unsigned int j = 0 ; j < ncols ; j++) {
                for(unsigned int k = 0 ; k < t0.ncols() ; k++) {
                    gf2x_mul(tmp.get(),
                            t0.part(i, k), b2w(t0.M.size),
                            t1.part(k, j), b2w(t1.M.size));
                    mpn_xor_n(t.part(i, j), t.part(i, j), tmp.get(),
                            runtime_numeric_cast<mp_size_t>(b2w(csize)));
                }
            }
        }
    }
}/*}}}*/

matpoly<true> matpoly<true>::mul(matpoly const & a, matpoly const & b)/*{{{*/
{
    size_t csize = a.size + b.size; csize -= (csize > 0);
    matpoly c(a.ab, a.m, b.n, csize);
    ASSERT_ALWAYS(a.n == b.m);
    c.zero();
    c.addmul(a, b);
    return c;
}/*}}}*/

void matpoly<true>::addmp(matpoly const & a, matpoly const & c)/*{{{*/
{
    // size_t fullsize = a.size + c.size; fullsize -= (fullsize > 0);
    size_t const nb = MAX(a.size, c.size) - MIN(a.size, c.size) + 1;
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

    matpoly<true>::addmp(*this, a, c);

    size = std::max(size, nb);
}/*}}}*/

void matpoly<true>::addmp(matpoly<true>::view_t t, matpoly<true>::const_view_t t0, matpoly<true>::const_view_t t1)/*{{{*/
{
    unsigned int const nrows = t0.nrows();
    unsigned int const ncols = t1.ncols();
    // size_t fullsize = t0.M.size + t1.M.size; fullsize -= (fullsize > 0);
    size_t const shift = MIN(t0.M.size, t1.M.size) - 1;
    size_t const nb = MAX(t0.M.size, t1.M.size) - MIN(t0.M.size, t1.M.size) + 1;
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
        /* It's possibly slightly more than fullsize */
        const std::unique_ptr<unsigned long[]> tmp(
                new unsigned long[b2w(t0.M.size) + b2w(t1.M.size)]);
#ifdef HAVE_OPENMP
#pragma omp for collapse(2)
#endif
        for(unsigned int i = 0 ; i < nrows ; i++) {
            for(unsigned int j = 0 ; j < ncols ; j++) {
                for(unsigned int k = 0 ; k < t0.ncols() ; k++) {
                    gf2x_mul(tmp.get(),
                            t0.part(i, k), b2w(t0.M.size),
                            t1.part(k, j), b2w(t1.M.size));
                    CopyBitsRsh(tmp.get(), tmp.get(), nb, shift);
                    mpn_xor_n(t.part(i, j), t.part(i, j), tmp.get(),
                            runtime_numeric_cast<mp_size_t>(b2w(nb)));
                }
            }
        }
    }
}/*}}}*/

matpoly<true> matpoly<true>::mp(matpoly const & a, matpoly const & c)/*{{{*/
{
    size_t const nb = std::max(a.size, c.size) - std::min(a.size, c.size) + 1;
    ASSERT_ALWAYS(a.n == c.m);
    matpoly b(a.ab, a.m, c.n, nb);
    b.zero();
    b.addmp(a, c);
    return b;
}/*}}}*/

int matpoly<true>::coeff_is_zero(size_t k) const
{
    size_t const kq = k / ULONG_BITS;
    size_t const kbit = 1UL << (k % ULONG_BITS);
    for(unsigned int i = 0 ; i < m ; i++)
        for(unsigned int j = 0; j < n; j++)
            if (part(i, j)[kq] & kbit)
                return 0;
    return 1;
}
void matpoly<true>::coeff_set_zero(size_t k)
{
    size_t const kq = k / ULONG_BITS;
    size_t const kbit = 1UL << (k % ULONG_BITS);
    for(unsigned int i = 0 ; i < m ; i++)
        for(unsigned int j = 0; j < n; j++)
            part(i, j)[kq] &= ~kbit;
}

matpoly<true> matpoly<true>::truncate_and_rshift(size_t truncated_size, size_t shiftcount)
{
    matpoly other(ab, m, n, size - shiftcount);
    other.rshift(*this, shiftcount);
    truncate(*this, truncated_size);
    shrink_to_fit();
    std::swap(*this, other);
    return other;
}

// NOLINTEND(readability-static-accessed-through-instance)
