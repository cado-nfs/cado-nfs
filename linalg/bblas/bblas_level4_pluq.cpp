#include "cado.h" // IWYU pragma: keep

#include <cstdint>
#include <cstring>

#include <algorithm>
#include <utility>

#include "bblas_level4.hpp"
#include "bblas_perm_matrix.hpp"
#include "bblas_simd.hpp"
#include "linalg/bblas/bblas_bitmat.hpp"   // for bitmat<>::vector_type
#include "linalg/bblas/bblas_level3a.hpp"  // for mat64_add
#include "linalg/bblas/bblas_level3b.hpp"  // for mul_6464_6464
#include "linalg/bblas/bblas_mat64.hpp"    // for mat64
#include "macros.h"                        // for ATTRIBUTE, ASSERT, ASSERT_...
#include "misc.h"      // cado_ctz64

// IWYU pragma: no_include <mmintrin.h>
// IWYU pragma: no_include <emmintrin.h>
// IWYU pragma: no_include <smmintrin.h>

/*  PLUQ stuff -- well we're not computing exactly PLUQ */

/* Compute matrices l and u such that l*a = u. phi[] is filled with the
 * column indices (shifted by col_offset) of he pivots in u: entry
 * (i,phi(i)-col_offset) is u is one. When row i in u has no important
 * non-zero coefficient, then phi[i] < 0.
 * In column phi[i]-col_offset of u, entries of row index >i are zero.
 */
static int PLUQ64_inner(int * phi, mat64 & l, mat64 & u, mat64 const & a, int col_offset)
{
    const int m = 64;
    const int n = 64;
    int phi0[64];
    if (phi == nullptr) {
        phi = phi0;
        for(int i = 0 ; i < 64 ; i++) phi[i]=-1;
    }

    u = a;
    l = 1;
    int rank = 0;
    uint64_t todo=~((uint64_t)0);
    for(int i = 0 ; i < m ; i++) {
        uint64_t const r=u[i];
        if (phi[i]>=0) continue;
        if (!(r&todo) || phi[i]>=0) continue;
        // this keeps only the least significant bit of r.
        uint64_t const v = r^(r&(r-1));
        uint64_t const j = cado_ctz64(r);
        phi[i] = col_offset + j;
#if defined(HAVE_SSE41) && !defined(VALGRIND) && !defined(__ICC)
        /* This code is fine, except that there's nowhere we convey the
         * information that we want mat64 objects 16-byte aligned. And
         * with icc, this indeed fails (for test_bitlinalg_matops). We
         * don't know exactly what to do here. 
         */
        int k = i+1;
        if (k&1) {      // alignment call
            uint64_t const w = -((u[k]&v)!=0);
            u[k]^=r&w;
            l[k]^=l[i]&w;
            k++;
        }
        /* ok, it's ugly, and requires sse 4.1.
         * but otoh is churns out data veeery fast */
        __m128i const vv = _cado_mm_set1_epi64(v);
        __m128i const pp = _cado_mm_set1_epi64(r);
        __m128i const ee = _cado_mm_set1_epi64(l[i]);
        auto * uu = (__m128i*) (u.data() + k);
        auto * ll = (__m128i*) (l.data() + k);
        for( ; k < n ; k+=2 ) {
            __m128i const ww = _mm_cmpeq_epi64(_mm_and_si128(*uu,vv),vv);
            *uu = _mm_xor_si128(*uu, _mm_and_si128(pp, ww));
            *ll = _mm_xor_si128(*ll, _mm_and_si128(ee, ww));
            uu++;
            ll++;
        }
#else
        uint64_t er = l[i];
        for(int k = i+1 ; k<n ; k++) {
            uint64_t w = -((u[k]&v)!=0);
            u[k]^=r&w;
            l[k]^=er&w;
        }
#endif
        todo^=v;
        rank++;
    }
#if defined(HAVE_SSE41) && !defined(VALGRIND) && !defined(__ICC)
    _mm_empty();
#endif
    return rank;
}

/* PLUQ -- well we're not computing exactly PLUQ 
 * PLUQ says: Any m*n matrix A with rank r , can be written A = P*L*U*Q
 * where P and Q are two permutation matrices, of dimension respectively
 * m*m and n*n, L is m*r unit lower triangular and U is r*n upper
 * triangular.
 *
 * Here we compute p,l,u,q such that p*l*a*transpose(q) = an upper
 * triangular matrix, whose diagonal has r one and n-r zeros.
 */
/* outer routine */
int PLUQ64(perm_matrix_ptr p, mat64 & l, mat64 & u, perm_matrix_ptr q, mat64 const & m)
{
    int phi[64];
    for(int i = 0 ; i < 64 ; i++) phi[i]=-1;
    int const r = PLUQ64_inner(phi,l,u,m,0);
    /* l*m = u */
    /* p*u*transpose(q) = diagonal.
     * p*l*m*transpose(q) = diagonal.
     */
    pqperms_from_phi(p,q,phi,64,64);
    return r;
}

int PLUQ64_n(int * phi, mat64 & l, mat64 * u, mat64 const * a, int n)
{
    const int m = 64;
    ASSERT_ALWAYS(n % 64 == 0);
    int const nb = n/m;
    l = 1;
    for(int i = 0 ; i < 64 ; i++) phi[i]=-1;
    int rank = 0;
    int b = 0;
    mat64::vector_type ls(nb);
    mat64 tl;
    for( ; b < nb && rank < m ; b++) {
        mat64 ta;
        mul_6464_6464(ta, l, a[b]);
        rank += PLUQ64_inner(phi, tl, u[b], ta, b*m);
        mul_6464_6464(l, tl, l);
        ls[b] = tl;
    }
    int const nspins = b;
    for(int c = b-2 ; c >= 0 ; c--) {
        mul_6464_6464(u[c], tl, u[c]);
        mul_6464_6464(tl, ls[c], tl);
    }
    for( ; b < nb ; b++)
        mul_6464_6464(u[b], l, a[b]);
    return nspins*m+b;
}

static inline void bli_64x64N_clobber(mat64 & h, mat64 * us, int const * phi, int nb)
{
    /* problem: we're modifying U here. So either we do a copy of U,
     * which can be probelmatic memory-wise, or we do an extraction ;
     * it's also possible to rebuild the original U from the extracted H
     * and U', merely with a product (_if ever_ we care about U, in
     * fact). However this latter option seems messy.
     */
    h = 1;
    const int m = 64;
    for(int i = 0 ; i < m ; i++) {
        int const j = phi[i];
        if (j<0) continue;
        uint64_t const m = ((uint64_t)1) << (j%64);
        int const d = j/64;
        /* TODO: use _mm_cmpeq_epi64 for this as well, of course */
        ASSERT(us[d][i]&m);
        int k = 0;
#if defined(HAVE_SSE41) && !defined(VALGRIND)
        __m128i const mm = _cado_mm_set1_epi64(m);
        auto * uu = (__m128i*) us[d].data();
        auto * hh = (__m128i*) h.data();
        __m128i const hi = _cado_mm_set1_epi64(h[i]);
        int const ii=i/2;
        for( ; k < ii ; k++) {
            __m128i const ww = _mm_cmpeq_epi64(_mm_and_si128(*uu++,mm),mm);
            for(int b = 0 ; b < nb ; b++) {
                // ((__m128i*)us[b])[k] ^= ww & _cado_mm_set1_epi64(us[b][i]);
                auto * z = ((__m128i*)us[b].data()) + k;
                *z = _mm_xor_si128(*z, _mm_and_si128(ww, _cado_mm_set1_epi64(us[b][i])));
            }
            hh[k] = _mm_xor_si128(hh[k], _mm_and_si128(ww, hi));
        }
        k*=2;
#endif
        for( ; k < i ; k++) {
            uint64_t const w = -((us[d][k]&m) != 0);
            for(int b = 0 ; b < nb ; b++) {
                us[b][k] ^= w & us[b][i];
            }
            h[k] ^= w & h[i];
        }
    }
#if defined(HAVE_SSE41) && !defined(VALGRIND)
    _mm_empty();
#endif
}

/* Given a 64x128 matrix u that is upper triangular up to some
 * permutation, written as a sequence of two 64x64 matrices,
 * and given a table phi such that either phi[i]<0, or entry (i,phi[i])
 * of u is non-zero, and the nonnegative values taken by phi are all
 * distinct, compute a matrix H such that H*U has exactly one non-zero
 * entry in each column whose index is a value taken by phi.
 */
static void bli_64x128(mat64 & h, mat64 * us, int * phi)
{
    mat64 uc[2] ATTRIBUTE((aligned(64)));
    std::ranges::copy(us, us + 2, std::begin(uc));
    bli_64x64N_clobber(h,uc,phi,2);
}

static void extract_cols_64_from_128(mat64 & t, mat64 const * m, int const * phi)
{
    // given the list of 64 integers phi, all in the range {-1} union
    // {0..127}, constitute a 64x64 matrix whose column of index j is
    // column of index phi[j] in the input matrix m. -1 means a zero
    // column.
    uint64_t s[2][64]={{0,},};
    for(int j = 0 ; j < 64 ; j++) {
        if (phi[j]<0) continue;
        s[phi[j]/64][j]=((uint64_t)1)<<(phi[j]%64);
    }
    t = 0;
    uint64_t mask = 1;
    for(int j = 0 ; j < 64 ; j++, mask<<=1) {
#if defined(HAVE_SSE41) && !defined(VALGRIND) && !defined(__ICC)
        /* This code is fine, except that there's nowhere we convey the
         * information that we want mat64 objects 16-byte aligned. And
         * with icc, this indeed fails (for test_bitlinalg_matops). We
         * don't know exactly what to do here. 
         */
        __m128i const ss[2] = {
            _cado_mm_set1_epi64(s[0][j]),
            _cado_mm_set1_epi64(s[1][j]) };
        __m128i * mm[2] = {(__m128i*)m[0].data(),(__m128i*)m[1].data()};
        __m128i * tt = (__m128i*)t.data();
        __m128i const mmk = _cado_mm_set1_epi64(mask);
        for(int i = 0 ; i < 64 ; i+=2) {
            // *tt ^= mmk & _mm_cmpeq_epi64((*mm[0]&ss[0])^(*mm[1]&ss[1]),ss[0]^ss[1]);
            *tt = _mm_xor_si128(*tt, _mm_and_si128(mmk,
                        _mm_cmpeq_epi64(
                            _mm_xor_si128(
                                _mm_and_si128(*mm[0], ss[0]),
                                _mm_and_si128(*mm[1], ss[1])
                                ),
                            _mm_xor_si128(ss[0], ss[1])
                            )
                        )
                    );
            mm[0]++,mm[1]++;
            tt++;
        }
#else
        for(int i = 0 ; i < 64 ; i++) {
            t[i] ^= mask & -(((m[0][i]&s[0][j]) ^ (m[1][i]&s[1][j])) != 0);
        }
#endif
    }
}

    /*
        __m128i vv = (__v2di) { v,v };
        __m128i pp = (__v2di) { r, r };
        __m128i ee = (__v2di) { l[i], l[i] };
    __m128i * tt = (__m128i*) (t);
    __m128i * mm[2] = { (__m128i*) m[0], (__m128i *) m[1] };
    __m128i * ss[2] = { (__m128i*) s[0], (__m128i *) s[1] };
    __v2di mask = (__v2di) {1,1};
    for(int i = 0 ; i < 64 ; i++) {
        *tt = mask & _mm_cmpeq_epi64(*mm&*ss,*ss);
        tt++,mm++,ss++,mask<<=1;
    }
    */

/* This code is here because someday, I vaguely had the idea of using it
 * as a building block for the binary lingen. In fact, the code fragments
 * here for PLUQ and such have never been put in production, so I'm
 * pretty sure they're quite fragile.
 */
int PLUQ128(perm_matrix_ptr p, mat64 * l, mat64 * u, perm_matrix_ptr q, mat64 const * m)
{
    /* This is really an outer routine. An inner routine will not have p
     * and q, but rather both merged as a phi argument, in the manner of
     * PLUQ64_inner (and of course the following lines would be changed a
     * bit).
     */
    int phi[128];
    for(int i = 0 ; i < 128 ; i++) phi[i]=-1;

    l[0] = 0;
    l[1] = 0;
    l[2] = 0;
    l[3] = 1;

    int r1 = PLUQ64_n(phi,l[0],u,m,128);
    r1 = r1 % 64;

    // andouille 7.65

    /* l[0] * m = u */

    mat64 h;
    bli_64x128(h, u, phi);
    /* h * u is "sort of" identity, at least up to permutation */

    mat64 l21;
    mat64 & S = l21;

    /* This is __very__ expensive w.r.t. what it really does :-(( */
    extract_cols_64_from_128(S, m+2, phi);

    /* Column i of S is column phi[i] in Mlow. Now by bli_64x128, in
     * column phi[i] of h*u, only the coefficient of row i is equal to 1,
     * so that column phi[i] of S*H*U is equal to column phi[i] of Mlow
     */

    // andouille 16.7 -- 17.5
    mul_6464_6464(l21, S, h);
    mul_6464_6464(l[2], l21, l[0]);

    // The matrix below has many zero columns (as many as the rank of
    // Mhigh).
    // Mlow+S*H*l[0]*Mhigh;

    /* The matrix [ l[0] 0 ] = L
     *            [ l[2] 1 ]
     * is equal to [ 1   0 ]   [ l[0]  0 ]
     *             [ l21 1 ] * [  0    1 ]
     *
     * Now based on l[0] * mhigh, compute t2 = (L*M)_low
     */
    mat64 t2[2] ATTRIBUTE((aligned(64)));
    mul_6464_6464(t2[0], l21, u[0]); mat64_add(t2[0], m[2], t2[0]);
    mul_6464_6464(t2[1], l21, u[1]); mat64_add(t2[1], m[3], t2[1]);

    /* And do pluq again on this low part. Most of it is zero. */
    int r2 = PLUQ64_n(phi + 64,l[3],u+2,t2,128);
    r2 = r2 % 64;

    /* need to adjust l[3] */
    mul_6464_6464(l[2], l[3], l[2]);

    pqperms_from_phi(p,q,phi,128,128);

    /* At this point P*L*M*Tranpose(Q) should be upper triangular with
     * unit diagonal */
    
    return r1 + r2;
}

/*  */
