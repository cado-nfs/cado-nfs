#include "cado.h"
#include <cstring>
#include "bblas_mat64.hpp"
#include "bblas_level5.hpp"
#include "bblas_level3a.hpp"  // for mat64_add
#include "bblas_level3b.hpp"  // for mul_6464_6464
#include "memory.h"      // malloc_aligned
#include "macros.h"                      // for ASSERT, ASSERT_ALWAYS

/* We reach the mpfq sources through the gf2x code base, and then these
 * are considered internal cantor-related stuff. We need to include the
 * gf2x config flags before including the mpfq sources.
 *
 * TODO: I don't understand this maze. Do we still need these crazy
 * includes ?
 */
#include "gf2x/gf2x-config-export.h" // IWYU pragma: keep
#include "gf2x/gf2x-impl-export.h" // IWYU pragma: keep
#if ULONG_BITS == 64
#include "mpfq/x86_64/mpfq_2_64.h"
#include "mpfq/x86_64/mpfq_2_128.h"
#elif ULONG_BITS == 32
#include "mpfq/i386/mpfq_2_64.h"
#include "mpfq/i386/mpfq_2_128.h"
#else
#error "neither 32 nor 64 ???"
#endif

/**********************************************************************/
/* level 5: polynomials with 64*64 matrix coefficients
 *      m64pol_add
 *      m64pol_mul
 *      m64pol_addmul
 *      m64pol_addmul_kara
 *      m64pol_mul_kara
 *
 *
 */

/* lengths of a1 and a2 are n1 and n2 */
void m64pol_add(mat64 * r, mat64 const * a1, mat64 const * a2, unsigned int n)/*{{{*/
{
    for(unsigned int i = 0 ; i < n ; i++) {
        mat64_add(r[i], a1[i], a2[i]);
    }
}/*}}}*/

void m64pol_mul(mat64 * r, mat64 const * a1, mat64 const * a2, unsigned int n1, unsigned int n2)/*{{{*/
{
    ASSERT_ALWAYS(r != a1 && r != a2);
    memset((void *) r, 0, (n1 + n2 - 1) * sizeof(mat64));
    m64pol_addmul(r, a1, a2, n1, n2);
}/*}}}*/

void m64pol_addmul(mat64 * r, mat64 const * a1, mat64 const * a2, unsigned int n1, unsigned int n2)/*{{{*/
{
    ASSERT(r != a1 && r != a2);
    for(unsigned int i = 0 ; i < n1 ; i++) {
        for(unsigned int j = 0 ; j < n2 ; j++) {
            mat64 x;
            mul_6464_6464(x, a1[i], a2[j]);
            mat64_add(r[i+j], r[i+j], x);
        }
    }
}/*}}}*/

void m64pol_mul_kara(mat64 * r, mat64 const * a1, mat64 const * a2, unsigned int n1, unsigned int n2)/*{{{*/
{
    ASSERT(r != a1 && r != a2);
    ASSERT(n1 == n2);
    ASSERT((n1 & (n1 - 1)) == 0);
    /* As is certainly not surprising, karatsuba wins as early on as one
     * can imagine */
    if (n1 == 1) {
        m64pol_mul(r, a1, a2, n1, n2);
        return;
    }
    unsigned int h = n1 >> 1;
    memset((void *) r, 0, (n1 + n2 - 1) * sizeof(mat64));

    m64pol_add(r, a1, a1 + h, h);
    m64pol_add(r + 2 * h, a2, a2 + h, h);

    mat64 * t = (mat64 *) malloc_aligned((2*h-1) * sizeof(mat64), 64);
    m64pol_mul_kara(t, r, r + 2 * h, h, h);

    m64pol_mul_kara(r, a1, a2, h, h);
    m64pol_mul_kara(r + 2 * h, a1 + h, a2 + h, h, h);
    m64pol_add(t, t, r, 2 * h - 1);
    m64pol_add(t, t, r + 2 * h, 2 * h - 1);
    m64pol_add(r + h, r + h, t, 2 * h - 1);
    free_aligned(t);
}/*}}}*/

void m64pol_addmul_kara(mat64 * r, mat64 const * a1, mat64 const * a2, unsigned int n1, unsigned int n2)/*{{{*/
{
    mat64 * t = (mat64 *) malloc_aligned((n1 + n2 -1) * sizeof(mat64), 64);
    m64pol_mul_kara(t, a1, a2, n1, n2);
    m64pol_add(r, r, t, n1 + n2 - 1);
    free_aligned(t);
}/*}}}*/

void m64pol_mul_gf2_64_bitslice(mat64 * r, mat64 const * a1, mat64 const * a2)/*{{{*/
{
    unsigned int n1 = 64;
    unsigned int n2 = 64;
    mat64 * t = (mat64 *) malloc_aligned((n1 + n2 -1) * sizeof(mat64), 64);
    m64pol_mul_kara(t, a1, a2, n1, n2);
    /* This reduces modulo the polynomial x^64+x^4+x^3+x+1 */
    for(unsigned int i = n1 + n2 - 2 ; i >= 64 ; i--) {
        mat64_add(t[i-64+4], t[i-64+4], t[i]);
        mat64_add(t[i-64+3], t[i-64+3], t[i]);
        mat64_add(t[i-64+1], t[i-64+1], t[i]);
        mat64_add(t[i-64  ], t[i-64  ], t[i]);
    }
    memcpy((void *) r, t, 64 * sizeof(mat64));
    free_aligned(t);
}/*}}}*/

void m64pol_scalmul_gf2_64_bitslice(mat64 * r, mat64 const * a, uint64_t * s)/*{{{*/
{
    unsigned int n1 = 64;
    unsigned int n2 = 64;
    mat64 * t = (mat64 *) malloc_aligned((n1 + n2 -1) * sizeof(mat64), 64);
    memset((void *) t, 0, (n1 + n2 -1) * sizeof(mat64));
    for(unsigned int i = 0 ; i < 64 ; i++) {
        if (((s[0]>>i)&UINT64_C(1))==0) continue;
        m64pol_add(t+i, t+i, a, 64);
        // for(unsigned int j = 0 ; j < 64 ; j++) { mat64_add(t[i+j], t[i+j], a[j]); }
    }
    /* This reduces modulo the polynomial x^64+x^4+x^3+x+1 */
    for(unsigned int i = n1 + n2 - 2 ; i >= 64 ; i--) {
        mat64_add(t[i-64+4], t[i-64+4], t[i]);
        mat64_add(t[i-64+3], t[i-64+3], t[i]);
        mat64_add(t[i-64+1], t[i-64+1], t[i]);
        mat64_add(t[i-64  ], t[i-64  ], t[i]);
    }
    memcpy((void *) r, t, 64 * sizeof(mat64));
    free_aligned(t);
}/*}}}*/

void m64pol_scalmul_gf2_64_bitslice2(mat64 * r, mat64 const * a, uint64_t * s)/*{{{*/
{
    /* Now try with precomputation of multiples. We'll do only four of
     * them to start with. */
    unsigned int n1 = 64;
    unsigned int n2 = 64;
    mat64 * t = (mat64 *) malloc_aligned((n1 + n2 -1) * sizeof(mat64), 64);
    memset((void *) t, 0, (n1 + n2 -1) * sizeof(mat64));

    /* Precompute multiples of a */
    /* The best value for NMULTS depends of course on the cache size. 2
     * and 4 ain't bad choices. Support for non-power-of-2 NMULTS should
     * be looked at (since presently 4 wins over 2).
     */
#define NMULTS  2
    mat64 * am_area = (mat64 *) malloc_aligned((1 << NMULTS) * (64 + NMULTS - 1) * sizeof(mat64), 64);
    mat64 * am[1 << NMULTS];

    for(unsigned int i = 0 ; i < (1 << NMULTS) ; i++) {
        am[i] = am_area + i * (64 + NMULTS - 1);
    }

    memset((void *) am_area, 0, (1 << NMULTS) * (64 + NMULTS - 1) * sizeof(mat64));
    memcpy((void *) am[1], a, 64 * sizeof(mat64));
    for(unsigned int j = 1 ; j < NMULTS ; j++) {
        /* Duplicate all stuff having msb set from level below */
        for(unsigned int i = (1u << (j-1)) ; i < (1u << j) ; i++) {
            memcpy((void *) (am[(i<<1)] + 1), am[i], (64 + j - 1) * sizeof(mat64));
            m64pol_add(am[(i<<1)+1], am[(i<<1)+1], am[1], 64);
        }
    }

    uint64_t v = *s;
    for(unsigned int i = 0 ; i < 64 ; i+=NMULTS, v>>=NMULTS) {
        m64pol_add(t + i, t + i, am[v & ((1<<NMULTS)-1)], 64+(NMULTS-1));
    }
    free_aligned(am_area);

    /* This reduces modulo the polynomial x^64+x^4+x^3+x+1 */
    for(unsigned int i = n1 + n2 - 2 ; i >= 64 ; i--) {
        mat64_add(t[i-64+4], t[i-64+4], t[i]);
        mat64_add(t[i-64+3], t[i-64+3], t[i]);
        mat64_add(t[i-64+1], t[i-64+1], t[i]);
        mat64_add(t[i-64  ], t[i-64  ], t[i]);
    }
    memcpy((void *) r, t, 64 * sizeof(mat64));
    free_aligned(t);
}/*}}}*/

void m64pol_mul_gf2_64_nobitslice(uint64_t * r, uint64_t * a1, uint64_t * a2)/*{{{*/
{
    /* WARNING: We are not considering the data in the same order as the
     * function above */
    for(unsigned int i = 0 ; i < 64 ; i++) {
        uint64_t * a1row = a1 + 64 * i;
        uint64_t * rrow = r + 64 * i;
        for(unsigned int j = 0 ; j < 64 ; j++) {
            uint64_t * a2col = a2 + j;
            uint64_t dst[2] = {0,};
            uint64_t sdst[2] = {0,};
            for(unsigned int k = 0 ; k < 64 ; k++) {
                mpfq_2_64_mul_ur(0, (unsigned long *) dst, (unsigned long *) (a1row + k), (unsigned long *) (a2col + 64*k));
                mpfq_2_64_elt_ur_add(0, (unsigned long *) sdst, (unsigned long *) sdst, (unsigned long *) dst);
            }
            mpfq_2_64_reduce(0, (unsigned long *) (rrow + j), (unsigned long *) sdst);
        }
    }
}/*}}}*/

void m64pol_scalmul_gf2_64_nobitslice(uint64_t * r, uint64_t * a, uint64_t * scalar)/*{{{*/
{
    /* WARNING: We are not considering the data in the same order as the
     * function above */
    for(unsigned int i = 0 ; i < 64 ; i++) {
        uint64_t * arow = a + 64 * i;
        uint64_t * rrow = r + 64 * i;
        for(unsigned int j = 0 ; j < 64 ; j++) {
            mpfq_2_64_mul(0, (unsigned long*) (rrow+j), (unsigned long*) (arow + j), (unsigned long*) scalar);
        }
    }
}/*}}}*/

void m64pol_mul_gf2_128_bitslice(mat64 * r, mat64 const * a1, mat64 const * a2)/*{{{*/
{
    unsigned int n1 = 128;
    unsigned int n2 = 128;
    mat64 * t = (mat64 *) malloc_aligned((n1 + n2 -1) * sizeof(mat64), 64);
    m64pol_mul_kara(t, a1, a2, n1, n2);
    /* This reduces modulo the polynomial x^128+x^7+x^2+x+1 */
    for(unsigned int i = n1 + n2 - 2 ; i >= 128; i--) {
        mat64_add(t[i-128+7], t[i-128+7], t[i]);
        mat64_add(t[i-128+2], t[i-128+2], t[i]);
        mat64_add(t[i-128+1], t[i-128+1], t[i]);
        mat64_add(t[i-128  ], t[i-128  ], t[i]);
    }
    memcpy((void *) r, t, 128 * sizeof(mat64));
    free_aligned(t);
}/*}}}*/

void m64pol_scalmul_gf2_128_bitslice(mat64 * r, mat64 const * a, uint64_t * s)/*{{{*/
{
    unsigned int n1 = 128;
    unsigned int n2 = 128;
    mat64 * t = (mat64 *) malloc_aligned((n1 + n2 -1) * sizeof(mat64), 64);
    memset((void *) t, 0, (n1 + n2 -1) * sizeof(mat64));
    for(unsigned int i = 0 ; i < 128 ; i++) {
        if (((s[i/64]>>(i&63))&UINT64_C(1))==0) continue;
        for(unsigned int j = 0 ; j < 64 ; j++) {
            mat64_add(t[i+j], t[i+j], a[j]);
        }
    }
    /* This reduces modulo the polynomial x^128+x^7+x^2+x+1 */
    for(unsigned int i = n1 + n2 - 2 ; i >= 128; i--) {
        mat64_add(t[i-128+7], t[i-128+7], t[i]);
        mat64_add(t[i-128+2], t[i-128+2], t[i]);
        mat64_add(t[i-128+1], t[i-128+1], t[i]);
        mat64_add(t[i-128  ], t[i-128  ], t[i]);
    }
    memcpy((void *) r, t, 128 * sizeof(mat64));
    free_aligned(t);
}/*}}}*/

void m64pol_mul_gf2_128_nobitslice(uint64_t * r, uint64_t * a1, uint64_t * a2)/*{{{*/
{
    /* WARNING: We are not considering the data in the same order as the
     * function above */
    for(unsigned int i = 0 ; i < 64 ; i++) {
        uint64_t * a1row = a1 + 64 * i;
        uint64_t * rrow = r + 64 * i;
        for(unsigned int j = 0 ; j < 64 ; j++) {
            uint64_t * a2col = a2 + j;
            uint64_t dst[4] = {0,};
            uint64_t sdst[4] = {0,};
            for(unsigned int k = 0 ; k < 64 ; k++) {
                mpfq_2_128_mul_ur(0, (unsigned long *) dst, (unsigned long *) (a1row + k), (unsigned long *) (a2col + 64*k));
                mpfq_2_128_elt_ur_add(0, (unsigned long *) sdst, (unsigned long *) sdst, (unsigned long *) dst);
            }
            mpfq_2_128_reduce(0, (unsigned long *) (rrow + j), (unsigned long *) sdst);
        }
    }
}/*}}}*/

void m64pol_scalmul_gf2_128_nobitslice(uint64_t * r, uint64_t * a, uint64_t * scalar)/*{{{*/
{
    /* WARNING: We are not considering the data in the same order as the
     * function above */
    for(unsigned int i = 0 ; i < 64 ; i++) {
        uint64_t * arow = a + 64 * i;
        uint64_t * rrow = r + 64 * i;
        for(unsigned int j = 0 ; j < 64 ; j++) {
            mpfq_2_128_mul(0, (unsigned long *) (rrow+j), (unsigned long *) (arow + j), (unsigned long *) scalar);
        }
    }
}/*}}}*/

void m64polblock_mul(mat64 * r, mat64 const * a1, mat64 const * a2, unsigned int n1, unsigned int n2, unsigned int K)/*{{{*/
{
    /* Same spirit, but treat multiplication of 64K by 64K matrices (of
     * polynomials).
     *
     * That is, we consider a K*K matrix of polynomials of 64*64
     * matrices.
     * 
     * We assume that there is no pointer aliasing, and that matrices are
     * stored row-major, with all polynomials contiguous (and of fixed
     * lengths n1, resp n2).
     */
    ASSERT(r != a1 && r != a2);
    memset((void *) r, 0, (n1 + n2 - 1) * K * K * sizeof(mat64));
    for(unsigned int i = 0 ; i < K ; i++) {
        mat64 const * ra1 = a1 + i * n1 * K;
        mat64 * rr = r + i * (n1 + n2 - 1) * K;
        for(unsigned int j = 0 ; j < K ; j++) {
            mat64 const * ca2 = a2 + j * n2;
            mat64 * pr = rr + j * (n1 + n2 - 1);
            for(unsigned int k = 0 ; k < K ; k++) {
                mat64 const * pa1 = ra1 + k * n1;
                mat64 const * pa2 = ca2 + k * n2 * K;
                m64pol_addmul(pr, pa1, pa2, n1, n2);
            }
        }
    }
}/*}}}*/

void m64polblock_mul_kara(mat64 * r, mat64 const * a1, mat64 const * a2, unsigned int n1, unsigned int n2, unsigned int K)/*{{{*/
{
    ASSERT(r != a1 && r != a2);
    memset((void *) r, 0, (n1 + n2 - 1) * K * K * sizeof(mat64));
    for(unsigned int i = 0 ; i < K ; i++) {
        mat64 const * ra1 = a1 + i * n1 * K;
        mat64 * rr = r + i * (n1 + n2 - 1) * K;
        for(unsigned int j = 0 ; j < K ; j++) {
            mat64 const * ca2 = a2 + j * n2;
            mat64 * pr = rr + j * (n1 + n2 - 1);
            for(unsigned int k = 0 ; k < K ; k++) {
                mat64 const * pa1 = ra1 + k * n1;
                mat64 const * pa2 = ca2 + k * n2 * K;
                m64pol_addmul_kara(pr, pa1, pa2, n1, n2);
            }
        }
    }
}/*}}}*/

