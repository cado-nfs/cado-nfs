#include "cado.h"
#include <cstring>
#include "bblas_mat64.hpp"
#include "bblas_level5.hpp"
#include "bblas_level3a.hpp"  // for mat64_add
#include "bblas_level3b.hpp"  // for mul_6464_6464
#include "memory.h"      // malloc_aligned
#include "macros.h"                      // for ASSERT, ASSERT_ALWAYS

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

