#include "cado.h"
#include <cstring>
#include "bblas_mat64.hpp"
#include "bblas_level5.hpp"
#include "bblas_level3a.hpp"  // for mat64_add
#include "bblas_level3b.hpp"  // for mul_6464_6464
#include "memory.h"      // malloc_aligned
#include "macros.h"                      // for ASSERT, ASSERT_ALWAYS

// #define WLEN ULONG_BITS
// #define GF2X_WORDSIZE ULONG_BITS
#define GF2X_MAYBE_UNUSED MAYBE_UNUSED
#define CANTOR_BASE_FIELD_SIZE 128
#include "gf2x-cantor-field-impl.h"


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
                Kmul_ur((unsigned long *) dst, (unsigned long *) (a1row + k), (unsigned long *) (a2col + 64*k));
                Kelt_ur_add((unsigned long *) sdst, (unsigned long *) sdst, (unsigned long *) dst);
            }
            Kreduce((unsigned long *) (rrow + j), (unsigned long *) sdst);
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
            Kmul((unsigned long *) (rrow+j), (unsigned long *) (arow + j), (unsigned long *) scalar);
        }
    }
}/*}}}*/

