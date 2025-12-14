#include "cado.h" // IWYU pragma: keep

#include <cstdint>
#include <cstring>

#include "bblas_mat64.hpp"
#include "bblas_level5.hpp"
#include "bblas_level3a.hpp"  // for mat64_add
#include "memory.h"      // malloc_aligned

#define GF2X_MAYBE_UNUSED MAYBE_UNUSED
#define GF2X_CANTOR_BASE_FIELD_SIZE 64

#include "gf2x-cantor-field-impl.h"

void m64pol_mul_gf2_64_bitslice(mat64 * r, mat64 const * a1, mat64 const * a2)/*{{{*/
{
    unsigned int const n1 = 64;
    unsigned int const n2 = 64;
    auto * t = (mat64 *) malloc_aligned((n1 + n2 -1) * sizeof(mat64), 64);
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
    unsigned int const n1 = 64;
    unsigned int const n2 = 64;
    auto * t = (mat64 *) malloc_aligned((n1 + n2 -1) * sizeof(mat64), 64);
    memset((void *) t, 0, (n1 + n2 -1) * sizeof(mat64));
    for(unsigned int i = 0 ; i < 64 ; i++) {
        if (((s[0]>>i)&uint64_t(1))==0) continue;
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
    unsigned int const n1 = 64;
    unsigned int const n2 = 64;
    auto * t = (mat64 *) malloc_aligned((n1 + n2 -1) * sizeof(mat64), 64);
    memset((void *) t, 0, (n1 + n2 -1) * sizeof(mat64));

    /* Precompute multiples of a */
    /* The best value for NMULTS depends of course on the cache size. 2
     * and 4 ain't bad choices. Support for non-power-of-2 NMULTS should
     * be looked at (since presently 4 wins over 2).
     */
#define NMULTS  2
    auto * am_area = (mat64 *) malloc_aligned((1 << NMULTS) * (64 + NMULTS - 1) * sizeof(mat64), 64);
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
            unsigned long dst[2 * sizeof(uint64_t) / sizeof(unsigned long)] = {0,};
            unsigned long sdst[2 * sizeof(uint64_t) / sizeof(unsigned long)] = {0,};
            for(unsigned int k = 0 ; k < 64 ; k++) {
                Kmul_ur(dst, (unsigned long *) (a1row + k), (unsigned long *) (a2col + 64*k));
                Kelt_ur_add(sdst, sdst, dst);
            }
            Kreduce((unsigned long *) (rrow + j), sdst);
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
            Kmul((unsigned long*) (rrow+j), (unsigned long*) (arow + j), (unsigned long*) scalar);
        }
    }
}/*}}}*/


