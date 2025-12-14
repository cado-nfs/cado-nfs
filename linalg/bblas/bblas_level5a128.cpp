#include "cado.h" // IWYU pragma: keep

#include <cstring>
#include <cstdint>

#include <memory>

#include "bblas_mat64.hpp"
#include "bblas_level5.hpp"
#include "bblas_level3a.hpp"  // for mat64_add
#include "memory.h"      // malloc_aligned
#include "macros.h"

#define GF2X_MAYBE_UNUSED MAYBE_UNUSED
#define GF2X_CANTOR_BASE_FIELD_SIZE 128

#include "gf2x-cantor-field-impl.h"

struct free_aligned_obj
{
    void operator()(void* x) { free_aligned(x); }
};

#if GNUC_VERSION_ATLEAST(6,1,0)
/* just because we're doing a unique_ptr on a type with an attribute.
 * Sigh.
 */
/* https://gcc.gnu.org/bugzilla/show_bug.cgi?id=69884 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
#endif

void m64pol_mul_gf2_128_bitslice(mat64 * r, mat64 const * a1, mat64 const * a2)/*{{{*/
{
    unsigned int const n1 = 128;
    unsigned int const n2 = 128;
    const std::unique_ptr<mat64, free_aligned_obj> t((mat64 *) malloc_aligned((n1 + n2 -1) * sizeof(mat64), 64));
    m64pol_mul_kara(t.get(), a1, a2, n1, n2);
    /* This reduces modulo the polynomial x^128+x^7+x^2+x+1 */
    for(unsigned int i = n1 + n2 - 2 ; i >= 128; i--) {
        mat64_add(t.get()[i-128+7], t.get()[i-128+7], t.get()[i]);
        mat64_add(t.get()[i-128+2], t.get()[i-128+2], t.get()[i]);
        mat64_add(t.get()[i-128+1], t.get()[i-128+1], t.get()[i]);
        mat64_add(t.get()[i-128  ], t.get()[i-128  ], t.get()[i]);
    }
    memcpy((void *) r, t.get(), 128 * sizeof(mat64));
}/*}}}*/

void m64pol_scalmul_gf2_128_bitslice(mat64 * r, mat64 const * a, uint64_t const * s)/*{{{*/
{
    unsigned int const n1 = 128;
    unsigned int const n2 = 128;
    const std::unique_ptr<mat64, free_aligned_obj> t((mat64 *) malloc_aligned((n1 + n2 -1) * sizeof(mat64), 64));
    memset((void *) t.get(), 0, (n1 + n2 -1) * sizeof(mat64));
    for(unsigned int i = 0 ; i < 128 ; i++) {
        if (((s[i/64]>>(i&63)) & uint64_t(1))==0) continue;
        for(unsigned int j = 0 ; j < 64 ; j++) {
            mat64_add(t.get()[i+j], t.get()[i+j], a[j]);
        }
    }
    /* This reduces modulo the polynomial x^128+x^7+x^2+x+1 */
    for(unsigned int i = n1 + n2 - 2 ; i >= 128; i--) {
        mat64_add(t.get()[i-128+7], t.get()[i-128+7], t.get()[i]);
        mat64_add(t.get()[i-128+2], t.get()[i-128+2], t.get()[i]);
        mat64_add(t.get()[i-128+1], t.get()[i-128+1], t.get()[i]);
        mat64_add(t.get()[i-128  ], t.get()[i-128  ], t.get()[i]);
    }
    memcpy((void *) r, t.get(), 128 * sizeof(mat64));
}/*}}}*/

#if GNUC_VERSION_ATLEAST(6,1,0)
#pragma GCC diagnostic pop
#endif

void m64pol_mul_gf2_128_nobitslice(uint64_t * r, uint64_t * a1, uint64_t * a2)/*{{{*/
{
    /* WARNING: We are not considering the data in the same order as the
     * function above */
    for(unsigned int i = 0 ; i < 64 ; i++) {
        uint64_t * a1row = a1 + 64 * i;
        uint64_t * rrow = r + 64 * i;
        for(unsigned int j = 0 ; j < 64 ; j++) {
            uint64_t * a2col = a2 + j;
            unsigned long dst[4 * sizeof(uint64_t) / sizeof(unsigned long)] = {0,};
            unsigned long sdst[4 * sizeof(uint64_t) / sizeof(unsigned long)] = {0,};
            for(unsigned int k = 0 ; k < 64 ; k++) {
                Kmul_ur(dst, (unsigned long *) (a1row + k), (unsigned long *) (a2col + 64*k));
                Kelt_ur_add(sdst, sdst, dst);
            }
            Kreduce((unsigned long *) (rrow + j), sdst);
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

