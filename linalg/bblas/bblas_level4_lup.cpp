#include "cado.h"
#include <cstdint>          // for uint64_t, UINT64_C
#include "bblas_mat64.hpp"  // for mat64
#include "bblas_level4.hpp"
#include "bblas_simd.hpp"
// the whole point of bblas_simd is to avoid including these files...
// IWYU pragma: no_include <emmintrin.h>
// IWYU pragma: no_include <smmintrin.h>
// IWYU pragma: no_include <immintrin.h>

/* Computes l,u,p, such that:
 *  - l is unit lower triangular
 *  - l*a=u
 *  - up=u*transpose(p) is upper triangular,
 *  -   up[k]==0 iff up[k,k] == 0
 *  -   rank(up)=rank(a)=# of diagonal 1's in up.
 */
int LUP64_imm(mat64 & l, mat64 & u, mat64 & p, mat64 const & a)
{
    u = a;
    uint64_t mask=1;
    uint64_t todo=~UINT64_C(0);
    int r = 0;
    for(int j = 0 ; j < 64 ; j++, mask<<=1) p[j]=l[j]=mask;
    mask=1;
    int store[64];
    int * ps = store;
    for(int j = 0 ; j < 64 ; j++, mask<<=1) {
#if 0
            mat64 t;
            mul_6464_6464(t,l,a); ASSERT_ALWAYS(mat64_eq(t,u));
#endif
        uint64_t pr=u[j];
        // ASSERT(!(pr&~todo));
        if (!(pr&todo)) { *ps++=j; continue; }
        // this keeps only the least significant bit of pr.
        uint64_t v = pr^(pr&(pr-1));
        p[j]=v;
#if defined(HAVE_SSE41) && !defined(VALGRIND) && !defined(__ICC)
        /* This code is fine, except that there's nowhere we convey the
         * information that we want mat64 objects 16-byte aligned. And
         * with icc, this indeed fails (for test_bitlinalg_matops). We
         * don't know exactly what to do here. 
         */
        int k = j+1;
        if (k&1) {      // alignment call
            uint64_t w = -((u[k]&v)!=0);
            u[k]^=pr&w;
            l[k]^=l[j]&w;
            k++;
        }
        /* ok, it's ugly, and requires sse 4.1.
         * but otoh is churns out data veeery fast */
        __m128i vv = _cado_mm_set1_epi64(v);
        __m128i pp = _cado_mm_set1_epi64(pr);
        __m128i ee = _cado_mm_set1_epi64(l[j]);
        __m128i * uu = (__m128i*) (u.data()+k);
        __m128i * ll = (__m128i*) (l.data()+k);
        for( ; k < 64 ; k+=2 ) {
            __m128i ww = _mm_cmpeq_epi64(_mm_and_si128(*uu,vv),vv);
            *uu = _mm_xor_si128(*uu, _mm_and_si128(pp, ww));
            *ll = _mm_xor_si128(*ll, _mm_and_si128(ee, ww));
            uu++;
            ll++;
        }
#else
        uint64_t er = l[j];
        for(int k = j+1 ; k<64 ; k++) {
            uint64_t w = -((u[k]&v)!=0);
            u[k]^=pr&w;
            l[k]^=er&w;
        }
#endif
        todo^=v;
        r++;
    }
    for(ps = store ; todo ; ) {
        uint64_t vv = todo^(todo&(todo-1));
        p[*ps++] = vv;
        todo^=vv;
    }
    return r;
}

