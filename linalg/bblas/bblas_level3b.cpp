#include "cado.h" // IWYU pragma: keep
#include <cstring>
#include <cstdint>                        // for uint64_t, UINT64_C
#ifdef HAVE_SSE2
#include <emmintrin.h>                     // for __m128i, _mm_and_si128
#endif
#include "bblas_level3b.hpp"
#include "bblas_level3c.hpp"  // for mul_6464lt_6464_lookup4
#include "bblas_mat64.hpp"    // for mat64
#include "bblas_simd.hpp"
#include "macros.h"                        // for MAYBE_UNUSED

/**********************************************************************/
/* level 3b: matrix multiplications
 *      mul_6464_6464   (several variants)
 *      addmul_6464_6464_fragment       TODO: missing test
 */

/* implemented here:
 *    - mul_6464_6464 ; self-explanatory
 *    - addmul_6464_6464_fragment: takes only a subset of rows of the
 *      first input, and a subset of columns of the second input.
 */

#if defined(HAVE_SSE2) && ULONG_BITS == 64
void mul_6464_6464_sse(mat64 & C, mat64 const & A, mat64 const & B)/*{{{*/
{
    int i;
    C = 0;
 
    __m128i *Cw = (__m128i *) C.data();
    __m128i *Aw = (__m128i *) A.data();

    for (int j = 0; j < 64; j += 2) {
	__m128i c = _mm_setzero_si128();
	__m128i a = *Aw++;

	__m128i one = _cado_mm_set1_epi64_c(1);
	for (i = 0; i < 64; i++) {
	    __m128i bw = _cado_mm_set1_epi64(B[i]);
	    // c ^= (bw & -(a & one));
            c = _mm_xor_si128(c, _mm_and_si128(bw, _mm_sub_epi64(_mm_setzero_si128(), _mm_and_si128(a, one))));
	    a = _mm_srli_epi64(a, 1);
	}
	*Cw++ = c;
    }
}/*}}}*/
#endif

void mul_6464_6464_v2(mat64 & C, mat64 const & A, mat64 const & B)/*{{{*/
{
    C = 0;

    for (int i = 0; i < 64; i++) {
        for (int j = 0; j < 64; j++) {
            if ((A[i]>>j)&1)
                C[i]^=B[j];
        }
    }
}/*}}}*/

void addmul_6464_6464_fragment_lookup4(mat64 & C,/*{{{*/
                   mat64 const & A,
                   mat64 const & B,
                   unsigned int i0,
                   unsigned int i1,
                   unsigned int yi0,
                   unsigned int yi1)
{
    uint64_t Bx[16][16];
    unsigned int j0 = yi0 / 4;
    unsigned int j1 = (yi1 + 3) / 4;
    for(unsigned int j = j0 ; j < j1 ; j++) {
        const uint64_t * bb = B.data() + 4 * j;
        uint64_t w = 0;
        Bx[j][0]  = w; w ^= bb[0];
        Bx[j][1]  = w; w ^= bb[1];
        Bx[j][3]  = w; w ^= bb[0];
        Bx[j][2]  = w; w ^= bb[2];
        Bx[j][6]  = w; w ^= bb[0];
        Bx[j][7]  = w; w ^= bb[1];
        Bx[j][5]  = w; w ^= bb[0];
        Bx[j][4]  = w; w ^= bb[3];
        Bx[j][12] = w; w ^= bb[0];
        Bx[j][13] = w; w ^= bb[1];
        Bx[j][15] = w; w ^= bb[0];
        Bx[j][14] = w; w ^= bb[2];
        Bx[j][10] = w; w ^= bb[0];
        Bx[j][11] = w; w ^= bb[1];
        Bx[j][9]  = w; w ^= bb[0];
        Bx[j][8]  = w;
    }
    uint64_t mask = (UINT64_C(1) << yi1) - (UINT64_C(1) << yi0);
    if (yi1 == 64)
        mask = - (UINT64_C(1) << yi0);
    for (size_t i = i0; i < i1; i++) {
        uint64_t aa = (A[i] & mask) >> (4 * j0);
        for(unsigned int j = j0 ; j < j1 ; j++) {
            C[i]^= Bx[j][aa & 15]; aa>>=4;
        }
    }
}/*}}}*/

/* {{{ final choices. These are static choices at this point, but it should
 * be the result of some tuning, ideally */
void mul_6464_6464(mat64 & C, mat64 const & A, mat64 const & B)
{
    mul_N64_6464_lookup4(C.data(), A.data(), B, 64);
}
void mul_6464lt_6464(mat64 & C, mat64 const & A, mat64 const & B)
{
    mul_6464lt_6464_lookup4(C.data(), A.data(), B);
}
/*}}}*/

void MAYBE_UNUSED addmul_6464_6464(mat64 & C,/*{{{*/
                   mat64 const & A,
                   mat64 const & B)
{
    addmul_6464_6464_fragment_lookup4(C, A, B, 0, 64, 0, 64);
}
/*}}}*/
void MAYBE_UNUSED addmul_6464_6464_fragment(mat64 & C,/*{{{*/
                   mat64 const & A,
                   mat64 const & B,
                   unsigned int i0,
                   unsigned int i1,
                   unsigned int yi0,
                   unsigned int yi1)
{
    addmul_6464_6464_fragment_lookup4(C, A, B, i0, i1, yi0, yi1);
}
/*}}}*/

