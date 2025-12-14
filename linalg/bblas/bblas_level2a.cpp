#include "cado.h" // IWYU pragma: keep

#include <cstdint>

#ifdef HAVE_SSE2
#include <emmintrin.h>                   // for _mm_setzero_si128, __m128i
#include <mmintrin.h>                    // for _mm_empty
#endif
#include "bblas_level2a.hpp"
#include "bblas_mat64.hpp"
#include "bblas_simd.hpp"

/**********************************************************************/
/* level 2a: rank-1 updates.
 *      addmul_To64_o64 (several variants)
 */

/* implemented here:
 *    - addmul_To64_o64 ; add to a 64*64 matrix the rank-1 product
 *      obtained my multiplying two vectors together.
 *
 */

/* compatible implementation for addmul_To64_o64 */
void addmul_To64_o64_lsb(mat64 & r, uint64_t a, uint64_t w)/*{{{*/
{
    /* One way... */
    mat64::datatype * rr = r.data();
    for (unsigned int i = 0; i < 64; i++) {
	*rr++ ^= w & -(a & 1);
	a >>= 1;
    }
}/*}}}*/

/* compatible implementation for addmul_To64_o64 */
void addmul_To64_o64_msb(mat64 & r, uint64_t a, uint64_t w)/*{{{*/
{
    /* ... and the other. weird enough, it seems a bit faster in some
     * cases.
     */
    mat64::datatype * rr = r.data();
    for (unsigned int i = 0; i < 64; i++) {
	*rr++ ^= w & (((int64_t) a) >> 63);
	a <<= 1;
    }
}/*}}}*/

/* compatible implementation for addmul_To64_o64 */
void addmul_To64_o64_lsb_packof2(mat64 & r, uint64_t a, uint64_t w)/*{{{*/
{
    /* À peu près comme la méthode 1, mais pas mieux */
    using mvec_t = uint64_t[2];
    mvec_t mb[4] = {
	{0, 0}, {w, 0}, {0, w}, {w, w},
    };
    mat64::datatype * rr = r.data();
    for (int i = 0; i < 64; i += 2) {
	const uint64_t *y = mb[a & 3];
	*rr++ ^= y[0];
	*rr++ ^= y[1];
	a >>= 2;
    }
}/*}}}*/

#if defined(HAVE_SSE2) && ULONG_BITS == 64
/* compatible implementation for addmul_To64_o64 */
void addmul_To64_o64_lsb_sse_v1(mat64 & r, uint64_t a, uint64_t w)/*{{{*/
{
    /* Using sse-2 */
    __m128i const mb[4] = {
	_mm_setzero_si128(),
	_cado_mm_setr_epi64(w, 0),
	_cado_mm_setr_epi64(0, w),
	_cado_mm_set1_epi64(w),
    };
    auto *sr = (__m128i *) r.data();
    for (int i = 0; i < 64; i += 2) {
	*sr = _mm_xor_si128(*sr, mb[a & 3]);
        sr++;
	a >>= 2;
    }
    _mm_empty();
}/*}}}*/
#endif

/* final choice. This is a static choice at this point, but it should be
 * the result of some tuning, ideally */
void addmul_To64_o64(mat64 & r, uint64_t a, uint64_t w)
{
#if defined(HAVE_SSE2) && ULONG_BITS == 64
    addmul_To64_o64_lsb_sse_v1(r,a,w);
#else
    addmul_To64_o64_lsb_packof2(r,a,w);
#endif
}
