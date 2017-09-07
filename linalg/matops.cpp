/* Various matrix operations.
 *  
 * Author: Emmanuel Thomé
 * 
 * Copyright 2010 INRIA
 * 
 * This file is part of the CADO project.
 * 
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2.1 of the License, or (at
 * your option) any later version.

 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.

 * You should have received a copy of the GNU Lesser General Public
 * License along with CADO-NFS; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

#include "cado.h"       /* HAVE_* macros ! */
#include "macros.h"

#include "gf2x.h"

#include "matops.h"
#include "utils/memory.h"
#include "utils/misc.h"

#if defined(HAVE_SSE2) && ULONG_BITS == 64
#include <emmintrin.h>
/* {{{ helper macros for sse-2. Copied from gf2x */
/* {{{ _mm_cvtsi64_m64 is not consistent across compiler versions... */
#if defined(__GNUC__) && __GNUC__ == 4 &&__GNUC_MINOR__ == 1
#define _cado_mm_cvtsi64_m64(u) _mm_cvtsi64x_m64((u))
#else
#define _cado_mm_cvtsi64_m64(u) _mm_cvtsi64_m64((u))
#endif
/* }}} */
/* {{{ _cado_mm_setr_epi64 _m128i from 2 int64_t's */
#define _cado_mm_setr_epi64(lo, hi)                      		\
    _mm_setr_epi64(                                      		\
            _cado_mm_cvtsi64_m64((int64_t) (lo)),       		\
            _cado_mm_cvtsi64_m64((int64_t) (hi))        		\
        )
/* }}} */
/* {{{ _cado_mm_set1_epi64 _m128i from 1 int64_t's */
#define _cado_mm_set1_epi64(u) _mm_set1_epi64( _cado_mm_cvtsi64_m64((int64_t) (u)))
/* }}} */
/* {{{ _cado_mm_setr_epi64_c _m128i from 2 int64_t CONSTANTS (and try to get suffix right) */
#define _cado_mm_setr_epi64_c(lo, hi)                    		\
    _mm_setr_epi64(                                      		\
            _cado_mm_cvtsi64_m64(INT64_C(lo)),          		\
            _cado_mm_cvtsi64_m64(INT64_C(hi))           		\
        )
/* }}} */
/* {{{ _cado_mm_set1_epi64_c _m128i from 1 int64_t CONSTANT (and try to get suffix right) */
#define _cado_mm_set1_epi64_c(u) _mm_set1_epi64( _cado_mm_cvtsi64_m64(INT64_C(u)))
/* }}} */
/* {{{ same for 32-bits (which, for some, have SSE-2) */
#define _cado_mm_setr_epi32(a0, a1, a2, a3)				\
    _mm_setr_epi32(                                      		\
            (int32_t) (a0),						\
            (int32_t) (a1),						\
            (int32_t) (a2),						\
            (int32_t) (a3)						\
            )
#define _cado_mm_set1_epi32(u) _mm_set1_epi32( (int32_t) (u))
#define _cado_mm_setr_epi32_c(a0, a1, a2, a3)				\
    _mm_setr_epi32(                                      		\
            (INT32_C(a0)),          					\
            (INT32_C(a1)),           					\
            (INT32_C(a2)),          					\
            (INT32_C(a3))           					\
        )
#define _cado_mm_set1_epi32_c(u) _mm_set1_epi32(INT32_C(u))
/* }}} */
/* }}} */
#endif

#ifdef  HAVE_PCLMUL
#include <wmmintrin.h>
#endif

#ifdef  HAVE_SSE41
#include <smmintrin.h>  // sse 4.1 _mm_cmpeq_epi64
#endif  /* HAVE_SSE41 */

/* We reach the mpfq sources through the gf2x code base, and then these
 * are considered internal cantor-related stuff. We need to include the
 * gf2x config flags before including the mpfq sources.
 */
#include "gf2x/gf2x-config-export.h"
#include "gf2x/gf2x-impl-export.h"
#if ULONG_BITS == 64
#include "mpfq/x86_64/mpfq_2_64.h"
#include "mpfq/x86_64/mpfq_2_128.h"
#elif ULONG_BITS == 32
#include "mpfq/i386/mpfq_2_64.h"
#include "mpfq/i386/mpfq_2_128.h"
#else
#error "neither 32 nor 64 ???"
#endif


/* Define this to check these functions against the m4ri library */
#define xxxHAVE_M4RI
#ifdef  HAVE_M4RI
#include "m4ri.h"
#endif  /* HAVE_M4RI */
#ifdef  HAVE_M4RIE
#include "finite_field_givaro.h"
#include "m4rie.h"
using namespace M4RIE;
#endif  /* HAVE_M4RIE */

/* Some older versions of m4ri had limbs in the wrong order. It's now
 * fixed, as of version 20111203 at least (possibly has been fixed
 * already for quite a while, no recent check of this was done).
 */
#define sometimes_bitrev(x)     (x)
#include "gauss.h"
#include "macros.h"

/* The following is **only** for 64 * 64 matrices */

#define WBITS   64
typedef uint64_t mat64[64] ATTRIBUTE((aligned(64)));
typedef uint64_t * mat64_ptr;
typedef const uint64_t * mat64_srcptr;

static inline uint64_t MAYBE_UNUSED bitrev(uint64_t a)/*{{{*/
{
    a = (a >> 32) ^ (a << 32);
    uint64_t m;
    m = UINT64_C(0x0000ffff0000ffff);
    a = ((a >> 16) & m) ^ ((a << 16) & ~m);
    m = UINT64_C(0x00ff00ff00ff00ff);
    a = ((a >> 8) & m) ^ ((a << 8) & ~m);
    m = UINT64_C(0x0f0f0f0f0f0f0f0f);
    a = ((a >> 4) & m) ^ ((a << 4) & ~m);
    m = UINT64_C(0x3333333333333333);
    a = ((a >> 2) & m) ^ ((a << 2) & ~m);
    m = UINT64_C(0x5555555555555555);
    a = ((a >> 1) & m) ^ ((a << 1) & ~m);
    return a;
}
/* like bitrev, but keep nibbles intact */
static inline uint64_t MAYBE_UNUSED nibrev(uint64_t a)
{
    a = (a >> 32) ^ (a << 32);
    uint64_t m;
    m = UINT64_C(0x0000ffff0000ffff);
    a = ((a >> 16) & m) ^ ((a << 16) & ~m);
    m = UINT64_C(0x00ff00ff00ff00ff);
    a = ((a >> 8) & m) ^ ((a << 8) & ~m);
    m = UINT64_C(0x0f0f0f0f0f0f0f0f);
    a = ((a >> 4) & m) ^ ((a << 4) & ~m);
    return a;
}/*}}}*/


/* level 1 */
#if defined(HAVE_SSE2) && ULONG_BITS == 64
void mul_6464_6464_sse(mat64_ptr C, mat64_srcptr A, mat64_srcptr B)
{
    int i;
    memset(C, 0, sizeof(mat64));
 
    __m128i *Cw = (__m128i *) C;
    __m128i *Aw = (__m128i *) A;

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
}
#endif

void mul_6464_6464_v2(mat64_ptr C, mat64_srcptr A, mat64_srcptr B)
{
    memset(C, 0, sizeof(mat64));
 
    for (int i = 0; i < 64; i++) {
    for (int j = 0; j < 64; j++) {
        if ((A[i]>>j)&1)
            C[i]^=B[j];
    }
    }
}

void add_6464_6464_C(mat64_ptr C, mat64_srcptr A, mat64_srcptr B)
{
    for (int j = 0; j < 64; j++) {
        C[j] = A[j] ^ B[j];
    }
}


void addmul_To64_o64_lsb(uint64_t * r, uint64_t a, uint64_t w)
{
    /* Dans un sens */
    for (unsigned int i = 0; i < 64; i++) {
	*r++ ^= w & -(a & 1);
	a >>= 1;
    }
}

void addmul_To64_o64_msb(uint64_t * r, uint64_t a, uint64_t w)
{
    /* Dans l'autre -- va un poil plus vite. */
    for (unsigned int i = 0; i < 64; i++) {
	*r++ ^= w & (((int64_t) a) >> 63);
	a <<= 1;
    }
}

void addmul_To64_o64_lsb_packof2(uint64_t * r, uint64_t a, uint64_t w)
{
    /* À peu près comme la méthode 1, mais pas mieux */
    typedef uint64_t mvec_t[2];
    mvec_t mb[4] = {
	{0, 0}, {w, 0}, {0, w}, {w, w},
    };
    for (int i = 0; i < 64; i += 2) {
	const uint64_t *y = mb[a & 3];
	*r++ ^= y[0];
	*r++ ^= y[1];
	a >>= 2;
    }
}

#if defined(HAVE_SSE2) && ULONG_BITS == 64
void addmul_To64_o64_lsb_sse_v1(uint64_t * r, uint64_t a, uint64_t w)
{
    /* Using sse-2 */
    __m128i mb[4] = {
	_mm_setzero_si128(),
	_cado_mm_setr_epi64(w, 0),
	_cado_mm_setr_epi64(0, w),
	_cado_mm_set1_epi64(w),
    };
    __m128i *sr = (__m128i *) r;
    for (int i = 0; i < 64; i += 2) {
	*sr = _mm_xor_si128(*sr, mb[a & 3]);
        sr++;
	a >>= 2;
    }
    _mm_empty();
}
#endif

/* implements mul_o64_6464 */
void mul_o64_6464_C_lsb(uint64_t * r, uint64_t a, mat64_srcptr w)
{
    uint64_t c = 0;
    for (unsigned int i = 0; i < WBITS; i++) {
	c ^= (w[i] & -(a & UINT64_C(1)));
	a >>= 1;
    }
    *r = c;
}

void mul_o64_6464_C_msb(uint64_t *r, uint64_t a, mat64_srcptr w)
{
    uint64_t c = 0;
    for (int i = 64 - 1; i >= 0; i--) {
        c ^= (w[i] & (((int64_t) a) >> (64 - 1)));
        a <<= 1;
    }
    *r = c;
}

void mul_o64_T6464_C_parity(uint64_t * w, uint64_t a, mat64_srcptr b)
{
    // Uses unoptimized __builtin_parityl function -- maybe better with gcc 4.3
    // note that popcnt is faster in asm than the more restricted parity
    // functions. So if it's available, it should be tested.
    uint64_t c = 0;
    for (unsigned int i = 0; i < WBITS; i++) {
        uint64_t p = cado_parity64(a & b[i]);
	c ^= p << i;
    }
    *w = c;
}

/* This is stolen from code by D. Harvey. (GPL, thus can't stay here) */
#define XMIX32(b, a) (((((a) << 32) ^ (a)) >> 32) + \
                     ((((b) >> 32) ^ (b)) << 32))
#define XMIX16(b, a) (((((a) >> 16) ^ (a)) & 0x0000FFFF0000FFFFll) +     \
                     ((((b) << 16) ^ (b)) & 0xFFFF0000FFFF0000ll));
#define XMIX8(b, a) (((((a) >> 8) ^ (a)) & 0x00FF00FF00FF00FFll) + \
                    ((((b) << 8) ^ (b)) & 0xFF00FF00FF00FF00ll));
#define XMIX4(b, a) (((((a) >> 4) ^ (a)) & 0x0F0F0F0F0F0F0F0Fll) + \
                    ((((b) << 4) ^ (b)) & 0xF0F0F0F0F0F0F0F0ll));
#define XMIX2(b, a) (((((a) >> 2) ^ (a)) & 0x3333333333333333ll) + \
                    ((((b) << 2) ^ (b)) & 0xCCCCCCCCCCCCCCCCll));
#define XMIX1(b, a) (((((a) >> 1) ^ (a)) & 0x5555555555555555ll) + \
                    ((((b) << 1) ^ (b)) & 0xAAAAAAAAAAAAAAAAll));
static inline uint64_t _parity64_helper2(const uint64_t* buf, uint64_t a)
{
   uint64_t a0, a1, b0, b1, c0, c1;
   a0 = XMIX32(buf[0x20] & a, buf[0x00] & a);
   a1 = XMIX32(buf[0x30] & a, buf[0x10] & a);
   b0 = XMIX16(a1, a0);
   a0 = XMIX32(buf[0x28] & a, buf[0x08] & a);
   a1 = XMIX32(buf[0x38] & a, buf[0x18] & a);
   b1 = XMIX16(a1, a0);
   c0 = XMIX8(b1, b0);
   a0 = XMIX32(buf[0x24] & a, buf[0x04] & a);
   a1 = XMIX32(buf[0x34] & a, buf[0x14] & a);
   b0 = XMIX16(a1, a0);
   a0 = XMIX32(buf[0x2C] & a, buf[0x0C] & a);
   a1 = XMIX32(buf[0x3C] & a, buf[0x1C] & a);
   b1 = XMIX16(a1, a0);
   c1 = XMIX8(b1, b0);
   return XMIX4(c1, c0);
}

void mul_o64_T6464_C_parity3(uint64_t * w, uint64_t a, mat64_srcptr b)
{
   uint64_t d0, d1, e0, e1;

   d0 = _parity64_helper2(b, a);
   d1 = _parity64_helper2(b + 2, a);
   e0 = XMIX2(d1, d0);

   d0 = _parity64_helper2(b + 1, a);
   d1 = _parity64_helper2(b + 3, a);
   e1 = XMIX2(d1, d0);

   *w = XMIX1(e1, e0);
}


void transp_6464(mat64_ptr dst, mat64_srcptr src)
{
    int i, j;
    ASSERT_ALWAYS(dst != src);
    for (i = 0; i < 64; i++) {
	dst[i] = 0;
	for (j = 0; j < 64; j++) {
	    dst[i] ^= ((src[j] >> i) & UINT64_C(1)) << j;
	}
    }
}

void copy_6464(mat64_ptr dst, mat64_srcptr src)
{
    memcpy(dst, src, sizeof(mat64));
}

/* level 2 */

static inline void MAYBE_UNUSED copy_N64(uint64_t * dst, const uint64_t * src, size_t m)
{
    memcpy(dst, src, m * sizeof(uint64_t));
}

static inline int MAYBE_UNUSED cmp_N64(const uint64_t * dst, const uint64_t * src, size_t m)
{
    return memcmp(dst, src, m * sizeof(uint64_t));
}


/* implements mul_N64_6464 */
/* This can work in place (C==A, or C==B, or both) */
void mul_N64_6464_lookup4(uint64_t *C,
                   const uint64_t *A,
                   const uint64_t *B, size_t m)
{
    uint64_t Bx[16][16];
    for(int j = 0 ; j < 16 ; j++) {
        const uint64_t * bb = B + 4 * j;
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
    /* We don't zero out C before the computation, but rather at the
     * moment we read A[i], so that A==C is supported */
    for (size_t i = 0; i < m; i++) {
        uint64_t aa = A[i];
        C[i] = Bx[0][aa & 15]; aa>>=4;
        C[i]^= Bx[1][aa & 15]; aa>>=4;
        C[i]^= Bx[2][aa & 15]; aa>>=4;
        C[i]^= Bx[3][aa & 15]; aa>>=4;
        C[i]^= Bx[4][aa & 15]; aa>>=4;
        C[i]^= Bx[5][aa & 15]; aa>>=4;
        C[i]^= Bx[6][aa & 15]; aa>>=4;
        C[i]^= Bx[7][aa & 15]; aa>>=4;
        C[i]^= Bx[8][aa & 15]; aa>>=4;
        C[i]^= Bx[9][aa & 15]; aa>>=4;
        C[i]^= Bx[10][aa & 15]; aa>>=4;
        C[i]^= Bx[11][aa & 15]; aa>>=4;
        C[i]^= Bx[12][aa & 15]; aa>>=4;
        C[i]^= Bx[13][aa & 15]; aa>>=4;
        C[i]^= Bx[14][aa & 15]; aa>>=4;
        C[i]^= Bx[15][aa];
    }
}
/* This can work in place (C==A, or C==B, or both) */
static inline void MAYBE_UNUSED addmul_N64_6464_lookup4(uint64_t *C,
                   const uint64_t *A,
                   const uint64_t *B, size_t m)
{
    uint64_t Bx[16][16];
    for(int j = 0 ; j < 16 ; j++) {
        const uint64_t * bb = B + 4 * j;
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
    for (size_t i = 0; i < m; i++) {
        uint64_t aa = A[i];
        C[i]^= Bx[0][aa & 15]; aa>>=4;
        C[i]^= Bx[1][aa & 15]; aa>>=4;
        C[i]^= Bx[2][aa & 15]; aa>>=4;
        C[i]^= Bx[3][aa & 15]; aa>>=4;
        C[i]^= Bx[4][aa & 15]; aa>>=4;
        C[i]^= Bx[5][aa & 15]; aa>>=4;
        C[i]^= Bx[6][aa & 15]; aa>>=4;
        C[i]^= Bx[7][aa & 15]; aa>>=4;
        C[i]^= Bx[8][aa & 15]; aa>>=4;
        C[i]^= Bx[9][aa & 15]; aa>>=4;
        C[i]^= Bx[10][aa & 15]; aa>>=4;
        C[i]^= Bx[11][aa & 15]; aa>>=4;
        C[i]^= Bx[12][aa & 15]; aa>>=4;
        C[i]^= Bx[13][aa & 15]; aa>>=4;
        C[i]^= Bx[14][aa & 15]; aa>>=4;
        C[i]^= Bx[15][aa];
    }
}

void mul_N64_6464_lookup8(uint64_t *C,
                   const uint64_t *A,
                   const uint64_t *B, size_t m)
{
    uint64_t Bx[8][256];
    for(int j = 0 ; j < 8 ; j++) {
        const uint64_t * bb = B + 8 * j;
        uint64_t w = 0;
        Bx[j][0] = w; w ^= bb[0];
        Bx[j][1] = w; w ^= bb[1];
        Bx[j][3] = w; w ^= bb[0];
        Bx[j][2] = w; w ^= bb[2];
        Bx[j][6] = w; w ^= bb[0];
        Bx[j][7] = w; w ^= bb[1];
        Bx[j][5] = w; w ^= bb[0];
        Bx[j][4] = w; w ^= bb[3];
        Bx[j][12] = w; w ^= bb[0];
        Bx[j][13] = w; w ^= bb[1];
        Bx[j][15] = w; w ^= bb[0];
        Bx[j][14] = w; w ^= bb[2];
        Bx[j][10] = w; w ^= bb[0];
        Bx[j][11] = w; w ^= bb[1];
        Bx[j][9] = w; w ^= bb[0];
        Bx[j][8] = w; w ^= bb[4];
        Bx[j][24] = w; w ^= bb[0];
        Bx[j][25] = w; w ^= bb[1];
        Bx[j][27] = w; w ^= bb[0];
        Bx[j][26] = w; w ^= bb[2];
        Bx[j][30] = w; w ^= bb[0];
        Bx[j][31] = w; w ^= bb[1];
        Bx[j][29] = w; w ^= bb[0];
        Bx[j][28] = w; w ^= bb[3];
        Bx[j][20] = w; w ^= bb[0];
        Bx[j][21] = w; w ^= bb[1];
        Bx[j][23] = w; w ^= bb[0];
        Bx[j][22] = w; w ^= bb[2];
        Bx[j][18] = w; w ^= bb[0];
        Bx[j][19] = w; w ^= bb[1];
        Bx[j][17] = w; w ^= bb[0];
        Bx[j][16] = w; w ^= bb[5];
        Bx[j][48] = w; w ^= bb[0];
        Bx[j][49] = w; w ^= bb[1];
        Bx[j][51] = w; w ^= bb[0];
        Bx[j][50] = w; w ^= bb[2];
        Bx[j][54] = w; w ^= bb[0];
        Bx[j][55] = w; w ^= bb[1];
        Bx[j][53] = w; w ^= bb[0];
        Bx[j][52] = w; w ^= bb[3];
        Bx[j][60] = w; w ^= bb[0];
        Bx[j][61] = w; w ^= bb[1];
        Bx[j][63] = w; w ^= bb[0];
        Bx[j][62] = w; w ^= bb[2];
        Bx[j][58] = w; w ^= bb[0];
        Bx[j][59] = w; w ^= bb[1];
        Bx[j][57] = w; w ^= bb[0];
        Bx[j][56] = w; w ^= bb[4];
        Bx[j][40] = w; w ^= bb[0];
        Bx[j][41] = w; w ^= bb[1];
        Bx[j][43] = w; w ^= bb[0];
        Bx[j][42] = w; w ^= bb[2];
        Bx[j][46] = w; w ^= bb[0];
        Bx[j][47] = w; w ^= bb[1];
        Bx[j][45] = w; w ^= bb[0];
        Bx[j][44] = w; w ^= bb[3];
        Bx[j][36] = w; w ^= bb[0];
        Bx[j][37] = w; w ^= bb[1];
        Bx[j][39] = w; w ^= bb[0];
        Bx[j][38] = w; w ^= bb[2];
        Bx[j][34] = w; w ^= bb[0];
        Bx[j][35] = w; w ^= bb[1];
        Bx[j][33] = w; w ^= bb[0];
        Bx[j][32] = w; w ^= bb[6];
        Bx[j][96] = w; w ^= bb[0];
        Bx[j][97] = w; w ^= bb[1];
        Bx[j][99] = w; w ^= bb[0];
        Bx[j][98] = w; w ^= bb[2];
        Bx[j][102] = w; w ^= bb[0];
        Bx[j][103] = w; w ^= bb[1];
        Bx[j][101] = w; w ^= bb[0];
        Bx[j][100] = w; w ^= bb[3];
        Bx[j][108] = w; w ^= bb[0];
        Bx[j][109] = w; w ^= bb[1];
        Bx[j][111] = w; w ^= bb[0];
        Bx[j][110] = w; w ^= bb[2];
        Bx[j][106] = w; w ^= bb[0];
        Bx[j][107] = w; w ^= bb[1];
        Bx[j][105] = w; w ^= bb[0];
        Bx[j][104] = w; w ^= bb[4];
        Bx[j][120] = w; w ^= bb[0];
        Bx[j][121] = w; w ^= bb[1];
        Bx[j][123] = w; w ^= bb[0];
        Bx[j][122] = w; w ^= bb[2];
        Bx[j][126] = w; w ^= bb[0];
        Bx[j][127] = w; w ^= bb[1];
        Bx[j][125] = w; w ^= bb[0];
        Bx[j][124] = w; w ^= bb[3];
        Bx[j][116] = w; w ^= bb[0];
        Bx[j][117] = w; w ^= bb[1];
        Bx[j][119] = w; w ^= bb[0];
        Bx[j][118] = w; w ^= bb[2];
        Bx[j][114] = w; w ^= bb[0];
        Bx[j][115] = w; w ^= bb[1];
        Bx[j][113] = w; w ^= bb[0];
        Bx[j][112] = w; w ^= bb[5];
        Bx[j][80] = w; w ^= bb[0];
        Bx[j][81] = w; w ^= bb[1];
        Bx[j][83] = w; w ^= bb[0];
        Bx[j][82] = w; w ^= bb[2];
        Bx[j][86] = w; w ^= bb[0];
        Bx[j][87] = w; w ^= bb[1];
        Bx[j][85] = w; w ^= bb[0];
        Bx[j][84] = w; w ^= bb[3];
        Bx[j][92] = w; w ^= bb[0];
        Bx[j][93] = w; w ^= bb[1];
        Bx[j][95] = w; w ^= bb[0];
        Bx[j][94] = w; w ^= bb[2];
        Bx[j][90] = w; w ^= bb[0];
        Bx[j][91] = w; w ^= bb[1];
        Bx[j][89] = w; w ^= bb[0];
        Bx[j][88] = w; w ^= bb[4];
        Bx[j][72] = w; w ^= bb[0];
        Bx[j][73] = w; w ^= bb[1];
        Bx[j][75] = w; w ^= bb[0];
        Bx[j][74] = w; w ^= bb[2];
        Bx[j][78] = w; w ^= bb[0];
        Bx[j][79] = w; w ^= bb[1];
        Bx[j][77] = w; w ^= bb[0];
        Bx[j][76] = w; w ^= bb[3];
        Bx[j][68] = w; w ^= bb[0];
        Bx[j][69] = w; w ^= bb[1];
        Bx[j][71] = w; w ^= bb[0];
        Bx[j][70] = w; w ^= bb[2];
        Bx[j][66] = w; w ^= bb[0];
        Bx[j][67] = w; w ^= bb[1];
        Bx[j][65] = w; w ^= bb[0];
        Bx[j][64] = w; w ^= bb[7];
        Bx[j][192] = w; w ^= bb[0];
        Bx[j][193] = w; w ^= bb[1];
        Bx[j][195] = w; w ^= bb[0];
        Bx[j][194] = w; w ^= bb[2];
        Bx[j][198] = w; w ^= bb[0];
        Bx[j][199] = w; w ^= bb[1];
        Bx[j][197] = w; w ^= bb[0];
        Bx[j][196] = w; w ^= bb[3];
        Bx[j][204] = w; w ^= bb[0];
        Bx[j][205] = w; w ^= bb[1];
        Bx[j][207] = w; w ^= bb[0];
        Bx[j][206] = w; w ^= bb[2];
        Bx[j][202] = w; w ^= bb[0];
        Bx[j][203] = w; w ^= bb[1];
        Bx[j][201] = w; w ^= bb[0];
        Bx[j][200] = w; w ^= bb[4];
        Bx[j][216] = w; w ^= bb[0];
        Bx[j][217] = w; w ^= bb[1];
        Bx[j][219] = w; w ^= bb[0];
        Bx[j][218] = w; w ^= bb[2];
        Bx[j][222] = w; w ^= bb[0];
        Bx[j][223] = w; w ^= bb[1];
        Bx[j][221] = w; w ^= bb[0];
        Bx[j][220] = w; w ^= bb[3];
        Bx[j][212] = w; w ^= bb[0];
        Bx[j][213] = w; w ^= bb[1];
        Bx[j][215] = w; w ^= bb[0];
        Bx[j][214] = w; w ^= bb[2];
        Bx[j][210] = w; w ^= bb[0];
        Bx[j][211] = w; w ^= bb[1];
        Bx[j][209] = w; w ^= bb[0];
        Bx[j][208] = w; w ^= bb[5];
        Bx[j][240] = w; w ^= bb[0];
        Bx[j][241] = w; w ^= bb[1];
        Bx[j][243] = w; w ^= bb[0];
        Bx[j][242] = w; w ^= bb[2];
        Bx[j][246] = w; w ^= bb[0];
        Bx[j][247] = w; w ^= bb[1];
        Bx[j][245] = w; w ^= bb[0];
        Bx[j][244] = w; w ^= bb[3];
        Bx[j][252] = w; w ^= bb[0];
        Bx[j][253] = w; w ^= bb[1];
        Bx[j][255] = w; w ^= bb[0];
        Bx[j][254] = w; w ^= bb[2];
        Bx[j][250] = w; w ^= bb[0];
        Bx[j][251] = w; w ^= bb[1];
        Bx[j][249] = w; w ^= bb[0];
        Bx[j][248] = w; w ^= bb[4];
        Bx[j][232] = w; w ^= bb[0];
        Bx[j][233] = w; w ^= bb[1];
        Bx[j][235] = w; w ^= bb[0];
        Bx[j][234] = w; w ^= bb[2];
        Bx[j][238] = w; w ^= bb[0];
        Bx[j][239] = w; w ^= bb[1];
        Bx[j][237] = w; w ^= bb[0];
        Bx[j][236] = w; w ^= bb[3];
        Bx[j][228] = w; w ^= bb[0];
        Bx[j][229] = w; w ^= bb[1];
        Bx[j][231] = w; w ^= bb[0];
        Bx[j][230] = w; w ^= bb[2];
        Bx[j][226] = w; w ^= bb[0];
        Bx[j][227] = w; w ^= bb[1];
        Bx[j][225] = w; w ^= bb[0];
        Bx[j][224] = w; w ^= bb[6];
        Bx[j][160] = w; w ^= bb[0];
        Bx[j][161] = w; w ^= bb[1];
        Bx[j][163] = w; w ^= bb[0];
        Bx[j][162] = w; w ^= bb[2];
        Bx[j][166] = w; w ^= bb[0];
        Bx[j][167] = w; w ^= bb[1];
        Bx[j][165] = w; w ^= bb[0];
        Bx[j][164] = w; w ^= bb[3];
        Bx[j][172] = w; w ^= bb[0];
        Bx[j][173] = w; w ^= bb[1];
        Bx[j][175] = w; w ^= bb[0];
        Bx[j][174] = w; w ^= bb[2];
        Bx[j][170] = w; w ^= bb[0];
        Bx[j][171] = w; w ^= bb[1];
        Bx[j][169] = w; w ^= bb[0];
        Bx[j][168] = w; w ^= bb[4];
        Bx[j][184] = w; w ^= bb[0];
        Bx[j][185] = w; w ^= bb[1];
        Bx[j][187] = w; w ^= bb[0];
        Bx[j][186] = w; w ^= bb[2];
        Bx[j][190] = w; w ^= bb[0];
        Bx[j][191] = w; w ^= bb[1];
        Bx[j][189] = w; w ^= bb[0];
        Bx[j][188] = w; w ^= bb[3];
        Bx[j][180] = w; w ^= bb[0];
        Bx[j][181] = w; w ^= bb[1];
        Bx[j][183] = w; w ^= bb[0];
        Bx[j][182] = w; w ^= bb[2];
        Bx[j][178] = w; w ^= bb[0];
        Bx[j][179] = w; w ^= bb[1];
        Bx[j][177] = w; w ^= bb[0];
        Bx[j][176] = w; w ^= bb[5];
        Bx[j][144] = w; w ^= bb[0];
        Bx[j][145] = w; w ^= bb[1];
        Bx[j][147] = w; w ^= bb[0];
        Bx[j][146] = w; w ^= bb[2];
        Bx[j][150] = w; w ^= bb[0];
        Bx[j][151] = w; w ^= bb[1];
        Bx[j][149] = w; w ^= bb[0];
        Bx[j][148] = w; w ^= bb[3];
        Bx[j][156] = w; w ^= bb[0];
        Bx[j][157] = w; w ^= bb[1];
        Bx[j][159] = w; w ^= bb[0];
        Bx[j][158] = w; w ^= bb[2];
        Bx[j][154] = w; w ^= bb[0];
        Bx[j][155] = w; w ^= bb[1];
        Bx[j][153] = w; w ^= bb[0];
        Bx[j][152] = w; w ^= bb[4];
        Bx[j][136] = w; w ^= bb[0];
        Bx[j][137] = w; w ^= bb[1];
        Bx[j][139] = w; w ^= bb[0];
        Bx[j][138] = w; w ^= bb[2];
        Bx[j][142] = w; w ^= bb[0];
        Bx[j][143] = w; w ^= bb[1];
        Bx[j][141] = w; w ^= bb[0];
        Bx[j][140] = w; w ^= bb[3];
        Bx[j][132] = w; w ^= bb[0];
        Bx[j][133] = w; w ^= bb[1];
        Bx[j][135] = w; w ^= bb[0];
        Bx[j][134] = w; w ^= bb[2];
        Bx[j][130] = w; w ^= bb[0];
        Bx[j][131] = w; w ^= bb[1];
        Bx[j][129] = w; w ^= bb[0];
        Bx[j][128] = w;
    }
    memset(C, 0, m * sizeof(uint64_t));
    for (size_t i = 0; i < m; i++) {
        uint64_t aa = A[i];
        C[i] = Bx[0][aa & 255]; aa>>=8;
        C[i]^= Bx[1][aa & 255]; aa>>=8;
        C[i]^= Bx[2][aa & 255]; aa>>=8;
        C[i]^= Bx[3][aa & 255]; aa>>=8;
        C[i]^= Bx[4][aa & 255]; aa>>=8;
        C[i]^= Bx[5][aa & 255]; aa>>=8;
        C[i]^= Bx[6][aa & 255]; aa>>=8;
        C[i]^= Bx[7][aa & 255];
    }
}

void mul_N64_6464_vec(uint64_t *C,
                   const uint64_t *A,
                   const uint64_t *B, size_t m)
{

    memset(C, 0, m * sizeof(uint64_t));
    for (size_t i = 0; i < m; i++) {
        mul_o64_6464(C++, *A++, B);
    }
}

void mul_N64_T6464_vec(uint64_t *C,
                   const uint64_t *A,
                   const uint64_t *B, size_t m)
{

    memset(C, 0, m * sizeof(uint64_t));
    for (size_t i = 0; i < m; i++) {
        mul_o64_T6464(C++, *A++, B);
    }
}

void mul_N64_6464_transB(uint64_t *C,
                   const uint64_t *A,
                   const uint64_t *B, size_t m)
{
    uint64_t *tb = (uint64_t *) malloc(64 * sizeof(uint64_t));
    transp_6464(tb, B);
    mul_N64_T6464(C, A, tb, m);
    free(tb);
}

void mul_N64_T6464_transB(uint64_t *C,
                   const uint64_t *A,
                   const uint64_t *B, size_t m)
{
    uint64_t *tb = (uint64_t *) malloc(64 * sizeof(uint64_t));
    transp_6464(tb, B);
    mul_N64_6464(C, A, tb, m);
    free(tb);
}

#if defined(HAVE_SSE2) && ULONG_BITS == 64
/* implements addmul_N64_6464 */
void addmul_N64_6464_sse(uint64_t *C,
		 const uint64_t *A,
		 const uint64_t *B, size_t m)
{
    size_t j;
    __m128i *Cw = (__m128i *) C;
    __m128i *Aw = (__m128i *) A;

    /* If m is odd, then we can't do sse all the way through because of
     * data width */
    for (j = 0; j < m - 1; j += 2) {
        __m128i c = _mm_setzero_si128();
	__m128i a = *Aw++;

        __m128i one = _cado_mm_set1_epi64_c(1);
	for (int i = 0; i < 64; i++) {
	    __m128i bw = _cado_mm_set1_epi64(B[i]);
	    // c ^= (bw & -(a & one));
            c = _mm_xor_si128(c, _mm_and_si128(bw, _mm_sub_epi64(_mm_setzero_si128(), _mm_and_si128(a, one))));
	    a = _mm_srli_epi64(a, 1);
	}
	*Cw = _mm_xor_si128(*Cw, c);
        Cw++;
    }
    C += j;
    A += j;
    for (; j < m; j++) {
	uint64_t c = UINT64_C(0);
	uint64_t a = *A++;
	for (int i = 0; i < 64; i++) {
	    c ^= (B[i] & -(a & UINT64_C(1)));
	    a >>= UINT64_C(1);
	}
	*C++ ^= c;
    }
}
/* can work in place, so not simply memset0 + addmul (the ^= have been
 * changed to =)
 */
void mul_N64_6464_sse(uint64_t *C,
		 const uint64_t *A,
		 const uint64_t *B, size_t m)
{
    size_t j;
    __m128i *Cw = (__m128i *) C;
    __m128i *Aw = (__m128i *) A;

    /* If m is odd, then we can't do sse all the way through because of
     * data width */
    for (j = 0; j < m - 1; j += 2) {
        __m128i c = _mm_setzero_si128();
	__m128i a = *Aw++;

        __m128i one = _cado_mm_set1_epi64_c(1);
	for (int i = 0; i < 64; i++) {
	    __m128i bw = _cado_mm_set1_epi64(B[i]);
            // c ^= (bw & -(a & one));
            c = _mm_xor_si128(c, _mm_and_si128(bw, _mm_sub_epi64(_mm_setzero_si128(), _mm_and_si128(a, one))));
	    a = _mm_srli_epi64(a, 1);
	}
	*Cw++ = c;
    }
    C += j;
    A += j;
    for (; j < m; j++) {
	uint64_t c = UINT64_C(0);
	uint64_t a = *A++;
	for (int i = 0; i < 64; i++) {
	    c ^= (B[i] & -(a & UINT64_C(1)));
	    a >>= UINT64_C(1);
	}
	*C++ = c;
    }
}
#endif

void mul_64N_N64_addmul(uint64_t *r, uint64_t *a, uint64_t *w, size_t n)
{
    memset(r, 0, 64 * sizeof(uint64_t));
    for (size_t i = 0; i < n; i++) {
        addmul_To64_o64(r, a[i], w[i]);
    }
}


void mul_TN32_N64_C(uint64_t * b, uint32_t * A, uint64_t * x, unsigned int ncol)
{
    uint32_t idx, i, rA;
    uint64_t rx;

    memset(b, 0, 32 * sizeof(uint64_t));
    for(idx = 0; idx < ncol; idx++) {
        rA = A[idx];
        rx = x[idx];
        for(i = 0; i < 32; i++) {
            b[i] ^= rx & -(rA & 1);
            rA >>= 1;
        }
    }
}

/* This takes, in row major order, an Nx64 matrix A (transpose of a 64xN
 * matrix), together with another Nx64 matrix B, and xors the output
 * matrix with the product transpose(A)*B -- this may as well be seen as
 * the block dot product of A and B.
 */
void addmul_TN64_N64_C(uint64_t * b, uint64_t * A, uint64_t * x, unsigned int ncol)
{
    uint64_t idx, i, rA;
    uint64_t rx;

    for(idx = 0; idx < ncol; idx++) {
        rA = A[idx];
        rx = x[idx];
        for(i = 0; i < 64; i++) {
            b[i] ^= rx & -(rA & 1);
            rA >>= 1;
        }
    }
}

void mul_TN64_N64_C(uint64_t * b, uint64_t * A, uint64_t * x, unsigned int ncol)
{
    memset(b, 0, 64 * sizeof(uint64_t));
    addmul_TN64_N64_C(b, A, x, ncol);
}

#if defined(HAVE_SSE2) && ULONG_BITS == 64
static inline void MAYBE_UNUSED mul_TN64K_N64_sse2(uint64_t * w, uint64_t * u, uint64_t * v, unsigned int n, unsigned int K)
{
    memset(w, 0, 64 * K * sizeof(uint64_t));
    for(unsigned int i = 0 ; i < n ; i++) {
        __m128i * w0 = (__m128i*) w;
        // TODO: It's possible to expand more, and use a __m128i
        // mb[4][2], or even [4]. This wouldn't change the code much
        // (see the u128 version), and is likely to speed things up a
        // wee bit maybe.
        __m128i mb[4] = {
            _mm_setzero_si128(),
            _cado_mm_setr_epi64(*v, 0 ),
            _cado_mm_setr_epi64(0,  *v),
            _cado_mm_set1_epi64(*v),
        };
        v++;
        __m128i *sw = w0;
        for(unsigned int k = 0 ; k < K ; k++) {
            uint64_t a = *u++;
            for (unsigned int j = 0; j < 64; j += 2) {
                *sw = _mm_xor_si128(*sw, mb[a & 3]);
                a >>= 2;
                sw ++;
            }
        }
    }
}
#endif

static inline void MAYBE_UNUSED mul_TN64K_N64_C(uint64_t * b, uint64_t * A, uint64_t * x, unsigned int ncol, unsigned int K)
{
    uint64_t idx, i, rA;
    uint64_t rx;

    memset(b, 0, 64 * K * sizeof(uint64_t));
    for(idx = 0; idx < ncol; idx++) {
        rx = x[idx];
        uint64_t* pb = b;
        for(unsigned int j = 0 ; j < K ; j++) {
            rA = *A++;
            for(i = 0; i < 64; i++) {
                *pb++ ^= rx & -(rA & 1);
                rA >>= 1;
            }
        }
    }
}
#if 0   /* haven't checked yet what the funny-named functions actually do... *//*{{{*/
static inline void TVUBit_v2(unsigned long m,
	       unsigned long n,
	       const unsigned long *A, const unsigned long *B,
	       unsigned long *C)
{
    unsigned long i, P, k;
    memset(C, 0, WBITS * sizeof(unsigned long));
    for (i = 0; i < m; i++) {
	P = *A++;
	for (k = 0; k < WBITS; k++) {
	    //if (P & UINT64_C(1)) C[k] ^= B[i];
	    C[k] ^= (B[i] & -(P & UINT64_C(1)));
	    P >>= UINT64_C(1);
	}
    }
}
#if 1
static inline void VUBit_v2(unsigned long m,
	      unsigned long n,
	      const unsigned long *A, const unsigned long *B,
	      unsigned long *C)
{
    unsigned long i, P, k;
    memset(C, 0, m * sizeof(unsigned long));
    for (i = 0; i < m; i++) {
	P = *A++;
	for (k = 0; k < n; k++) {
	    C[i] ^= (B[k] & -(P & UINT64_C(1)));
	    P >>= UINT64_C(1);
	}
    }
}
#endif
#if 0
static inline void VUBit(unsigned long m,
	      unsigned long n,
	      const unsigned long *A, const unsigned long *B,
	      unsigned long *C)
{
    unsigned long i, P, k, j;
    memset(C, 0, m * sizeof(unsigned long));
    for (i = 0; i < m; i++) {
	P = *A++;
	j = C[i];
	for (k = 0; k < n; k++) {
	    j ^= (B[k] & -(P & UINT64_C(1)));
	    P >>= UINT64_C(1);
	}
	C[i] = j;
    }
}
#endif
#endif/*}}}*/


/* polynomials */

/* lengths of a1 and a2 are n1 and n2 */
void m64pol_addmul(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2, unsigned int n1, unsigned int n2)
{
    assert(r != a1 && r != a2);
    for(unsigned int i = 0 ; i < n1 ; i++) {
        for(unsigned int j = 0 ; j < n2 ; j++) {
            mat64 x;
            mul_6464_6464(x, a1[i], a2[j]);
            add_6464_6464(r[i+j], r[i+j], x);
        }
    }
}

void m64pol_add(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2, unsigned int n)
{
    for(unsigned int i = 0 ; i < n ; i++) {
        add_6464_6464(r[i], a1[i], a2[i]);
    }
}

void m64pol_mul(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2, unsigned int n1, unsigned int n2)
{
    memset(r, 0, (n1 + n2 - 1) * sizeof(mat64));
    m64pol_addmul(r, a1, a2, n1, n2);
}

void m64pol_mul_kara(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2, unsigned int n1, unsigned int n2)
{
    assert(r != a1 && r != a2);
    assert(n1 == n2);
    assert((n1 & (n1 - 1)) == 0);
    /* As is certainly not surprising, karatsuba wins as early on as one
     * can imagine */
    if (n1 == 1) {
        m64pol_mul(r, a1, a2, n1, n2);
        return;
    }
    unsigned int h = n1 >> 1;
    memset(r, 0, (n1 + n2 - 1) * sizeof(mat64));

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
}

void m64pol_addmul_kara(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2, unsigned int n1, unsigned int n2)
{
    mat64 * t = (mat64 *) malloc_aligned((n1 + n2 -1) * sizeof(mat64), 64);
    m64pol_mul_kara(t, a1, a2, n1, n2);
    m64pol_add(r, r, t, n1 + n2 - 1);
    free_aligned(t);
}

void m64pol_mul_gf2_64_bitslice(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2)
{
    unsigned int n1 = 64;
    unsigned int n2 = 64;
    mat64 * t = (mat64 *) malloc_aligned((n1 + n2 -1) * sizeof(mat64), 64);
    m64pol_mul_kara(t, a1, a2, n1, n2);
    /* This reduces modulo the polynomial x^64+x^4+x^3+x+1 */
    for(unsigned int i = n1 + n2 - 2 ; i >= 64 ; i--) {
        add_6464_6464(t[i-64+4], t[i-64+4], t[i]);
        add_6464_6464(t[i-64+3], t[i-64+3], t[i]);
        add_6464_6464(t[i-64+1], t[i-64+1], t[i]);
        add_6464_6464(t[i-64  ], t[i-64  ], t[i]);
    }
    memcpy(r, t, 64 * sizeof(mat64));
    free_aligned(t);
}

void m64pol_scalmul_gf2_64_bitslice(m64pol_ptr r, m64pol_srcptr a, uint64_t * s)
{
    unsigned int n1 = 64;
    unsigned int n2 = 64;
    mat64 * t = (mat64 *) malloc_aligned((n1 + n2 -1) * sizeof(mat64), 64);
    memset(t, 0, (n1 + n2 -1) * sizeof(mat64));
    for(unsigned int i = 0 ; i < 64 ; i++) {
        if (((s[0]>>i)&UINT64_C(1))==0) continue;
        m64pol_add(t+i, t+i, a, 64);
        // for(unsigned int j = 0 ; j < 64 ; j++) { add_6464_6464(t[i+j], t[i+j], a[j]); }
    }
    /* This reduces modulo the polynomial x^64+x^4+x^3+x+1 */
    for(unsigned int i = n1 + n2 - 2 ; i >= 64 ; i--) {
        add_6464_6464(t[i-64+4], t[i-64+4], t[i]);
        add_6464_6464(t[i-64+3], t[i-64+3], t[i]);
        add_6464_6464(t[i-64+1], t[i-64+1], t[i]);
        add_6464_6464(t[i-64  ], t[i-64  ], t[i]);
    }
    memcpy(r, t, 64 * sizeof(mat64));
    free_aligned(t);
}

void m64pol_scalmul_gf2_64_bitslice2(m64pol_ptr r, m64pol_srcptr a, uint64_t * s)
{
    /* Now try with precomputation of multiples. We'll do only four of
     * them to start with. */
    unsigned int n1 = 64;
    unsigned int n2 = 64;
    mat64 * t = (mat64 *) malloc_aligned((n1 + n2 -1) * sizeof(mat64), 64);
    memset(t, 0, (n1 + n2 -1) * sizeof(mat64));

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

    memset(am_area, 0, (1 << NMULTS) * (64 + NMULTS - 1) * sizeof(mat64));
    memcpy(am[1], a, 64 * sizeof(mat64));
    for(unsigned int j = 1 ; j < NMULTS ; j++) {
        /* Duplicate all stuff having msb set from level below */
        for(unsigned int i = (1u << (j-1)) ; i < (1u << j) ; i++) {
            memcpy(am[(i<<1)] + 1, am[i], (64 + j - 1) * sizeof(mat64));
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
        add_6464_6464(t[i-64+4], t[i-64+4], t[i]);
        add_6464_6464(t[i-64+3], t[i-64+3], t[i]);
        add_6464_6464(t[i-64+1], t[i-64+1], t[i]);
        add_6464_6464(t[i-64  ], t[i-64  ], t[i]);
    }
    memcpy(r, t, 64 * sizeof(mat64));
    free_aligned(t);
}

void m64pol_mul_gf2_64_nobitslice(uint64_t * r, uint64_t * a1, uint64_t * a2)
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
}

void m64pol_scalmul_gf2_64_nobitslice(uint64_t * r, uint64_t * a, uint64_t * scalar)
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
}


void m64pol_mul_gf2_128_bitslice(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2)
{
    unsigned int n1 = 128;
    unsigned int n2 = 128;
    mat64 * t = (mat64 *) malloc_aligned((n1 + n2 -1) * sizeof(mat64), 64);
    m64pol_mul_kara(t, a1, a2, n1, n2);
    /* This reduces modulo the polynomial x^128+x^7+x^2+x+1 */
    for(unsigned int i = n1 + n2 - 2 ; i >= 128; i--) {
        add_6464_6464(t[i-128+7], t[i-128+7], t[i]);
        add_6464_6464(t[i-128+2], t[i-128+2], t[i]);
        add_6464_6464(t[i-128+1], t[i-128+1], t[i]);
        add_6464_6464(t[i-128  ], t[i-128  ], t[i]);
    }
    memcpy(r, t, 128 * sizeof(mat64));
    free_aligned(t);
}

void m64pol_scalmul_gf2_128_bitslice(m64pol_ptr r, m64pol_srcptr a, uint64_t * s)
{
    unsigned int n1 = 128;
    unsigned int n2 = 128;
    mat64 * t = (mat64 *) malloc_aligned((n1 + n2 -1) * sizeof(mat64), 64);
    memset(t, 0, (n1 + n2 -1) * sizeof(mat64));
    for(unsigned int i = 0 ; i < 128 ; i++) {
        if (((s[i/64]>>(i&63))&UINT64_C(1))==0) continue;
        for(unsigned int j = 0 ; j < 64 ; j++) {
            add_6464_6464(t[i+j], t[i+j], a[j]);
        }
    }
    /* This reduces modulo the polynomial x^128+x^7+x^2+x+1 */
    for(unsigned int i = n1 + n2 - 2 ; i >= 128; i--) {
        add_6464_6464(t[i-128+7], t[i-128+7], t[i]);
        add_6464_6464(t[i-128+2], t[i-128+2], t[i]);
        add_6464_6464(t[i-128+1], t[i-128+1], t[i]);
        add_6464_6464(t[i-128  ], t[i-128  ], t[i]);
    }
    memcpy(r, t, 128 * sizeof(mat64));
    free_aligned(t);
}


void m64pol_mul_gf2_128_nobitslice(uint64_t * r, uint64_t * a1, uint64_t * a2)
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
}

void m64pol_scalmul_gf2_128_nobitslice(uint64_t * r, uint64_t * a, uint64_t * scalar)
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
}



// 2.11 -- 512 512 based on basecase @ 512
// 1.49 -- 512 512 based on basecase @ 256
// 1.07 -- 512 512 based on basecase @ 128
// 0.76 -- 512 512 based on basecase @ 64
// 0.56 -- 512 512 based on basecase @ 32
// 0.41 -- 512 512 based on basecase @ 16
// 0.31 -- 512 512 based on basecase @ 8
// 0.24 -- 512 512 based on basecase @ 4
// 0.19 -- 512 512 based on basecase @ 2
// 0.14 -- 512 512 based on basecase @ 1

/* Same spirit, but treat multiplication of 64K by 64K matrices (of
 * polynomials).
 *
 * That is, we consider a K*K matrix of polynomials of 64*64 matrices.
 * 
 * We assume that there is no pointer aliasing, and that matrices are
 * stored row-major, with all polynomials contiguous (and of fixed
 * lengths n1, resp n2).
 */
void m64polblock_mul(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2, unsigned int n1, unsigned int n2, unsigned int K)
{
    assert(r != a1 && r != a2);
    memset(r, 0, (n1 + n2 - 1) * K * K * sizeof(mat64));
    for(unsigned int i = 0 ; i < K ; i++) {
        m64pol_srcptr ra1 = a1 + i * n1 * K;
        m64pol_srcptr rr = r + i * (n1 + n2 - 1) * K;
    for(unsigned int j = 0 ; j < K ; j++) {
        m64pol_srcptr ca2 = a2 + j * n2;
        m64pol_srcptr pr = rr + j * (n1 + n2 - 1);
    for(unsigned int k = 0 ; k < K ; k++) {
        m64pol_srcptr pa1 = ra1 + k * n1;
        m64pol_srcptr pa2 = ca2 + k * n2 * K;
        m64pol_addmul(pr, pa1, pa2, n1, n2);
    }
    }
    }
}

void m64polblock_mul_kara(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2, unsigned int n1, unsigned int n2, unsigned int K)
{
    assert(r != a1 && r != a2);
    memset(r, 0, (n1 + n2 - 1) * K * K * sizeof(mat64));
    for(unsigned int i = 0 ; i < K ; i++) {
        m64pol_srcptr ra1 = a1 + i * n1 * K;
        m64pol_srcptr rr = r + i * (n1 + n2 - 1) * K;
    for(unsigned int j = 0 ; j < K ; j++) {
        m64pol_srcptr ca2 = a2 + j * n2;
        m64pol_srcptr pr = rr + j * (n1 + n2 - 1);
    for(unsigned int k = 0 ; k < K ; k++) {
        m64pol_srcptr pa1 = ra1 + k * n1;
        m64pol_srcptr pa2 = ca2 + k * n2 * K;
        m64pol_addmul_kara(pr, pa1, pa2, n1, n2);
    }
    }
    }
}

#ifdef  HAVE_M4RI
void my_mzd_randomize(mzd_t * A)
{
    for(size_t i = 0 ; i < A->nrows ; i++) {
        uint64_t * ptr = (uint64_t *) A->rows[i];
        for(size_t j = 0 ; j < A->width ; j++) {
            ptr[j] = rand64();
        }
    }
}

void mzd_set_mem(mzd_t * M, const uint64_t * s, unsigned int n)
{
    assert(M->width == 1);
    assert(M->nrows == n);
    for(size_t i = 0 ; i < M->nrows ; i++) {
        uint64_t * ptr = (uint64_t *) M->rows[i];
        ptr[0] = sometimes_bitrev(s[i]);
    }
}

void mzd_set_memT(mzd_t * M, const uint64_t * s, unsigned int n)
{
    mzd_set_mem(M, s, n);
    mzd_transpose(M, M);
}

void mzd_check_mem(mzd_t * M, uint64_t * s, unsigned int n)
{
    assert(M->width == 1);
    assert(M->nrows == n);
    for(size_t i = 0 ; i < M->nrows ; i++) {
        uint64_t * ptr = (uint64_t *) M->rows[i];
        if (ptr[0] != sometimes_bitrev(s[i])) {
            fprintf(stderr, "Rows %zu differ: %016" PRIx64 " != %016" PRIx64 "\n",
                    i, sometimes_bitrev(ptr[0]), s[i]);
            abort();
        }
    }
}
void mzd_check_memT(mzd_t * M, uint64_t * s, unsigned int n)
{
    mzd_t * tmp = mzd_transpose(NULL, M);
    mzd_check_mem(M, s, n);
    mzd_free(tmp);
}
#endif  /* HAVE_M4RI */


/* level 3 -- linear systems */

int gauss_6464_C(mat64 mm, mat64 e, mat64 m)
{
    memcpy(mm,m,sizeof(mat64));
    uint64_t * ee[64];
    for(int j = 0 ; j < 64 ; j++) ee[j] = &(e[j]);
    int r = kernel((mp_limb_t *) mm, (mp_limb_t **) ee, 64, 64, 64/ULONG_BITS, 64/ULONG_BITS);
    return r;
}

int gauss_6464_imm(mat64 mm, mat64 e, mat64 m)
{
    memcpy(mm,m,sizeof(mat64));
    uint64_t mask=1;
    uint64_t taken=0;
    uint64_t cancelled_cols=0;
    int r = 0;
    for(int j = 0 ; j < 64 ; j++, mask<<=1) e[j]=mask;
    mask=1;
    for(int j = 0 ; j < 64 ; j++, mask<<=1) {
        int k = 0;
        uint64_t z = UINT64_C(1);
        uint64_t pr;
        for(k = 0 ; z && !(((pr=mm[k])&mask) && !(taken&z)); k++, z<<=1) ;
        if (!z) continue;
        taken|=z;
        r++;
        cancelled_cols|=mask;
        uint64_t er = e[k];
        for(k++ ; k < 64 ; k++) {
            uint64_t w = -((mm[k]&mask)!=0);
            mm[k]^=pr&w;
            e[k]^=er&w;
        }
    }
    return r;
}

/* Computes l,u,p, such that:
 *  - l is unit lower triangular
 *  - l*a=u
 *  - up=u*transpose(p) is upper triangular,
 *  -   up[k]==0 iff up[k,k] == 0
 *  -   rank(up)=rank(a)=# of diagonal 1's in up.
 */
int LUP64_imm(mat64 l, mat64 u, mat64 p, mat64 a)
{
    memcpy(u,a,sizeof(mat64));
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
        __m128i * uu = (__m128i*) (u+k);
        __m128i * ll = (__m128i*) (l+k);
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

/* Computes e,mm such that mm=e*m is in row echelon form */
int full_echelon_6464_imm(mat64 mm, mat64 e, mat64 m)
{
    memcpy(mm,m,sizeof(mat64));
    uint64_t mask=1;
    uint64_t cancelled_cols=0;
    int r = 0;
    for(int j = 0 ; j < 64 ; j++, mask<<=1) e[j]=mask;
    mask=1;
    for(int j = 0 ; j < 64 ; j++, mask<<=1) {
        int k = 0;
        uint64_t z = UINT64_C(1);
        uint64_t pr;
        for(k = 0 ; z ; k++, z<<=1) {
            pr=mm[k];
            if ((pr&mask)&&!(pr&(mask-1)))
                break;
        }
        if (!z) continue;
        z=mm[k];mm[k]=mm[j];mm[j]=z;
        z=e[k];e[k]=e[j];e[j]=z;
        r++;
        cancelled_cols|=mask;
        uint64_t er = e[j];
        for(k = 0 ; k < 64 ; k++) {
            if (k==j) continue;
            uint64_t w = -((mm[k]&mask)!=0);
            mm[k]^=pr&w;
            e[k]^=er&w;
        }
    }
    return r;
}

int gauss_128128_C(uint64_t * m)
{
    mat64 mm[4] ATTRIBUTE((aligned(64))); /* handy, even though it does not properly reflect how data is used */
    memcpy(mm,m,4*sizeof(mat64));
    int r = kernel((mp_limb_t*)mm, NULL, 128, 128, 128/ULONG_BITS, 128/ULONG_BITS);
    return r;
}

#if 0
int gauss_128128_imm(uint64_t * m)
{
    mat64 mm[4] ATTRIBUTE((aligned(64)));
    uint64_t * pm = m;
    for(int j = 0 ; j < 64 ; j++, pm+=2) {
        mm[0][j] = pm[0];
        mm[1][j] = pm[1];
    }
    for(int j = 0 ; j < 64 ; j++, pm+=2) {
        mm[2][j] = pm[0];
        mm[3][j] = pm[1];
    }


    mat64 e;
    memcpy(mm,m,sizeof(mat64));
    uint64_t mask=1;
    uint64_t taken=0;
    uint64_t cancelled_cols=0;
    int r = 0;
    for(int j = 0 ; j < 64 ; j++, mask<<=1) e[j]=mask;
    mask=1;
    for(int j = 0 ; j < 64 ; j++, mask<<=1) {
        int k = 0;
        uint64_t z = UINT64_C(1);
        uint64_t pr;
        for(k = 0 ; z && !(((pr=mm[k])&mask) && !(taken&z)); k++, z<<=1) ;
        if (!z) continue;
        taken|=z;
        r++;
        cancelled_cols|=mask;
        uint64_t er = e[k];
        int k0=k;
#define TRIANGULAR_ONLY /* speeds up things by 20 to 25% */
#ifndef  TRIANGULAR_ONLY
        k = 0;
#endif
        for( ; k < 64 ; k++) {
            if (k==k0) continue;
            uint64_t w = -((mm[k]&mask)!=0);
            mm[k]^=pr&w;
            e[k]^=er&w;
        }
    }
    return r;
}
#endif

void pmat_6464(mat64 m)
{
    for(int i = 0; i < 64 ; i++) {
        uint64_t mask=1;
        for(int j = 0 ; j < 64 ; j++, mask<<=1) {
            putchar((m[i]&mask) ?'1':'0');
        }
        putchar('\n');
    }
}

void pmat_mn(mat64 * m, int rb, int cb)
{
    ASSERT_ALWAYS(rb%64 == 0); rb /= 64;
    ASSERT_ALWAYS(cb%64 == 0); cb /= 64;
    for(int ib = 0 ; ib < rb ; ib++) {
        mat64 * mr = m + ib * cb;
        for(int i = 0; i < 64 ; i++) {
            for(int jb = 0 ; jb < cb ; jb++) {
                uint64_t mask=1;
                for(int j = 0 ; j < 64 ; j++, mask<<=1) {
                    putchar((mr[jb][i]&mask) ?'1':'0');
                }
            }
            putchar('\n');
        }
    }
}


/*** table of best functions ***/

void mul_6464_6464(mat64 C, mat64 A, mat64 B)
{
    mul_N64_6464_lookup4(C,A,B,64);
}
void add_6464_6464(mat64_ptr C, mat64_srcptr A, mat64_srcptr B)
{
    add_6464_6464_C(C,A,B);
}

/* Given an N*64 matrix A (N uint64_t's) and a 64*64 matrix B, compute
 * the product.
 *
 * With respect to endianness, we match the column of (A[i]&1)'s with
 * B[0].
 */
void mul_N64_6464(uint64_t *C,
		 const uint64_t *A,
		 const uint64_t *B, size_t m)
{
/* The chosen function is optimal (among the ones here) for N about
 * 20000. At N=2000000, a twice faster version can be obtained. However,
 * it's not critical for cado, so we stick with the slower version.
 */
#if defined(HAVE_SSE2) && ULONG_BITS == 64
    mul_N64_6464_sse(C,A,B,m);
#else
    mul_N64_6464_lookup4(C,A,B,m);
#endif
}
void addmul_N64_6464(uint64_t *C,
		 const uint64_t *A,
		 const uint64_t *B, size_t m)
{
#if defined(HAVE_SSE2) && ULONG_BITS == 64
    addmul_N64_6464_sse(C,A,B,m);
#else
    addmul_N64_6464_lookup4(C,A,B,m);
#endif
}

void mul_N64_T6464(uint64_t *C,
                   const uint64_t *A,
                   const uint64_t *B, size_t m)
{
    mul_N64_T6464_transB(C,A,B,m);
}
void addmul_To64_o64(uint64_t * r, uint64_t a, uint64_t w)
{
#if defined(HAVE_SSE2) && ULONG_BITS == 64
    addmul_To64_o64_lsb_sse_v1(r,a,w);
#else
    addmul_To64_o64_lsb_packof2(r,a,w);
#endif
}
/* Given a 64-bit vector a, and a 64*64 matrix w, compute the product.
 *
 * With respect to endianness, we match a&1 with w[0]
 */
void mul_o64_6464(uint64_t * r, uint64_t a, mat64_srcptr w)
{
    mul_o64_6464_C_lsb(r,a,w);
}
void mul_o64_T6464(uint64_t * w, uint64_t a, mat64_srcptr b)
{
    mul_o64_T6464_C_parity(w,a,b);
}

void mul_TN64_N64(uint64_t * b, uint64_t * A, uint64_t * x, unsigned int ncol)
{
    mul_TN64_N64_C(b, A, x, ncol);
}

void addmul_TN64_N64(uint64_t * b, uint64_t * A, uint64_t * x, unsigned int ncol)
{
    addmul_TN64_N64_C(b, A, x, ncol);
}


/*************/


int mat64_is_uppertriangular(mat64_srcptr u)
{
    uint64_t mask = 1;
    for(int k =0 ; k < 64 ; k++,mask<<=1) {
        if (u[k]&(mask-1)) return 0;
    }
    return 1;
}

int mat64_is_lowertriangular(mat64_srcptr u)
{
    uint64_t mask = -2;
    for(int k =0 ; k < 64 ; k++,mask<<=1) {
        if (u[k]&mask) return 0;
    }
    return 1;
}

int mat64_triangular_is_unit(mat64_srcptr u)
{
    uint64_t mask = 1;
    for(int k =0 ; k < 64 ; k++,mask<<=1) {
        if (!(u[k]&mask)) return 0;
    }
    return 1;
}

int mat64_eq(mat64_srcptr a, mat64_srcptr b)
{
    for(int k =0 ; k < 64 ; k++) {
        if (a[k]!=b[k]) return 0;
    }
    return 1;
}

#ifdef  HAVE_M4RI
static inline void mzd_mypluq(mzd_t * LU, mzd_t * M, mzp_t * P, mzp_t * Q, int c)
{
    mzd_copy(LU, M);
    mzd_pluq(LU,P,Q,c);
}
static inline void mzd_myechelonize_m4ri(mzd_t * E, mzd_t * M, int full, int k)
{
    mzd_copy(E, M);
    mzd_echelonize_m4ri(E,full,k);
}
static inline void mzd_myechelonize_pluq(mzd_t * E, mzd_t * M, int full)
{
    mzd_copy(E, M);
    mzd_echelonize_pluq(E,full);
}
#endif

/*{{{ pmat stuff */
void pmat_init(pmat_ptr x, int n)
{
    x->n=n;
    x->v=(int*) malloc(n*sizeof(int));
    for(int i = 0 ; i < n ; i++) x->v[i]=i;
}
void pmat_clear(pmat_ptr x)
{
    free(x->v); x->v = NULL; x->n = 0;
}
void pmat_transpose(pmat_ptr x, pmat_srcptr y)
{
    if (x == y) {
        pmat t;
        pmat_init(t, y->n);
        pmat_transpose(x, t);
        pmat_clear(t);
        return;
    }
    ASSERT_ALWAYS(x->n == y->n);
    for(int i = 0 ; i < x->n ; i++)
        x->v[i]=-1;
    for(int i = 0 ; i < x->n ; i++) {
        if (y->v[i] >= 0)
            x->v[y->v[i]]=i;
    }
}

void pmat_get_matrix(mat64 * qm, pmat_ptr qp)
{
    int * phi = qp->v;
    int n = qp->n;
    ASSERT_ALWAYS((n%64)==0);
    int nb = n/64;
    memset(qm, 0, n*nb*sizeof(uint64_t));
    uint64_t * qq = (uint64_t*) qm;
    for(int k = 0 ; k < n ; ) {
        for(int jq = 0 ; jq < 64 ; jq++, k++, qq++) {
            int v = phi[k];
            ASSERT_ALWAYS(v >= 0);
            qq[v&~63]=((uint64_t)1)<<(v%64);
        }
        qq+=(nb-1)*64;
    }
}


/* phi_p must point to a zone where at least as many limbs as the number
 * of set bits in the bits[] array
 */
void pmattab_complete(int * phi, uint64_t * bits, int nbits)
{
    ASSERT_ALWAYS(nbits % 64 == 0);
    for(int offset=0 ; offset < nbits ; offset+=64) {
        for(uint64_t w = bits[offset/64], z ; w ; w^=z) {
            z = w^(w&(w-1));
            uint64_t j = cado_ctz64(w);
            *phi++ = offset + j;
        }
    }
}

/* Given phi as in PLUQ64_inner below, return two permutation matrices
 * such that p*u*transpose(q) is a diagonal matrix.
 */
void pqperms_from_phi(pmat_ptr p, pmat_ptr q, int * phi, int m, int n)
{
    int ip = 0;ASSERT_ALWAYS((m%64)==0);ASSERT_ALWAYS(p->n==m);
    int jq = 0;ASSERT_ALWAYS((n%64)==0);ASSERT_ALWAYS(q->n==n);
    uint64_t ms[m/64]; for(int i = 0 ; i < m/64 ; i++) ms[i]=~((uint64_t)0);
    uint64_t ns[n/64]; for(int j = 0 ; j < n/64 ; j++) ns[j]=~((uint64_t)0);
    uint64_t w;
    for(int k = 0 ; k < m ; k++) {
        int v = phi[k];
        if (v < 0) continue;
        p->v[ip++] = k; w=((uint64_t)1)<<(k%64); ms[k/64]^=w;
        q->v[jq++] = v; w=((uint64_t)1)<<(v%64); ns[v/64]^=w;
    }
    pmattab_complete(p->v + ip, ms, m);
    pmattab_complete(q->v + jq, ns, n);
}
/*}}}*/

void mat64_set_identity(mat64_ptr m)
{
    uint64_t mask = 1;
    for(int j = 0 ; j < 64 ; j++, mask<<=1) m[j]=mask;
}
void mat64_copy(mat64_ptr b, mat64_srcptr a)
{
    memcpy(b,a,sizeof(mat64));
}
void mat64_set_zero(mat64_ptr m)
{
    memset(m,0,sizeof(mat64));
}

/* {{{ PLUQ stuff -- well we're not computing exactly PLUQ */

/* Compute matrices l and u such that l*a = u. phi[] is filled with the
 * column indices (shifted by col_offset) of he pivots in u: entry
 * (i,phi(i)-col_offset) is u is one. When row i in u has no important
 * non-zero coefficient, then phi[i] < 0.
 * In column phi[i]-col_offset of u, entries of row index >i are zero.
 */
int PLUQ64_inner(int * phi, mat64 l, mat64 u, mat64 a, int col_offset)
{
    const int m = 64;
    const int n = 64;
    int phi0[64];
    if (phi == NULL) {
        phi = phi0;
        for(int i = 0 ; i < 64 ; i++) phi[i]=-1;
    }

    mat64_copy(u, a);
    mat64_set_identity(l);
    int rank = 0;
    uint64_t todo=~((uint64_t)0);
    for(int i = 0 ; i < m ; i++) {
        uint64_t r=u[i];
        if (phi[i]>=0) continue;
        if (!(r&todo) || phi[i]>=0) continue;
        // this keeps only the least significant bit of r.
        uint64_t v = r^(r&(r-1));
        uint64_t j = cado_ctz64(r);
        phi[i] = col_offset + j;
#if defined(HAVE_SSE41) && !defined(VALGRIND) && !defined(__ICC)
        /* This code is fine, except that there's nowhere we convey the
         * information that we want mat64 objects 16-byte aligned. And
         * with icc, this indeed fails (for test_bitlinalg_matops). We
         * don't know exactly what to do here. 
         */
        int k = i+1;
        if (k&1) {      // alignment call
            uint64_t w = -((u[k]&v)!=0);
            u[k]^=r&w;
            l[k]^=l[i]&w;
            k++;
        }
        /* ok, it's ugly, and requires sse 4.1.
         * but otoh is churns out data veeery fast */
        __m128i vv = _cado_mm_set1_epi64(v);
        __m128i pp = _cado_mm_set1_epi64(r);
        __m128i ee = _cado_mm_set1_epi64(l[i]);
        __m128i * uu = (__m128i*) (u+k);
        __m128i * ll = (__m128i*) (l+k);
        for( ; k < n ; k+=2 ) {
            __m128i ww = _mm_cmpeq_epi64(_mm_and_si128(*uu,vv),vv);
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
int PLUQ64(pmat_ptr p, mat64 * l, mat64 * u, pmat_ptr q, mat64 * m)
{
    int phi[64];
    for(int i = 0 ; i < 64 ; i++) phi[i]=-1;
    int r = PLUQ64_inner(phi,l[0],u[0],m[0],0);
    /* l*m = u */
    /* p*u*transpose(q) = diagonal.
     * p*l*m*transpose(q) = diagonal.
     */
    pqperms_from_phi(p,q,phi,64,64);
    return r;
}

int PLUQ64_n(int * phi, mat64 l, mat64 * u, mat64 * a, int n)
{
    const int m = 64;
    ASSERT_ALWAYS(n % 64 == 0);
    int nb = n/m;
    mat64_set_identity(l);
    for(int i = 0 ; i < 64 ; i++) phi[i]=-1;
    int rank = 0;
    int b = 0;
#ifdef ALLOC_LS
    mat64 ** ls = (mat64**) malloc(nb * sizeof(mat64*));
#else
    mat64 ls[nb] ATTRIBUTE((aligned(64)));
#endif
    mat64 tl;
    for( ; b < nb && rank < m ; b++) {
        mat64 ta;
        mul_6464_6464(ta, l, a[b]);
        rank += PLUQ64_inner(phi, tl, u[b], ta, b*m);
        mul_6464_6464(l, tl, l);
#ifdef  ALLOC_LS
        ls[b]=(mat64*) malloc_aligned(sizeof(mat64), 64);
        mat64_copy(*ls[b], tl);
#else
        mat64_copy(ls[b], tl);
#endif
    }
    int nspins = b;
#ifdef  ALLOC_LS
    free(ls[b-1]);
    for(int c = b-2 ; c >= 0 ; c--) {
        mul_6464_6464(u[c], tl, u[c]);
        mul_6464_6464(tl, *ls[c], tl);
        free_aligned(ls[c]);
    }
    free(ls);
#else
    for(int c = b-2 ; c >= 0 ; c--) {
        mul_6464_6464(u[c], tl, u[c]);
        mul_6464_6464(tl, ls[c], tl);
    }
#endif
    for( ; b < nb ; b++)
        mul_6464_6464(u[b], l, a[b]);
    return nspins*m+b;
}

static inline void bli_64x64N_clobber(mat64 h, mat64 * us, int * phi, int nb)
{
    /* problem: we're modifying U here. So either we do a copy of U,
     * which can be probelmatic memory-wise, or we do an extraction ;
     * it's also possible to rebuild the original U from the extracted H
     * and U', merely with a product (_if ever_ we care about U, in
     * fact). However this latter option seems messy.
     */
    mat64_set_identity(h);
    const int m = 64;
    for(int i = 0 ; i < m ; i++) {
        int j = phi[i];
        if (j<0) continue;
        uint64_t m = ((uint64_t)1) << (j%64);
        int d = j/64;
        /* TODO: use _mm_cmpeq_epi64 for this as well, of course */
        ASSERT(us[d][i]&m);
        int k = 0;
#if defined(HAVE_SSE41) && !defined(VALGRIND)
        __m128i mm = _cado_mm_set1_epi64(m);
        __m128i * uu = (__m128i*) us[d];
        __m128i * hh = (__m128i*) h;
        __m128i hi = _cado_mm_set1_epi64(h[i]);
        int ii=i/2;
        for( ; k < ii ; k++) {
            __m128i ww = _mm_cmpeq_epi64(_mm_and_si128(*uu++,mm),mm);
            for(int b = 0 ; b < nb ; b++) {
                // ((__m128i*)us[b])[k] ^= ww & _cado_mm_set1_epi64(us[b][i]);
                __m128i * z = ((__m128i*)us[b]) + k;
                *z = _mm_xor_si128(*z, _mm_and_si128(ww, _cado_mm_set1_epi64(us[b][i])));
            }
            hh[k] = _mm_xor_si128(hh[k], _mm_and_si128(ww, hi));
        }
        k*=2;
#endif
        for( ; k < i ; k++) {
            uint64_t w = -((us[d][k]&m) != 0);
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
void bli_64x128(mat64 h, mat64 * us, int * phi)
{
    mat64 uc[2] ATTRIBUTE((aligned(64)));
    memcpy(uc,us,2*sizeof(mat64));
    bli_64x64N_clobber(h,uc,phi,2);
}

void extract_cols_64_from_128(mat64 t, mat64 * m, int * phi)
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
    memset(t, 0, sizeof(mat64));
    uint64_t mask = 1;
    for(int j = 0 ; j < 64 ; j++, mask<<=1) {
#if defined(HAVE_SSE41) && !defined(VALGRIND) && !defined(__ICC)
        /* This code is fine, except that there's nowhere we convey the
         * information that we want mat64 objects 16-byte aligned. And
         * with icc, this indeed fails (for test_bitlinalg_matops). We
         * don't know exactly what to do here. 
         */
        __m128i ss[2] = {
            _cado_mm_set1_epi64(s[0][j]),
            _cado_mm_set1_epi64(s[1][j]) };
        __m128i * mm[2] = {(__m128i*)m[0],(__m128i*)m[1]};
        __m128i * tt = (__m128i*)t;
        __m128i mmk = _cado_mm_set1_epi64(mask);
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
int PLUQ128(pmat_ptr p, mat64 * l, mat64 * u, pmat_ptr q, mat64 * m)
{
    /* This is really an outer routine. An inner routine will not have p
     * and q, but rather both merged as a phi argument, in the manner of
     * PLUQ64_inner (and of course the following lines would be changed a
     * bit).
     */
    int phi[128];
    for(int i = 0 ; i < 128 ; i++) phi[i]=-1;

    mat64_set_zero(l[0]);
    mat64_set_zero(l[1]);
    mat64_set_zero(l[2]);
    mat64_set_identity(l[3]);

    int r1 = PLUQ64_n(phi,l[0],u,m,128);
    r1 = r1 % 64;

    // andouille 7.65

    /* l[0] * m = u */

    mat64 h;
    bli_64x128(h, u, phi);
    /* h * u is "sort of" identity, at least up to permutation */

    mat64 l21;
    mat64_ptr S = l21;

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
    mul_6464_6464(t2[0], l21, u[0]); add_6464_6464(t2[0], m[2], t2[0]);
    mul_6464_6464(t2[1], l21, u[1]); add_6464_6464(t2[1], m[3], t2[1]);

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

/* }}} */

