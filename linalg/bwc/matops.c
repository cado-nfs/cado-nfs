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

/*
 * Warning: this code is specific to x86_64 architecture.
 */


#include "cado.h"       /* HAVE_* macros ! */

#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <emmintrin.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <inttypes.h>


#ifdef  HAVE_SSE41
#include <smmintrin.h>  // sse 4.1 _mm_cmpeq_epi64
#endif  /* HAVE_SSE41 */


/* Define this to check these functions against the m4ri library */
#define xxxHAVE_M4RI
#ifdef  HAVE_M4RI
#include "m4ri/m4ri.h"
#endif  /* HAVE_M4RI */
#include "gauss.h"
#include "macros.h"

/* This is **only** for 64 * 64 matrices */

#define WBITS   64
typedef uint64_t mat64[64];
typedef uint64_t * mat64_ptr;
typedef const uint64_t * mat64_srcptr;

#if 0/*{{{*/
static int urandom_fd = -1;

static inline uint64_t rand64()
{
    if (urandom_fd == -1)
        urandom_fd = open("/dev/urandom", O_RDONLY);
    if (urandom_fd < 0) abort();
    uint64_t r;
    size_t s = read(urandom_fd, &r, sizeof(uint64_t));
    if (s != sizeof(uint64_t)) abort();
    return r;
}

static inline void rand64_mem(uint64_t * ptr, size_t z)
{
    if (urandom_fd == -1)
        urandom_fd = open("/dev/urandom", O_RDONLY);
    if (urandom_fd < 0) abort();
    size_t s = read(urandom_fd, ptr, z * sizeof(uint64_t));
    if (s != z) abort();
}
#else
static inline uint64_t rand64()
{
    uint64_t r0 = rand();
    uint64_t r1 = rand(); r1 <<= 22;
    uint64_t r2 = rand(); r2 <<= 44;
    return r0 ^r1 ^ r2;
}
static inline void rand64_mem(uint64_t * ptr, size_t z)
{
    for( ; z-- ; ptr[z] = rand64());
}
#endif/*}}}*/

static inline uint64_t bitrev(uint64_t a)/*{{{*/
{
    a = (a >> 32) ^ (a << 32);
    uint64_t m;
    m = 0x0000ffff0000ffffUL;
    a = ((a >> 16) & m) ^ ((a << 16) & ~m);
    m = 0x00ff00ff00ff00ffUL;
    a = ((a >> 8) & m) ^ ((a << 8) & ~m);
    m = 0x0f0f0f0f0f0f0f0fUL;
    a = ((a >> 4) & m) ^ ((a << 4) & ~m);
    m = 0x3333333333333333UL;
    a = ((a >> 2) & m) ^ ((a << 2) & ~m);
    m = 0x5555555555555555UL;
    a = ((a >> 1) & m) ^ ((a << 1) & ~m);
    return a;
}
/* like bitrev, but keep nibbles intact */
static inline uint64_t nibrev(uint64_t a)
{
    a = (a >> 32) ^ (a << 32);
    uint64_t m;
    m = 0x0000ffff0000ffffUL;
    a = ((a >> 16) & m) ^ ((a << 16) & ~m);
    m = 0x00ff00ff00ff00ffUL;
    a = ((a >> 8) & m) ^ ((a << 8) & ~m);
    m = 0x0f0f0f0f0f0f0f0fUL;
    a = ((a >> 4) & m) ^ ((a << 4) & ~m);
    return a;
}/*}}}*/

/* prototypes for the functions defined */
static inline void mul_6464_6464(mat64 C, mat64 A, mat64 B);
static inline void add_6464_6464(mat64_ptr C, mat64_srcptr A, mat64_srcptr B);
static inline void mul_o64_6464(uint64_t *r, uint64_t a, mat64_srcptr w);
static inline void mul_N64_6464(uint64_t *C, const uint64_t *A,
		 const uint64_t *B, unsigned long m);
static inline void mul_N64_T6464(uint64_t *C, const uint64_t *A,
                   const uint64_t *B, unsigned long m);
static inline void addmul_To64_o64(uint64_t * r, uint64_t a, uint64_t w);
static inline void mul_o64_6464(uint64_t * r, uint64_t a, mat64_srcptr w);
static inline void mul_o64_T6464(uint64_t * w, uint64_t a, mat64_srcptr b);


/* level 1 */
#ifdef  HAVE_SSE2
static inline void mul_6464_6464_sse(mat64_ptr C, mat64_srcptr A, mat64_srcptr B)
{
    int i;
    memset(C, 0, sizeof(mat64));
 
    /* As per emmintrin.h, only the __m128i type is declared with
     * attribute __may_alias__, meaning that accesses through this
     * pointer may alias other types (e.g. uint64_t's, in our cases). For
     * this reason, we _must_ use this pointer type, and not __v2di, for
     * our pointer types. Reading from a __m128 * into a __v2di, or the
     * converse, are legal operations.
     */
    __m128i *Cw = (__m128i *) C;
    __m128i *Aw = (__m128i *) A;

    for (int j = 0; j < 64; j += 2) {
	__v2di c = { 0, 0 };
	__v2di a = *Aw++;

	__v2di one = { 1, 1, };
#define SHR(x,r) _mm_srli_epi64((x),(r))
	for (i = 0; i < 64; i++) {
	    __v2di bw = { B[i], B[i], };

	    c ^= (bw & -(a & one));
	    a = SHR(a, 1);
	}
#undef  SHR
	*Cw++ = c;
    }
}
#endif

static inline void mul_6464_6464_v2(mat64_ptr C, mat64_srcptr A, mat64_srcptr B)
{
    memset(C, 0, sizeof(mat64));
 
    for (int i = 0; i < 64; i++) {
    for (int j = 0; j < 64; j++) {
        if ((A[i]>>j)&1)
            C[i]^=B[j];
    }
    }
}

static inline void add_6464_6464_C(mat64_ptr C, mat64_srcptr A, mat64_srcptr B)
{
    for (int j = 0; j < 64; j++) {
        C[j] = A[j] ^ B[j];
    }
}


static inline void addmul_To64_o64_lsb(uint64_t * r, uint64_t a, uint64_t w)
{
    /* Dans un sens */
    for (unsigned int i = 0; i < 64; i++) {
	*r++ ^= w & -(a & 1);
	a >>= 1;
    }
}

static inline void addmul_To64_o64_msb(uint64_t * r, uint64_t a, uint64_t w)
{
    /* Dans l'autre -- va un poil plus vite. */
    for (unsigned int i = 0; i < 64; i++) {
	*r++ ^= w & (((int64_t) a) >> 63);
	a <<= 1;
    }
}

static inline void addmul_To64_o64_lsb_packof2(uint64_t * r, uint64_t a, uint64_t w)
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
#ifdef  HAVE_SSE2
static inline void addmul_To64_o64_lsb_sse_v1(uint64_t * r, uint64_t a, uint64_t w)
{
    /* Avec des sse-2 */
    __v2di mb[4] = {
	(__v2di) {0, 0},
	(__v2di) {w, 0},
	(__v2di) {0, w},
	(__v2di) {w, w},
    };
    __m128i *sr = (__m128i *) r;
    for (int i = 0; i < 64; i += 2) {
	*sr++ ^= mb[a & 3];
	a >>= 2;
    }
}
#endif

static inline void mul_o64_6464_C_lsb(uint64_t * r, uint64_t a, mat64_srcptr w)
{
    uint64_t c = 0;
    for (unsigned int i = 0; i < WBITS; i++) {
	c ^= (w[i] & -(a & 1UL));
	a >>= 1UL;
    }
    *r = c;
}

static inline void mul_o64_6464_C_msb(uint64_t *r,
                   uint64_t a,
                   mat64_srcptr w)
{
    uint64_t c = 0UL;
    for (int i = 64 - 1; i >= 0; i--) {
        c ^= (w[i] & (((int64_t) a) >> (64 - 1)));
        a <<= 1UL;
    }
    *r = c;
}


static inline void mul_o64_T6464_C_parity(uint64_t * w, uint64_t a, mat64_srcptr b)
{
    // Uses unoptimized __builtin_parityl function -- maybe better with gcc 4.3
    // note that popcnt is faster in asm than the more restricted parity
    // functions. So if it's available, it should be tested.
    uint64_t c = 0;
    for (unsigned int i = 0; i < WBITS; i++) {
        uint64_t p = __builtin_parityl(a & b[i]);
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

static inline void mul_o64_T6464_C_parity3(uint64_t * w, uint64_t a, mat64_srcptr b)
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


static inline void transp_6464(mat64_ptr dst, mat64_srcptr src)
{
    int i, j;
    for (i = 0; i < 64; i++) {
	dst[i] = 0;
	for (j = 0; j < 64; j++) {
	    dst[i] ^= ((src[j] >> i) & 1UL) << j;
	}
    }
}

static inline void copy_6464(mat64_ptr dst, mat64_srcptr src)
{
    memcpy(dst, src, sizeof(mat64));
}


/* level 2 */

static inline void copy_N64(uint64_t * dst, const uint64_t * src, unsigned long m)
{
    memcpy(dst, src, m * sizeof(uint64_t));
}

static inline int cmp_N64(const uint64_t * dst, const uint64_t * src, unsigned long m)
{
    return memcmp(dst, src, m * sizeof(uint64_t));
}


/* This can work in place (C==A, or C==B, or both) */
static inline void mul_N64_6464_lookup4(uint64_t *C,
                   const uint64_t *A,
                   const uint64_t *B, unsigned long m)
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
    for (unsigned long i = 0; i < m; i++) {
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
static inline void addmul_N64_6464_lookup4(uint64_t *C,
                   const uint64_t *A,
                   const uint64_t *B, unsigned long m)
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
    memset(C, 0, m * sizeof(uint64_t));
    for (unsigned long i = 0; i < m; i++) {
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

static inline void mul_N64_6464_lookup8(uint64_t *C,
                   const uint64_t *A,
                   const uint64_t *B, unsigned long m)
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
    memset(C, 0, m * sizeof(unsigned long));
    for (unsigned long i = 0; i < m; i++) {
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

static inline void mul_N64_6464_vec(uint64_t *C,
                   const uint64_t *A,
                   const uint64_t *B, unsigned long m)
{

    memset(C, 0, m * sizeof(uint64_t));
    for (unsigned long i = 0; i < m; i++) {
        mul_o64_6464(C++, *A++, B);
    }
}

static inline void mul_N64_T6464_vec(uint64_t *C,
                   const uint64_t *A,
                   const uint64_t *B, unsigned long m)
{

    memset(C, 0, m * sizeof(uint64_t));
    for (unsigned long i = 0; i < m; i++) {
        mul_o64_T6464(C++, *A++, B);
    }
}

static inline void mul_N64_6464_transB(uint64_t *C,
                   const uint64_t *A,
                   const uint64_t *B, unsigned long m)
{
    uint64_t *tb = malloc(64 * sizeof(uint64_t));
    transp_6464(tb, B);
    mul_N64_T6464(C, A, tb, m);
    free(tb);
}

static inline void mul_N64_T6464_transB(uint64_t *C,
                   const uint64_t *A,
                   const uint64_t *B, unsigned long m)
{
    uint64_t *tb = malloc(64 * sizeof(uint64_t));
    transp_6464(tb, B);
    mul_N64_6464(C, A, tb, m);
    free(tb);
}

#ifdef  HAVE_SSE2
static inline void mul_N64_6464_sse(uint64_t *C,
		 const uint64_t *A,
		 const uint64_t *B, unsigned long m)
{
    unsigned long j;
    memset(C, 0, m * sizeof(uint64_t));

    __m128i *Cw = (__m128i *) C;
    __m128i *Aw = (__m128i *) A;

    for (j = 0; j < m; j += 2) {
	__v2di c = { 0, 0 };
	__v2di a = *Aw++;

	__v2di one = { 1, 1, };
#define SHR(x,r) _mm_srli_epi64((x),(r))
	for (int i = 0; i < 64; i++) {
	    __v2di bw = { B[i], B[i], };

	    c ^= (bw & -(a & one));
	    a = SHR(a, 1);
	}
#undef  SHR
	*Cw++ = c;
    }
    C += j;
    A += j;
    for (; j < m; j++) {
	uint64_t c = 0UL;
	uint64_t a = *A++;
	for (int i = 0; i < 64; i++) {
	    c ^= (B[i] & -(a & 1UL));
	    a >>= 1UL;
	}
	*C++ = c;
    }
}
#endif

static inline void mul_64N_N64_addmul(uint64_t *r, uint64_t *a, uint64_t *w, unsigned long n)
{
    memset(r, 0, 64 * sizeof(uint64_t));
    for (unsigned long i = 0; i < n; i++) {
        addmul_To64_o64(r, a[i], w[i]);
    }
}


static inline void mul_TN32_N64_C(uint64_t * b, uint32_t * A, uint64_t * x, unsigned int ncol)
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

static inline void mul_TN64_N64_C(uint64_t * b, uint64_t * A, uint64_t * x, unsigned int ncol)
{
    uint64_t idx, i, rA;
    uint64_t rx;

    memset(b, 0, 64 * sizeof(uint64_t));
    for(idx = 0; idx < ncol; idx++) {
        rA = A[idx];
        rx = x[idx];
        for(i = 0; i < 64; i++) {
            b[i] ^= rx & -(rA & 1);
            rA >>= 1;
        }
    }
}

#ifdef  HAVE_SSE2
static inline void mul_TN64K_N64_sse2(uint64_t * w, uint64_t * u, uint64_t * v, unsigned int n, unsigned int K)
{
    memset(w, 0, 64 * K * sizeof(uint64_t));
    for(unsigned int i = 0 ; i < n ; i++) {
        __v2di * w0 = (__v2di*) w;
        // TODO: It's possible to expand more, and use a __v2di
        // mb[4][2], or even [4]. This wouldn't change the code much
        // (see the u128 version), and is likely to speed things up a
        // wee bit maybe.
        __v2di mb[4] = {
            (__v2di) {0, 0},
            (__v2di) {*v, 0},
            (__v2di) {0, *v},
            (__v2di) {*v, *v},
        };
        v++;
        __v2di *sw = w0;
        for(unsigned int k = 0 ; k < K ; k++) {
            uint64_t a = *u++;
            for (unsigned int j = 0; j < 64; j += 2) {
                *sw ^= mb[a & 3];
                a >>= 2;
                sw ++;
            }
        }
    }
}
#endif

static inline void mul_TN64K_N64_C(uint64_t * b, uint64_t * A, uint64_t * x, unsigned int ncol, unsigned int K)
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
	    //if (P & 1UL) C[k] ^= B[i];
	    C[k] ^= (B[i] & -(P & 1UL));
	    P >>= 1UL;
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
	    C[i] ^= (B[k] & -(P & 1UL));
	    P >>= 1UL;
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
	    j ^= (B[k] & -(P & 1UL));
	    P >>= 1UL;
	}
	C[i] = j;
    }
}
#endif
#endif/*}}}*/


/* polynomials */

/* lengths of a1 and a2 are n1 and n2 */
typedef uint64_t (*m64pol_ptr)[64];
typedef uint64_t (*const m64pol_srcptr)[64];

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
    if (n1 == 1) {
        m64pol_mul(r, a1, a2, n1, n2);
        return;
    }
    unsigned int h = n1 >> 1;
    memset(r, 0, (n1 + n2 - 1) * sizeof(mat64));

    m64pol_add(r, a1, a1 + h, h);
    m64pol_add(r + 2 * h, a2, a2 + h, h);

    mat64 * t = (mat64 *) malloc((2*h-1) * sizeof(mat64));
    m64pol_mul_kara(t, r, r + 2 * h, h, h);

    m64pol_mul_kara(r, a1, a2, h, h);
    m64pol_mul_kara(r + 2 * h, a1 + h, a2 + h, h, h);
    m64pol_add(t, t, r, 2 * h - 1);
    m64pol_add(t, t, r + 2 * h, 2 * h - 1);
    m64pol_add(r + h, r + h, t, 2 * h - 1);
    free(t);
}

void m64pol_addmul_kara(m64pol_ptr r, m64pol_srcptr a1, m64pol_srcptr a2, unsigned int n1, unsigned int n2)
{
    mat64 * t = (mat64 *) malloc((n1 + n2 -1) * sizeof(mat64));
    m64pol_mul_kara(t, a1, a2, n1, n2);
    m64pol_add(r, r, t, n1 + n2 - 1);
    free(t);
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

#define t_and_unit_from_clock__bare(t, unit, t1, j)                     \
    double t = t1;							\
    t /= CLOCKS_PER_SEC;						\
    if (j) t /= j; else t = 0;						\
    char * unit = "s";							\
    if (t < 1.0e-7) { unit = "ns"; t *= 1.0e9;			        \
    } else if (1 || t < 1.0e-4) { unit = "micros"; t *= 1.0e6;		\
    } else if (t < 1.0e-1) { unit = "ms"; t *= 1.0e3; }                 \
    do { } while (0)

#define TIME1__bare(maxtime, what, args) 		        	\
    clock_t measuring_time = maxtime * CLOCKS_PER_SEC;			\
    clock_t t0, t1;							\
    int j;								\
    t0 = clock();							\
    for (j = 0; ; j++) {						\
        what args;							\
        t1 = clock() - t0;						\
        if (j && t1 > measuring_time)					\
            break;							\
    }									\
    t_and_unit_from_clock__bare(t, unit, t1, j);

#define TIME1(maxtime, what, args) do {			        	\
    TIME1__bare(maxtime, what, args)                                    \
    printf(#what " \t%d times in %.4f %s each\n",       		\
            j, t, unit);		                        	\
} while (0)

#define TIME1N(maxtime, what, args) do {		        	\
    TIME1__bare(maxtime, what, args)                                    \
    printf(#what "(n=%d) \t%d times in %.4f %s each\n", n,     	        \
            j, t, unit);		                        	\
} while (0)

#define TIME1N_spins(rexpr, maxtime, what, args, spinexpr, spinmax) do {        \
    clock_t ts[spinmax];						\
    int ns[spinmax];                                                    \
    for(int s = 0 ; s < spinmax ; s++) ts[s] = ns[s] = 0;		\
    clock_t t0, t1;							\
    int j;								\
    t0 = clock();							\
    clock_t t = t0;                                                     \
    clock_t fence = t0 + maxtime * CLOCKS_PER_SEC;			\
    for (j = 0; ; j++) {						\
        rexpr;                                                          \
        int ret = what args;						\
        int s = spinexpr;                                               \
        ts[s] += (t1 = clock()) - t;                                    \
        ns[s] ++;                                                       \
        t = t1;                                                         \
        if (j && t1 > fence)					        \
            break;							\
    }									\
    int nch=0;                                                          \
    for(int s = 0 ; s < spinmax ; s++) {				\
        if (s == 0) nch=printf(#what"(n=%d)", n);			\
        else for(int k = nch ; k-- ; putchar(' '));		        \
        t_and_unit_from_clock__bare(t, unit, ts[s], ns[s]);		\
        printf(" \t[%d] %d times in %.4fs %s each\n",s,ns[s],t,unit);	\
    }									\
} while (0)



/* Same spirit, but treat multiplication of 64K by 64K matrices (of
 * polynomials).
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
        ptr[0] = bitrev(s[i]);
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
        if (ptr[0] != bitrev(s[i])) {
            fprintf(stderr, "Rows %zu differ: %016"PRIx64" != %016"PRIx64"\n",
                    i, bitrev(ptr[0]), s[i]);
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
    int r = kernel(mm, ee, 64, 64, 1, 1);
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
        uint64_t z = 1UL;
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
    uint64_t todo=~((uint64_t)0);
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
#if defined(HAVE_SSE41) && !defined(VALGRIND)
        int k = j+1;
        if (k&1) {      // alignment call
            uint64_t w = -((u[k]&v)!=0);
            u[k]^=pr&w;
            l[k]^=l[j]&w;
            k++;
        }
        /* ok, it's ugly, and requires sse 4.1.
         * but otoh is churns out data veeery fast */
        __m128i vv = (__v2di) { v,v };
        __m128i pp = (__v2di) { pr, pr };
        __m128i ee = (__v2di) { l[j], l[j] };
        __m128i * uu = (__m128i*) (u+k);
        __m128i * ll = (__m128i*) (l+k);
        for( ; k < 64 ; k+=2 ) {
            __v2di ww = _mm_cmpeq_epi64(*uu&vv,vv);
            *uu++ ^= pp & ww;
            *ll++ ^= ee & ww;
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
        uint64_t z = 1UL;
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
    mat64 mm[4]; /* handy, even though it does not properly reflect how data is used */
    memcpy(mm,m,4*sizeof(mat64));
    int r = kernel((uint64_t*)mm, NULL, 128, 128, 2, 2);
    return r;
}

#if 0
int gauss_128128_imm(uint64_t * m)
{
    mat64 mm[4];
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
        uint64_t z = 1UL;
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

static inline void mul_6464_6464(mat64 C, mat64 A, mat64 B)
{
    mul_N64_6464_lookup4(C,A,B,64);
}
static inline void add_6464_6464(mat64_ptr C, mat64_srcptr A, mat64_srcptr B)
{
    add_6464_6464_C(C,A,B);
}
static inline void mul_N64_6464(uint64_t *C,
		 const uint64_t *A,
		 const uint64_t *B, unsigned long m)
{
#ifdef  HAVE_SSE2
    mul_N64_6464_sse(C,A,B,m);
#else
    mul_N64_6464_lookup4(C,A,B,m);
#endif
}
static inline void mul_N64_T6464(uint64_t *C,
                   const uint64_t *A,
                   const uint64_t *B, unsigned long m)
{
    mul_N64_T6464_transB(C,A,B,m);
}
static inline void addmul_To64_o64(uint64_t * r, uint64_t a, uint64_t w)
{
#ifdef HAVE_SSE2
    addmul_To64_o64_lsb_sse_v1(r,a,w);
#else
    addmul_To64_o64_lsb_packof2(r,a,w);
#endif
}
static inline void mul_o64_6464(uint64_t * r, uint64_t a, mat64_srcptr w)
{
    mul_o64_6464_C_lsb(r,a,w);
}
static inline void mul_o64_T6464(uint64_t * w, uint64_t a, mat64_srcptr b)
{
    mul_o64_T6464_C_parity(w,a,b);
}


/* {{{ test routines */
struct l1_data_s {
    uint64_t * xr;
    uint64_t * r;
    uint64_t * a;
    uint64_t * w;
    uint64_t * wt;
#ifdef  HAVE_M4RI
    mzd_t *R;
    mzd_t *A;
    mzd_t *W;
    mzd_t *WT;
#endif  /* HAVE_M4RI */
    unsigned int n;
};

typedef struct l1_data_s l1_data[1];
typedef struct l1_data_s * l1_data_ptr;
typedef const struct l1_data_s * l1_data_srcptr;

void l1_data_init_set(l1_data_ptr D, unsigned int n)
{
    D->n = n;
    D->xr = (uint64_t *) malloc(n * sizeof(uint64_t));
    D->r = (uint64_t *) malloc(n * sizeof(uint64_t));
    D->a = (uint64_t *) malloc(n * sizeof(uint64_t));
    D->w = (uint64_t *) malloc(64 * sizeof(uint64_t));
    D->wt = (uint64_t *) malloc(64 * sizeof(uint64_t));
#ifdef  HAVE_M4RI
    D->R = mzd_init(n, 64);
    D->A = mzd_init(n, 64);
    D->W = mzd_init(64, 64);
    D->WT = mzd_init(64, 64);
#endif  /* HAVE_M4RI */

    rand64_mem(D->a, n);
    rand64_mem(D->w, 64);
    transp_6464(D->wt, D->w);
#ifdef  HAVE_M4RI
    mzd_set_mem(D->A, D->a, n);
    mzd_set_mem(D->W, D->w, 64);
    mzd_set_mem(D->WT, D->wt, 64);
#endif  /* HAVE_M4RI */
}

void l1_data_clear(l1_data_ptr D)
{
    free(D->xr);
    free(D->r);
    free(D->a);
    free(D->w);
    free(D->wt);
#ifdef  HAVE_M4RI
    mzd_free(D->R);
    mzd_free(D->A);
    mzd_free(D->W);
    mzd_free(D->WT);
#endif  /* HAVE_M4RI */
}

#define SCOPE_L1_DATA_MEMBERS_STANDARD(D)				\
    unsigned int n __attribute__((unused)) = D->n;			\
    uint64_t * xr __attribute__((unused)) = D->xr;      		\
    uint64_t * r __attribute__((unused)) = D->r;			\
    const uint64_t * a __attribute__((unused)) = D->a;			\
    const uint64_t * w __attribute__((unused)) = D->w;			\
    const uint64_t * wt __attribute__((unused)) = D->wt

#define SCOPE_L1_DATA_MEMBERS_M4RI(D)   				\
    mzd_t * R __attribute__((unused)) = D->R;				\
    mzd_t * A __attribute__((unused)) = D->A;				\
    mzd_t * W __attribute__((unused)) = D->W;				\
    mzd_t * WT __attribute__((unused)) = D->WT

#ifdef  HAVE_M4RI
#define SCOPE_L1_DATA_MEMBERS(D)                \
        SCOPE_L1_DATA_MEMBERS_STANDARD(D);      \
        SCOPE_L1_DATA_MEMBERS_M4RI(D)
#else  /* HAVE_M4RI */
#define SCOPE_L1_DATA_MEMBERS(D)                \
        SCOPE_L1_DATA_MEMBERS_STANDARD(D)
#endif  /* HAVE_M4RI */

void level1_basic_tests()
{
    l1_data D;
    l1_data_init_set(D, 64);
    SCOPE_L1_DATA_MEMBERS(D);

    /* transposition */
    transp_6464(r, a); memcpy(xr, r, 64 * sizeof(uint64_t));

    transp_6464(r, a);
    if (memcmp(xr, r, 64 * sizeof(uint64_t))) abort();
    TIME1(1, transp_6464, (r,a));

#ifdef  HAVE_M4RI
    mzd_transpose(R, A);
    mzd_check_mem(R, xr, 64);
    TIME1(1, mzd_transpose, (R, A));
#endif  /* HAVE_M4RI */

    /* copy */
    copy_6464(r, a); memcpy(xr, r, 64 * sizeof(uint64_t));

    copy_6464(r, a);
    if (memcmp(xr, r, 64 * sizeof(uint64_t))) abort();
    TIME1(1, copy_6464, (r,a));

#ifdef  HAVE_M4RI
    mzd_copy(R, A);
    mzd_check_mem(R, xr, 64);
    TIME1(1, mzd_copy, (R, A));
#endif  /* HAVE_M4RI */

    /* add */
    add_6464_6464(r, a, w); memcpy(xr, r, 64 * sizeof(uint64_t));

    add_6464_6464_C(r, a, w);
    if (memcmp(xr, r, 64 * sizeof(uint64_t))) abort();
    TIME1(1, add_6464_6464_C, (r,a,w));

#ifdef  HAVE_M4RI
    mzd_add(R, A, W);
    mzd_check_mem(R, xr, 64);
    TIME1(1, mzd_add, (R, A, W));
#endif  /* HAVE_M4RI */

    /*******/

    l1_data_clear(D);
}

void level1_mul_tests_N_list(l1_data_ptr D)
{
    SCOPE_L1_DATA_MEMBERS(D);

    mul_N64_6464_vec(r, a, w, n);
    if (memcmp(xr, r, n * sizeof(uint64_t))) abort();
    TIME1N(2, mul_N64_6464_vec, (r,a,w,n));

    mul_N64_6464_transB(r, a, w, n);
    if (memcmp(xr, r, n * sizeof(uint64_t))) abort();
    TIME1N(2, mul_N64_6464_transB, (r,a,w,n));

#ifdef  HAVE_SSE2
    mul_N64_6464_sse(r, a, w, n);
    if (memcmp(xr, r, n * sizeof(uint64_t))) abort();
    TIME1N(2, mul_N64_6464_sse, (r,a,w,n));
#endif

    mul_N64_6464_lookup4(r, a, w, n);
    if (memcmp(xr, r, n * sizeof(uint64_t))) abort();
    TIME1N(2, mul_N64_6464_lookup4, (r,a,w,n));

    mul_N64_6464_lookup8(r, a, w, n);
    if (memcmp(xr, r, n * sizeof(uint64_t))) abort();
    TIME1N(2, mul_N64_6464_lookup8, (r,a,w,n));

    mul_N64_T6464_vec(r, a, wt, n);
    if (memcmp(xr, r, n * sizeof(uint64_t))) abort();
    TIME1N(2, mul_N64_T6464_vec, (r,a,w,n));

    mul_N64_T6464_transB(r, a, wt, n);
    if (memcmp(xr, r, n * sizeof(uint64_t))) abort();
    TIME1N(2, mul_N64_T6464_transB, (r,a,w,n));

#ifdef  HAVE_M4RI
    mzd_mul_naive(R, A, W);
    mzd_check_mem(R, xr, n);
    TIME1(1, mzd_mul_naive, (R, A, W));

    _mzd_mul_naive(R, A, WT, 1);
    mzd_check_mem(R, xr, n);
    TIME1(1, _mzd_mul_naive, (R, A, WT, 1));

    mzd_mul_m4rm(R, A, W, 0);
    mzd_check_mem(R, xr, n);
    TIME1(1, mzd_mul_m4rm, (R, A, W, 0));

    _mzd_mul_m4rm(R, A, W, 0, 1);
    mzd_check_mem(R, xr, n);
    TIME1(1, _mzd_mul_m4rm, (R, A, W, 0, 1));
#endif  /* HAVE_M4RI */
}


void level1_mul_tests_N(unsigned int n)
{
    l1_data D;
    l1_data_init_set(D, n);

    /* multiplication of a by a matrix */
    mul_N64_6464_vec(D->r, D->a, D->w, D->n);
    memcpy(D->xr, D->r, D->n * sizeof(uint64_t));

    level1_mul_tests_N_list(D);

    l1_data_clear(D);
}

void level1_mul_tests_1_list(l1_data_ptr D)
{
    SCOPE_L1_DATA_MEMBERS(D);
    assert(n == 1);

    mul_o64_6464_C_lsb(r, *a, w);
    if (memcmp(xr, r, n * sizeof(uint64_t))) abort();
    TIME1(1, mul_o64_6464_C_lsb, (r, *a, w));

    mul_o64_6464_C_msb(r, *a, w);
    if (memcmp(xr, r, n * sizeof(uint64_t))) abort();
    TIME1(1, mul_o64_6464_C_msb, (r, *a, w));

    mul_o64_T6464_C_parity(r, *a, wt);
    if (memcmp(xr, r, n * sizeof(uint64_t))) abort();
    TIME1(1, mul_o64_T6464_C_parity, (r, *a, wt));

    mul_o64_T6464_C_parity3(r, *a, wt);
    if (memcmp(xr, r, n * sizeof(uint64_t))) abort();
    TIME1(1, mul_o64_T6464_C_parity3, (r, *a, wt));
}

void level1_mul_tests_1()
{
    unsigned int n = 1;
    l1_data D;
    l1_data_init_set(D, n);

    /* multiplication of a by a matrix */
    mul_o64_6464(D->r, *D->a, D->w);
    memcpy(D->xr, D->r, D->n * sizeof(uint64_t));

    level1_mul_tests_1_list(D);
    /* Functions which can do any n can also do n=1 */
    level1_mul_tests_N_list(D);

    l1_data_clear(D);
}
void level1_mul_tests_64_list(l1_data_ptr D)
{
    SCOPE_L1_DATA_MEMBERS(D);
    assert(n == 64);

#ifdef  HAVE_SSE2
    mul_6464_6464_sse(r, a, w);
    if (memcmp(xr, r, n * sizeof(uint64_t))) abort();
    TIME1(1, mul_6464_6464_sse, (r, a, w));
#endif

    mul_6464_6464_v2(r, a, w);
    if (memcmp(xr, r, n * sizeof(uint64_t))) abort();
    TIME1(1, mul_6464_6464_v2, (r, a, w));
}

void level1_mul_tests_64()
{
    unsigned int n = 64;
    l1_data D;
    l1_data_init_set(D, n);

    /* multiplicate of two 64x64 matrices */
    mul_6464_6464(D->r, D->a, D->w);
    memcpy(D->xr, D->r, D->n * sizeof(uint64_t));
    
    level1_mul_tests_64_list(D);
    /* Functions which can do any n can also do n=64 */
    level1_mul_tests_N_list(D);

    l1_data_clear(D);
}
/* }}} */

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

struct pmat_s {/*{{{*/
    int * v;
    int n;
};
typedef struct pmat_s pmat[1];
typedef struct pmat_s * pmat_ptr;
typedef const struct pmat_s * pmat_srcptr;

void pmat_init(pmat_ptr x, int n)
{
    x->n=n;
    x->v=malloc(n*sizeof(int));
    for(int i = 0 ; i < n ; i++) x->v[i]=i;
}
void pmat_clear(pmat_ptr x)
{
    free(x->v); x->v = NULL; x->n = 0;
}
static inline int pmat_get(pmat_srcptr x, int k) { return x->v[k]; }
static inline void pmat_set(pmat_srcptr x, int k, int w) { x->v[k]=w; }

static inline void pmat_transpose(pmat_ptr x, pmat_srcptr y)
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
static void pmattab_complete(int * phi, uint64_t * bits, int nbits)
{
    ASSERT_ALWAYS(nbits % 64 == 0);
    for(int offset=0 ; offset < nbits ; offset+=64) {
        for(uint64_t w = bits[offset/64], z ; w ; w^=z) {
            z = w^(w&(w-1));
            uint64_t j;
            __asm__ ("bsfq %1,%0" : "=r" (j) : "rm" ((uint64_t)w));
            *phi++ = offset + j;
        }
    }
}

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

/* {{{ PLUQ stuff */
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
        uint64_t j;
        __asm__ ("bsfq %1,%0" : "=r" (j) : "rm" ((uint64_t)r));
        phi[i] = col_offset + j;
#if defined(HAVE_SSE41) && !defined(VALGRIND)
        int k = i+1;
        if (k&1) {      // alignment call
            uint64_t w = -((u[k]&v)!=0);
            u[k]^=r&w;
            l[k]^=l[i]&w;
            k++;
        }
        /* ok, it's ugly, and requires sse 4.1.
         * but otoh is churns out data veeery fast */
        __m128i vv = (__v2di) { v,v };
        __m128i pp = (__v2di) { r, r };
        __m128i ee = (__v2di) { l[i], l[i] };
        __m128i * uu = (__m128i*) (u+k);
        __m128i * ll = (__m128i*) (l+k);
        for( ; k < n ; k+=2 ) {
            __v2di ww = _mm_cmpeq_epi64(*uu&vv,vv);
            *uu++ ^= pp & ww;
            *ll++ ^= ee & ww;
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
    return rank;
}

/* outer routine */
int PLUQ64(pmat_ptr p, mat64 * l, mat64 * u, pmat_ptr q, mat64 * m)
{
    int phi[64];
    for(int i = 0 ; i < 64 ; i++) phi[i]=-1;
    int r = PLUQ64_inner(phi,l[0],u[0],m[0],0);
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
    mat64 ** ls = malloc(nb * sizeof(mat64*));
#else
    mat64 ls[nb];
#endif
    mat64 tl;
    for( ; b < nb && rank < m ; b++) {
        mat64 ta;
        mul_6464_6464(ta, l, a[b]);
        rank += PLUQ64_inner(phi, tl, u[b], ta, b*m);
        mul_6464_6464(l, tl, l);
#ifdef  ALLOC_LS
        ls[b]=malloc(sizeof(mat64));
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
        free(ls[c]);
    }
    free(ls);
#else
    for(int c = b-2 ; c >= 0 ; c--) {
        mul_6464_6464(u[c], tl, u[c]);
        mul_6464_6464(tl, ls[c], tl);
    }
#endif
    for( ; b < nb ; b++)
        mul_6464_6464(u[b], l, u[b]);
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
        __v2di mm = (__v2di) { m,m};
        __m128i * uu = (__m128i*) us[d];
        __m128i * hh = (__m128i*) h;
        __m128i hi = (__v2di) { h[i], h[i] };
        int ii=i/2;
        for( ; k < ii ; k++) {
            __v2di ww = _mm_cmpeq_epi64(*uu++&mm,mm);
            for(int b = 0 ; b < nb ; b++) {
                ((__m128i*)us[b])[k] ^= ww & (__v2di) {us[b][i],us[b][i]};
            }
            hh[k] ^= ww & hi;
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
}

void bli_64x128(mat64 h, mat64 * us, int * phi)
{
    mat64 uc[2];
    memcpy(uc,us,2*sizeof(mat64));
    bli_64x64N_clobber(h,uc,phi,2);
}

void extract_cols_64_from_128(mat64 t, mat64 * m, int * phi)
{
    // given the list of 64 integers phi, all in the rage {-1} union
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
#if defined(HAVE_SSE41) && !defined(VALGRIND)
        __v2di ss[2] = {(__v2di) {s[0][j],s[0][j]}, (__v2di) {s[1][j],s[1][j]}};
        __m128i * mm[2] = {(__m128i*)m[0],(__m128i*)m[1]};
        __m128i * tt = (__m128i*)t;
        __v2di mmk = (__v2di) {mask,mask};
        for(int i = 0 ; i < 64 ; i+=2) {
            *tt ^= mmk & _mm_cmpeq_epi64((*mm[0]&ss[0])^(*mm[1]&ss[1]),ss[0]^ss[1]);
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

int PLUQ128(pmat_ptr p, mat64 * l, mat64 * u, pmat_ptr q, mat64 * m)
{
    /* This is really an outer routine. An inner routine will not have p
     * and q, but rather both merged as a phi argument, in the manner of
     * PLUQ64_inner (and of course the following lines would be changed a
     * bit).
     */
    int phi[128];
    for(int i = 0 ; i < 128 ; i++) phi[i]=-1;

    int r1 = PLUQ64_n(phi,l[0],u,m,128);
    r1 = r1 % 64;

    // andouille 7.65

    mat64 h;
    bli_64x128(h, u, phi);
    mat64 l21;

    /* This is __very__ expensive w.r.t. what it really does :-(( */
    extract_cols_64_from_128(l21, m+2, phi);

    // andouille 16.7 -- 17.5
    mul_6464_6464(l21, l21, h);
    mul_6464_6464(l[2], l21, l[0]);

    mat64 t2[2];
    mul_6464_6464(t2[0], l21, u[0]); add_6464_6464(t2[0], m[2], t2[0]);
    mul_6464_6464(t2[1], l21, u[1]); add_6464_6464(t2[1], m[3], t2[1]);

    int r2 = PLUQ64_n(phi + 64,l[3],u+2,t2,128);
    r2 = r2 % 64;
    mul_6464_6464(l[2], l[3], l[2]);

    pqperms_from_phi(p,q,phi,128,128);
    
    return r1 + r2;
}

void check_pluq(pmat_ptr p, mat64 * l, mat64 * u, pmat_ptr q, mat64 * m, int n)
{
    mat64 pm[(n/64)*(n/64)];
    pmat_get_matrix(pm, p);

    pmat qt;
    pmat_init(qt, n);
    pmat_transpose(qt, q);

    mat64 qmt[(n/64)*(n/64)];
    pmat_get_matrix(qmt, qt);

    /* compute p*u*q^-1 */
    mat64 pu[(n/64)*(n/64)];
    memset(pu, 0, (n/64)*(n/64)*sizeof(mat64));

    for(int i = 0 ; i < (n/64) ; i++ )
    for(int j = 0 ; j < (n/64) ; j++ )
    for(int k = 0 ; k < (n/64) ; k++ ) {
        mat64 tmp;
        mul_6464_6464(tmp, pm[i*(n/64)+k], u[k*(n/64)+j]);
        add_6464_6464(pu[i*(n/64)+j], pu[i*(n/64)+j], tmp);
    }

    mat64 puq[(n/64)*(n/64)];
    memset(puq, 0, (n/64)*(n/64)*sizeof(mat64));

    for(int i = 0 ; i < (n/64) ; i++ )
    for(int j = 0 ; j < (n/64) ; j++ )
    for(int k = 0 ; k < (n/64) ; k++ ) {
        mat64 tmp;
        mul_6464_6464(tmp, pu[i*(n/64)+k], qmt[k*(n/64)+j]);
        add_6464_6464(puq[i*(n/64)+j], puq[i*(n/64)+j], tmp);
    }

    mat64 lm[(n/64)*(n/64)];
    memset(lm, 0, (n/64)*(n/64)*sizeof(mat64));

    for(int i = 0 ; i < (n/64) ; i++ )
    for(int j = 0 ; j < (n/64) ; j++ )
    for(int k = 0 ; k <= i ; k++ ) {
        mat64 tmp;
        mul_6464_6464(tmp, l[i*(n/64)+k], m[k*(n/64)+j]);
        add_6464_6464(lm[i*(n/64)+j], lm[i*(n/64)+j], tmp);
    }

    for(int i = 0 ; i < (n/64) ; i++ ) {
        ASSERT_ALWAYS(mat64_is_lowertriangular(l[i*(n/64)+i]));
        ASSERT_ALWAYS(mat64_triangular_is_unit(l[i*(n/64)+i]));
        for(int j = 0 ; j < (n/64) ; j++ ) {
            ASSERT_ALWAYS(mat64_eq(lm[i*(n/64)+j], u[i*(n/64)+j]));
        }
    }
    for(int i = 0 ; i < (n/64) ; i++ ) {
        ASSERT_ALWAYS(mat64_is_uppertriangular(puq[i*(n/64)+i]));
    }
}
/* }}} */

void level3_gauss_tests_N(int n __attribute__((unused)))
{
#ifdef  HAVE_M4RI
    mzd_t * M;
    mzd_t * LU;
    mzp_t * P, *Q;
    M = mzd_init(n, n);
#if 0
    mzd_set_mem(M, m, n);
    uint64_t m = malloc((n*n/64)*sizeof(uint64_t));
    rand64_mem(m, 64);
    free(m);
#else
    my_mzd_randomize(M);
#endif
    LU = mzd_init(n,n);
    P = mzp_init(n);
    Q = mzp_init(n);
    TIME1N(2, mzd_mypluq, (LU, M, P, Q, 0));
    TIME1N(2, mzd_myechelonize_m4ri, (LU, M, 0, 0));
    TIME1N(2, mzd_myechelonize_pluq, (LU, M, 0));
    mzd_free(M);
    mzd_free(LU);
    mzp_free(P);
    mzp_free(Q);
#endif
}

int main()
{
  unsigned int n = 2 * 1000 * 1000;

    if (0) {
        uint64_t * r = (uint64_t *) malloc(64 * sizeof(uint64_t));
        uint64_t a = rand64();
        uint64_t w = rand64();

        printf("-- level-1 benches --\n");
        level1_basic_tests();
        level1_mul_tests_1();
        printf("-- level-1 benches, n=64 --\n");
        level1_mul_tests_64();
        printf("-- level-1 benches, misc --\n");
        TIME1(1, addmul_To64_o64_lsb, (r, a, w));
        TIME1(1, addmul_To64_o64_msb, (r, a, w));
        TIME1(1, addmul_To64_o64_lsb_packof2, (r, a, w));
#ifdef  HAVE_SSE2
        TIME1(1, addmul_To64_o64_lsb_sse_v1, (r, a, w));
#endif
        TIME1(1, addmul_To64_o64, (r, a, w));
        free(r);
    }

    if (1) {
        uint64_t * r = (uint64_t *) malloc(64 * sizeof(uint64_t));
        uint64_t * a = (uint64_t *) malloc(n * sizeof(uint64_t));
        uint64_t * w = (uint64_t *) malloc(n * sizeof(uint64_t));
        rand64_mem(a, n);
        rand64_mem(w, n);

        printf("-- level-2 benches (N=%u) --\n", n);
        level1_mul_tests_N(n);
        TIME1N(1, mul_64N_N64_addmul, (r,a,w,n));
        TIME1N(5, mul_TN32_N64_C, (r,(uint32_t*)a,w,n));
        TIME1N(5, mul_TN64_N64_C, (r,a,w,n));

        free(r); free(a); free(w);
    }

    if (1) {
        mat64 m[4],l[4],u[4];
        srand(1728);

        {
            pmat p,q;
            pmat_init(p, 64);
            pmat_init(q, 64);
            rand64_mem((uint64_t*)m, 64);
            PLUQ64(p,l,u,q,m);
            check_pluq(p,l,u,q,m,64);
            pmat_clear(p);
            pmat_clear(q);
        }

        {
            pmat p,q;
            pmat_init(p, 128);
            pmat_init(q, 128);
            rand64_mem((uint64_t*)m, 4*64);
            PLUQ128(p,l,u,q,m);
            check_pluq(p,l,u,q,m,128);
            pmat_clear(p);
            pmat_clear(q);
        }

        printf("PLUQ128 stub executed ok\n");


            // pmat_6464(m[0]); printf("\n");
            // pmat_6464(u); printf("\n");
            // pmat_6464(l); printf("\n");
            // pmat_6464(p); printf("\n");
            // pmat_6464(t); printf("\n");
            

            // mul_N64_T6464(t,u,p,64);
            // pmat_6464(t); printf("\n");

            
#if 0
            memset(e,0,sizeof(e));
            rand64_mem(m[0], 64);
            rand64_mem(m[1], 64);
            rand64_mem(m[2], 64);
            rand64_mem(m[3], 64);
            // pmat_6464(m[0]); printf("\n");
            full_echelon_6464_imm(mm,e[0],m[0]);

            /* 
            mul_6464

            pmat_6464(mm); printf("\n");
            pmat_6464(e); printf("\n"); printf("\n");
            */
#endif
    }

    if (1) {
        mat64 m;
        mat64 e;
        mat64 mm;
        mat64 l,u,p;
        rand64_mem(m, 64);
        rand64_mem(e, 64);
        rand64_mem(mm, 64);
        mat64 m4[4];
        mat64 u4[4];
        rand64_mem((uint64_t*)m4, 256);
        // printf("-- for reference: best matrix mult, 64x64 --\n");
        // TIME1(2, mul_6464_6464, (mm,e,m));
        // TIME1(2, mul_N64_T6464, (mm,e,m,64));
        printf("-- level-3 (reduction) benches, n=64 --\n");
        // TIME1(2, gauss_6464_C, (mm,e,m));
        // TIME1(2, gauss_6464_imm, (mm,e,m));
        // TIME1(2, PLUQ64_inner, (NULL,l,u,m,0));
        int phi[128];
        {
            pmat p,q;
            mat64 m[4],l[2],u[4];
            rand64_mem((uint64_t*)m, 4*64);
            pmat_init(p, 128);
            pmat_init(q, 128);
            TIME1(2, PLUQ128, (p,l,u,q,m));
        }
        int n=2;
        TIME1N(2, rand64_mem((uint64_t*)m4, n*64), );
        TIME1N_spins(, 2, PLUQ64_n, (phi,l,u4,m4,64*n), ret/64,n+1);
        TIME1N_spins(rand64_mem((uint64_t*)m4, n*64), 2, PLUQ64_n, (phi,l,u4,m4,64*n), ret/64,n+1);
        TIME1(2, LUP64_imm, (l,u,p,m));
        TIME1(2, full_echelon_6464_imm, (mm,e,m));
        TIME1(2, gauss_128128_C, (m));
        level3_gauss_tests_N(64);
        level3_gauss_tests_N(128);
        level3_gauss_tests_N(256);
        level3_gauss_tests_N(512);
        level3_gauss_tests_N(1024);
    }

    if (0) {
        size_t n = 512;
        mat64 * A = (mat64 *) malloc(n * sizeof(mat64));
        mat64 * B = (mat64 *) malloc(n * sizeof(mat64));
        mat64 * C = (mat64 *) malloc(2 *n * sizeof(mat64));
        printf("-- polynomials (N=%zu) --\n", n);
        TIME1(5, m64pol_mul, (C,A,B,n,n));
        TIME1(5, m64pol_mul_kara, (C,A,B,n,n));
        free(A);
        free(B);
        free(C);
    }

    if (0) {
        size_t n = 512;
        mat64 * A = (mat64 *) malloc(2 * 2 * n * sizeof(mat64));
        mat64 * B = (mat64 *) malloc(2 * 2 * n * sizeof(mat64));
        mat64 * C = (mat64 *) malloc(2 * 2 * 2 *n * sizeof(mat64));
        printf("-- polynomials, larger matrices (N=%zu) --\n", n);
        TIME1(5, m64polblock_mul, (C,A,B,n,n,2));
        TIME1(5, m64polblock_mul_kara, (C,A,B,n,n,2));
        free(A);
        free(B);
        free(C);
    }

    return 0;
}
