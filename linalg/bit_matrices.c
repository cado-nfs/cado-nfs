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

/* This builds on top of 64x64 matrices, for which we hope to have some fast
 * arithmetic. Note that we do not expect the code here to be fast on 32-bit
 * hardware.
 */

#include <string.h>
#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif
#include "bit_matrices.h"

/* copy-pasted from mul_TN64_N64_C in linalg/bwc/matops.c
 */
void addmul_TN64_N64(mat64 b, uint64_t * A, uint64_t * x, unsigned int ncol)
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

void transp_6464(mat64 dst, mat64 src)
{
  int i, j;

  for (i = 0; i < 64; i++) {
    dst[i] = 0;
    for (j = 0; j < 64; j++) {
      dst[i] ^= ((src[j] >> i) & 1UL) << j;
    }
  }
}

/* With a table of 256 precomputed multiples, we get a nice speed
 * improvement for m ~ 10^6 or above. However this not used in critical
 * areas presently.
 */
void addmul_N64_6464(uint64_t *C,
                   const uint64_t *A,
                   mat64 B, unsigned long m)
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

#ifdef  HAVE_SSE2
void mul_6464_6464(mat64_ptr C, mat64_srcptr A, mat64_srcptr B)
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
#else
/* slow fallback */
void mul_6464_6464(mat64_ptr C, mat64_srcptr A, mat64_srcptr B)
{
    memset(C, 0, sizeof(mat64));
 
    for (int i = 0; i < 64; i++) {
    for (int j = 0; j < 64; j++) {
        if ((A[i]>>j)&1)
            C[i]^=B[j];
    }
    }
}
#endif
