/* This file is part of the gf2x library.

   Copyright 2010, 2013, 2015, 2016
   Richard Brent, Pierrick Gaudry, Emmanuel Thome', Paul Zimmermann

   This program is free software; you can redistribute it and/or modify it
   under the terms of either:
    - If the archive contains a file named toom-gpl.c (not a trivial
    placeholder), the GNU General Public License as published by the Free
    Software Foundation; either version 3 of the License, or (at your
    option) any later version.
    - If the archive contains a file named toom-gpl.c which is a trivial
    placeholder, the GNU Lesser General Public License as published by
    the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.
   
   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.  See the license text for more details.
   
   You should have received a copy of the GNU General Public License as
   well as the GNU Lesser General Public License along with this program;
   see the files COPYING and COPYING.LIB.  If not, write to the Free
   Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
   02110-1301, USA.
*/

#ifndef GF2X_MUL5_H_
#define GF2X_MUL5_H_

#include "gf2x.h"
/* All gf2x source files for lowlevel functions must include gf2x-small.h
 * This is mandatory for the tuning mechanism. */
#include "gf2x/gf2x-small.h"

#if GF2X_WORDSIZE != 64
#error "This code is for 64-bit only"
#endif

#ifndef GF2X_HAVE_PCLMUL_SUPPORT
#error "This code needs pclmul support"
#endif

/* TODO: if somebody comes up with a neat way to improve the interface so
 * as to remove the false dependency on pclmul, that would be nice.
 */
/* This specialized version avoids loads, and relies on the destination
 * being aligned, so that aligned stores are possible */
#define PXOR(lop, rop) _mm_xor_si128((lop), (rop))
#define PXOR3(op1, op2, op3) PXOR(op1, PXOR(op2, op3))
#define PXOR4(op1, op2, op3, op4) PXOR(op1, PXOR3(op2, op3, op4))
#define PZERO    _mm_setzero_si128()

static inline __m128i
GF2X_FUNC(mul5clk2_mul1)(__m128i ss1, __m128i ss2)
{
    return _mm_clmulepi64_si128(ss1, ss2, 0);
}

static inline void
GF2X_FUNC(mul5clk2_mul2)(__m128i * t, __m128i ss1, __m128i ss2)
{
    __m128i t00 = _mm_clmulepi64_si128(ss1, ss2, 0);
    __m128i t11 = _mm_clmulepi64_si128(ss1, ss2, 0x11);
    ss1 = PXOR(ss1, _mm_shuffle_epi32(ss1, _MM_SHUFFLE(1,0,3,2)));
    ss2 = PXOR(ss2, _mm_shuffle_epi32(ss2, _MM_SHUFFLE(1,0,3,2)));
    __m128i tk = PXOR(PXOR(t00, t11), _mm_clmulepi64_si128(ss1, ss2, 0));
    t[0] = PXOR(t00, _mm_unpacklo_epi64(PZERO, tk));
    t[1] = PXOR(t11, _mm_unpackhi_epi64(tk, PZERO));
}


/* variant with 5 calls to mul2 and 1 to mul1, i.e., 16 multiplications */
GF2X_STORAGE_CLASS_mul5
void gf2x_mul5 (unsigned long *c, const unsigned long *a, const unsigned long *b)
{
    __m128i aa[3], bb[3];
    __m128i p0[2], p1[2], p2[1];
    __m128i pp0[2], pp1[2], pp2[2];
    __m128i a0 = _mm_loadu_si128((__m128i*)(a));
    __m128i a1 = _mm_loadu_si128((__m128i*)(a+2));
    /* load lower half only, using movq */
    __m128i a2 = _mm_loadl_epi64((__m128i*)(a+4));
    __m128i b0 = _mm_loadu_si128((__m128i*)(b));
    __m128i b1 = _mm_loadu_si128((__m128i*)(b+2));
    /* load lower half only, using movq */
    __m128i b2 = _mm_loadl_epi64((__m128i*)(b+4));
    aa[0] = PXOR(a1, a2);
    aa[1] = PXOR(a0, a2);
    aa[2] = PXOR(a0, a1);
    bb[0] = PXOR(b1, b2);
    bb[1] = PXOR(b0, b2);
    bb[2] = PXOR(b0, b1);
    GF2X_FUNC(mul5clk2_mul2)(p0, a0, b0);
    GF2X_FUNC(mul5clk2_mul2)(p1, a1, b1);
    p2[0] = GF2X_FUNC(mul5clk2_mul1)(a2, b2);
    GF2X_FUNC(mul5clk2_mul2)(pp0, aa[0], bb[0]);
    GF2X_FUNC(mul5clk2_mul2)(pp1, aa[1], bb[1]);
    GF2X_FUNC(mul5clk2_mul2)(pp2, aa[2], bb[2]);
    _mm_storeu_si128((__m128i*)(c + 0), p0[0]);
    _mm_storeu_si128((__m128i*)(c + 2),
            PXOR(PXOR3(p0[0], p1[0], pp2[0])       , p0[1]));
    _mm_storeu_si128((__m128i*)(c + 4),
            PXOR(PXOR4(p0[0], p1[0], p2[0], pp1[0]), PXOR3(p0[1], p1[1], pp2[1])));
    _mm_storeu_si128((__m128i*)(c + 6),
            PXOR(PXOR3(pp0[0], p1[0], p2[0])       , PXOR3(p0[1], p1[1], pp1[1])));
    _mm_storeu_si128((__m128i*)(c + 8),
            PXOR(p2[0]                             , PXOR(pp0[1], p1[1])));
}

#undef PXOR
#undef PXOR3
#undef PXOR4
#undef PZERO
#endif  /* GF2X_MUL5_H_ */
