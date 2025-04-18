/* This file is part of the gf2x library.

   Copyright 2010, 2012, 2013, 2015, 2016
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

#ifndef GF2X_MUL9_H_
#define GF2X_MUL9_H_

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
static inline __m128i
GF2X_FUNC(mul9clk2_mul1) (unsigned long a, unsigned long b)
{
    __m128i aa = _gf2x_mm_setr_epi64(a, 0);
    __m128i bb = _gf2x_mm_setr_epi64(b, 0);
    return _mm_clmulepi64_si128(aa, bb, 0);
}

#define PXOR(lop, rop) _mm_xor_si128((lop), (rop))
#define PXOR3(op1, op2, op3) PXOR(op1, PXOR(op2, op3))
#define PXOR4(op1, op2, op3, op4) PXOR(op1, PXOR3(op2, op3, op4))
#define PXOR5(op1, op2, op3, op4, op5) PXOR(op1, PXOR4(op2, op3, op4, op5))
#define PXOR6(op1, op2, op3, op4, op5, op6) PXOR(op1, PXOR5(op2, op3, op4, op5, op6))
#define PXOR7(op1, op2, op3, op4, op5, op6, op7) PXOR(op1, PXOR6(op2, op3, op4, op5, op6, op7))
#define PZERO    _mm_setzero_si128()

/* uses the variant of Karatsuba with 6 multiplications */
static void GF2X_FUNC(mul9clk2_mul3) (__m128i *ce, __m128i *co, const unsigned long *a, const unsigned long *b)
{
  unsigned long aa[3], bb[3];
  __m128i p0, p1, p2;
  __m128i pp0, pp1, pp2;
  aa[0] = a[1]^a[2];
  aa[1] = a[0]^a[2];
  aa[2] = a[0]^a[1];
  bb[0] = b[1]^b[2];
  bb[1] = b[0]^b[2];
  bb[2] = b[0]^b[1];
  p0  = GF2X_FUNC(mul9clk2_mul1)(a[0], b[0]);
  p1  = GF2X_FUNC(mul9clk2_mul1)(a[1], b[1]);
  p2  = GF2X_FUNC(mul9clk2_mul1)(a[2], b[2]);
  pp0 = GF2X_FUNC(mul9clk2_mul1)(aa[0], bb[0]);
  pp1 = GF2X_FUNC(mul9clk2_mul1)(aa[1], bb[1]);
  pp2 = GF2X_FUNC(mul9clk2_mul1)(aa[2], bb[2]);

  ce[0] = p0;
  ce[1] = PXOR4(p0, p1, p2, pp1);
  ce[2] = p2;

  co[0] = PXOR3(p0, p1, pp2);
  co[1] = PXOR3(pp0, p1, p2);
}

/* recursive application of the Karatsuba-3 algorithm with 6 multiplies,
   i.e., with a total of 6*6=36 multiplies */
GF2X_STORAGE_CLASS_mul9
void gf2x_mul9 (unsigned long *c, const unsigned long *a, const unsigned long *b)
{
  unsigned long aa[9], bb[9];

  /*
   * a = al(t) + t^3 * am(t) + t^6 * ah(t)
   * b = bl(t) + t^3 * bm(t) + t^6 * bh(t)
   *
   * aa0 = am + ah
   * aa1 = al + ah
   * aa2 = al + am
   * 
   * p0e + z^64 * p0o = al * bl
   * p1e + z^64 * p1o = am * bm
   * p2e + z^64 * p2o = ah * bh
   *
   * aa0 * bb0   = (am*bm + ah*bh + am*bh + ah*bm)
   * aa1 * bb1   = (al*bl + ah*bh + al*bh + ah*bl)
   * aa2 * bb2   = (al*bl + am*bm + al*bm + am*bl)
   *
   * a*b =  al*bl * t^0
   *      + (aa2*bb2 + al*bl + am*bm) * t^3
   *      + (am*bm + al*bl + ah*bh + aa1 * bb1) * t^6
   *      + (aa0*bb0 + ah*bh + am*bm) * t^9
   *      + ah*bh * t^12
   */
  aa[0] = a[3]^a[6];
  aa[1] = a[4]^a[7];
  aa[2] = a[5]^a[8];

  aa[3] = a[0]^a[6];
  aa[4] = a[1]^a[7];
  aa[5] = a[2]^a[8];

  aa[6] = a[0]^a[3];
  aa[7] = a[1]^a[4];
  aa[8] = a[2]^a[5];

  bb[0] = b[3]^b[6];
  bb[1] = b[4]^b[7];
  bb[2] = b[5]^b[8];

  bb[3] = b[0]^b[6];
  bb[4] = b[1]^b[7];
  bb[5] = b[2]^b[8];

  bb[6] = b[0]^b[3];
  bb[7] = b[1]^b[4];
  bb[8] = b[2]^b[5];


  __m128i p0e[3], p0o[2];
  __m128i p1e[3], p1o[2];
  __m128i p2e[3], p2o[2];
  __m128i q0e[3], q0o[2];
  __m128i q1e[3], q1o[2];
  __m128i q2e[3], q2o[2];
  GF2X_FUNC(mul9clk2_mul3) (p0e, p0o, a+0, b+0);
  GF2X_FUNC(mul9clk2_mul3) (p1e, p1o, a+3, b+3);
  GF2X_FUNC(mul9clk2_mul3) (p2e, p2o, a+6, b+6);
  GF2X_FUNC(mul9clk2_mul3) (q0e, q0o, aa+0, bb+0);
  GF2X_FUNC(mul9clk2_mul3) (q1e, q1o, aa+3, bb+3);
  GF2X_FUNC(mul9clk2_mul3) (q2e, q2o, aa+6, bb+6);

  __m128i e,h,l;
  e = p0e[0];
  l = p0o[0];
  _mm_storeu_si128((__m128i*)(c), PXOR(e, _mm_slli_si128(l, 8)));

  e = p0e[1];                      h=l; l= PXOR4(p0o[1], p1e[0], q2e[0], p0e[0]);
  _mm_storeu_si128((__m128i*)(c+2), PXOR3(e, _mm_slli_si128(l, 8), _mm_srli_si128(h, 8)));

  e = PXOR4(p0e[2], p0o[0], p1o[0], q2o[0]); h=l; l= PXOR3(p1e[1], q2e[1], p0e[1]);
  _mm_storeu_si128((__m128i*)(c+4), PXOR3(e, _mm_slli_si128(l, 8), _mm_srli_si128(h, 8)));

  e = PXOR7(p0e[0], q1e[0], p2e[0], p1e[0], p0o[1], p1o[1], q2o[1]); h=l; l= PXOR7(p0o[0], q1o[0], p2o[0], p1o[0], p1e[2], q2e[2], p0e[2]);
  _mm_storeu_si128((__m128i*)(c+6), PXOR3(e, _mm_slli_si128(l, 8), _mm_srli_si128(h, 8)));

  e = PXOR4(p0e[1], q1e[1], p2e[1], p1e[1]); h=l; l= PXOR7(p0o[1], q1o[1], p2o[1], p1o[1], q0e[0], p1e[0], p2e[0]);
  _mm_storeu_si128((__m128i*)(c+8), PXOR3(e, _mm_slli_si128(l, 8), _mm_srli_si128(h, 8)));

  e = PXOR7(p0e[2], q1e[2], p2e[2], p1e[2], q0o[0], p1o[0], p2o[0]); h=l; l= PXOR3(q0e[1], p1e[1], p2e[1]);
  _mm_storeu_si128((__m128i*)(c+10), PXOR3(e, _mm_slli_si128(l, 8), _mm_srli_si128(h, 8)));

  e = PXOR4(p2e[0], q0o[1], p1o[1], p2o[1]); h=l; l= PXOR4(p2o[0], q0e[2], p1e[2], p2e[2]);
  _mm_storeu_si128((__m128i*)(c+12), PXOR3(e, _mm_slli_si128(l, 8), _mm_srli_si128(h, 8)));

  e = p2e[1]; h=l; l=p2o[1];
  _mm_storeu_si128((__m128i*)(c+14), PXOR3(e, _mm_slli_si128(l, 8), _mm_srli_si128(h, 8)));

  e = p2e[2]; h=l;
  _mm_storeu_si128((__m128i*)(c+16), PXOR(e, _mm_srli_si128(h, 8)));

  /*
   * a*b =  al*bl * t^0
   *      + (aa2*bb2 + al*bl + am*bm) * t^3
   *      + (am*bm + al*bl + ah*bh + aa1 * bb1) * t^6
   *      + (aa0*bb0 + ah*bh + am*bm) * t^9
   *      + ah*bh * t^12
   *
   * this is the 64-bit breakdown
  c[0]  =                   p0e[0]                   ^                                           ;
  c[1]  =                   p0e[1]                   ^                   p0o[0]                  ;
  c[2]  =                   p0e[2]                   ^                   p0o[1]                  ;
  c[3]  =                   p0e[3]                   ^                   p0o[2]p1e[0]q2e[0]p0e[0];
  c[4]  =                   p0e[4]p0o[0]p1o[0]q2o[0] ^                   p0o[3]p1e[1]q2e[1]p0e[1];
  c[5]  =                   p0e[5]p0o[1]p1o[1]q2o[1] ^                         p1e[2]q2e[2]p0e[2];
  c[6]  = p0e[0]q1e[0]p2e[0]p1e[0]p0o[2]p1o[2]q2o[2] ^                         p1e[3]q2e[3]p0e[3];
  c[7]  = p0e[1]q1e[1]p2e[1]p1e[1]p0o[3]p1o[3]q2o[3] ^ p0o[0]q1o[0]p2o[0]p1o[0]p1e[4]q2e[4]p0e[4];
  c[8]  = p0e[2]q1e[2]p2e[2]p1e[2]                   ^ p0o[1]q1o[1]p2o[1]p1o[1]p1e[5]q2e[5]p0e[5];
  c[9]  = p0e[3]q1e[3]p2e[3]p1e[3]                   ^ p0o[2]q1o[2]p2o[2]p1o[2]q0e[0]p1e[0]p2e[0];
  c[10] = p0e[4]q1e[4]p2e[4]p1e[4]q0o[0]p1o[0]p2o[0] ^ p0o[3]q1o[3]p2o[3]p1o[3]q0e[1]p1e[1]p2e[1];
  c[11] = p0e[5]q1e[5]p2e[5]p1e[5]q0o[1]p1o[1]p2o[1] ^                         q0e[2]p1e[2]p2e[2];
  c[12] = p2e[0]                  q0o[2]p1o[2]p2o[2] ^                         q0e[3]p1e[3]p2e[3];
  c[13] = p2e[1]                  q0o[3]p1o[3]p2o[3] ^                   p2o[0]q0e[4]p1e[4]p2e[4];
  c[14] = p2e[2]                                     ^                   p2o[1]q0e[5]p1e[5]p2e[5];
  c[15] = p2e[3]                                     ^                   p2o[2]                  ;
  c[16] = p2e[4]                                     ^                   p2o[3]                  ;
  c[17] = p2e[5]                                     ^                                           ;
  */
}

#undef PXOR
#undef PXOR3
#undef PXOR4
#undef PXOR5
#undef PXOR6
#undef PXOR7
#undef PZERO
#endif  /* GF2X_MUL9_H_ */
