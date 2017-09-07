/* This file is part of the gf2x library.

   Copyright 2007, 2008, 2009, 2010, 2013, 2015
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

/* Implements 128x128 -> 256 bit product using SSE2 instructions. */

#ifndef GF2X_MUL4_H_
#define GF2X_MUL4_H_

#include "gf2x.h"
/* All gf2x source files for lowlevel functions must include gf2x-small.h
 * This is mandatory for the tuning mechanism. */
#include "gf2x/gf2x-small.h"

#if GF2X_WORDSIZE != 32
#error "This code is for 32-bit only"
#endif

#ifndef GF2X_HAVE_SSE2_SUPPORT
#error "This code needs sse-2 support"
#endif

#if defined(__GNUC__) && __GNUC__ == 4 && __GNUC_MINOR__ == 3 && (__GNUC_PATCHLEVEL__ == 0 || __GNUC_PATCHLEVEL__ == 1)
#warning "Your GCC version is buggy. Binary fields may fail randomly"
/* Gcc bug reports 37101 and 37340 -- the only convenient fix is to
 * upgrade to 4.3.2 */
#endif

/* This has been edited without testing */

GF2X_STORAGE_CLASS_mul4
void gf2x_mul4(unsigned long *t, unsigned long const *s1,
               unsigned long const *s2)
{
#define SHL(x, r) _mm_slli_epi64((x), (r))
#define SHR(x, r) _mm_srli_epi64((x), (r))
#define SHLD(x, r) _mm_slli_si128((x), (r) >> 3)
#define SHRD(x, r) _mm_srli_si128((x), (r) >> 3)
#define PZERO    _mm_setzero_si128()
#define PADD(x, y)      _mm_add_epi64((x), (y))
#define PSUB(x, y)      _mm_sub_epi64((x), (y))
#define PNEG(x)      PSUB(PZERO, (x))
#define PXOR(lop, rop) _mm_xor_si128((lop), (rop))
#define XOREQ(lop, rop) lop = _mm_xor_si128((lop), (rop))
#define PAND(lop, rop) _mm_and_si128((lop), (rop))

    __m128i u;
    __m128i t0;
    __m128i t1;
    __m128i t2;

    __m128i g[16];
    /* sequence update walk */
    g[0] = PZERO;
    g[1] = _gf2x_mm_setr_epi32(s2[0], s2[1], s2[2], s2[3]);
    g[2] = SHL(g[1], 1);
    g[3] = PXOR(g[2], g[1]);
    g[4] = SHL(g[2], 1);
    g[5] = PXOR(g[4], g[1]);
    g[6] = SHL(g[3], 1);
    g[7] = PXOR(g[6], g[1]);
    g[8] = SHL(g[4], 1);
    g[9] = PXOR(g[8], g[1]);
    g[10] = SHL(g[5], 1);
    g[11] = PXOR(g[10], g[1]);
    g[12] = SHL(g[6], 1);
    g[13] = PXOR(g[12], g[1]);
    g[14] = SHL(g[7], 1);
    g[15] = PXOR(g[14], g[1]);

    /* round 0 */
    u = g[s1[0]       & 15];
    t0  = u;
    u = g[s1[0] >>  4 & 15];
    XOREQ(t0, SHL(u,  4)); t1  = SHR(u, 60);
    u = g[s1[0] >>  8 & 15];
    XOREQ(t0, SHL(u,  8)); XOREQ(t1, SHR(u, 56));
    u = g[s1[0] >> 12 & 15];
    XOREQ(t0, SHL(u, 12)); XOREQ(t1, SHR(u, 52));
    u = g[s1[0] >> 16 & 15];
    XOREQ(t0, SHL(u, 16)); XOREQ(t1, SHR(u, 48));
    u = g[s1[0] >> 20 & 15];
    XOREQ(t0, SHL(u, 20)); XOREQ(t1, SHR(u, 44));
    u = g[s1[0] >> 24 & 15];
    XOREQ(t0, SHL(u, 24)); XOREQ(t1, SHR(u, 40));
    u = g[s1[0] >> 28 & 15];
    XOREQ(t0, SHL(u, 28)); XOREQ(t1, SHR(u, 36));
    u = g[s1[1]       & 15];
    XOREQ(t0, SHL(u, 32)); XOREQ(t1, SHR(u, 32));
    u = g[s1[1] >>  4 & 15];
    XOREQ(t0, SHL(u, 36)); XOREQ(t1, SHR(u, 28));
    u = g[s1[1] >>  8 & 15];
    XOREQ(t0, SHL(u, 40)); XOREQ(t1, SHR(u, 24));
    u = g[s1[1] >> 12 & 15];
    XOREQ(t0, SHL(u, 44)); XOREQ(t1, SHR(u, 20));
    u = g[s1[1] >> 16 & 15];
    XOREQ(t0, SHL(u, 48)); XOREQ(t1, SHR(u, 16));
    u = g[s1[1] >> 20 & 15];
    XOREQ(t0, SHL(u, 52)); XOREQ(t1, SHR(u, 12));
    u = g[s1[1] >> 24 & 15];
    XOREQ(t0, SHL(u, 56)); XOREQ(t1, SHR(u,  8));
    u = g[s1[1] >> 28 & 15];
    XOREQ(t0, SHL(u, 60)); XOREQ(t1, SHR(u,  4));

    /* round 1 */
    u = g[s1[2]       & 15];
    XOREQ(t1, u);
    u = g[s1[2] >>  4 & 15];
    XOREQ(t1, SHL(u,  4)); t2  = SHR(u, 60);
    u = g[s1[2] >>  8 & 15];
    XOREQ(t1, SHL(u,  8)); XOREQ(t2, SHR(u, 56));
    u = g[s1[2] >> 12 & 15];
    XOREQ(t1, SHL(u, 12)); XOREQ(t2, SHR(u, 52));
    u = g[s1[2] >> 16 & 15];
    XOREQ(t1, SHL(u, 16)); XOREQ(t2, SHR(u, 48));
    u = g[s1[2] >> 20 & 15];
    XOREQ(t1, SHL(u, 20)); XOREQ(t2, SHR(u, 44));
    u = g[s1[2] >> 24 & 15];
    XOREQ(t1, SHL(u, 24)); XOREQ(t2, SHR(u, 40));
    u = g[s1[2] >> 28 & 15];
    XOREQ(t1, SHL(u, 28)); XOREQ(t2, SHR(u, 36));
    u = g[s1[3]       & 15];
    XOREQ(t1, SHL(u, 32)); XOREQ(t2, SHR(u, 32));
    u = g[s1[3] >>  4 & 15];
    XOREQ(t1, SHL(u, 36)); XOREQ(t2, SHR(u, 28));
    u = g[s1[3] >>  8 & 15];
    XOREQ(t1, SHL(u, 40)); XOREQ(t2, SHR(u, 24));
    u = g[s1[3] >> 12 & 15];
    XOREQ(t1, SHL(u, 44)); XOREQ(t2, SHR(u, 20));
    u = g[s1[3] >> 16 & 15];
    XOREQ(t1, SHL(u, 48)); XOREQ(t2, SHR(u, 16));
    u = g[s1[3] >> 20 & 15];
    XOREQ(t1, SHL(u, 52)); XOREQ(t2, SHR(u, 12));
    u = g[s1[3] >> 24 & 15];
    XOREQ(t1, SHL(u, 56)); XOREQ(t2, SHR(u,  8));
    u = g[s1[3] >> 28 & 15];
    XOREQ(t1, SHL(u, 60)); XOREQ(t2, SHR(u,  4));
    /* end */

    /* repair steps */
    /* repair section 200711-200803 */
    __m128i v1 = SHR(_gf2x_mm_setr_epi32(s1[0], s1[1], s1[0], s1[1]), 1);
    __m128i v2 = SHR(_gf2x_mm_setr_epi32(s1[2], s1[3], s1[2], s1[3]), 1);
    __m128i m = _gf2x_mm_set1_epi32_c(0x77777777);
    __m128i w = PNEG(SHR(g[1],63));
    v1 = PAND(v1, m);
    XOREQ(t1, PAND(v1, w));
    v2 = PAND(v2, m);
    XOREQ(t2, PAND(v2, w));
    w = PNEG(SHR(g[2],63));
    v1 = SHR(v1, 1) & m;
    XOREQ(t1, PAND(v1, w));
    v2 = SHR(v2, 1) & m;
    XOREQ(t2, PAND(v2, w));
    w = PNEG(SHR(g[4],63));
    v1 = SHR(v1, 1) & m;
    XOREQ(t1, PAND(v1, w));
    v2 = SHR(v2, 1) & m;
    XOREQ(t2, PAND(v2, w));

    /* store result */
    _mm_storeu_si128((__m128i *)(t),  PXOR(t0, _mm_unpacklo_epi64(PZERO, t1)));
    _mm_storeu_si128((__m128i *)(t+4),PXOR(t2, _mm_unpackhi_epi64(t1, PZERO)));
#undef PAND
#undef XOREQ
#undef PXOR
#undef PZERO
#undef PNEG
#undef PADD
#undef PSUB
#undef SHL
#undef SHR
#undef SHLD
#undef SHRD
}
#endif  /* GF2X_MUL4_H_ */
