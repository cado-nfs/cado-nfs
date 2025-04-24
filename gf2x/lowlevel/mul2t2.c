/* This file is part of the gf2x library.

   Copyright 2007, 2008, 2009, 2010, 2012, 2013, 2015, 2017
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

#ifndef GF2X_MUL2_H_
#define GF2X_MUL2_H_

#include "gf2x.h"
/* All gf2x source files for lowlevel functions must include gf2x-small.h
 * This is mandatory for the tuning mechanism. */
#include "gf2x/gf2x-small.h"

#if GF2X_WORDSIZE != 64
#error "This code is for 64-bit only"
#endif

#ifndef GF2X_HAVE_SSE2_SUPPORT
#error "This code needs sse-2 support"
#endif

#if defined(__GNUC__) && __GNUC__ == 4 && __GNUC_MINOR__ == 3 && (__GNUC_PATCHLEVEL__ == 0 || __GNUC_PATCHLEVEL__ == 1)
#warning "Your GCC version is buggy. Binary fields may fail randomly"
/* Gcc bug reports 37101 and 37340 -- the only convenient fix is to
 * upgrade to 4.3.2 */
#endif

/* This code has been modified in comparison to the previously existing
 * one, in order to use Intel intrinsics exclusively, and not rely on gcc
 * syntax.
 *
 * The Intel C++ Intrinsics Reference is a good read. However I can't
 * find the official link right now. An old version can be found in the wild:
 *
 * http://www.info.univ-angers.fr/~richer/ens/l3info/ao/intel_intrinsics.pdf
 *
 * I don't even know if this document evolved. Probably not.
 *
 * It says: ``Document Number: 312482-003US''
 */
GF2X_STORAGE_CLASS_mul2
void gf2x_mul2(unsigned long * t, unsigned long const * s1,
        unsigned long const * s2)
{
#define SHL(x, r) _mm_slli_epi64((x), (r))
#define SHR(x, r) _mm_srli_epi64((x), (r))
#define SHLD(x, r) _mm_slli_si128((x), (r) >> 3)
#define SHRD(x, r) _mm_srli_si128((x), (r) >> 3)
#define XOREQ(lop, rop) lop = _mm_xor_si128((lop), (rop))
#define PXOR(lop, rop) _mm_xor_si128((lop), (rop))
#define PAND(lop, rop) _mm_and_si128((lop), (rop))
#define PZERO    _mm_setzero_si128()
#define PSUB(x, y)      _mm_sub_epi64((x), (y))
#define PNEG(x)      PSUB(PZERO, (x))

    __m128i u;
    __m128i t0;
    __m128i t1;
    __m128i t2;

    __m128i g[16];
    __m128i w;
    // __m128i m = _mm_set1_epi32(0xeeeeeeee);
    __m128i m = _gf2x_mm_set1_epi64_c(0xeeeeeeeeeeeeeeee);
    /* sequence update walk */
    __m128i b0 = _mm_loadu_si128((__m128i*) s2);
    g[ 0] = _mm_setzero_si128();
    g[ 1] = b0;
    __m128i v1 = _mm_loadu_si128((__m128i*) s1);
    w = PNEG(SHR(b0,63));
    __m128i v2 = _mm_unpackhi_epi64(v1, v1);
    v1 = _mm_unpacklo_epi64(v1, v1);
    v1 = SHR(PAND(v1, m), 1); t1 = PAND(v1, w);
    g[ 2] = SHL(b0, 1); g[ 3] = PXOR(g[ 2], b0);
    v2 = SHR(PAND(v2, m), 1); t2 = PAND(v2, w);
    g[ 4] = SHL(g[ 2], 1); g[ 5] = PXOR(g[ 4], b0);
    w = PNEG(SHR(g[ 2],63));
    g[ 6] = SHL(g[ 3], 1); g[ 7] = PXOR(g[ 6], b0);
    v1 = SHR(PAND(v1, m), 1); XOREQ(t1, PAND(v1, w));
    g[ 8] = SHL(g[ 4], 1); g[ 9] = PXOR(g[ 8], b0);
    v2 = SHR(PAND(v2, m), 1); XOREQ(t2, PAND(v2, w));
    g[10] = SHL(g[ 5], 1); g[11] = PXOR(g[10], b0);
    w = PNEG(SHR(g[4],63));
    g[12] = SHL(g[ 6], 1); g[13] = PXOR(g[12], b0);
    v1 = SHR(PAND(v1, m), 1); XOREQ(t1, PAND(v1, w));
    g[14] = SHL(g[ 7], 1); g[15] = PXOR(g[14], b0);
    v2 = SHR(PAND(v2, m), 1); XOREQ(t2, PAND(v2, w));


    /* round 0 */
    u = g[s1[0]       & 15]; t0  = u;
    u = g[s1[0] >>  4 & 15]; XOREQ(t0, SHL(u,  4)); XOREQ(t1, SHR(u, 60));
    u = g[s1[0] >>  8 & 15]; XOREQ(t0, SHL(u,  8)); XOREQ(t1, SHR(u, 56));
    u = g[s1[0] >> 12 & 15]; XOREQ(t0, SHL(u, 12)); XOREQ(t1, SHR(u, 52));
    u = g[s1[0] >> 16 & 15]; XOREQ(t0, SHL(u, 16)); XOREQ(t1, SHR(u, 48));
    u = g[s1[0] >> 20 & 15]; XOREQ(t0, SHL(u, 20)); XOREQ(t1, SHR(u, 44));
    u = g[s1[0] >> 24 & 15]; XOREQ(t0, SHL(u, 24)); XOREQ(t1, SHR(u, 40));
    u = g[s1[0] >> 28 & 15]; XOREQ(t0, SHL(u, 28)); XOREQ(t1, SHR(u, 36));
    u = g[s1[0] >> 32 & 15]; XOREQ(t0, SHL(u, 32)); XOREQ(t1, SHR(u, 32));
    u = g[s1[0] >> 36 & 15]; XOREQ(t0, SHL(u, 36)); XOREQ(t1, SHR(u, 28));
    u = g[s1[0] >> 40 & 15]; XOREQ(t0, SHL(u, 40)); XOREQ(t1, SHR(u, 24));
    u = g[s1[0] >> 44 & 15]; XOREQ(t0, SHL(u, 44)); XOREQ(t1, SHR(u, 20));
    u = g[s1[0] >> 48 & 15]; XOREQ(t0, SHL(u, 48)); XOREQ(t1, SHR(u, 16));
    u = g[s1[0] >> 52 & 15]; XOREQ(t0, SHL(u, 52)); XOREQ(t1, SHR(u, 12));
    u = g[s1[0] >> 56 & 15]; XOREQ(t0, SHL(u, 56)); XOREQ(t1, SHR(u,  8));
    u = g[s1[0] >> 60 & 15]; XOREQ(t0, SHL(u, 60)); XOREQ(t1, SHR(u,  4));

    /* round 1 */
    u = g[s1[1]       & 15]; XOREQ(t1, u);
    u = g[s1[1] >>  4 & 15]; XOREQ(t1, SHL(u,  4)); XOREQ(t2, SHR(u, 60));
    u = g[s1[1] >>  8 & 15]; XOREQ(t1, SHL(u,  8)); XOREQ(t2, SHR(u, 56));
    u = g[s1[1] >> 12 & 15]; XOREQ(t1, SHL(u, 12)); XOREQ(t2, SHR(u, 52));
    u = g[s1[1] >> 16 & 15]; XOREQ(t1, SHL(u, 16)); XOREQ(t2, SHR(u, 48));
    u = g[s1[1] >> 20 & 15]; XOREQ(t1, SHL(u, 20)); XOREQ(t2, SHR(u, 44));
    u = g[s1[1] >> 24 & 15]; XOREQ(t1, SHL(u, 24)); XOREQ(t2, SHR(u, 40));
    u = g[s1[1] >> 28 & 15]; XOREQ(t1, SHL(u, 28)); XOREQ(t2, SHR(u, 36));
    u = g[s1[1] >> 32 & 15]; XOREQ(t1, SHL(u, 32)); XOREQ(t2, SHR(u, 32));
    u = g[s1[1] >> 36 & 15]; XOREQ(t1, SHL(u, 36)); XOREQ(t2, SHR(u, 28));
    u = g[s1[1] >> 40 & 15]; XOREQ(t1, SHL(u, 40)); XOREQ(t2, SHR(u, 24));
    u = g[s1[1] >> 44 & 15]; XOREQ(t1, SHL(u, 44)); XOREQ(t2, SHR(u, 20));
    u = g[s1[1] >> 48 & 15]; XOREQ(t1, SHL(u, 48)); XOREQ(t2, SHR(u, 16));
    u = g[s1[1] >> 52 & 15]; XOREQ(t1, SHL(u, 52)); XOREQ(t2, SHR(u, 12));
    u = g[s1[1] >> 56 & 15]; XOREQ(t1, SHL(u, 56)); XOREQ(t2, SHR(u,  8));
    u = g[s1[1] >> 60 & 15]; XOREQ(t1, SHL(u, 60)); XOREQ(t2, SHR(u,  4));
    /* end */

    /* store result */
    _mm_storeu_si128((__m128i*)t, PXOR(t0, SHLD(t1, 64)));
    _mm_storeu_si128((__m128i*)(t+2), PXOR(t2, SHRD(t1, 64)));

#undef PNEG
#undef PSUB
#undef PZERO
#undef SHL
#undef SHR
#undef SHLD
#undef SHRD
#undef XOREQ
#undef PXOR
#undef PAND
}

#endif  /* GF2X_MUL2_H_ */
