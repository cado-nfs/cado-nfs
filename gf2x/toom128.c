/* This file is part of the gf2x library.

   Copyright 2007-2016
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

/* Variants of Toom-Cook using SSE2 128-bit XOR, copied from toom.c.
   Those routines assume that "unsigned long" is a 64-bit type. */

#include <stdio.h>
#include <string.h> /* for memcpy() */

#include "gf2x.h"
#include "gf2x/gf2x-impl.h"

#ifdef HAVE_KARAX
#include <emmintrin.h>

/* Let spx(n) be the space requirement (in number of 128-bit words) for stk
   in gf2x_mul_karax_internal(n), and sp(n) the space requirement (in number
   of 64-bit words) for stk in the gf2x_mul_kara() routine:
   (1) if 2*n < GF2X_MUL_KARA_THRESHOLD then spx(n) <= ceil(sp(2*n)/2)
   (2) otherwise spx(n) <= 3*ceil(n/2) + spx(ceil(n/2)).

   The parameter 'odd' can be 0 or 1:
   * if odd=0, then a and b have 2n words of 64 bits, the result c has 4n
               words of 64 bits, i.e., 2n words of 128 bits
   * if odd=1, then a and b have 2n-1 words of 64 bits, the result c has 4n-2
               words of 64 bits, i.e., 2n-1 words of 128 bits

   Assumes stk is 128-bit aligned.

   FIXME: write a 256-bit variant using AVX2:
   VPXOR: __m256i _mm256_xor_si256 ( __m256i a, __m256i b)
 */
#define PXOR(a,b) _mm_xor_si128(a,b)
static void
gf2x_mul_karax_internal (__m128i *c, const __m128i *a,
                         const __m128i *b, long n, __m128i * stk, int odd)
{
    __m128i *aa, *bb, *cc;
    long j, d, n2;

    /* since this routine is usually faster than gf2x_mul_kara(), it calls
       itself recursively instead of calling gf2x_mul_toom(), thus we have
       to stop the recursion by calling gf2x_mul_basecase() */
    if (2 * n - odd < GF2X_MUL_KARA_THRESHOLD)
      {
        gf2x_mul_basecase ((unsigned long*) c, (unsigned long*) a, 2 * n - odd,
                           (unsigned long*) b, 2 * n - odd);
        return;
      }

    n2 = (n + 1) / 2;		/* ceil(n/2) */
    d = n & 1;			/* 2*n2 - n = 1 if n odd, 0 if n even */
    aa = stk;			/* Size n2   */
    bb = aa + n2;		/* Size n2   */
    cc = bb + n2;		/* Size n2   */

    stk = cc + n2;		/* sp(n) = 3*ceil(n/2)) + sp(ceil(n/2)) */

    const __m128i *a1 = a + n2;	/* a[n2] */
    const __m128i *b1 = b + n2;	/* b[n2] */
    __m128i *c1 = c + n2;	/* c[n2]   */
    __m128i *c2 = c1 + n2;	/* c[2*n2] */
    __m128i *c3 = c2 + n2;	/* c[3*n2] */

    gf2x_mul_karax_internal (c, a, b, n2, stk, 0);	/* Low */

    /* now {c, 2*n2} contains A0 * B0 */

    gf2x_mul_karax_internal (c2, a1, b1, n2 - d, stk, odd); /* High */

    /* {c2, 2*(n2-d)} contains A1 * B1 */

    /* {c2, 2*(n2-d)} contains 2*(n2-d)-odd/2 words of 128 bits, i.e.,
       4*(n2-d)-odd words of 64 bits */

    /* now compute in {aa, n2} the sum A0+A1, and in {bb, n2} the sum B0+B1,
       and in parallel in {cc, n2} the n2 low 12-bit words of A0*B0 + A1*B1 */

    /* {a, n2} has n2 words of 128 bits, i.e., 2*n2 words of 64 bits
       {a1, n2-d} has n2-d-odd/2 words of 128 bits, i.e., 2*(n2-d)-odd
       words of 64 bits */

    for (j = 0; j < n2 - d - odd; j++)
      {
        // aa[j] = PXOR (a[j], a1[j]);
        // aa is 128-bit aligned, thus we can use aa[j] = ...
        aa[j] = PXOR (_mm_loadu_si128 (a + j), _mm_loadu_si128 (a1 + j));
        // bb[j] = PXOR (b[j], b1[j]);
        bb[j] = PXOR (_mm_loadu_si128 (b + j), _mm_loadu_si128 (b1 + j));
        // cc[j] = PXOR (c1[j], c2[j]);
        cc[j] = PXOR (_mm_loadu_si128 (c1 + j), _mm_loadu_si128 (c2 + j));
      }
    for (; j < n2 - d; j++) /* one loop only, and only when odd=1 */
      {
        /* zero the upper 64 bits of a1[j] and b1[j] */
        // aa[j] = PXOR (a[j], LOW(a1[j]));
        aa[j] = PXOR (_mm_loadu_si128 (a + j), _mm_loadl_epi64 (a1 + j));
        // bb[j] = PXOR (b[j], LOW(b1[j]));
        bb[j] = PXOR (_mm_loadu_si128 (b + j), _mm_loadl_epi64 (b1 + j));
        // cc[j] = PXOR (c1[j], c2[j]);
        cc[j] = PXOR (_mm_loadu_si128 (c1 + j), _mm_loadu_si128 (c2 + j));
      }
    for (; j < n2; j++)
      {	/* Only when n odd */
	// aa[j] = a[j];
        aa[j] = _mm_loadu_si128 (a + j);
	// bb[j] = b[j];
        bb[j] = _mm_loadu_si128 (b + j);
	// cc[j] = PXOR (c1[j], c2[j]);
        cc[j] = PXOR (_mm_loadu_si128 (c1 + j), _mm_loadu_si128 (c2 + j));
    }

    /* now we have:
       {c, n2} = low(A0*B0)
       {c3, n2-2*d} = high(A1*B1)
       {cc, n2} = high(A0*B0) + low(A1*B1) */

    gf2x_mul_karax_internal (c1, aa, bb, n2, stk, 0);	/* Middle */

    /* now we have:
       {c, n2} = low(A0*B0)
       {c1, 2*n2} = (A0+B0)*(A1+B1)
       {c3, n2-2*d} = high(A1*B1)
       {cc, n2} = high(A0*B0) + low(A1*B1)
       It remains to add {cc, n2} + {c, n2} to {c1, n2},
       and {cc, n2} + {c3, n2-2*d-odd} to {c2, n2} */

    for (j = 0; j < n2 - 2 * d - odd; j++)
      {
	// c1[j] = PXOR (c1[j], PXOR (cc[j], c[j]));
        _mm_storeu_si128 (c1 + j, PXOR (_mm_loadu_si128 (c1 + j), PXOR (cc[j], _mm_loadu_si128 (c + j))));
	// c2[j] = PXOR (c2[j], PXOR (cc[j], c3[j]));
        _mm_storeu_si128 (c2 + j, PXOR (_mm_loadu_si128 (c2 + j), PXOR (cc[j], _mm_loadu_si128 (c3 + j))));
      }
    for (; j < n2; j++)
      {	/* Only when n odd */
	// c1[j] = PXOR (c1[j], PXOR (cc[j], c[j]));
        _mm_storeu_si128 (c1 + j, PXOR (_mm_loadu_si128 (c1 + j), PXOR (cc[j], _mm_loadu_si128 (c + j))));
        // c2[j] = PXOR (c2[j], cc[j]);
	_mm_storeu_si128 (c2 + j, PXOR (_mm_loadu_si128 (c2 + j), cc[j]));
      }
}

void
gf2x_mul_karax (unsigned long *c, const unsigned long *a,
                const unsigned long *b, long n, unsigned long *stk)
{
  if (((uintptr_t) stk) % 16)
    stk ++; /* ensure stk is 128-bit aligned */
  gf2x_mul_karax_internal ((__m128i*) c, (__m128i*) a,
                           (__m128i*) b, (n + 1) >> 1, (__m128i*) stk, n & 1);
}

#endif /* HAVE_KARAX */
