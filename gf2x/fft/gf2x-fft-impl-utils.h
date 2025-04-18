/* This file is part of the gf2x library.

   Copyright 2019, 2023
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
#ifndef GF2X_FFT_IMPL_UTILS_H_
#define GF2X_FFT_IMPL_UTILS_H_

/* Assume wordlength  GF2X_WORDSIZE is 32 or 64 */

#include <stddef.h>

/* Don't define MIN, MAX, ABS as inlines, as they're already quite
 *
 * customarily defined as macros */
#ifndef MAX
#define MAX(h,i) ((h) > (i) ? (h) : (i))
#endif
#ifndef MIN
#define MIN(h,i) ((h) < (i) ? (h) : (i))
#endif
#ifndef ABS
#define ABS(h) ((h) < 0 ? -(h) : (h))
#endif


/** Some support functions. **/

/* CEIL(a,b) = ceiling(a/b) */
static inline size_t CEIL(size_t a, size_t b)
{
    return ((a)+(b)-1)/(b);
}

/* W(b) is the number of words needed to store b bits */
static inline size_t W(size_t b)
{
    return CEIL(b, GF2X_WORDSIZE);
}

/* I(b) is the index word of bit b, assuming bits 0..GF2X_WORDSIZE-1
   have index 0 */
static inline size_t I(size_t b)
{
    return b / GF2X_WORDSIZE;
}

static inline size_t R(size_t b)
{
    return b % GF2X_WORDSIZE;
}

static inline size_t R2(size_t b)       /* remaining bits */
{
    return (-b) % GF2X_WORDSIZE;
}

static inline unsigned long MASK(size_t x)
{
    ASSERT(x < GF2X_WORDSIZE);
    return ((1UL << (x)) - 1UL);
}

/* GETBIT(a,i)   gets the i-th bit of the bit-array starting at a[0],
   XORBIT(a,i,x) xors this bit with the bit x, where x = 0 or 1.   */
static inline unsigned long GETBIT(unsigned long *a, size_t i)
{
    return (a[I(i)] >> R(i)) & 1UL;
}

static inline void XORBIT(unsigned long *a, size_t i, unsigned long x)
{
    ASSERT((x & ~1UL) == 0);
    a[I(i)] ^= x << R(i);
}

static inline void Copy(unsigned long *a, const unsigned long *b, size_t n)
{
    memcpy(a, b, n * sizeof(unsigned long));
}

static inline void Zero(unsigned long *a, size_t n)
{
    memset(a, 0, n * sizeof(unsigned long));
}

static inline void Clear(unsigned long *a, size_t low, size_t high)
{
    if (high > low)
       memset (a + low, 0, (high - low) * sizeof(unsigned long));
}

/** Now the specific things */

/* a <- b + c */

static inline void AddMod(unsigned long *a, unsigned long *b, unsigned long *c,
		   size_t n)
{
    for (size_t i = 0; i < n; i++)
       a[i] = b[i] ^ c[i];
}

/* a <- b + c + d */
static inline void
AddMod3 (unsigned long *a, unsigned long *b, unsigned long *c,
         unsigned long *d, size_t n)
{
  for (size_t i = 0; i < n; i++)
    a[i] = b[i] ^ c[i] ^ d[i];
}

/* c <- a * x^k, return carry out, 0 <= k < GF2X_WORDSIZE */
static inline unsigned long
Lsh (unsigned long *c, unsigned long *a, size_t n, size_t k)
{
    if (k == 0) {
	if (c != a)
	    Copy(c, a, n);
	return 0;
    }

    /* {c, n} and {a, n} should not overlap */
    ASSERT(c <= a || a + n <= c);
    ASSERT(k > 0);

    unsigned long t, cy = 0;
    for (size_t i = 0; i < n; i++) {
	t = (a[i] << k) | cy;
	cy = a[i] >> (GF2X_WORDSIZE - k);
	c[i] = t;
    }
    return cy;
}

/* c <- c + a * x^k, return carry out, 0 <= k < GF2X_WORDSIZE */

static inline unsigned long AddLsh(unsigned long *c, unsigned long *a, size_t n,
			    size_t k)
{
    unsigned long t, cy = 0;

    if (k == 0) {
	AddMod(c, c, a, n);
	return 0;
    }

    /* {c, n} and {a, n} should not overlap */
    ASSERT(c <= a || a + n <= c);

    ASSERT(k > 0);

    for (size_t i = 0; i < n; i++) {
	t = (a[i] << k) | cy;
	cy = a[i] >> (GF2X_WORDSIZE - k);
	c[i] ^= t;
    }
    return cy;
}

/* c <- a / x^k, return carry out, 0 <= k < GF2X_WORDSIZE */
static inline unsigned long
Rsh (unsigned long *c, const unsigned long *a, size_t n, size_t k)
{
    if (k == 0) {
	if (c != a)
	    Copy(c, a, n);
	return 0;
    }

    ASSERT(k > 0);

    unsigned long t, cy = 0;
    for (size_t i = n; i-- ; ) {
	t = (a[i] >> k) | cy;
	cy = a[i] << (GF2X_WORDSIZE - k);
	c[i] = t;
    }
    return cy;
}

/* c <- c + a / x^k, return carry out, 0 <= k < GF2X_WORDSIZE */

static inline unsigned long AddRsh(unsigned long *c, unsigned long *a, size_t n,
			    size_t k)
{
    unsigned long t, cy = 0;

    if (k == 0) {
	AddMod(c, c, a, n);
	return 0;
    }

    ASSERT(k > 0);

    for (size_t i = n ; i-- ; ) {
	t = (a[i] >> k) | cy;
	cy = a[i] << (GF2X_WORDSIZE - k);
	c[i] ^= t;
    }
    return cy;
}

/* Copy bits_c bit from c1, starting from position shift. High bits of
 * the last word of c are cleared.
 */
static inline void CopyBitsRsh(unsigned long * c, const unsigned long * c1, size_t bits_c, size_t shift)
{
    if (shift == 0) {
        Copy(c, c1, W(bits_c));
    } else {
        size_t t = W(bits_c);
        size_t words_full = W(bits_c + shift);
        size_t pick = I(shift);
        size_t cnt = R(shift);
        size_t tnc = GF2X_WORDSIZE - cnt;
        Rsh(c, c1 + pick, t, cnt);
        /* words_full - pick - t is either 0 or 1 */
        if (words_full - pick == t + 1)
            c[t - 1] |= c1[pick + t] << tnc;
    }
    if (R(bits_c))
        c[I(bits_c)] &= MASK(R(bits_c));
}

#endif	/* GF2X_FFT_IMPL_UTILS_H_ */
