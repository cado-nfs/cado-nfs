/* This file is part of the gf2x library.

   Copyright 2007, 2008, 2009
   Richard Brent, Pierrick Gaudry, Emmanuel Thome', Paul Zimmermann

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 2.1 of the License, or (at
   your option) any later version.
   
   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
   License for more details.
   
   You should have received a copy of the GNU Lesser General Public
   License along with CADO-NFS; see the file COPYING.  If not, write to
   the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
   Boston, MA 02110-1301, USA.
*/
/* Multiplication over GF(2)[x] using Fast Fourier Transform.
 * Bit-aligned version. */

#include <stdio.h>
#include <stdint.h> /* for uint64_t in Lshift */
#include <stdlib.h>
#include <string.h>

#include <stdarg.h> /* for the fft interface */

#include "macros.h"
#include "gf2x.h"
#include "gf2x/gf2x-impl.h"
#include "gf2x/gf2x-tfft.h"

// #define DEBUG
// #define DEBUG_LSHIFT
// #define DEBUG_MULMOD
// #define VERBOSE
// #define TIMING

/* The check in tunefft finds that gf2x_mul_fft does not work correctly for
   balanced operands of size <= 13 words.  Hence this threshhold. */


#ifndef MUL_FFT_THRESHOLD
#define MUL_FFT_THRESHOLD 28
#endif

#if (MUL_FFT_THRESHOLD < 28)
#error "MUL_FFT_THRESHOLD too small, should be at least 28"
#endif

/* Assume wordlength  WLEN is 32 or 64 */

# define WLEN GF2X_WORDSIZE

/** Some support functions. **/

/* CEIL(a,b) = ceiling(a/b) */
static inline size_t CEIL(size_t a, size_t b)
{
    return ((a)+(b)-1)/(b);
}

/* W(b) is the number of words needed to store b bits */
static inline size_t W(size_t b)
{
    return CEIL(b, WLEN);
}

/* I(b) is the index word of bit b, assuming bits 0..WLEN-1
   have index 0 */
static inline size_t I(size_t b)
{
    return b / WLEN;
}

static inline size_t R(size_t b)
{
    return b % WLEN;
}

static inline size_t R2(size_t b)       /* remaining bits */
{
    return (-b) % WLEN;
}

static inline unsigned long MASK(size_t x)
{
    ASSERT(x < WLEN);
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

/* Don't define MIN, MAX, ABS as inlines, as they're already quite
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

static void *malloc_or_die(size_t size)
{
    void *res = malloc(size);
    if (res == NULL)
	abort();
    return res;
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
	memset(a + low, 0, (high - low) * sizeof(unsigned long));
}

/** Now the specific things */

/* a <- b + c */

static void AddMod(unsigned long *a, unsigned long *b, unsigned long *c,
		   size_t n)
{
    for (size_t i = 0; i < n; i++)
	a[i] = b[i] ^ c[i];
}

/* c <- a * x^k, return carry out, 0 <= k < WLEN */

static unsigned long Lsh(unsigned long *c, unsigned long *a, size_t n, size_t k)
{
    unsigned long t, cy = 0;

    if (k == 0) {
	if (c != a)
	    Copy(c, a, n);
	return 0;
    }

    /* {c, n} and {a, n} should not overlap */
    ASSERT(c <= a || a + n <= c);

    ASSERT(k > 0);

    for (size_t i = 0; i < n; i++) {
	t = (a[i] << k) | cy;
	cy = a[i] >> (WLEN - k);
	c[i] = t;
    }
    return cy;
}

/* c <- c + a * x^k, return carry out, 0 <= k < WLEN */

static unsigned long AddLsh(unsigned long *c, unsigned long *a, size_t n,
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
	cy = a[i] >> (WLEN - k);
	c[i] ^= t;
    }
    return cy;
}

/* c <- a / x^k, return carry out, 0 <= k < WLEN */

static unsigned long Rsh(unsigned long *c, const unsigned long *a, size_t n,
			 size_t k)
{
    unsigned long t, cy = 0;

    if (k == 0) {
	if (c != a)
	    Copy(c, a, n);
	return 0;
    }

    ASSERT(k > 0);

    for (size_t i = n; i-- ; ) {
	t = (a[i] >> k) | cy;
	cy = a[i] << (WLEN - k);
	c[i] = t;
    }
    return cy;
}

/* c <- c + a / x^k, return carry out, 0 <= k < WLEN */

static unsigned long AddRsh(unsigned long *c, unsigned long *a, size_t n,
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
	cy = a[i] << (WLEN - k);
	c[i] ^= t;
    }
    return cy;
}

#if 0
// unused
/* c <- x^j * a / x^k */

static unsigned long LshRsh(unsigned long *c, unsigned long *a, size_t n,
			    size_t j, size_t k)
{
    ASSERT(n >= 0);
    if (j > k)
	Lsh(c, a, n, j - k);
    else
	Rsh(c, a, n, k - j);
}

#endif

#if (defined(DEBUG) || defined(DEBUG_LSHIFT) || defined(DEBUG_MULMOD))

static void dump(const unsigned long *a, size_t n)
{
    for (size_t i = 0; i < n; i++) {
	printf("+%lu*X^%lu", a[i], i);
	if ((i + 1) % 3 == 0)
	    printf("\n");
    }
    printf(":\n");
}
#endif

/* a <- b * x^k mod x^(2*N)+x^N+1.
   Assume a and b do not overlap.
*/
static void Lshift(unsigned long *a, unsigned long *b, uint64_t k, size_t N)
{
    size_t r, h, l, n, ih, il;
    unsigned long s1=0, s2=0; // hush gcc

    n = W(2 * N);

#ifdef DEBUG_LSHIFT
    printf("R:=x^%u+x^%u+1: k:=%ld:\nb:=", 2 * N, N, k);
    dump(b, n);
#endif

    k = k % (3 * N);

    ASSERT(k < 3 * N);

    if (k == 0) {
	if (a != b)
	    Copy(a, b, n);
    } else if (k <= N) {
	/*  ------------------------------------------
	   |  L0  |      L1      |  L2  |      H      |
	   ------------------------------------------
	   ------------------------------------------
	   |      H      |  L0  |    L1+H      |  L2  |
	   ------------------------------------------
	   L0 has l bits, L1 has h bits, L2 has l bits, H has h bits
	 */
	h = k;			/* 1 <= h <= N */
        ASSERT(h <= N);
        ASSERT(1 <= h);
	l = N - h;		/* 0 <= l < N */
	/* A <- 0:L0:L1:L2 */
	s1 = Lsh(a + I(h), b, W(N + l), R(h));
	/* b[W(N + l)-1] has R2(N+l) bits from H, thus s1 has no bit from L
	   if R(h) <= R2(N+l), and R(h)-R2(N+l) bits from L otherwise */
	if (R(h) > R2(N + l))
	    a[I(h) + W(N + l)] = s1 & MASK(R(h) - R2(N + l));
	if (n - I(N + l) > I(h))
	    s1 = a[I(h)];
	if (n - I(N + l) > I(h) + 1)
	    s2 = a[I(h) + 1];
	ASSERT(n - I(N + l) <= I(h) + 2);
	/* copy H in low words: will clobber a[0]..a[n-I(N+l)-1] */
	Rsh(a, b + I(N + l), n - I(N + l), R(N + l));
	if (n - I(N + l) > I(h))
	    a[I(h)] ^= s1;
	if (n - I(N + l) > I(h) + 1)
	    a[I(h) + 1] ^= s2;
	/* we now have a = H:L0:L1:L2 */
	/* add H in high words */
	if (R(N + l) <= R(N)) {
	    s1 = AddLsh(a + I(N), b + I(N + l), n - I(N + l),
			R(N) - R(N + l));
	    /* we have to remove the high R(N+l) bits from L2 */
	    a[I(N)] ^= (b[I(N + l)] & MASK(R(N + l))) << (R(N) - R(N + l));
	    if (s1)
		a[I(N) + n - I(N + l)] ^= s1;
	} else {
	    /* the low word from H in b contains R(N + l) bits from L2, and the
	       low word from L1+H in a contains R(N) bits from L0. After the
	       shift, it remains R(N) bits from L2 */
	    AddRsh(a + I(N), b + I(N + l), n - I(N + l), R(N + l) - R(N));
	    if (R(N))
		a[I(N)] ^=
		    (b[I(N + l)] & MASK(R(N + l))) >> (R(N + l) - R(N));
	}
    } else if (k <= 2 * N) {
	/*  ------------------------------------------
	   |  L   |      H0      |  H1  |     H2      |
	   ------------------------------------------
	   ------------------------------------------
	   |     H0+H2   |  H1  |     H0      |  L+H1 |
	   ------------------------------------------
	   L has l bits, H0 has h bits, H1 has l bits, H2 has h bits
	 */
	h = k - N;		/* 1 <= h <= N */
	l = N - h;		/* 0 <= l < N */
	r = R(l);
	Rsh(a, b + I(l), W(N + l) - I(l), r);
	/* mask high bits */
	if (R(N))
	    a[I(N)] &= MASK(R(N));
	/* now we have a = H0:H1 */
	if (R(N) > 0)
	    s1 = a[I(N)];
	// ASSERT(I(N) >= 0);
	s2 = Lsh(a + I(N), a, I(N), R(N));
	if (R(N) > 0) {
	    a[I(N)] ^= s1;	/* restore high R(N) bits */
	    a[2 * I(N)] = s2 ^ (s1 << R(N));
	    if (2 * R(N) > WLEN)
		a[2 * I(N) + 1] = s1 >> (WLEN - R(N));
	}
	/* now we have a = H0:H1:H0:H1 */
	AddRsh(a, b + I(N + l), W(2 * N) - I(N + l), R(N + l));
	/* now we have a = H0+H2:H1:H0:H1 */
	s1 = AddLsh(a + I(N + h), b, W(l), R(N + h));
	if (r > 0)		/* mask shifted low bits from H0 */
	    a[I(N + h) + I(l)] ^= (b[I(l)] & ~MASK(R(l))) << R(N + h);
	/* b[W(l)-1] contains R2(l) bits from H0, thus if R(N+h) > R2(l),
	   then s1 contains R(N+h)-R2(l) bits from L */
	if (R(N + h) > R2(l))
	    a[I(N + h) + W(l)] ^= s1 & MASK(R(N + h) - R2(l));
    } else {			/* 2*N < k < 3*N */

	/*  ------------------------------------------
	   |  L   |      H0      |  H1  |     H2      |
	   ------------------------------------------
	   ------------------------------------------
	   |     H0+H2   |  H1  |     H0      |  L+H1 | after k-N
	   ------------------------------------------
	   ------------------------------------------
	   |    H0       | L+H1 |     H2      |   L   |
	   ------------------------------------------
	   L has l bits, H0 has h bits, H1 has l bits, H2 has h bits
	 */
	h = k - 2 * N;		/* 1 <= h < N */
	l = N - h;		/* 0 < l < N */
	il = W(2 * N) - I(l);
	Rsh(a, b + I(l), il, R(l));
	/* we now have a = H0:H1:H2:0 */
	s1 = AddLsh(a + I(h), b, W(l), R(h));	/* H0:H1+L:H2:0 */
	if (R(h) > R2(l))
	    /* add least R(h) - R2(l) bits from s1 */
	    a[I(h) + W(l)] ^= s1 & MASK(R(h) - R2(l));
	else if (R2(l) > R(h))
	    /* mask low bits from H0 */
	    a[I(h) + W(l) - 1] ^= (b[W(l) - 1] & ~MASK(R(l))) << R(h);
	/* we now have a = H0:H1+L:H2:0 */
	ih = I(N + h);		/* index of lowest word from L in a */
	if (R(N + h) > 0)	/* H2 and L share a common word a[ih] */
	    s2 = a[ih];
	// ASSERT(W(l) >= 0);
	s1 = Lsh(a + ih, b, W(l), R(N + h));	/* H0:H1+L:H2:L */
	/* b[W(l)-1] contains R2(l) bits of H0, thus if R(N+h) > R2(l),
	   s1 contains R(N+h) - R2(l) bits of L */
	if (R(N + h) > R2(l))
	    a[ih + W(l)] = s1 & MASK(R(N + h) - R2(l));
	if (R(N + h) > 0)
	    a[ih] ^= s2;
    }
    /* mask high bits of result */
    r = R(2 * N);
    if (r > 0)
	a[n - 1] &= MASK(r);
#ifdef DEBUG_LSHIFT
    printf("a:=");
    dump(a, n);
#endif
}

/* a <- b * c mod x^(2*N)+x^(N)+1.
   Assumes t has space for 2n words, and u for gf2x_toomspace(n) words,
   where n = ceil(2N/WLEN).
   a and b may be equal.
   a must have space for n words.
*/

static void MulMod(unsigned long *a, const unsigned long *b, const unsigned long *c,
		   size_t N, unsigned long *t, unsigned long *u)
{
    size_t n = W(2 * N), sh, l;

    /* FIXME: in practice N is divisible by a multiple of 3, thus if WLEN is
       a power of two, R(N) > 0 and W(N) = I(N) + 1, thus we can avoid a few
       tests below. */

#ifdef DEBUG_MULMOD
    printf("b:=");
    dump(b, n);
    printf("c:=");
    dump(c, n);
#endif
    gf2x_mul_toom(t, b, c, n, u);	/* t = | L0:N | L1:N | H0:N | H1:N | */
    /* t mod x^(2*N)+x^(N)+1 = | L0+H0+H1:N | L1+H0:N | */
    l = W(N);
    Rsh(a, t + I(2 * N), W(3 * N) - I(2 * N), R(2 * N));	/* H0 */
    if (R(N) > 0)		/* mask low bits from H1 if necessary */
	a[I(N)] &= MASK(R(N));
    if (I(N) < l)		/* necessarily l = I(N)+1 */
	u[0] = a[l - 1];
    // ASSERT(I(N) >= 0);
    u[1] = Lsh(a + I(N), a, I(N), R(N));	/* H0:H0 */
    /* u[0] and a[1] have R(N) bits */
    if (I(N) < l) {
	ASSERT(2 * l - 2 < n);
	a[2 * l - 2] = u[1] ^ (u[0] << R(N));	/* 2R(N) bits */
	if (2 * R(N) > WLEN) {
	    ASSERT(2 * l - 1 < n);
	    a[2 * l - 1] = u[0] >> (WLEN - R(N));
	}
	a[l - 1] ^= u[0];	/* restore low R(N) bits */
    }
    AddMod(a, a, t, W(2 * N));	/* H0+L0:H0+L1 */
    AddRsh(a, t + I(3 * N), W(4 * N) - I(3 * N), R(3 * N));	/* L0+H0+H1:H0+L1 */
    /* mask high bits */
    sh = R(2 * N);
    if (sh > 0)
	a[n - 1] &= MASK(sh);
#ifdef DEBUG_MULMOD
    printf("a:=");
    dump(a, W(2 * N));
#endif
}

/* starting from i=0, j=0, K and Z=1, enumerates all pairs (i,j)
   for 0 <= i < K such that j = bitrev(i).
*/

static void bitrev(size_t i, size_t j, size_t K, size_t Z, size_t *perm)
{
    if (K == 1)
	perm[i] = j;
    else {
	bitrev(i, j, K / 3, 3 * Z, perm);
	bitrev(i + K / 3, j + Z, K / 3, 3 * Z, perm);
	bitrev(i + 2 * K / 3, j + 2 * Z, K / 3, 3 * Z, perm);
    }
}

/* performs an FFT of length K on A[0], A[stride], A[(K-1)*stride] with
   omega=x^j as root of unity, where all computations are
   done modulo x^(2Np) + x^Np + 1.
   Each A[i] has space for 2np words where np=ceil(Np/WLEN).
   Assume omega^K = 1 mod x^(2Np) + x^Np + 1, i.e., mod x^(3Np)+1.
   t1, t2, t3 are buffers of 2np words.
   Assumes 0 <= j < 3*Np.
   For 0 <= i < K, p[i] = tritrev(i), i.e., the integer j such that the
   base-3 representation of i and j are reverse one from each other.
*/

static void fft(unsigned long **A, uint64_t K, uint64_t j, size_t Np, size_t stride,
		unsigned long *t1, unsigned long *t2, unsigned long *t3,
		size_t *p)
{
    ASSERT(j < 3 * Np);

    if (K == 1)
	return;

    size_t i, k = K / 3, twonp = W(2 * Np);
    uint64_t ii;

    fft(A, k, (3 * j) % (3 * Np), Np, 3 * stride, t1, t2, t3, p);
    fft(A + stride, k, (3 * j) % (3 * Np), Np, 3 * stride, t1, t2, t3, p);
    fft(A + 2 * stride, k, (3 * j) % (3 * Np), Np, 3 * stride, t1, t2, t3, p);

#define a A[3*i*stride]
#define b A[(3*i+1)*stride]
#define c A[(3*i+2)*stride]
    for (i = 0; i < k; i++) {
	ii = p[3 * stride * i];	/* bitrev(i,K/3) = bitrev(3*stride*i,K) */
	Lshift(t1, b, ii * j, Np);	/* t1 = w^ii*b */
	Lshift(t2, b, (ii + 2 * K / 3) * j, Np);	/* t2 = w^(ii+2K/3)*b */
	Lshift(t3, b, (ii + K / 3) * j, Np);
	AddMod(b, a, t3, twonp);
	Lshift(t3, c, (2 * ii + 2 * K / 3) * j, Np);
	AddMod(b, b, t3, twonp);
	Lshift(t3, c, (2 * ii + 4 * K / 3) * j, Np);
	AddMod(t2, t2, t3, twonp);
	Lshift(t3, c, 2 * ii * j, Np);	/* t3 = w^(2ii)*c */
	AddMod(c, a, t2, twonp);
	AddMod(a, a, t1, twonp);
	AddMod(a, a, t3, twonp);
    }
#undef a
#undef b
#undef c
}

/* allocate A[0]...A[K-1], and put there {a, an} cut into K chunks of M bits;
   return pointer to block of memory containing A[0]...A[K-1] (to be freed
   by the calling routine) */

static void decompose(unsigned long **A, const unsigned long *a,
				size_t an, size_t M, size_t K, size_t np)
{
    size_t i, j, k, l, sh = R(M);

    for (i = 0, j = 0, k = 0; i < K; i++) {
	/* invariant: we have already used k bits of a[j] */
	l = W(k + M);		/* number of words for k+M bits */
	if (j + l > an)
	    l = (an > j) ? an - j : 0;	/* remains l words a[j]..a[j+l-1] */
	if (l > 0) {
	    ASSERT(j + l <= an);
	    Rsh(A[i], a + j, l, k);
	    /* mask last bits */
	    if (sh) {
		ASSERT(I(M) < 2 * np);
		A[i][I(M)] &= MASK(sh);
	    }
	    if (l > W(M))
		l = W(M);
	}
	ASSERT(l < 2 * np);
	Zero(A[i] + l, 2 * np - l);
	k += M;
	j += k / WLEN;
	k %= WLEN;
    }
}

static void recompose(unsigned long * c, size_t cn, unsigned long **C, size_t K, size_t M, size_t Np)
{
    // size_t np = W(Np);	       	// Words to store Np bits

    /* reconstruct C = sum(C[i]*X^(M*i) mod x^N+1.
       We first compute sum(C[i]*X^(M*i), then reduce it using the wrap function.
       Since we know the result has at most cn words, any value exceeding cn
       words is necessarily zero.
       Each C[i] has 2*Np bits, thus the full C has (K-1)*M+2*Np
       = N - M + 2*Np >= N + Np >= 2*n*WLEN + Np bits.
       Thus exactly 2*Np-M bits wrap around mod x^N+1.
     */

    size_t l = 2 * Np - M;/* number of overlapping bits with previous C[i] */
    size_t i, j, k;
    size_t j1, k1;
    size_t z;

    for (i = 0, j = 0, k = 0, j1 = I(l), k1 = R(l); i < K; i++) {
	// unsigned long cy;
	/* invariants:
	   - first bit of C[i] is bit k of c[j]
	   - first bit of C[i] non overlapping with C[i-1] is bit k1 of c[j1] */
	/* add C[i] shifted by i*M to c, where K*M = N, and Np >= M */
	if (i == 0)		/* FIXME: we could set C[0] to c to avoid this copy */
	    Copy(c, C[i], W(2 * Np) < cn ? W(2 * Np) : cn);
	else {
	    /* we have already set bit k of c[j] up to bit k1 of c[j1]
	       (excluded), i.e., words c[j] up to c[j1 - (k1 == 0)] */
	    z = j1 + (k1 != 0) - j;	/* number of overlapping words */
	    /* first treat the high (non overlapping) words of C[i], i.e.,
	       {C[i] + z, W(2*Np) - z} */
	    if (j + W(2 * Np) < cn) {
		if (z < W(2 * Np)) {
		    // ASSERT(W(2 * Np) - z >= 0);
		    ASSERT((j + W(2 * Np) < cn));
		    c[j + W(2 * Np)] =
			Lsh(c + j + z, C[i] + z, W(2 * Np) - z, k);
		} else if (z == W(2 * Np)) {	/* all words overlap with C[i-1] */
		    ASSERT((j + W(2 * Np) < cn));
		    c[j + W(2 * Np)] = 0UL;
		}
	    } else if (j + z < cn) {
		Lsh(c + j + z, C[i] + z, cn - j - z, k);
	    }
	    /* then deal with the low bits of C[i], overlapping with C[i-1] */
	    if (j + z < cn)
		c[j + z] ^= AddLsh(c + j, C[i], z, k);
	    else if (j < cn)
		AddLsh(c + j, C[i], cn - j, k);
	}

	k += M;
	j += k / WLEN;
	k %= WLEN;
	k1 += M;
	j1 += k1 / WLEN;
	k1 %= WLEN;
    }
}

static inline size_t compute_Np(size_t M, size_t K)
{
    size_t Mp = CEIL(M, K / 3);	// ceil(M/(K/3))
    size_t Np = Mp * (K / 3);	// Np >= M, Np multiple of K/3
    return Np;
}

static inline size_t compute_np(size_t M, size_t K)
{
    size_t np = W(compute_Np(M,K));// Words to store Np bits
    return np;
}

// Wraps the polynomial represented by c mod x^N + 1
// Assumes wraps at most once, i.e. deg(c) < 2N.
// The high part of c (bits N to WLEN*cn) are cleared.
// RPB 20070429

static void wrap(unsigned long *c, size_t bits_c, size_t N)
{
    size_t i;
    size_t Nw = I(N);
    size_t Nb = R(N);
    size_t Nbc = WLEN - Nb;
    size_t cn = W(bits_c);

    // Perhaps most of this could be done by a call to AddLsh ?

    if (N >= bits_c)
        return;

    // xor bits N .. WLEN*cn of c to c[0...]
    if (Nb == 0) {
        for (i = 0; i < cn - Nw - 1; i++)
            c[i] ^= c[i + Nw];
    } else {
        for (i = 0; i < cn - Nw - 1; i++)
            c[i] ^= (c[i + Nw] >> Nb) | (c[i + Nw + 1] << Nbc);
    }
    ASSERT(cn > Nw);
    c[cn - Nw - 1] ^= (c[cn - 1] >> Nb);

    // Now clear remaining bits of c

    c[Nw] &= MASK(Nb);	// Clear high Nbc bits of c[Nw]
    Clear(c, Nw + 1, cn);
}

static void split_reconstruct(unsigned long * c, unsigned long * c1, unsigned long * c2, size_t cn, size_t K, size_t m1)
{
    size_t n = WLEN * cn;	// Max bit-size of full product
    size_t delta = K;		// delta = n1 - n2;
    size_t jw, jn1w, jn1b, jn1bc, jdw, jdb, jdbc;
    unsigned long t, next;

// Now extract the result. First do a partial word bit-by-bit.

    size_t m2 = m1 - 1;

    size_t n2 = K * m2;		// n2 smallest possible multiple of K
    size_t n1 = K * m1;		// next possible multiple of K
    size_t j;

    for (j = n - n1 - 1; (j % WLEN) != (WLEN - 1); j--) {
	t = GETBIT(c1, j + delta) ^ GETBIT(c2, j + delta);
	XORBIT(c1, j + n1, t);	// XOR assumes high part of c1 was zero
	XORBIT(c1, j, t);
    }

// Now do the rest using full-word operations.

    j -= WLEN - 1;
    jdb = R(j + delta);
    jdbc = WLEN - 1 - jdb;
    jn1b = R(j + n1);
    jw = I(j);
    jdw = I(j + delta);
    jn1w = I(j + n1);
    next = c1[jdw + 1] ^ c2[jdw + 1];

    if (jn1b == 0) {		// Unusual case
	for (jw++ ; jw-- ; jdw--, jn1w--) {
	    t = (next << 1) << jdbc;
	    next = c1[jdw] ^ c2[jdw];
	    t ^= next >> jdb;
	    c1[jw] ^= t;
	    c1[jn1w] = t;
	}
    } else {			// Usual case
	for (jn1bc = WLEN - jn1b, jw++; jw-- ; jdw--, jn1w--) {
	    t = (next << 1) << jdbc;
	    next = c1[jdw] ^ c2[jdw];
	    t ^= next >> jdb;
	    c1[jw] ^= t;
	    c1[jn1w] ^= t << jn1b;
	    c1[jn1w + 1] ^= t >> jn1bc;
	}
    }

// Do a consistency check. This is cheap and detects most errors.
// If DEBUG defined we check the first delta bits, otherwise we only
// check the first WLEN bits.

#ifdef DEBUG
    for (j = 0; j < delta; j++) {
	if ((GETBIT(c2, j) ^ GETBIT(c1, j) ^ GETBIT(c1, j + n2)) != 0) {
	    printf("Consistency check failed in gf2x_mul_fft2, bit %ld\n", j);
	    exit(1);
	}
    }
#endif

    t = c2[0] ^ c1[0] ^ (c1[n2 / WLEN] >> n2 % WLEN) ^
	((c1[n2 / WLEN + 1] << 1) << (WLEN - 1 - n2 % WLEN));
    if (t != 0) {
	printf("Consistency check failed in gf2x_mul_fft2, low word %lx\n", t);
	exit(1);
    }

    Copy(c, c1, cn);	// Copy result
}


/** now the external calls for the gf2x_tfft interface **/

size_t gf2x_tfft_size(gf2x_tfft_info_srcptr o)
{
    size_t K = o->K;
    if (K == 0) {
        /* special fall-back case */
        return W(o->bits_a) + W(o->bits_b);
    } else if (!o->split) {
        return 2 * K * compute_np(o->M, K);
    } else {    /* K < 0 : FFT split in two. */
        /* m2 is just m1-1, so to make things simpler, we claim that the
         * size is made only of blocks corresponding to size m1. But yes,
         * half with m1 and half with m2 would make sense. But wouldn't
         * change the picture so much anyway.
         */
        return 4 * o->K * compute_np(o->M, o->K);
    }
}

void gf2x_tfft_zero(gf2x_tfft_info_srcptr o, gf2x_tfft_ptr ptr, size_t n)
{
    memset(ptr, 0, n * gf2x_tfft_size(o) * sizeof(gf2x_tfft_t));
}

void gf2x_tfft_cpy(gf2x_tfft_info_srcptr o, gf2x_tfft_ptr y, gf2x_tfft_srcptr x)
{
    memcpy(y, x, gf2x_tfft_size(o) * sizeof(gf2x_tfft_t));
}

gf2x_tfft_ptr gf2x_tfft_get(gf2x_tfft_info_srcptr o, gf2x_tfft_ptr ptr, size_t k)
{
    return ptr + k * gf2x_tfft_size(o);
}

gf2x_tfft_ptr gf2x_tfft_alloc(gf2x_tfft_info_srcptr o, size_t n)
{
    return malloc_or_die(n * gf2x_tfft_size(o) * sizeof(gf2x_tfft_t));
}

void gf2x_tfft_free(gf2x_tfft_info_srcptr o MAYBE_UNUSED, gf2x_tfft_ptr ptr, size_t n MAYBE_UNUSED)
{
    free(ptr);
}

static void gf2x_tfft_dft_inner(gf2x_tfft_info_srcptr o, gf2x_tfft_ptr tr, const unsigned long * a, size_t bits_a, size_t M)
{
    size_t K = o->K;
    size_t Mp = CEIL(M, K / 3);	// ceil(M/(K/3))
    size_t Np = Mp * (K / 3);	// Np >= M, Np multiple of K/3
    size_t np = W(Np);	       	// Words to store Np bits

    // allocate the array of pointers. It's just temporary stuff.
    unsigned long ** A = malloc_or_die(K * sizeof(unsigned long *));
    for (size_t i = 0; i < K; i++) A[i] = tr + 2 * i * np;
    decompose(A, a, W(bits_a), M, K, np);
    unsigned long * tmp1, * tmp2, * tmp3;
    tmp1 = o->tmp;
    tmp2 = o->tmp + 2 * np;
    tmp3 = o->tmp + 4 * np; /* max(2np,gf2x_toomspace(2np)) words */
    fft(A, K, Mp, Np, 1, tmp1, tmp2, tmp3, o->perm);
    free(A);
}

static void gf2x_tfft_dft_inner_split(gf2x_tfft_info_srcptr o, gf2x_tfft_ptr tr, const unsigned long * a, size_t bits_a, size_t M, unsigned long * buf, size_t bufsize)
{
    size_t K = o->K;
    size_t N = K * M;

    // ASSERT(K >= WLEN);  // why ? I've seen this comment somewhere, but why ?
    // ASSERT(2 * K * M >= W(bits_a + bits_b));

    // FIXME: This wrapping, and use of extra buffer space, should be
    // merged into decompose().
    Copy(buf, a, W(bits_a));
    Clear(buf, W(bits_a), bufsize);		// Clear upper part of a
    wrap(buf, bits_a, N);
    gf2x_tfft_dft_inner(o, tr, buf, MIN(N, bits_a), M);
}

/* bits_a is a number of BITS */
void gf2x_tfft_dft(gf2x_tfft_info_srcptr o, gf2x_tfft_ptr tr, const unsigned long * a, size_t bits_a)
{
    // bits_a is a number of BITS.
    if (o->K == 0) {
        Copy(tr, a, W(bits_a));
        /* zeroing out the bits isn't really needed. */
        Clear(tr, W(bits_a), W(o->bits_a) + W(o->bits_b));
    } else if (!o->split) {
        gf2x_tfft_dft_inner(o, tr, a, bits_a, o->M);
    } else {
        size_t m1 = o->M;
        size_t m2 = o->M - 1;
        size_t K = o->K;
        /* there's some work to be done prior to doing the decomposition:
         * wrapping.  We need some temporary space for that.  */

        size_t bufsize = MAX(W(bits_a), W((size_t) m1));

        unsigned long * buf = malloc_or_die(bufsize * sizeof(unsigned long));

        gf2x_tfft_dft_inner_split(o, tr, a, bits_a, m1, buf, bufsize);
        tr += 2 * K * compute_np(m1, K);
        gf2x_tfft_dft_inner_split(o, tr, a, bits_a, m2, buf, bufsize);

        free(buf);
    }
}

static void gf2x_tfft_compose_inner(gf2x_tfft_info_srcptr o, gf2x_tfft_ptr tc, gf2x_tfft_srcptr ta, gf2x_tfft_srcptr tb, size_t M)
{
    size_t K = o->K;
    // tc, ta, tb may happily alias each other.
    size_t Np = compute_Np(M, K);
    size_t np = W(Np);
    unsigned long * tmp1 = o->tmp;
    // tmp2 = o->tmp + 2 * np;
    /* max(2np,gf2x_toomspace(2np)) words */
    unsigned long * tmp3 = o->tmp + 4 * np;

    for (size_t i = 0; i < K; i++) {
	MulMod(tc, ta, tb, Np, tmp1, tmp3);
        ta += 2 * np;
        tb += 2 * np;
        tc += 2 * np;
    }
}

void gf2x_tfft_compose(gf2x_tfft_info_srcptr o, gf2x_tfft_ptr tc, gf2x_tfft_srcptr ta, gf2x_tfft_srcptr tb)
{
    if (o->K == 0) {
        gf2x_mul(tc, ta, W(o->bits_a), tb, W(o->bits_b));
    } else if (!o->split){
        gf2x_tfft_compose_inner(o, tc, ta, tb, o->M);
    } else {
        size_t m1 = o->M;
        size_t m2 = o->M - 1;
        size_t K = o->K;
        gf2x_tfft_compose_inner(o, tc, ta, tb, m1);
        size_t offset = 2 * K * compute_np(m1, K);
        tc += offset;
        ta += offset;
        tb += offset;
        gf2x_tfft_compose_inner(o, tc, ta, tb, m2);
    }
}

void gf2x_tfft_add(gf2x_tfft_info_srcptr o, gf2x_tfft_ptr tc, gf2x_tfft_srcptr ta, gf2x_tfft_srcptr tb)
{
    for (size_t i = 0; i < gf2x_tfft_size(o) ; i++) {
        tc[i] = ta[i] ^ tb[i];
    }
}

void gf2x_tfft_addcompose(gf2x_tfft_info_srcptr o, gf2x_tfft_ptr tc, gf2x_tfft_srcptr ta, gf2x_tfft_srcptr tb)
{
    gf2x_tfft_t * t = gf2x_tfft_alloc(o, 1);
    gf2x_tfft_compose(o, t, ta, tb);
    gf2x_tfft_add(o, tc, tc, t);
    gf2x_tfft_free(o, t, 1);
}

void gf2x_tfft_ift_inner(gf2x_tfft_info_srcptr o, unsigned long * a, size_t bits_a, gf2x_tfft_ptr tr, size_t M)
{
    size_t K = o->K;
    size_t Mp = CEIL(M, K / 3);	// ceil(M/(K/3))
    size_t Np = Mp * (K / 3);	// Np >= M, Np multiple of K/3
    size_t np = W(Np);	       	// Words to store Np bits
    size_t i;

    unsigned long * tmp1, * tmp2, * tmp3;
    tmp1 = o->tmp;
    tmp2 = o->tmp + 2 * np;
    tmp3 = o->tmp + 4 * np; /* max(2np,gf2x_toomspace(2np)) words */

    // allocate the array of pointers. It's just temporary stuff.
    unsigned long ** A = malloc_or_die(K * sizeof(unsigned long *));
    for (i = 0; i < K; i++) A[i] = tr + 2 * i * np;
    unsigned long ** Ap = malloc_or_die(K * sizeof(unsigned long *));

    for (i = 0; i < K; i++) Ap[i] = A[o->perm[i]]; 
    fft(Ap, K, 3 * Np - Mp, Np, 1, tmp1, tmp2, tmp3, o->perm);
    for (i = 0; i < K; i++) ASSERT(A[i] == Ap[o->perm[i]]);

    free(Ap);
    recompose(a, W(bits_a), A, K, M, Np);
    free(A);
}

void gf2x_tfft_ift(gf2x_tfft_info_srcptr o, unsigned long * c, size_t bits_c, gf2x_tfft_ptr tr)
{
    if (o->K == 0) {
        Copy(c, tr, W(bits_c));
    } else if (!o->split) {
        gf2x_tfft_ift_inner(o, c, bits_c, tr, o->M);
    } else {
        size_t K = o->K;
        size_t m1 = o->M;
        size_t m2 = m1 - 1;
        size_t cn = W(2 * K * m1);
        size_t cn0 = W(o->bits_a) + W(o->bits_b);
        ASSERT(cn0 <= cn);

        size_t cn1 = W(MIN(K*m1,o->bits_a)) + W(MIN(K*m1,o->bits_b));
        unsigned long * c1 = malloc_or_die(cn * sizeof(unsigned long));
        Clear(c1, I(K * m1), cn);
        gf2x_tfft_ift_inner(o, c1, cn * WLEN, tr, m1);
        wrap(c1, cn1 * WLEN, K * m1);

        tr += 2 * K * compute_np(m1, K);

        size_t cn2 = W(MIN(K*m2,o->bits_a)) + W(MIN(K*m2,o->bits_b));
        unsigned long * c2 = malloc_or_die(cn * sizeof(unsigned long));
        Clear(c2, I(K * m2), cn);
        gf2x_tfft_ift_inner(o, c2, cn * WLEN, tr, m2);
        wrap(c2, cn2 * WLEN, K * m2);

        split_reconstruct(c, c1, c2, cn0, K, m1);
        free(c1);
        free(c2);
    }
}


/* multiplies {a, an} by {b, bn} using an FFT of length K,
   and stores the result into {c, an+bn}. If an+bn is too small
   then Toom-Cook is used.  */

// arrange so that we multiply polynomials having respectively n1 and n2
// _BITS_ ; which means degree+1 (thus can be unsigned).
//
// because this algorithm needs to know about the K value, which in turns
// depends on proper tuning, we ask for it to be provided by the caller.
// Negative values of K mean to use FFT2.
void gf2x_tfft_init(gf2x_tfft_info_ptr o, size_t bits_a, size_t bits_b, ...)
{
    o->bits_a = bits_a;
    o->bits_b = bits_b;

    size_t nwa = W(bits_a);
    size_t nwb = W(bits_b);

    va_list ap;
    va_start(ap, bits_b);
    long K = va_arg(ap, long);

    /* Since we have a dangerous interface with a variable argument which
     * is a possibly _signed_ long, better try to catch the expected
     * mistakes */
    for(int i = 2*(K>0)-1 ; K/i > 1 ; i*=3) {
        if ((K/i)%3) {
            fprintf(stderr, "extra argument to gf2x_tfft_init (of type long) must be a power of 3 (got %ld)\n", K);
            abort();
        }
    }

    size_t M;
    if (K > 0) {
        M = CEIL((nwa + nwb) * WLEN, K);	// ceil(bits(product)/K)
        o->K = K;
        o->M = M;
        o->split = 0;
    } else {
        ASSERT(-K >= WLEN);
        size_t cn2 = CEIL(nwa + nwb, 2);	// Space for half product
        size_t m2 = CEIL(cn2 * WLEN, -K);	// m2 = ceil(cn2*WLEN/K)
        size_t m1 = m2 + 1;		// next possible M
        M = m1;
        o->K = -K;
        o->M = M;
        o->split = 1;
    }

    if (nwa + nwb < MUL_FFT_THRESHOLD) {
        // make this special.
        o->K = 0;
        o->M = 0;
        o->tmp = NULL;
        o->perm = NULL;
        return;
    }

    /* We also have to allocate the temporary space used by this FFT */
    size_t np = compute_np(M, o->K);

    size_t i = gf2x_toomspace(2 * np);
    if (i < 2 * np)
	i = 2 * np;
    size_t ltmp = 4 * np + i;
    o->tmp = (unsigned long *) malloc_or_die(ltmp * sizeof(unsigned long));
    o->perm = (size_t *) malloc_or_die(o->K * sizeof(size_t));
    bitrev(0, 0, o->K, 1, o->perm);
    va_end(ap);
}

void gf2x_tfft_init_similar(gf2x_tfft_info_ptr o, size_t bits_a, size_t bits_b, gf2x_tfft_info_srcptr other)
{
    gf2x_tfft_init(o, bits_a, bits_b, other->K);
}

int gf2x_tfft_compatible(gf2x_tfft_info_srcptr o1, gf2x_tfft_info_srcptr o2)
{
    return o1->K == o2->K && o1->M == o2->M;
}

void gf2x_tfft_clear(gf2x_tfft_info_ptr o)
{
    if (o->K) {
        free(o->tmp);
        free(o->perm);
    }
    memset(o, 0, sizeof(gf2x_tfft_info_t));
}

/** gf2x_mul_fft merely wraps around the calls above **/

/* multiplies {a, an} by {b, bn} using an FFT of length K,
   and stores the result into {c, an+bn}.
   The result is computed mod (x^N+1) where N = K*M.
   Thus for a full product K*M should be >= (an+bn)*WLEN,
   the size of product in bits. For the result mod (x^N+1)
   it is only required that 2*K*M >= (an+bn)*WLEN */

// here an and bn denote numbers of WORDS, while the gf2x_tfft_* routines
// are interested in number of BITS.
void gf2x_mul_fft(unsigned long *c, const unsigned long *a, size_t an,
	    const unsigned long *b, size_t bn, long K)
{
    gf2x_tfft_info_t o;
    gf2x_tfft_init(o, an * WLEN, bn * WLEN, K);
    
    if (o->K == 0) {
	printf("gf2x_mul_fft: arguments (%zu, %zu) too small\n", an, bn);
        /* Note that actually the routines below do work, because they're
         * specified for working. However, this contradicts the fact that
         * via this entry point, we have explicitly asked for _not_
         * falling back to standard gf2x routines. So it's a caller bug
         */
        abort();
    }
    gf2x_tfft_ptr ta = gf2x_tfft_alloc(o, 1);
    gf2x_tfft_ptr tb = gf2x_tfft_alloc(o, 1);
    gf2x_tfft_ptr tc = gf2x_tfft_alloc(o, 1);

    gf2x_tfft_dft(o, ta, a, an * WLEN);
    gf2x_tfft_dft(o, tb, b, bn * WLEN);

    gf2x_tfft_compose(o, tc, ta, tb);

    gf2x_tfft_ift(o, c, (an+bn)*WLEN, tc);

    gf2x_tfft_free(o, ta, 1);
    gf2x_tfft_free(o, tb, 1);
    gf2x_tfft_free(o, tc, 1);

    gf2x_tfft_clear(o);
}
