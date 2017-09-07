#ifndef MPFQ_H_
#define MPFQ_H_

/* This header contains common declarations used by mpfq modules */

/* we always include stdio.h, otherwise our inclusion of gmp.h might
 * prevent gmp's I/O functions to ever be exposed... */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef GF2X_WORDSIZE
#error "This file has been modified for use within gf2x. Please arrange so that gf2x-small.h is included before this file"
#endif

#ifdef __cplusplus
extern "C" {
#endif



/*** Constants for field_specify ***/

#define MPFQ_DONE 0             /* At the end of the variadic option functions */
// #define MPFQ_PRIME_MPN 1        /* mp_limb_t *, size depending on implementation. Prefer MPFQ_PRIME_MPZ */
#define MPFQ_POLYNOMIAL 2       /* this expects an mpfq polynomial */
#define MPFQ_DEGREE 3           /* int */
#define MPFQ_IO_TYPE 4          /* for setopt */
#define MPFQ_GROUPSIZE 5        /* int (SIMD group size) */
// #define MPFQ_PRIME_MPZ 6        /* mpz_t */
#define MPFQ_MANDATORY_TAG 7    /* force the tag to be this one ; this is
                                 * of course pointless for the low-level
                                 * implementation, but
                                 * mpfq_vbase_oo_field_init_byfeatures
                                 * uses it. */

#define LEXGE2(X,Y,A,B) (X>A || (X == A && Y >= B))
#define LEXGE3(X,Y,Z,A,B,C) (X>A || (X == A && LEXGE2(Y,Z,B,C)))
#define LEXLE2(X,Y,A,B) LEXGE2(A,B,X,Y)
#define LEXLE3(X,Y,Z,A,B,C) LEXGE3(A,B,C,X,Y,Z)

#ifndef GNUC_VERSION
#ifndef __GNUC__
#define GNUC_VERSION(X,Y,Z) 0
#else
#define GNUC_VERSION(X,Y,Z)     \
(__GNUC__ == X && __GNUC_MINOR__ == Y && __GNUC_PATCHLEVEL__ == Z)
#endif
#endif

#ifndef GNUC_VERSION_ATLEAST
#ifndef __GNUC__
#define GNUC_VERSION_ATLEAST(X,Y,Z) 0
#else
#define GNUC_VERSION_ATLEAST(X,Y,Z)     \
LEXGE3(__GNUC__,__GNUC_MINOR__,__GNUC_PATCHLEVEL__,X,Y,Z)
#endif
#endif

#ifndef GNUC_VERSION_ATMOST
#ifndef __GNUC__
#define GNUC_VERSION_ATMOST(X,Y,Z) 0
#else
#define GNUC_VERSION_ATMOST(X,Y,Z)     \
LEXLE3(__GNUC__,__GNUC_MINOR__,__GNUC_PATCHLEVEL__,X,Y,Z)
#endif
#endif

/* typedef unsigned long ulong; */

#ifndef MAX
#define MAX(h,i) ((h) > (i) ? (h) : (i))
#endif

#ifndef MIN
#define MIN(l,o) ((l) < (o) ? (l) : (o))
#endif

#ifndef __cplusplus
#ifndef max
#define max(a,b) ((a) > (b) ? (a) : (b))
#endif

#ifndef min
#define min(a,b) ((a) < (b) ? (a) : (b))
#endif
#endif	/* __cplusplus */


#ifndef ABS
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#endif

#if GNUC_VERSION_ATLEAST(3,4,0)
/* according to
 * http://gcc.gnu.org/onlinedocs/gcc-3.1.1/gcc/Variable-Attributes.html#Variable%20Attributes
 * the 'unused' attribute already existed in 3.1.1 ; however the rules
 * for its usage remained quirky until 3.4.0, so we prefer to stick to
 * the more modern way of using the unused attribute, and recommend
 * setting the -Wno-unused flag for pre-3.4 versions of gcc
 */
#ifndef	MAYBE_UNUSED
#define	MAYBE_UNUSED	__attribute__((unused))
#endif
#endif

#if GNUC_VERSION_ATLEAST(3,1,0) /* apparently */
#ifndef MPFQ_EXPECT
#define MPFQ_EXPECT(x,val)   __builtin_expect(x,val)
#endif
#endif

#if GNUC_VERSION_ATLEAST(3,4,0)
#define mpfq_clzl(x)         __builtin_clzl(x)
#define mpfq_ctzl(x)         __builtin_ctzl(x)
#define mpfq_parityl(x)      __builtin_parityl(x)
#else
/* provide slow fallbacks */
static inline int mpfq_clzl(unsigned long x)
{
        static const int t[4] = { 2, 1, 0, 0 };
        int a = 0;
        int res;
#if (GF2X_WORDSIZE == 64)
        if (x >> 32) { a += 32; x >>= 32; }
#endif  
        if (x >> 16) { a += 16; x >>= 16; }
        if (x >>  8) { a +=  8; x >>=  8; }
        if (x >>  4) { a +=  4; x >>=  4; }
        if (x >>  2) { a +=  2; x >>=  2; }
        res = GF2X_WORDSIZE - 2 - a + t[x];
        return res;
}
static inline int mpfq_ctzl(unsigned long x)
{
	return GF2X_WORDSIZE - mpfq_clzl(x & - x);
}
static inline int mpfq_parityl(unsigned long x)
{
	static const int t[4] = { 0, 1, 1, 0, };
#if (GF2X_WORDSIZE == 64)
	x ^= (x >> 32);
#endif
	x ^= (x >> 16);
	x ^= (x >>  8);
	x ^= (x >>  4);
	x ^= (x >>  2);
	return t[x & 3UL];
}
#endif

#ifndef	MAYBE_UNUSED
#define	MAYBE_UNUSED	/**/
#endif

#ifndef MPFQ_EXPECT
#define MPFQ_EXPECT(x,val)   (x)
#endif

#ifndef	MPFQ_UNLIKELY
#define	MPFQ_UNLIKELY(x)	MPFQ_EXPECT(x, 0)
#endif
#ifndef	MPFQ_LIKELY
#define	MPFQ_LIKELY(x)	MPFQ_EXPECT(x, 1)
#endif



static inline int mpfq_clzlx(unsigned long * x, int n)
{
	int r = 0;
	for( ; n > 0 && MPFQ_UNLIKELY(!x[n-1]) ; --n) r+=GF2X_WORDSIZE;
	if (n == 0) return r;
	r += mpfq_clzl(x[n-1]);
	return r;
}

static inline int mpfq_ctzlx(unsigned long * x, int n)
{
	int r = 0;
	for( ; n > 0 && MPFQ_UNLIKELY(!*x) ; --n,++x) r+=GF2X_WORDSIZE;
	if (n == 0) return r;
	r += mpfq_ctzl(*x);
	return r;
}

/***  Some useful macros ***/

/* use these only for the large malloc()s, please */
static inline void * mpfq_malloc_check(size_t s) {
    void * r = malloc(s);
#ifdef  MPFQ_TRACK_MALLOC
    if (s>>28) { fprintf(stderr, "MALLOC(%.1f)\n", s/1048576.); }
#endif
    if (!r) {
        fprintf(stderr, "malloc(%zu) failed\n", s);
        abort();
    }
    return r;
}

static inline void * mpfq_realloc_check(void * p, size_t os, size_t s) {
    void * r = realloc(p, s);
#ifdef  MPFQ_TRACK_MALLOC
    if (s>>28) { fprintf(stderr, "REALLOC(%.1f, %.1f)\n", os/1048576., s/1048576.); }
#endif
    if (s && !r) {
        fprintf(stderr, "realloc(%zu, %zu) failed\n", os, s);
        abort();
    }
    return r;
}

static inline void mpfq_free(void * p, size_t s MAYBE_UNUSED)
{
#ifdef  MPFQ_TRACK_MALLOC
    if (s>>28) { fprintf(stderr, "FREE(%.1f)\n", s/1048576.); }
#endif
    free(p);
}

static inline void malloc_failed() {
    fprintf(stderr, "malloc() failed\n");
    abort();
}


#if 0

/* Given the fact that copies are always very small, we're probably
 * better off giving the compiler the opportunity to optimize all this
 * away.
 */

/* dst and src known to not overlap, except possibly if dst == src */
static inline void mpfq_copy(mp_limb_t * dst, const mp_limb_t * src, mp_size_t n) {
    // if (dst != src) mpn_copyi(dst, src, n);
    for( ; n-- ; ) *dst++ = *src++;
}

/* dst and src possibly overlap, copy increasingly so that src >= dst is ok */
static inline void mpfq_copyi(mp_limb_t * dst, const mp_limb_t * src, mp_size_t n) {
    // mpn_copyi(dst, src, n);
    for( ; n-- ; ) *dst++ = *src++;
}

/* dst and src possibly overlap, copy decreasingly so that src <= dst is ok */
static inline void mpfq_copyd(mp_limb_t * dst, const mp_limb_t * src, mp_size_t n) {
    // mpn_copyd(dst, src, n);
    for(dst += n, src += n ; n-- ; ) *--dst = *--src;
}

static inline void mpfq_zero(mp_limb_t * dst, mp_size_t n) {
    // mpn_zero(dst, 0, n);
    for( ; n-- ; ) *dst++ = 0;
}

#endif


#ifdef __cplusplus
}
#endif

#endif	/* MPFQ_H_ */
