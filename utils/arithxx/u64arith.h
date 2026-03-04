/* Some commonly used assembly helper macros for arithmetic on 
   unsigned 64-bit integers. 
   Defining U64ARITH_VERBOSE_ASM puts in the asm output a line for each
   function call that shows what registers/memory locations the operands 
   are in.
   Defining U64ARITH_NO_ASM avoids asm macros and uses the C fallback code
   where available. This can be used for testing the fallback code.
*/

#ifndef CADO_U64ARITH_H
#define CADO_U64ARITH_H

#include <stdint.h>
#include <stdio.h>
#include "macros.h"

#if defined(HAVE_INT128)
typedef union {uint64_t x[2]; unsigned __int128 y;} _u64arith_union2_t;
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** Test whether a1:a2 > b1:b2
 * 
 * Let a = a1 + 2^64 * a2, b = b1 + 2^64 * b2. Return 1 if a > b,
 * and 0 otherwise. 
 */
static inline int
u64arith_gt_2_2(uint64_t a1, uint64_t a2, uint64_t b1, uint64_t b2) ATTRIBUTE((const));
static inline int
u64arith_gt_2_2(const uint64_t a1, const uint64_t a2,
                const uint64_t b1, const uint64_t b2)
{
#if defined(GENUINE_GCC) && __GNUC__ >= 8 && defined(HAVE_INT128)
    _u64arith_union2_t a = {{a1, a2}}, b = {{b1, b2}};
    return a > b;
#else
    return a2 > b2 || (a2 == b2 && a1 > b1);
#endif
}

/** Test whether a1:a2 >= b1:b2
 * 
 * Let a = a1 + 2^64 * a2, b = b1 + 2^64 * b2. Return 1 if a >= b,
 * and 0 otherwise. 
 */
static inline int
u64arith_ge_2_2(uint64_t a1, uint64_t a2, uint64_t b1, uint64_t b2) ATTRIBUTE((const));
static inline int
u64arith_ge_2_2(const uint64_t a1, const uint64_t a2,
                const uint64_t b1, const uint64_t b2)
{
#if defined(GENUINE_GCC) && __GNUC__ >= 8 && defined(HAVE_INT128)
    _u64arith_union2_t a = {{a1, a2}, b = {b1, b2}};
    return a >= b;
#else
  return a2 > b2 || (a2 == b2 && a1 >= b1);
#endif
    
}

/** Test whether a1:a2 < b1:b2
 * 
 * Let a = a1 + 2^64 * a2, b = b1 + 2^64 * b2. Return 1 if a < b,
 * and 0 otherwise. 
 */
static inline int
u64arith_lt_2_2(uint64_t a1, uint64_t a2, uint64_t b1, uint64_t b2) ATTRIBUTE((const));
static inline int
u64arith_lt_2_2(const uint64_t a1, const uint64_t a2,
                const uint64_t b1, const uint64_t b2)
{
#if defined(GENUINE_GCC) && __GNUC__ >= 8 && defined(HAVE_INT128)
    _u64arith_union2_t a = {{a1, a2}}, b = {{b1, b2}};
    return a < b;
#else
  return a2 < b2 || (a2 == b2 && a1 < b1);
#endif
    
}

/** Test whether a1:a2 <= b1:b2
 * 
 * Let a = a1 + 2^64 * a2, b = b1 + 2^64 * b2. Return 1 if a <= b,
 * and 0 otherwise. 
 */
static inline int
u64arith_le_2_2(uint64_t a1, uint64_t a2, uint64_t b1, uint64_t b2) ATTRIBUTE((const));
static inline int
u64arith_le_2_2(const uint64_t a1, const uint64_t a2,
                const uint64_t b1, const uint64_t b2)
{
#if defined(GENUINE_GCC) && __GNUC__ >= 8 && defined(HAVE_INT128)
    _u64arith_union2_t a = {{a1, a2}}, b = {{b1, b2}};
    return a <= b;
#else
  return a2 < b2 || (a2 == b2 && a1 <= b1);
#endif
    
}

/** Add r1:r2 += a1:a2 with carry propagation
 * 
 * Let r = r1 + 2^64 * r2, a = a1 + 2^64 * a2. Set r1:r2 := (r+a) mod 2^128.
 */
static inline void
u64arith_add_2_2 (uint64_t *r1, uint64_t *r2, 
		  const uint64_t a1, const uint64_t a2)
{
#if !defined (U64ARITH_NO_ASM) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM) \
    && ( defined(GENUINE_GNUC) && __GNUC__ <= 7 || defined(__INTEL_COMPILER) )
#ifdef U64ARITH_VERBOSE_ASM
  __asm__ ("# u64arith_add_2_2 (%0, %1, %2, %3)\n" : : 
           "X" (*r1), "X" (*r2), "X" (a1), "X" (a2));
#endif
  __asm__ __VOLATILE (
    "addq %2, %0\n\t"
    "adcq %3, %1\n"
    : "+&r" (*r1), "+r" (*r2)
    : "rme" (a1), "rme" (a2)
    : "cc");
#else
  /* On x86_64, Clang 6.0, 7.0, 8.0, 9.0, gcc 8.0 and 9.0 produce good code
   * from this. Intel icc however produces terrible code. */
  *r1 += a1;
  *r2 += a2 + (*r1 < a1);
#endif
}

/** Add r1:r2 += a with carry propagation
 * 
 * Let r = r1 + 2^64 * r2. Set r1:r2 := (r+a) mod 2^128.
 */
static inline void
u64arith_add_1_2 (uint64_t *r1, uint64_t *r2, const uint64_t a)
{
  /* We have assembly only for x86_64 and on that architecture, this function
     would generate the same code as u64arith_add_2_2() with an immediate
     zero value for the high word, so we fall back to that. */
  u64arith_add_2_2(r1, r2, a, 0);
}


/* Adds two uint64_t from two uint64_t with carry propagation 
   from low word (r1) to high word (r2). Returns 1 if there was a carry out 
   from high word, otherwise returns 0. */

static inline char
u64arith_add_2_2_cy (uint64_t *r1, uint64_t *r2, 
		     const uint64_t a1, const uint64_t a2)
{
  char cy;
#ifdef U64ARITH_VERBOSE_ASM
  __asm__ ("# u64arith_add_2_2_cy (%0, %1, %2, %3)\n" : : 
           "X" (*r1), "X" (*r2), "X" (a1), "X" (a2));
#endif
#if !defined (U64ARITH_NO_ASM) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
  __asm__ __VOLATILE (
    "addq %3, %0\n\t"
    "adcq %4, %1\n\t"
    "setc %2\n"
    : "+&r" (*r1), "+r" (*r2), "=r" (cy)
    : "rme" (a1), "rme" (a2)
    : "cc");
#else
  uint64_t u1 = *r1 + a1,
           u2 = *r2 + a2 + (u1 < a1);
  /* Overflow occurred iff the sum is smaller than one of the summands */
  cy = u64arith_gt_2_2(a1, a2, u1, u2);
  *r1 = u1;
  *r2 = u2;
#endif
  return cy;
}


/** Compute r = (a + b) % m
 * 
 * Requires a < m and b <= m, then r == a+b (mod m) and r < m.
 */
static inline void
u64arith_addmod_1_1 (uint64_t *r, const uint64_t a,
                     const uint64_t b, const uint64_t m)
{
  ASSERT_EXPENSIVE (a < m && b <= m);

#if !defined (U64ARITH_NO_ASM) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
  {
    uint64_t t = a + b, tr = a - m;

    __asm__ __VOLATILE (
      "addq %2, %0\n\t"   /* tr += b */
      "cmovnc %1, %0\n\t"  /* if (!cy) tr = t */
      : "+&r" (tr)
      : "rm" (t), "rme" (b)
      : "cc"
    );
    ASSERT_EXPENSIVE (tr == ((a >= m - b) ? (a - (m - b)) : (a + b)));
    r[0] = tr;
  }
#else
  r[0] = (b >= m - a) ? (b - (m - a)) : (a + b);
#endif

  ASSERT_EXPENSIVE (r[0] < m);
}


/** Add r1:r2 -= a1:a2 with borrow propagation
 * 
 * Let r = r1 + 2^64 * r2, a = a1 + 2^64 * a2. Set r1:r2 := (r-a) mod 2^128.
 */
static inline void
u64arith_sub_2_2 (uint64_t *r1, uint64_t *r2, 
		  const uint64_t a1, const uint64_t a2)
{
#ifdef U64ARITH_VERBOSE_ASM
  __asm__ ("# u64arith_sub_2_2 (%0, %1, %2, %3)\n" : : 
           "X" (*r1), "X" (*r2), "X" (a1), "X" (a2));
#endif
#if !defined (U64ARITH_NO_ASM) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
  __asm__ __VOLATILE (
    "subq %2, %0\n\t"
    "sbbq %3, %1\n"
    : "+&r" (*r1), "+r" (*r2)
    : "rme" (a1), "rme" (a2)
    : "cc");
#else
  uint64_t u = *r1;
  *r1 -= a1;
  *r2 -= a2 + (*r1 > u);
#endif
}


/** Add r1:r2 -= a with borrow propagation
 * 
 * Let r = r1 + 2^64 * r2. Set r1:r2 := (r-a) mod 2^128.
 */
static inline void
u64arith_sub_1_2 (uint64_t *r1, uint64_t *r2, 
                  const uint64_t a)
{
    u64arith_sub_2_2(r1, r2, a, 0);
}


/* Subtract two uint64_t from two uint64_t with borrow propagation 
   from low word (r1) to high word (r2). Returns 1 if there was a borrow out 
   from high word, otherwise returns 0. */

static inline char
u64arith_sub_2_2_cy (uint64_t *r1, uint64_t *r2, 
                     const uint64_t a1, const uint64_t a2)
{
  char cy;
#ifdef U64ARITH_VERBOSE_ASM
  __asm__ ("# u64arith_sub_2_2_cy (%0, %1, %2, %3)\n" : : 
           "X" (*r1), "X" (*r2), "X" (a1), "X" (a2));
#endif
#if !defined (U64ARITH_NO_ASM) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
  __asm__ __VOLATILE (
    "subq %3, %0\n\t"
    "sbbq %4, %1\n\t"
    "setc %2\n"
    : "+&r" (*r1), "+r" (*r2), "=r" (cy)
    : "rme" (a1), "rme" (a2)
    : "cc");
#else
  uint64_t u1 = *r1 - a1, u2 = *r2 - a2;
  if (a1 > *r1)
    u2--;
  cy = u64arith_gt_2_2(a1, a2, *r1, *r2);
  *r1 = u1;
  *r2 = u2;
#endif
  return cy;
}


/* Subtract only if result is non-negative */
static inline void
u64arith_sub_1_1_ge (uint64_t *r, const uint64_t a)
{
  uint64_t t = *r;
#ifdef U64ARITH_VERBOSE_ASM
  __asm__ ("# u64arith_sub_1_1_ge (%0, %1, %2, %3)\n" : : 
           "X" (*r1), "X" (*r2), "X" (a1), "X" (a2));
#endif
#if !defined (U64ARITH_NO_ASM) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
  __asm__ __VOLATILE (
    "subq %2, %0\n\t" /* r -= a */
    "cmovc %1, %0\n\t" /* If there's a borrow, restore r from t */
    : "+&r" (*r)
    : "r" (t), "rme" (a)
    : "cc");
#else
  t -= a;
  if (*r >= a)
    *r = t;
#endif
}

/* Subtract only if result is non-negative */
static inline void
u64arith_sub_2_2_ge (uint64_t *r1, uint64_t *r2, 
                     const uint64_t a1, const uint64_t a2)
{
  uint64_t t1 = *r1, t2 = *r2;
#ifdef U64ARITH_VERBOSE_ASM
  __asm__ ("# u64arith_sub_2_2_ge (%0, %1, %2, %3)\n" : : 
           "X" (*r1), "X" (*r2), "X" (a1), "X" (a2));
#endif
#if !defined (U64ARITH_NO_ASM) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
  __asm__ __VOLATILE (
    "subq %4, %0\n\t" /* r1 -= a1 */
    "sbbq %5, %1\n\t" /* r2 -= a2 + cy */
    "cmovc %2, %0\n\t" /* If there's a borrow, restore r1 from t1 */
    "cmovc %3, %1\n\t" /* and r2 from t2 */
    : "+&r" (*r1), "+&r" (*r2)
    : "r" (t1), "r" (t2), "rme" (a1), "rme" (a2)
    : "cc");
#else
  if (!u64arith_gt_2_2(a1, a2, *r1, *r2))
    {
      *r1 = t1 - a1;
      *r2 = t2 - a2 - (a1 > t1);
    }
#endif
}


/** Compute r = (a - b) % m
 * 
 * Requires a < m and b < m, then r == a+b (mod m) and r < m.
 */
static inline void
u64arith_submod_1_1 (uint64_t *r, const uint64_t a,
                     const uint64_t b, const uint64_t m)
{
    ASSERT_EXPENSIVE (a < m && b < m);
#if !defined (U64ARITH_NO_ASM) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
    /* Could do tr = a+m-b; t = a-b; if (!carry) tr = t; */
    uint64_t tr, t = a;
    __asm__ __VOLATILE (
      "sub %2, %1\n\t"  /* t -= b ( = a - b) */
      "lea (%1,%3,1), %0\n\t" /* tr = t + m ( = a - b + m) */
      "cmovnc %1, %0\n\t" /* if (a >= b) tr = t */
      : "=&r" (tr), "+&r" (t)
      : "rme" (b), "r" (m)
      : "cc"
    );
    r[0] = tr;
#elif 1
  /* Seems to be faster than the one below */
    uint64_t t = 0, tr;
    if ((tr = a - b) > a)
      t = m;
    r[0] = tr + t;
#else
    r[0] = (a < b) ? (a - b + m) : (a - b);
#endif

    ASSERT_EXPENSIVE (r[0] < m);
}


/* Multiply two uint64_t "a" and "b" and put the result as 
   r2:r1 (r2 being the high word) */

static inline void
u64arith_mul_1_1_2 (uint64_t *r1, uint64_t *r2, 
		    const uint64_t a, const uint64_t b)
{
#ifdef U64ARITH_VERBOSE_ASM
  __asm__ ("# u64arith_mul_1_1_2 (%0, %1, %2, %3)\n" : : 
           "X" (*r1), "X" (*r2), "X" (a), "X" (b));
#endif
#if !defined (U64ARITH_NO_ASM) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
  __asm__ __VOLATILE (
    "mulq %3"
    : "=a" (*r1), "=d" (*r2)
    : "%0" (a), "rm" (b)
    : "cc");
#elif defined(HAVE_INT128)
    unsigned __int128 r = (unsigned __int128) a * b;
    *r1 = r;
    *r2 = r >> 64;
#else
  const uint64_t mask = ((uint64_t)1 << 32) - 1;
  const uint64_t ah = a >> 32, al = a & mask, bh = b >> 32, bl = b & mask;
  uint64_t t1, t2, p1, p2;

  t1 = al * bl;
  t2 = 0;
  p1 = ah * bl;
  p2 = p1 >> 32;
  p1 = (p1 & mask) << 32;
  u64arith_add_2_2 (&t1, &t2, p1, p2);
  p1 = al * bh;
  p2 = p1 >> 32;
  p1 = (p1 & mask) << 32;
  u64arith_add_2_2 (&t1, &t2, p1, p2);
  t2 += ah * bh;
  *r1 = t1; 
  *r2 = t2;
#endif
}


static inline void
u64arith_sqr_1_2 (uint64_t *r1, uint64_t *r2, 
		  const uint64_t a)
{
#ifdef U64ARITH_VERBOSE_ASM
  __asm__ ("# u64arith_sqr_1_2 (%0, %1, %2)\n" : : 
           "X" (*r1), "X" (*r2), "X" (a));
#endif
#if !defined (U64ARITH_NO_ASM) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
  __asm__ __VOLATILE (
    "mulq %%rax"
    : "=a" (*r1), "=d" (*r2)
    : "0" (a)
    : "cc");
#elif defined(HAVE_INT128)
    unsigned __int128 r = (unsigned __int128) a * a;
    *r1 = r;
    *r2 = r >> 64;
#else
  const uint64_t mask = ((uint64_t)1 << 32) - 1;
  const uint64_t ah = a >> 32, al = a & mask;
  uint64_t t1, t2, p1, p2;

  t1 = al * al;
  t2 = 0;
  p1 = ah * al;
  p2 = p1 >> 32;
  p1 = (p1 & mask) << 32;
  u64arith_add_2_2 (&t1, &t2, p1, p2);
  u64arith_add_2_2 (&t1, &t2, p1, p2);
  t2 += ah * ah;
  *r1 = t1; 
  *r2 = t2;
#endif
}


/* Set *r to lo shifted right by i bits, filling in the low bits from hi into the high
   bits of *r. I.e., *r = (hi * 2^64 + lo) / 2^i. Assumes 0 <= i < 64. */
static inline void
u64arith_shrd (uint64_t *r, const uint64_t hi, const uint64_t lo,
              const unsigned char i)
{
  ASSERT_EXPENSIVE (i < 64);
#ifdef U64ARITH_VERBOSE_ASM
/* Disable the "uninitialized" warning here, as *r is only written to and
   does not need to be initialized, but we need to write (*r) here so the
   "X" constraint can be resolved even when r does not have an address, e.g.,
   when it is passed around in a register. It seems that "X" is assumed by
   gcc as possibly referring to an input, and since "X" matches anything,
   that's probably a neccessary assumtion to make. */
#if GNUC_VERSION_ATLEAST(4,4,0)
#if GNUC_VERSION_ATLEAST(4,6,0)
#pragma GCC diagnostic push
#endif
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
  __asm__ ("# u64arith_shrd (*r=%0, hi=%1, lo=%2, i=%3)\n" : : 
           "X" (*r), "X" (hi), "X" (lo), "X" (i));
#if GNUC_VERSION_ATLEAST(4,6,0)
#pragma GCC diagnostic pop
#endif
#endif

#if !defined (U64ARITH_NO_ASM) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
  __asm__ __VOLATILE (
    "shrdq %3, %1, %0\n"
    : "=rm" (*r)
    : "r" (hi), "0" (lo), "cJ" (i) /* i can be in %cl or a literal constant < 64 */
    : "cc");
#else
  if (i > 0) /* shl by 64 is no-op on x86! */
    *r = (lo >> i) | (hi << (64 - i));
  else
    *r = lo;
#endif
}


/* Shift the 128-bit integer in lo:hi right by i bits, 0 <= i < 64 */

static inline void
u64arith_shr_2 (uint64_t *lo, uint64_t *hi, const unsigned char i)
{
  ASSERT_EXPENSIVE (i < 64);
  u64arith_shrd(lo, *hi, *lo, i);
  *hi >>= i;
}

/* Set *r to hi shifted left by i bits, filling in the high bits from lo into the low
   bits of *r. I.e., *r = (hi + lo*2^-64) * 2^i. Assumes 0 <= i < 64. */
static inline void
u64arith_shld (uint64_t *r, const uint64_t lo, const uint64_t hi,
              const unsigned char i)
{
  ASSERT_EXPENSIVE (i < 64);
#ifdef U64ARITH_VERBOSE_ASM
#if GNUC_VERSION_ATLEAST(4,4,0)
#if GNUC_VERSION_ATLEAST(4,6,0)
#pragma GCC diagnostic push
#endif
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
  __asm__ ("# u64arith_shld (*r=%0, lo=%1, hi=%2, i=%3)\n" : : 
           "X" (*r), "X" (lo), "X" (hi), "X" (i));
#if GNUC_VERSION_ATLEAST(4,6,0)
#pragma GCC diagnostic pop
#endif
#endif

#if !defined (U64ARITH_NO_ASM) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
  __asm__ __VOLATILE (
    "shldq %3, %1, %0\n"
    : "=rm" (*r)
    : "r" (lo), "0" (hi), "cJ" (i) /* i can be in %cl or a literal constant < 64 */
    : "cc");
#else
  if (i > 0) /* shr by 64 is no-op on x86! */
    *r = (hi << i) | (lo >> (64 - i));
  else
    *r = hi;
#endif
}


/* Shift the 128-bit integer in lo:hi left by i bitss, 0 <= i < 64 */

static inline void
u64arith_shl_2 (uint64_t *lo, uint64_t *hi, const unsigned char i)
{
  ASSERT_EXPENSIVE (i < 64);
  u64arith_shld(hi, *lo, *hi, i);
  *lo <<= i;
}

/* Returns number of trailing zeros in a. a must not be zero */
static inline constexpr int
u64arith_ctz (const uint64_t a)
{
#if !defined (U64ARITH_NO_ASM) && GNUC_VERSION_ATLEAST(3, 4, 0)
  ASSERT_EXPENSIVE (a != 0);
  return __builtin_ctzll(a);
#else
  static const unsigned char trailing_zeros[256] =
    {8,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     7,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
     5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0};
  char lsh, t = 0;
  uint64_t y = a;
  ASSERT_EXPENSIVE (a != 0);
  do {
    lsh = trailing_zeros [(unsigned char) y];
    y >>= lsh;
    t += lsh;
  } while (lsh == 8);
  return (int) t;
#endif
}

/* Returns number of leading zeros in a. a must not be zero */
static inline constexpr int
u64arith_clz (const uint64_t a)
{
#if !defined (U64ARITH_NO_ASM) && GNUC_VERSION_ATLEAST(3, 4, 0)
  ASSERT_EXPENSIVE (a != 0);
  STATIC_ASSERT(sizeof(uint64_t) == sizeof(unsigned long) || sizeof(uint64_t) == sizeof(unsigned long long), please_implement_me);
  if (sizeof(uint64_t) == sizeof(unsigned long))
    return __builtin_clzl(a);
  else
    return __builtin_clzll(a);
#else
  uint64_t t = (uint64_t)1 << (64 - 1);
  int i;
  ASSERT_EXPENSIVE (a != 0);
  for (i = 0; (a & t) == 0; i++)
    t >>= 1;
  return i;
#endif
}


/* Let A = a1 + a2*2^64. Returns q = floor(A / b), r = A mod b. Requires
   a2 < b. Slow, simple binary grammar-school division. Used as reference
   for testing faster code on machines without hardware division. */

static inline void
u64arith_divqr_2_1_1_slow (uint64_t *q, uint64_t *r,
		      const uint64_t a1, const uint64_t a2,
		      const uint64_t b)
{
  uint64_t R1 = a1, R2 = a2, D1 = 0, D2 = b, M = (uint64_t)1 << 63, Q = 0;

  ASSERT(a2 < b); /* Or there will be quotient overflow */

  for (int i = 0; i < 64; i++) {
    u64arith_shr_2(&D2, &D1, 1);
    /* if R >= D */
    if (!u64arith_gt_2_2(D1, D2, R1, R2)) {
      u64arith_sub_2_2(&R1, &R2, D1, D2); /* R := R - D */
      Q += M;
    }
    M >>= 1;
  }

#ifdef WANT_ASSERT_EXPENSIVE
  ASSERT_EXPENSIVE(R2 == 0 && R1 < b);
  uint64_t P1, P2;
  u64arith_mul_1_1_2(&P1, &P2, Q, b);
  u64arith_add_2_2(&P1, &P2, R1, R2);
  ASSERT_EXPENSIVE(P1 == a1 && P2 == a2);
#endif
  *q = Q;
  *r = R1;
}


/* Computes (2^128-1) / d - 1 for 2^63 <= d < 2^64. */

static inline uint64_t
u64arith_reciprocal_for_div(const uint64_t d)
{
  /* Assumes 2^63 <= d <= 2^64-1 */
  const uint64_t one = 1;
  ASSERT(d >= (one << 63));
  /* lut[i] = (2^19 - 3 * 256) / (i + 256) */
  static const uint16_t lut[256] = {2045, 2037, 2029, 2021, 2013, 2005, 1998,
    1990, 1983, 1975, 1968, 1960, 1953, 1946, 1938, 1931, 1924, 1917, 1910,
    1903, 1896, 1889, 1883, 1876, 1869, 1863, 1856, 1849, 1843, 1836, 1830,
    1824, 1817, 1811, 1805, 1799, 1792, 1786, 1780, 1774, 1768, 1762, 1756,
    1750, 1745, 1739, 1733, 1727, 1722, 1716, 1710, 1705, 1699, 1694, 1688,
    1683, 1677, 1672, 1667, 1661, 1656, 1651, 1646, 1641, 1636, 1630, 1625,
    1620, 1615, 1610, 1605, 1600, 1596, 1591, 1586, 1581, 1576, 1572, 1567,
    1562, 1558, 1553, 1548, 1544, 1539, 1535, 1530, 1526, 1521, 1517, 1513,
    1508, 1504, 1500, 1495, 1491, 1487, 1483, 1478, 1474, 1470, 1466, 1462,
    1458, 1454, 1450, 1446, 1442, 1438, 1434, 1430, 1426, 1422, 1418, 1414,
    1411, 1407, 1403, 1399, 1396, 1392, 1388, 1384, 1381, 1377, 1374, 1370,
    1366, 1363, 1359, 1356, 1352, 1349, 1345, 1342, 1338, 1335, 1332, 1328,
    1325, 1322, 1318, 1315, 1312, 1308, 1305, 1302, 1299, 1295, 1292, 1289,
    1286, 1283, 1280, 1276, 1273, 1270, 1267, 1264, 1261, 1258, 1255, 1252,
    1249, 1246, 1243, 1240, 1237, 1234, 1231, 1228, 1226, 1223, 1220, 1217,
    1214, 1211, 1209, 1206, 1203, 1200, 1197, 1195, 1192, 1189, 1187, 1184,
    1181, 1179, 1176, 1173, 1171, 1168, 1165, 1163, 1160, 1158, 1155, 1153,
    1150, 1148, 1145, 1143, 1140, 1138, 1135, 1133, 1130, 1128, 1125, 1123,
    1121, 1118, 1116, 1113, 1111, 1109, 1106, 1104, 1102, 1099, 1097, 1095,
    1092, 1090, 1088, 1086, 1083, 1081, 1079, 1077, 1074, 1072, 1070, 1068,
    1066, 1064, 1061, 1059, 1057, 1055, 1053, 1051, 1049, 1047, 1044, 1042,
    1040, 1038, 1036, 1034, 1032, 1030, 1028, 1026, 1024};

  const uint64_t
      d0 = d % 2,
      d9 = d >> 55, /* 2^8 <= d9 < 2^9 */
      d40 = (d >> 24) + 1, /* 2^39 + 1 <= d40 < 2^40+1 */
      d63 = d / 2 + d0; /* ceil(d/2), 2^62 <= d63 <= 2^63 */
  
  uint64_t
    v0 = lut[d9 - 256],
    v1 = (1<<11)*v0 - (v0*v0*d40 >> 40) - 1,
    v2 = (v1 << 13) + (v1 * ((one << 60) - v1 * d40) >> 47),
    p1, p2;
  
  u64arith_mul_1_1_2(&p1, &p2, v2, ((v2 / 2) & (-d0)) - v2 * d63);
  uint64_t v3 = (v2 << 31) + p2 / 2;
  
  u64arith_mul_1_1_2(&p1, &p2, v3, d);
  u64arith_add_2_2(&p1, &p2, d, d);
  uint64_t v4 = v3 - p2;

  return v4;
}

/* Given D=(d0,d1), compute (v, shift) such that
 * 2^shift * D in [2^127, 2^128) and 
 * v + 2^64 = (2^192-1) / (2^shift*D)
 */
static inline uint64_t
u64arith_reciprocal_for_div_3by2(const uint64_t d0, const uint64_t d1)
{
    uint64_t v = u64arith_reciprocal_for_div(d1);
    uint64_t p = d1*v + d0;
    if (p < d0) {
        v--;
        if (p >= d1) {
            v--;
            p -= d1;
        }
        p -= d1;
    }
    uint64_t t0, t1;
    u64arith_mul_1_1_2(&t0, &t1, v, d0);
    p += t1;
    if (p < t1) {
        v--;
        if (u64arith_ge_2_2(t0, p, d0, d1)) {
            v--;
        }
    }
    return v;
}


static inline void
u64arith_divqr_2_1_1_recip_precomp (uint64_t *q, uint64_t *r,
		      const uint64_t a1, const uint64_t a2,
		      const uint64_t d, const uint64_t v,
		      const int s)
{
  uint64_t u0 = a1, u1 = a2;
  /* Adjust dividend to match divisor */
  ASSERT_EXPENSIVE(0 <= s && s < 64);
  ASSERT_EXPENSIVE((u1 & ~(UINT64_MAX >> s)) == 0);
  u64arith_shl_2(&u0, &u1, s);

  uint64_t q0, q1;
  u64arith_mul_1_1_2(&q0, &q1, v, u1);
  u64arith_add_2_2(&q0, &q1, u0, u1);
  q1++;
  uint64_t r0 = u0 - q1*d;
  if (r0 > q0) {
    q1--;
    r0 += d;
  }
  if (r0 >= d) {
    q1++;
    r0 -= d;
  }
  r0 >>= s;

  *q = q1;
  *r = r0;
}

/* Integer division of two uint64_t values a2:a1 by a uint64_t divisor b.
   Returns quotient and remainder. Uses the algorithm described in
   Niels Möller and Torbjörn Granlund, Improved Division by Invariant
   Integers, IEEE Transactions on Computers (Volume: 60, Issue: 2,
   Feb. 2011) */

static inline void
u64arith_divqr_2_1_1_recip (uint64_t *q, uint64_t *r,
		      const uint64_t a1, const uint64_t a2,
		      const uint64_t b)
{
  if (a2 == 0) {
    *q = a1 / b;
    *r = a1 % b;
    return;
  }

  const int s = u64arith_clz(b);
  uint64_t d = b << s;
  /* Left-adjust divisor */
  const uint64_t v = u64arith_reciprocal_for_div(d);
  u64arith_divqr_2_1_1_recip_precomp(q, r, a1, a2, d, v, s);

#ifdef WANT_ASSERT_EXPENSIVE
  uint64_t P1, P2;
  u64arith_mul_1_1_2(&P1, &P2, *q, b);
  unsigned char cy = u64arith_add_2_2_cy(&P1, &P2, *r, 0);
  // printf("a1=%lu, a2=%lu, b=%lu, s=%d, q1=%lu, r0=%lu\n", a1, a2, b, s, q1, r0);
  ASSERT_EXPENSIVE(P1 == a1 && P2 == a2 && cy == 0);
#endif
}


/* Integer division of two uint64_t values a2:a1 by a uint64_t divisor b.
   Returns quotient and remainder. */

static inline void
u64arith_divqr_2_1_1 (uint64_t *q, uint64_t *r,
		      const uint64_t a1, const uint64_t a2,
		      const uint64_t b)
{
  ASSERT(a2 < b); /* Or there will be quotient overflow */
#ifdef U64ARITH_VERBOSE_ASM
  __asm__ ("# u64arith_divqr_2_1_1 (%0, %1, %2, %3, %4)\n" : : 
           "X" (*q), "X" (*r), "X" (a1), "X" (a2), "X" (b));
#endif
#if !defined (U64ARITH_NO_ASM) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
  __asm__ __VOLATILE (
    "divq %4"
    : "=a" (*q), "=d" (*r)
    : "0" (a1), "1" (a2), "rm" (b)
    : "cc");
#else 
  u64arith_divqr_2_1_1_recip(q, r, a1, a2, b);
#endif
}


/* Integer division of two uint64_t values by a uint64_t divisor. Returns
   only remainder. */

static inline void
u64arith_divr_2_1_1 (uint64_t *r, uint64_t a1, const uint64_t a2,
                     const uint64_t b)
{
  ASSERT(a2 < b); /* Or there will be quotient overflow */
#ifdef U64ARITH_VERBOSE_ASM
  __asm__ ("# u64arith_divr_2_1_1 (%0, %1, %2, %3)\n" : : 
           "X" (*r), "X" (a1), "X" (a2), "X" (b));
#endif
#if !defined (U64ARITH_NO_ASM) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
  __asm__ __VOLATILE (
    "divq %3"
    : "+a" (a1), "=d" (*r)
    : "1" (a2), "rm" (b)
    : "cc");
#else
  uint64_t dummy;
  u64arith_divqr_2_1_1(&dummy, r, a1, a2, b);
#endif
}

/* divide U=(u0,u1,u2) by D=(d0,d1), compute quotient in q[0] and
 * remainder in [*R0,*R1]. U must be < 2^64*D.
 *
 * Helper information is (v, shift) as computed by
 * u64arith_reciprocal_for_div_3by2
 */
static inline void
u64arith_divqr_3_2_1_recip_precomp (
    uint64_t *q, uint64_t *R0, uint64_t *R1,
    uint64_t u0, uint64_t u1, uint64_t u2,
    const uint64_t d0, const uint64_t d1, const uint64_t v,
    const int s)
{
    uint64_t q0, q1, t0, t1, r0, r1;
    ASSERT((d1 & (UINT64_C(1) << 63)) != 0);
    ASSERT(u64arith_lt_2_2(u1, u2, d0, d1));
    ASSERT(0 <= s && s < 64);
    ASSERT((u2 & ~(UINT64_MAX >> s)) == 0);

    u64arith_shld(&u2, u1, u2, s);
    u64arith_shld(&u1, u0, u1, s);
    u0 <<= s;

    u64arith_mul_1_1_2(&q0, &q1, u2, v);
    u64arith_add_2_2(&q0, &q1, u1, u2);
    r1 = u1 - q1*d1;
    u64arith_mul_1_1_2(&t0, &t1, d0, q1);
    r0 = u0;
    u64arith_sub_2_2(&r0, &r1, t0, t1);
    u64arith_sub_2_2(&r0, &r1, d0, d1);
    q1++;
    if (r1 >= q0) {
        q1--;
        u64arith_add_2_2(&r0, &r1, d0, d1);
    }
    if (u64arith_ge_2_2(r0, r1, d0, d1)) {
        q1++;
        u64arith_sub_2_2(&r0, &r1, d0, d1);
    }
    u64arith_shr_2(&r0, &r1, s);
    *R0 = r0;
    *R1 = r1;
    *q = q1;
}

static inline void
u64arith_divqr_3_2_1 (uint64_t *q, uint64_t *R0, uint64_t *R1,
    uint64_t u0, uint64_t u1, uint64_t u2,
    uint64_t d0, uint64_t d1)
{
    ASSERT(d1 != 0);
    ASSERT(u64arith_lt_2_2(u1, u2, d0, d1));
    int s = u64arith_clz(d1);
    u64arith_shl_2(&d0, &d1, s);
    const uint64_t v = u64arith_reciprocal_for_div_3by2(d0, d1);
    u64arith_divqr_3_2_1_recip_precomp(q, R0, R1, u0, u1, u2, d0, d1, v, s);
}

/* Compute 1/n (mod 2^wordsize) */
static inline uint64_t
u64arith_invmod (const uint64_t n)
{
  uint64_t r;

  ASSERT (n % 2 != 0);
  
  /* Suggestion from PLM: initing the inverse to (3*n) XOR 2 gives the
     correct inverse modulo 32, then 3 (for 32 bit) or 4 (for 64 bit) 
     Newton iterations are enough. */
  r = (3 * n) ^ 2;
  /* Newton iteration */
  r = 2 * r - (uint32_t) r * (uint32_t) r * (uint32_t) n;
  r = 2 * r - (uint32_t) r * (uint32_t) r * (uint32_t) n;
  r = 2 * r - (uint32_t) r * (uint32_t) r * (uint32_t) n;
#if 0
  r = 2 * r - r * r * n;
#else
  /*
    r*n = k*2^32 + 1
        
    r' = 2 * r - r * r * n
    r' = 2 * r - r * (k*2^32 + 1)
    r' = 2 * r - r * k*2^32 - r
    r' = r - r * k*2^32
    r' = r - ((r * k) % 2^32) * 2^32
  */
  uint32_t k = (uint32_t)(r * n >> 32);
  k *= (uint32_t) r;
  r = r - ((uint64_t)k << 32);
#endif

  ASSERT_EXPENSIVE(n * r == 1);
  return r;
}

/* Compute n/2 (mod m), where m must be odd. */
static inline uint64_t
u64arith_div2mod (const uint64_t n, const uint64_t m)
{
#if !defined (U64ARITH_NO_ASM) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
    uint64_t N = n, M = m/2, t = 0;
    ASSERT_EXPENSIVE (m % 2 != 0);
    __asm__ __VOLATILE(
        "shr $1, %1\n\t" /* N /= 2 */
        "cmovc %0, %2\n\t" /* if (cy) {t = M;} */
        "adc %2, %1\n\t" /* N += t + cy */
        : "+&r" (M), "+&r" (N), "+&r" (t)
        : : "cc"
    );
  return N;
#else
  ASSERT_EXPENSIVE (m % 2 != 0);
  if (n % 2 == 0)
    return n / 2;
  else
    return n / 2 + m / 2 + 1;
#endif
}


/* Integer (truncated) square root of n */
static inline uint64_t
u64arith_sqrt (const uint64_t n)
{
    /* TODO: use Hensel lifting instead */
  unsigned int i;
  uint64_t xs, c, d, s2;
  const unsigned int l = 63 - (unsigned int)__builtin_clzl(n);

  d = n; /* d = n - x^2 */
  xs = 0;
  s2 = (uint64_t)1 << (l - l % 2);

  for (i = l / 2; i != 0; i--)
    {
      /* Here, s2 = 1 << (2*i) */
      /* xs = x << (i + 1), the value of x shifted left i+1 bits */

      c = xs + s2; /* c = (x + 2^i) ^ 2 - x^2 = 2^(i+1) * x + 2^(2*i) */
      xs >>= 1; /* Now xs is shifted only i positions */
      if (d >= c)
        {
          d -= c;
          xs |= s2; /* x |= 1 << i <=> xs |= 1 << (2*i) */
        }
      s2 >>= 2;
    }

  c = xs + s2;
  xs >>= 1;   
  if (d >= c) 
    xs |= s2;
  
  return xs;
}

/* Given r = -rem/p (mod den), we want num/(den*2^k) (mod p) ==
   (ratio + rem/den)/2^k (mod p).
   Using (a variant of) Bezout's identity, we have, for some non-negative
   integer t,
   r * p - t * den = -rem, or
   r * p + rem = t * den,
   thus den | (r * p + rem), and thus
   t = (r * p + rem) / den is an integer and satisfies
   t = rem/den (mod p).

   We have 0 <= r <= den-1 and rem <= den-1, and thus
   0 <= t = p * r/den + rem/den <=
   p (1 - 1/den) + 1 - 1/den =
   p + 1 - (p + 1)/den < p + 1.
   Thus t is almost a properly reduced residue for rem/den (mod p).
   As p fits in uint64_t, so does t, and we can compute t modulo
   2^64; since den is odd, we can multiply by den^{-1} mod 2^64
   to effect division by den.

   Finally we compute (t + ratio)/2^k mod p = num/(den*2^k) mod p.  */

static inline uint64_t
u64arith_post_process_inverse(const uint64_t r, const uint64_t p,
  const uint64_t rem, const uint64_t den_inv,
  const uint64_t ratio, const int k)
{
  uint64_t t = (r * p + rem) * den_inv;
  const uint64_t ratio_p = (ratio >= p) ? ratio % p : ratio;
  ASSERT_ALWAYS(t <= p); /* Cheap and fairly strong test */
  /* u64arith_addmod_1_1() accepts third operand == p and still produces
     a properly reduced sum mod p. */
  u64arith_addmod_1_1 (&t, ratio_p, t, p);

  ASSERT_EXPENSIVE(t < p);
  ASSERT_EXPENSIVE(k == 0 || p % 2 == 1);
  for (int j = 0; j < k; j++) {
    t = u64arith_div2mod(t, p);
  }
  return t;
}


/* Computes r = ((phigh * 2^64 + plow) / 2^64) % m
   Requires phigh < m and invm = -1/m (mod 2^64). */

static inline void
u64arith_redc(uint64_t *r, const uint64_t plow,
             const uint64_t phigh, const uint64_t m,
             const uint64_t invm)
{
  uint64_t t = phigh;

  ASSERT_EXPENSIVE (phigh < m);

#if !defined (U64ARITH_NO_ASM) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)

  /* TODO: are the register constraints watertight?
     %rax gets modified but putting tlow as an output constraint with "+"
     will keep r from getting allocated in %rax, which is a shame
     since we often want the result in %rax for the next multiply. */

  __asm__ __VOLATILE (
    "imulq %[invm], %%rax\n\t"
    "cmpq $1, %%rax \n\t"                /* if plow != 0, increase t */
    "sbbq $-1, %[t]\n\t"
    "mulq %[m]\n\t"
    "lea (%[t],%%rdx,1), %[r]\n\t"  /* compute (rdx + thigh) mod m */
    "subq %[m], %[t]\n\t"
    "addq %%rdx, %[t]\n\t"
    "cmovcq %[t], %[r]\n\t"
    : [t] "+&r" (t), [r] "=&r" (r[0])
    : [invm] "rm" (invm), [m] "rm" (m), "a" (plow)
    : "%rdx", "cc"
  );
#else
  uint64_t tlow, thigh;
  tlow = plow * invm;
  u64arith_mul_1_1_2 (&tlow, &thigh, tlow, m);
  /* Let w = 2^wordsize. We know (phigh * w + plow) + (thigh * w + tlow)
     == 0 (mod w) so either plow == tlow == 0, or plow !=0 and tlow != 0.
     In the former case we want phigh + thigh + 1, in the latter
     phigh + thigh. Since t = phigh < m, and modredcul_add can handle the
     case where the second operand is equal to m, adding 1 is safe */

  t += (plow != 0) ? 1 : 0; /* Does not depend on the mul */

  u64arith_addmod_1_1(r, t, thigh, m);
#endif
  ASSERT_EXPENSIVE (r[0] < m);
}


#ifdef __cplusplus
}
#endif


#endif /* ifndef CADO_U64ARITH_H */
