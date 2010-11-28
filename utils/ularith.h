/* Some commonly used assembly helper macros for arithmetic on 
   unsigned long. 
   Defining ULARITH_VERBOSE_ASM puts in the asm output a line for each
   function call that shows what registers/memory locations the operands 
   are in.
   Defining ULARITH_NO_ASM avoids asm macros and uses the C fallback code
   where available.
*/

#ifndef UL_ARITH_H__

#define UL_ARITH_H__

#include <assert.h>
#include <limits.h>
#include <gmp.h>
#include "macros.h"

/* <limits.h> defines LONG_BIT only with _XOPEN_SOURCE defined, but if 
   another header (such as <stdio.h>) already included <features.h> before 
   _XOPEN_SOURCE was set to 1, future includes of <features.h> are
   short-circuited and _XOPEN_SOURCE is ignored. */

#ifndef LONG_BIT
#ifdef LONG_MAX /* ISO C99 requires LONG_MAX in limits.h */
#if LONG_MAX == 2147483647L
#define LONG_BIT 32
#else
#define LONG_BIT 64
#endif /* if LONG_MAX == 2147483647L */
#elif defined __LONG_MAX__
#if __LONG_MAX__ == 2147483647L
#define LONG_BIT 32
#else
#define LONG_BIT 64
#endif /* if __LONG_MAX__ == 2147483647L */
#else /* elif defined __LONG_MAX__ */
#error Cannot guess LONG_BIT, please define
#endif /* elif defined __LONG_MAX__ else */
#endif /* ifndef LONG_BIT */

#ifndef ASSERT
#define ASSERT(x)	assert(x)
#endif

#ifdef WANT_ASSERT_EXPENSIVE
#define ASSERT_EXPENSIVE(x) ASSERT(x)
#else
#define ASSERT_EXPENSIVE(x)
#endif

/* On 32 bit x86, the general constrains for, e.g., the source operand
   of add is "g". For x86_64, it is "rme", since immediate constants
   must be 32 bit. */
#if defined(__i386__) && defined(__GNUC__)
#define ULARITH_CONSTRAINT_G "g"
#elif defined(__x86_64__) && defined(__GNUC__)
#define ULARITH_CONSTRAINT_G "rme"
#endif


/* Increases r if a != 0 */
static inline void
ularith_inc_nz (unsigned long *r, const unsigned long a)
{
#ifdef ULARITH_VERBOSE_ASM
  __asm__ ("# ularith_inc_nz (%0, %1)\n" : : "X" (*r), "X" (a));
#endif
#if !defined (ULARITH_NO_ASM) && defined(__x86_64__) && defined(__GNUC__)
  __asm__ ( "cmpq $1, %1\n\t"
            "sbbq $-1, %0\n\t"
    : "+r" (*r)
    : "rm" (a)
    : "cc");
#elif !defined (ULARITH_NO_ASM) && defined(__i386__) && defined(__GNUC__)
  __asm__ ( "cmpl $1, %1\n\t"
            "sbbl $-1, %0\n\t"
    : "+r" (*r)
    : "rm" (a)
    : "cc");
#else
  if (a != 0)
    r += 1;
#endif
}

/* Add an unsigned long to two unsigned longs with carry propagation from 
   low word (r1) to high word (r2). Any carry out from high word is lost. */

static inline void
ularith_add_ul_2ul (unsigned long *r1, unsigned long *r2, 
			const unsigned long a)
{
#ifdef ULARITH_VERBOSE_ASM
  __asm__ ("# ularith_add_ul_2ul (%0, %1, %2)\n" : : 
           "X" (*r1), "X" (*r2), "X" (a));
#endif
#if !defined (ULARITH_NO_ASM) && defined(__x86_64__) && defined(__GNUC__)
  __asm__ ( "addq %2, %0\n\t"
            "adcq $0, %1\n"
            : "+&r" (*r1), "+r" (*r2) 
            : "rme" (a)
            : "cc"); /* TODO: add commutativity and alternative for add to 
                        memory */
#elif !defined (ULARITH_NO_ASM) && defined(__i386__) && defined(__GNUC__)
  __asm__ ( "addl %2, %0\n\t"
            "adcl $0, %1\n"
            : "+&r" (*r1), "+r" (*r2)
            : "g" (a)
            : "cc");
#else
  *r1 += a;
  if (*r1 < a)
    (*r2)++;
#endif
}


/* Add two unsigned longs to two unsigned longs with carry propagation from 
   low word (r1) to high word (r2). Any carry out from high word is lost. */

static inline void
ularith_add_2ul_2ul (unsigned long *r1, unsigned long *r2, 
			 const unsigned long a1, const unsigned long a2)
{
#ifdef ULARITH_VERBOSE_ASM
  __asm__ ("# ularith_add_2ul_2ul (%0, %1, %2, %3)\n" : : 
           "X" (*r1), "X" (*r2), "X" (a1), "X" (a2));
#endif
#if !defined (ULARITH_NO_ASM) && defined(__x86_64__) && defined(__GNUC__)
  __asm__ ( "addq %2, %0\n\t"
            "adcq %3, %1\n"
            : "+&r" (*r1), "+r" (*r2)
            : "rme" (a1), "rme" (a2)
            : "cc");
#elif !defined (ULARITH_NO_ASM) && defined(__i386__) && defined(__GNUC__)
  __asm__ ( "addl %2, %0\n\t"
            "adcl %3, %1\n"
            : "+&r" (*r1), "+r" (*r2)
            : "g" (a1), "g" (a2)
            : "cc");
#else
  *r1 += a1;
  *r2 += a2 + (*r1 < a1);
#endif
}

/* Adds two unsigned longs from two unsigned longs with carry propagation 
   from low word (r1) to high word (r2). Returns 1 if there was a carry out 
   from high word, otherwise returns 0. */

static inline char
ularith_add_2ul_2ul_cy (unsigned long *r1, unsigned long *r2, 
			const unsigned long a1, const unsigned long a2)
{
  char cy;
#ifdef ULARITH_VERBOSE_ASM
  __asm__ ("# ularith_add_2ul_2ul_cy (%0, %1, %2, %3)\n" : : 
           "X" (*r1), "X" (*r2), "X" (a1), "X" (a2));
#endif
#if !defined (ULARITH_NO_ASM) && defined(__x86_64__) && defined(__GNUC__)
  __asm__ ( "addq %3, %0\n\t"
            "adcq %4, %1\n\t"
	    "setc %2\n"
            : "+&r" (*r1), "+r" (*r2), "=r" (cy)
            : "rme" (a1), "rme" (a2)
            : "cc");
#elif !defined (ULARITH_NO_ASM) && defined(__i386__) && defined(__GNUC__)
  __asm__ ( "addl %3, %0\n\t"
            "adcl %4, %1\n"
	    "setc %2\n"
            : "+&r" (*r1), "+r" (*r2), "=r" (cy)
            : "g" (a1), "g" (a2)
            : "cc");
#else
  unsigned long u1 = *r1 + a1, u2 = *r2 + a2;
  if (u1 < *r1)
    u2++;
  cy = (u2 < *r2) ? 1 : 0;
  *r1 = u1;
  *r2 = u2;
#endif
  return cy;
}

/* Subtract an unsigned long from two unsigned longs with borrow propagation 
   from low word (r1) to high word (r2). Any borrow out from high word is 
   lost. */

static inline void
ularith_sub_ul_2ul (unsigned long *r1, unsigned long *r2, 
			const unsigned long a)
{
#ifdef ULARITH_VERBOSE_ASM
  __asm__ ("# ularith_sub_ul_2ul  (%0, %1, %2)\n" : : 
           "X" (*r1), "X" (*r2), "X" (a));
#endif
#if !defined (ULARITH_NO_ASM) && defined(__x86_64__) && defined(__GNUC__)
  __asm__ ( "subq %2, %0\n\t"
            "sbbq $0, %1\n"
            : "+&r" (*r1), "+r" (*r2)
            : "rme" (a)
            : "cc");
#elif !defined (ULARITH_NO_ASM) && defined(__i386__) && defined(__GNUC__)
  __asm__ ( "subl %2, %0\n\t"
            "sbbl $0, %1\n"
            : "+&r" (*r1), "+r" (*r2)
            : "g" (a)
            : "cc");
#else
  unsigned long u = *r1;
  *r1 -= a;
  if (*r1 > u)
    (*r2)--;
#endif
}


/* Subtract two unsigned longs from two unsigned longs with borrow propagation 
   from low word (r1) to high word (r2). Any borrow out from high word is 
   lost. */

static inline void
ularith_sub_2ul_2ul (unsigned long *r1, unsigned long *r2, 
			 const unsigned long a1, const unsigned long a2)
{
#ifdef ULARITH_VERBOSE_ASM
  __asm__ ("# ularith_sub_2ul_2ul (%0, %1, %2, %3)\n" : : 
           "X" (*r1), "X" (*r2), "X" (a1), "X" (a2));
#endif
#if !defined (ULARITH_NO_ASM) && defined(__x86_64__) && defined(__GNUC__)
  __asm__ ( "subq %2, %0\n\t"
            "sbbq %3, %1\n"
            : "+&r" (*r1), "+r" (*r2)
            : "rme" (a1), "rme" (a2)
            : "cc");
#elif !defined (ULARITH_NO_ASM) && defined(__i386__) && defined(__GNUC__)
  __asm__ ( "subl %2, %0\n\t"
            "sbbl %3, %1\n"
            : "+&r" (*r1), "+r" (*r2)
            : "g" (a1), "g" (a2)
            : "cc");
#else
  unsigned long u = *r1;
  *r1 -= a1;
  *r2 -= a2;
  if (*r1 > u)
    (*r2)--;
#endif
}

/* Subtract two unsigned longs from two unsigned longs with borrow propagation 
   from low word (r1) to high word (r2). Returns 1 if there was a borrow out 
   from high word, otherwise returns 0. */

static inline char
ularith_sub_2ul_2ul_cy (unsigned long *r1, unsigned long *r2, 
			const unsigned long a1, const unsigned long a2)
{
  char cy;
#ifdef ULARITH_VERBOSE_ASM
  __asm__ ("# ularith_sub_2ul_2ul_cy (%0, %1, %2, %3)\n" : : 
           "X" (*r1), "X" (*r2), "X" (a1), "X" (a2));
#endif
#if !defined (ULARITH_NO_ASM) && defined(__x86_64__) && defined(__GNUC__)
  __asm__ ( "subq %3, %0\n\t"
            "sbbq %4, %1\n\t"
	    "setc %2\n"
            : "+&r" (*r1), "+r" (*r2), "=r" (cy)
            : "rme" (a1), "rme" (a2)
            : "cc");
#elif !defined (ULARITH_NO_ASM) && defined(__i386__) && defined(__GNUC__)
  __asm__ ( "subl %3, %0\n\t"
            "sbbl %4, %1\n"
	    "setc %2\n"
            : "+&r" (*r1), "+r" (*r2), "=q" (cy)
            : "g" (a1), "g" (a2)
            : "cc");
#else
  unsigned long u1 = *r1 - a1, u2 = *r2 - a2;
  if (u1 > *r1)
    u2--;
  cy = (u2 > *r2) ? 1 : 0;
  *r1 = u1;
  *r2 = u2;
#endif
  return cy;
}

/* Subtract only if result is non-negative */

static inline void
ularith_sub_2ul_2ul_ge (unsigned long *r1, unsigned long *r2, 
			const unsigned long a1, const unsigned long a2)
{
  unsigned long t1 = *r1, t2 = *r2;
#ifdef ULARITH_VERBOSE_ASM
  __asm__ ("# ularith_sub_2ul_2ul_ge (%0, %1, %2, %3)\n" : : 
           "X" (*r1), "X" (*r2), "X" (a1), "X" (a2));
#endif
#if !defined (ULARITH_NO_ASM) && defined(__x86_64__) && defined(__GNUC__)
  __asm__ ( "subq %4, %0\n\t" /* r1 -= a1 */
            "sbbq %5, %1\n\t" /* r2 -= a2 + cy */
	    "cmovc %2, %0\n\t" /* If there's a borrow, restore r1 from t1 */
	    "cmovc %3, %1\n\t" /* and r2 from t2 */
            : "+&r" (*r1), "+&r" (*r2)
            : "r" (t1), "r" (t2), "rme" (a1), "rme" (a2)
            : "cc");
#elif !defined (ULARITH_NO_ASM) && defined(__i386__) && defined(__GNUC__)
  __asm__ ( "subl %4, %0\n\t"
            "sbbl %5, %1\n\t"
	    "cmovc %2, %0\n\t"
	    "cmovc %3, %1\n\t"
            : "+&r" (*r1), "+&r" (*r2)
            : "r" (t1), "r" (t2), "g" (a1), "g" (a2)
            : "cc");
#else
  t1 -= a1;
  t2 -= a2 + (t1 > *r1);
  if (t2 <= *r2)
    {
      *r1 = t1;
      *r2 = t2;
    }
#endif
}


/* Multiply two unsigned long "a" and "b" and put the result as 
   r2:r1 (r2 being the high word) */

static inline void
ularith_mul_ul_ul_2ul (unsigned long *r1, unsigned long *r2, 
			   const unsigned long a, const unsigned long b)
{
#ifdef ULARITH_VERBOSE_ASM
  __asm__ ("# ularith_mul_ul_ul_2ul (%0, %1, %2, %3)\n" : : 
           "X" (*r1), "X" (*r2), "X" (a), "X" (b));
#endif
#if !defined (ULARITH_NO_ASM) && defined(__x86_64__) && defined(__GNUC__)
  __asm__ ( "mulq %3"
	    : "=a" (*r1), "=d" (*r2)
	    : "%0" (a), "rm" (b)
	    : "cc");
#elif !defined (ULARITH_NO_ASM) && defined(__i386__) && defined(__GNUC__)
  __asm__ ( "mull %3"
	    : "=a" (*r1), "=d" (*r2)
	    : "%0" (a), "rm" (b)
	    : "cc");
#else
  const int half = LONG_BIT / 2;
  const unsigned long mask = (1UL << half) - 1UL;
  unsigned long t1, t2, p1, p2;

  t1 = (a & mask) * (b & mask);
  t2 = 0UL;
  p1 = (a >> half) * (b & mask);
  p2 = p1 >> half;
  p1 = (p1 & mask) << half;
  ularith_add_2ul_2ul (&t1, &t2, p1, p2);
  p1 = (a & mask) * (b >> half);
  p2 = p1 >> half;
  p1 = (p1 & mask) << half;
  ularith_add_2ul_2ul (&t1, &t2, p1, p2);
  t2 += (a >> half) * (b >> half);
  *r1 = t1; 
  *r2 = t2;
#endif
}


static inline void
ularith_sqr_ul_2ul (unsigned long *r1, unsigned long *r2, 
		    const unsigned long a)
{
#ifdef ULARITH_VERBOSE_ASM
  __asm__ ("# ularith_mul_ul_ul_2ul (%0, %1, %2)\n" : : 
           "X" (*r1), "X" (*r2), "X" (a));
#endif
#if !defined (ULARITH_NO_ASM) && defined(__x86_64__) && defined(__GNUC__)
  __asm__ ( "mulq %%rax"
	    : "=a" (*r1), "=d" (*r2)
	    : "0" (a)
	    : "cc");
#elif !defined (ULARITH_NO_ASM) && defined(__i386__) && defined(__GNUC__)
  __asm__ ( "mull %%eax"
	    : "=a" (*r1), "=d" (*r2)
	    : "0" (a)
	    : "cc");
#else
  const int half = LONG_BIT / 2;
  const unsigned long mask = (1UL << half) - 1UL;
  unsigned long t1, t2, p1, p2;

  t1 = (a & mask) * (a & mask);
  t2 = 0UL;
  p1 = (a >> half) * (a & mask);
  p2 = p1 >> half;
  p1 = (p1 & mask) << half;
  ularith_add_2ul_2ul (&t1, &t2, p1, p2);
  ularith_add_2ul_2ul (&t1, &t2, p1, p2);
  t2 += (a >> half) * (a >> half);
  *r1 = t1; 
  *r2 = t2;
#endif
}


/* Integer division of a two ulong value a2:a1 by a ulong divisor. Returns
   quotient and remainder. */

static inline void
ularith_div_2ul_ul_ul (unsigned long *q, unsigned long *r, 
			   const unsigned long a1, const unsigned long a2, 
			   const unsigned long b)
{
  ASSERT(a2 < b); /* Or there will be quotient overflow */
#ifdef ULARITH_VERBOSE_ASM
  __asm__ ("# ularith_div_2ul_ul_ul (%0, %1, %2, %3, %4)\n" : : 
           "X" (*q), "X" (*r), "X" (a1), "X" (a2), "X" (b));
#endif
#if defined(__x86_64__) && defined(__GNUC__)
  __asm__ ( "divq %4"
            : "=a" (*q), "=d" (*r)
            : "0" (a1), "1" (a2), "rm" (b)
            : "cc");
#elif defined(__i386__) && defined(__GNUC__)
  __asm__ ( "divl %4"
            : "=a" (*q), "=d" (*r)
            : "0" (a1), "1" (a2), "rm" (b)
            : "cc");
#else
  mp_limb_t A[2] = {a1, a2};
  ASSERT(sizeof(unsigned long) == sizeof(mp_limb_t));
  r[0] = mpn_divmod_1 (A, A, 2, b);
  q[0] = A[0];
#endif
}


/* Integer division of a two longint value by a longint divisor. Returns
   only remainder. */

static inline void
ularith_div_2ul_ul_ul_r (unsigned long *r, unsigned long a1,
                 const unsigned long a2, const unsigned long b)
{
  ASSERT(a2 < b); /* Or there will be quotient overflow */
#ifdef ULARITH_VERBOSE_ASM
  __asm__ ("# ularith_div_2ul_ul_ul_r (%0, %1, %2, %3)\n" : : 
           "X" (*r), "X" (a1), "X" (a2), "X" (b));
#endif
#if defined(__x86_64__) && defined(__GNUC__)
  __asm__ ( "divq %3"
            : "+a" (a1), "=d" (*r)
            : "1" (a2), "rm" (b)
            : "cc");
#elif defined(__i386__) && defined(__GNUC__)
  __asm__ ( "divl %3"
            : "+a" (a1), "=d" (*r)
            : "1" (a2), "rm" (b)
            : "cc");
#else
  mp_limb_t A[2] = {a1, a2};
  ASSERT(sizeof(unsigned long) == sizeof(mp_limb_t));
  r[0] = mpn_divmod_1 (A, A, 2, b);
#endif
}


/* Shift *r right by i bits, filling in the low bits from a into the high
   bits of *r. Assumes i < LONG_BIT */
MAYBE_UNUSED
static inline void
ularith_shrd (unsigned long *r, const unsigned long a, const int i)
{
  ASSERT_EXPENSIVE (i < LONG_BIT);
#ifdef ULARITH_VERBOSE_ASM
  __asm__ ("# ularith_shrd (%0, %1, %2)\n" : : 
           "X" (*r), "X" (a), "X" (i));
#endif
#if !defined (ULARITH_NO_ASM) && defined(__x86_64__) && defined(__GNUC__)
  __asm__ ("shrdq %b2, %1, %0\n": 
  /* the b modifier makes gcc print the byte part of the register (%cl) */
           "+rm" (*r) :
           "r" (a), "cJ" (i) : /* i can be in %cx or a constant < 64 */
           "cc"
          );
#elif !defined (ULARITH_NO_ASM) && defined(__i386__) && defined(__GNUC__)
  __asm__ ("shrdl %b2, %1, %0\n": 
           "+rm" (*r) :
           "r" (a), "cI" (i) : /* i can be in %cx or a constant < 32 */
           "cc"
          );
#else
  if (i > 0) /* shl by LONG_BIT is no-op on x86! */
    *r = (*r >> i) | (a << (LONG_BIT - i));
#endif
}

/* Shift *r left by i bits, filling in the high bits from a into the low
   bits of *r. Assumes i < LONG_BIT */
MAYBE_UNUSED
static inline void
ularith_shld (unsigned long *r, const unsigned long a, const int i)
{
  ASSERT_EXPENSIVE (i < LONG_BIT);
#ifdef ULARITH_VERBOSE_ASM
  __asm__ ("# ularith_shld (%0, %1, %2)\n" : : 
           "X" (*r), "X" (a), "X" (i));
#endif
#if !defined (ULARITH_NO_ASM) && defined(__x86_64__) && defined(__GNUC__)
  __asm__ ("shldq %b2, %1, %0\n": 
  /* the b modifier makes gcc print the byte part of the register (%cl) */
           "+rm" (*r) :
           "r" (a), "cJ" (i) : /* i can be in %cx or a constant < 64 */
           "cc"
          );
#elif !defined (ULARITH_NO_ASM) && defined(__i386__) && defined(__GNUC__)
  __asm__ ("shldl %b2, %1, %0\n": 
           "+rm" (*r) :
           "r" (a), "cI" (i) : /* i can be in %cx or a constant < 32 */
           "cc"
          );
#else
  if (i > 0) /* shr by LONG_BIT is no-op on x86! */
  *r = (*r << i) | (a >> (LONG_BIT - i));
#endif
}

/* Returns number of trailing zeros in a. a must not be zero */
MAYBE_UNUSED
static inline int
ularith_ctz (const unsigned long a)
{
#if !defined (ULARITH_NO_ASM) && defined(__GNUC__) && \
    (__GNUC__ >= 4 || __GNUC__ >= 3 && __GNUC_MINOR__ >= 4)
  ASSERT_EXPENSIVE (a != 0UL);
  return __builtin_ctzl(a);
#else
  unsigned long t = a;
  int i;
  ASSERT_EXPENSIVE (a != 0UL);
  for (i = 0; (t & 1UL) == 0UL; i++)
    t >>= 1;
  return i;
#endif
}

/* Returns number of leading zeros in a. a must not be zero */
MAYBE_UNUSED
static inline int
ularith_clz (const unsigned long a)
{
#if !defined (ULARITH_NO_ASM) && defined(__GNUC__) && \
    (__GNUC__ >= 4 || __GNUC__ >= 3 && __GNUC_MINOR__ >= 4)
  ASSERT_EXPENSIVE (a != 0UL);
  return __builtin_clzl(a);
#else
  unsigned long t = 1UL << (LONG_BIT - 1);
  int i;
  ASSERT_EXPENSIVE (a != 0UL);
  for (i = 0; (a & t) == 0UL; i++)
    t >>= 1;
  return i;
#endif
}


/* Compute 1/n (mod 2^wordsize) */
MAYBE_UNUSED
static inline unsigned long
ularith_invmod (const unsigned long n)
{
  unsigned long r;

  ASSERT (n % 2UL != 0UL);
  
  /* Suggestion from PLM: initing the inverse to (3*n) XOR 2 gives the
     correct inverse modulo 32, then 3 (for 32 bit) or 4 (for 64 bit) 
     Newton iterations are enough. */
  r = (3UL * n) ^ 2UL;
  /* Newton iteration */
  r = 2UL * r - (unsigned int) r * (unsigned int) r * (unsigned int) n;
  r = 2UL * r - (unsigned int) r * (unsigned int) r * (unsigned int) n;
  if (sizeof (unsigned long) == 4)
    {
      r = 2UL * r - r * r * n;
    }
  else
    {
      r = 2UL * r - (unsigned int) r * (unsigned int) r * (unsigned int) n;
      r = 2UL * r - r * r * n;
    }

  return r;
}

#endif /* ifndef UL_ARITH_H__ */
