/* Some functions for modular arithmetic with residues and modulus in
   unsigned long variables. Due to inlining, this file must be included
   in the caller's source code with #include */

/* Naming convention: all function start with modul, for 
   MODulus Unsigned Long, followed by underscore, functionality of function
  (add, mul, etc), and possibly underscore and specification of what argument
  types the function takes (_ul, etc). */

#ifndef MOD_UL_H
#define MOD_UL_H

/**********************************************************************/
#include <assert.h>
#include <limits.h>
#include "macros.h"
#include "ularith.h"

#ifndef ASSERT
#define ASSERT(x)	assert(x)
#endif

/* Even simple assertions are relatively expensive in very simple functions.
   If we want them anyway to hunt a bug, define WANT_ASSERT_EXPENSIVE */
#ifdef WANT_ASSERT_EXPENSIVE
#define ASSERT_EXPENSIVE(x) ASSERT(x)
#else
#define ASSERT_EXPENSIVE(x)
#endif

/*********************************************************************/
/* Helper macros, see also ularith.h */

/* A macro for function renaming. All functions here start with modul_ */
#define MODUL_RENAME(x) modul_##x

#define MODUL_SIZE 1
#define MODUL_MAXBITS LONG_BIT

typedef unsigned long residueul_t[MODUL_SIZE];
typedef unsigned long modintul_t[MODUL_SIZE];
typedef unsigned long modulusul_t[MODUL_SIZE];

/* ================= Functions that are part of the API ================= */

/* Some functions for integers of the same width as the modulus */

MAYBE_UNUSED
static inline void
modul_intset (modintul_t r, const modintul_t s)
{
  r[0] = s[0];
}


MAYBE_UNUSED
static inline void
modul_intset_ul (modintul_t r, const unsigned long s)
{
  r[0] = s;
}


MAYBE_UNUSED
static inline int
modul_intequal (const modintul_t a, const modintul_t b)
{
  return (a[0] == b[0]);
}


MAYBE_UNUSED
static inline int
modul_intequal_ul (const modintul_t a, const unsigned long b)
{
  return (a[0] == b);
}


MAYBE_UNUSED
static inline int
modul_intcmp (const modintul_t a, const modintul_t b)
{
  return (a[0] < b[0]) ? -1 : (a[0] == b[0]) ? 0 : 1;
}


MAYBE_UNUSED
static inline int
modul_intcmp_ul (const modintul_t a, const unsigned long b)
{
  return (a[0] < b) ? -1 : (a[0] == b) ? 0 : 1;
}


MAYBE_UNUSED
static inline int
modul_intfits_ul (const modintul_t a MAYBE_UNUSED)
{
  return 1;
}


MAYBE_UNUSED
static inline void
modul_intadd (modintul_t r, const modintul_t a, const modintul_t b)
{
  r[0] = a[0] + b[0];
}


MAYBE_UNUSED
static inline void
modul_intsub (modintul_t r, const modintul_t a, const modintul_t b)
{
  r[0] = a[0] - b[0];
}


/* Returns the number of bits in a, that is, floor(log_2(n))+1. 
   For n==0 returns 0. */
MAYBE_UNUSED
static inline int
modul_intbits (const modintul_t a)
{
  int bits = 0;
  unsigned long n = a[0];
  while (n > 0UL) /* TODO: use clzl */
    {
      bits++;
      n >>= 1;
    }
  return bits;
}


/* r = n/d. We require d|n */
MAYBE_UNUSED
static inline void
modul_intdivexact (modintul_t r, const modintul_t n, const modintul_t d)
{
  r[0] = n[0] / d[0]; 
}


/* r = n%d */
MAYBE_UNUSED
static inline void
modul_intmod (modintul_t r, const modintul_t n, 
              const modintul_t d)
{
  r[0] = n[0] % d[0]; 
}


/* Functions for the modulus */

MAYBE_UNUSED
static inline void
modul_initmod_ul (modulusul_t m, const unsigned long s)
{
  m[0] = s;
}


MAYBE_UNUSED
static inline void
modul_initmod_uls (modulusul_t m, const modintul_t s)
{
  m[0] = s[0];
}


MAYBE_UNUSED
static inline unsigned long
modul_getmod_ul (const modulusul_t m)
{
  return m[0];
}


MAYBE_UNUSED
static inline void
modul_getmod_uls (modintul_t r, const modulusul_t m)
{
  r[0] = m[0];
}


MAYBE_UNUSED
static inline void
modul_clearmod (modulusul_t m MAYBE_UNUSED)
{
  return;
}


/* Functions for residues */

/* Initialises a residue_t type and sets it to zero */
MAYBE_UNUSED
static inline void
modul_init (residueul_t r, const modulusul_t m MAYBE_UNUSED)
{
  r[0] = 0UL;
  return;
}


/* Initialises a residue_t type, but does not set it to zero. For fixed length
   residue_t types, that leaves nothing to do at all. */
MAYBE_UNUSED
static inline void
modul_init_noset0 (residueul_t r MAYBE_UNUSED, 
		   const modulusul_t m MAYBE_UNUSED)
{
  return;
}


MAYBE_UNUSED
static inline void
modul_clear (residueul_t r MAYBE_UNUSED, 
             const modulusul_t m MAYBE_UNUSED)
{
  return;
}


MAYBE_UNUSED
static inline void
modul_set (residueul_t r, const residueul_t s, 
           const modulusul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (s[0] < m[0]);
  r[0] = s[0];
}


MAYBE_UNUSED
static inline void
modul_set_ul (residueul_t r, const unsigned long s, const modulusul_t m)
{
  r[0] = s % m[0];
}


/* Sets the residue_t to the class represented by the integer s. Assumes that
   s is reduced (mod m), i.e. 0 <= s < m */

MAYBE_UNUSED
static inline void
modul_set_ul_reduced (residueul_t r, const unsigned long s, 
                      const modulusul_t m MAYBE_UNUSED)
{
  ASSERT (s < m[0]);
  r[0] = s;
}


MAYBE_UNUSED
static inline void
modul_set_uls (residueul_t r, const modintul_t s, 
		   const modulusul_t m)
{
  r[0] = s[0] % m[0];
}


MAYBE_UNUSED
static inline void
modul_set_uls_reduced (residueul_t r, const modintul_t s, 
		       const modulusul_t m MAYBE_UNUSED)
{
  ASSERT (s[0] < m[0]);
  r[0] = s[0];
}


/* These two are so trivial that we don't really require m in the
 * interface. For 1 we might, as the internal representation might 
 * not use "1" for 1 (e.g. when using Montgomery's REDC.)
 * For interface homogeneity we make even modul_set0 take the m parameter.
 */
MAYBE_UNUSED 
static inline void 
modul_set0 (residueul_t r, const modulusul_t m MAYBE_UNUSED) 
{ 
  r[0] = 0UL;
}


MAYBE_UNUSED 
static inline void 
modul_set1 (residueul_t r, const modulusul_t m MAYBE_UNUSED) 
{ 
  if (m[0] == 1UL)
    r[0] = 0UL;
  else
    r[0] = 1UL;
}


/* Exchanges the values of the two arguments */

MAYBE_UNUSED
static inline void
modul_swap (residueul_t a, residueul_t b, 
            const modulusul_t m MAYBE_UNUSED)
{
  unsigned long t;
  ASSERT_EXPENSIVE (a[0] < m[0] && b[0] < m[0]);
  t = a[0];
  a[0] = b[0];
  b[0] = t;
}
                          

MAYBE_UNUSED
static inline unsigned long
modul_get_ul (const residueul_t s, const modulusul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (s[0] < m[0]);
  return s[0];
}


MAYBE_UNUSED
static inline void
modul_get_uls (modintul_t r, const residueul_t s, 
		   const modulusul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (s[0] < m[0]);
  r[0] = s[0];
}


MAYBE_UNUSED
static inline int
modul_equal (const residueul_t a, const residueul_t b, 
             const modulusul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (a[0] < m[0] && b[0] < m[0]);
  return (a[0] == b[0]);
}


MAYBE_UNUSED
static inline int
modul_is0 (const residueul_t a, const modulusul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (a[0] < m[0]);
  return (a[0] == 0UL);
}


MAYBE_UNUSED
static inline int
modul_is1 (const residueul_t a, const modulusul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (a[0] < m[0]);
  return (a[0] == 1UL);
}


MAYBE_UNUSED
static inline void
modul_add (residueul_t r, const residueul_t a, const residueul_t b, 
	   const modulusul_t m)
{
  ASSERT_EXPENSIVE (a[0] < m[0] && b[0] < m[0]);
#ifdef MODTRACE
  printf ("modul_add: a = %lu, b = %lu", a[0], b[0]);
#endif

#if (defined(__i386__) || defined(__x86_64__)) && defined(__GNUC__)
  {
    unsigned long t = a[0] - m[0], tr = a[0] + b[0];
    
    __asm__ (
      "add %2, %1\n\t"   /* t += b */
      "cmovc %1, %0\n\t"  /* if (cy) tr = t */
      : "+r" (tr), "+&r" (t)
      : "g" (b[0])
      : "cc"
    );
    ASSERT_EXPENSIVE (tr == ((a[0] >= m[0] - b[0]) ? (a[0] - (m[0] - b[0])) : (a[0] + b[0])));
    r[0] = tr;
  }
#else
  r[0] = (a[0] >= m[0] - b[0]) ? (a[0] - (m[0] - b[0])) : (a[0] + b[0]);
#endif

#ifdef MODTRACE
  printf (", r = %lu\n", r[0]);
#endif
  ASSERT_EXPENSIVE (r[0] < m[0]);
}


MAYBE_UNUSED
static inline void
modul_add_ul (residueul_t r, const residueul_t a, const unsigned long b, 
	      const modulusul_t m)
{
  unsigned long b2 = b % m[0];
  ASSERT_EXPENSIVE (a[0] < m[0]);
  r[0] = (a[0] >= m[0] - b2) ? (a[0] - (m[0] - b2)) : (a[0] + b2);
}

MAYBE_UNUSED
static inline void
modul_sub (residueul_t r, const residueul_t a, const residueul_t b, 
	   const modulusul_t m)
{
  ASSERT_EXPENSIVE (a[0] < m[0] && b[0] < m[0]);
#ifdef MODTRACE
  printf ("submod_ul: a = %lu, b = %lu", a[0], b[0]);
#endif

#if (defined(__i386__) || defined(__x86_64__)) && defined(__GNUC__)
  {
    unsigned long tr, t = a[0];
    __asm__ (
      "sub %2, %1\n\t"  /* t -= b ( = a - b) */
      "lea (%1,%3,1), %0\n\t" /* tr = t + m ( = a - b + m) */
      "cmovnc %1, %0\n\t" /* if (a >= b) tr = t */
      : "=&r" (tr), "+&r" (t)
      : "g" (b[0]), "r" (m[0])
      : "cc"
    );
    r[0] = tr;
  }
#elif 1
  /* Seems to be faster than the one below */
  {
    unsigned long t = 0, tr;
    tr = a[0] - b[0];
    if (tr > a[0])
      t = m[0];
    r[0] = tr + t;
  }
#else
  r[0] = (a[0] < b[0]) ? (a[0] - b[0] + m[0]) : (a[0] - b[0]);
#endif

#ifdef MODTRACE
  printf (", r = %lu\n", r[0]);
#endif
  ASSERT_EXPENSIVE (r[0] < m[0]);
}


MAYBE_UNUSED
static inline void
modul_sub_ul (residueul_t r, const residueul_t a, const unsigned long b, 
	      const modulusul_t m)
{
  unsigned long b2 = b % m[0];
  ASSERT_EXPENSIVE (a[0] < m[0]);
  r[0] = (a[0] >= b2) ? (a[0] - b2) : (a[0] + (m[0] - b2));
}


MAYBE_UNUSED
static inline void
modul_neg (residueul_t r, const residueul_t a, const modulusul_t m)
{
  ASSERT_EXPENSIVE (a[0] < m[0]);
  if (a[0] == 0UL)
    r[0] = a[0];
  else
    r[0] = m[0] - a[0];
}


MAYBE_UNUSED
static inline void
modul_mul (residueul_t r, const residueul_t a, const residueul_t b, 
           const modulusul_t m)
{
  unsigned long t1, t2;
#if defined(MODTRACE)
  printf ("mulmod_ul: a = %lu, b = %lu", a[0], b[0]);
#endif
  
  ASSERT_EXPENSIVE (a[0] < m[0] && b[0] < m[0]);
  ularith_mul_ul_ul_2ul (&t1, &t2, a[0], b[0]);
  ularith_div_2ul_ul_ul_r (r, t1, t2, m[0]);
  
#if defined(MODTRACE)
  printf (", r = %lu\n", r);
#endif
}

MAYBE_UNUSED
static inline void
modul_sqr (residueul_t r, const residueul_t a, const modulusul_t m)
{
  unsigned long t1, t2;
  
  ASSERT_EXPENSIVE (a[0] < m[0]);
  ularith_mul_ul_ul_2ul (&t1, &t2, a[0], a[0]);
  ularith_div_2ul_ul_ul_r (r, t1, t2, m[0]);
}

/* Computes (a * 2^wordsize) % m */
MAYBE_UNUSED
static inline void
modul_tomontgomery (residueul_t r, const residueul_t a, const modulusul_t m)
{
  ASSERT_EXPENSIVE (a[0] < m[0]);
  ularith_div_2ul_ul_ul_r (r, 0UL, a[0], m[0]);
}


/* Computes (a / 2^wordsize) % m */
MAYBE_UNUSED
static inline void
modul_frommontgomery (residueul_t r, const residueul_t a, 
                      const unsigned long invm, const modulusul_t m)
{
  unsigned long tlow, thigh;
  ularith_mul_ul_ul_2ul (&tlow, &thigh, a[0], invm); /* FIXME: why not IMUL? */
  ularith_mul_ul_ul_2ul (&tlow, &thigh, tlow, m[0]);
  r[0] = thigh + (a[0] != 0UL ? 1UL : 0UL);
}


/* Computes (a / 2^wordsize) % m, but result can be r = m. 
   Input a must not be equal 0 */
MAYBE_UNUSED
static inline void
modul_redcsemi_ul_not0 (residueul_t r, const unsigned long a, 
                        const unsigned long invm, const modulusul_t m)
{
  unsigned long tlow, thigh;

  ASSERT (a != 0);

  tlow = a * invm; /* tlow <= 2^w-1 */
  ularith_mul_ul_ul_2ul (&tlow, &thigh, tlow, m[0]);
  /* thigh:tlow <= (2^w-1) * m */
  r[0] = thigh + 1UL; 
  /* (thigh+1):tlow <= 2^w + (2^w-1) * m  <= 2^w + 2^w*m - m 
                    <= 2^w * (m + 1) - m */
  /* r <= floor ((2^w * (m + 1) - m) / 2^w) <= floor((m + 1) - m/2^w)
       <= m */
}


/* Computes ((a + b) / 2^wordsize) % m. a <= m is permissible */
MAYBE_UNUSED
static inline void
modul_addredc_ul (residueul_t r, const residueul_t a, const unsigned long b, 
                  const unsigned long invm, const modulusul_t m)
{
  unsigned long slow, shigh, tlow, thigh;
  
  ASSERT_EXPENSIVE (a[0] <= m[0]);
  slow = b;
  shigh = 0UL;
  ularith_add_ul_2ul (&slow, &shigh, a[0]);
  
  tlow = slow * invm;
  ularith_mul_ul_ul_2ul (&tlow, &thigh, tlow, m[0]);
  ASSERT_EXPENSIVE (slow + tlow == 0UL);
  r[0] = thigh + shigh + (slow != 0UL ? 1UL : 0UL);
  
  /* r = ((a+b) + (((a+b)%2^w * invm) % 2^w) * m) / 2^w  Use a<=m-1, b<=2^w-1
     r <= (m + 2^w - 1 + (2^w - 1) * m) / 2^w
        = (m - 1 + 2^w + m*2^w - m) / 2^w
        = (- 1 + 2^w + m2^w) / 2^w
        = m + 1 - 1/2^w
     r <= m, since r is an integer
  */
  if (r[0] == m[0])
    r[0] = 0UL;
}


/* Computes ((a + b) / 2^wordsize) % m, but result can be == m.
   a <= m is permissible */
MAYBE_UNUSED
static inline void
modul_addredcsemi_ul (residueul_t r, const residueul_t a, 
                      const unsigned long b, const unsigned long invm, 
                      const modulusul_t m)
{
  unsigned long slow, shigh, tlow;
  
  ASSERT_EXPENSIVE(a[0] <= m[0]);
  slow = b;
#if defined(__x86_64__) && defined(__GNUC__)
  {
    unsigned char sb;
    __asm__ ( "addq %2, %0\n\t" /* cy * 2^w + slow = a + b */
	      "setne %1\n\t"     /* if (slow != 0) sb = 1 */
	      "adcb $0, %1\n"    /* sb += cy */
	      : "+&r" (slow), "=qm" (sb)
	      : "rm" (a[0])
	      : "cc");
    shigh = sb;
  }
#elif defined(__i386__) && defined(__GNUC__)
  {
    unsigned char sb;
    __asm__ ( "addl %2, %0\n\t"
	      "setne %1\n\t"
	      "adcb $0, %1\n"
	      : "+&r" (slow), "=qm" (sb)
	      : "rm" (a[0])
	      : "cc");
    shigh = sb;
  }
#else
  shigh = 0UL;
  ularith_add_ul_2ul (&slow, &shigh, a[0]);
  shigh += (slow != 0UL ? 1UL : 0UL);
#endif

  tlow = slow * invm;
  ularith_mul_ul_ul_2ul (&tlow, r, tlow, m[0]);
  ASSERT_EXPENSIVE (slow + tlow == 0UL);
  r[0] += shigh;

  /* r = ((a+b) + (((a+b)%2^w * invm) % 2^w) * m) / 2^w
     r <= ((a+b) + (2^w - 1) * m) / 2^w
     r <= (m + 2^w-1 + m*2^w - m) / 2^w
     r <= (2^w -1 + p2^w) / 2^w
     r <= p + 1 - 1/2^w
     r <= p
  */
}


MAYBE_UNUSED
static inline void
modul_mulredc (residueul_t r, const residueul_t a, const residueul_t b,
               const unsigned long invm, const modulusul_t m)
{
  unsigned long plow, phigh, tlow, thigh;

  ASSERT_EXPENSIVE (m[0] % 2 != 0);
  ASSERT_EXPENSIVE (a[0] < m[0] && b[0] < m[0]);
#if defined(MODTRACE)
  printf ("(%lu * %lu / 2^%ld) %% %lu", 
          a[0], b[0], 8 * sizeof(unsigned long), m[0]);
#endif
  
  ularith_mul_ul_ul_2ul (&plow, &phigh, a[0], b[0]);
  tlow = plow * invm;
  ularith_mul_ul_ul_2ul (&tlow, &thigh, tlow, m[0]);
  /* Let w = 2^wordsize. We know (phigh * w + plow) + (thigh * w + tlow) 
     == 0 (mod w) so either plow == tlow == 0, or plow !=0 and tlow != 0. 
     In the former case we want phigh + thigh + 1, in the latter 
     phigh + thigh */
  /* Since a <= p-1 and b <= p-1, and p <= w-1, a*b <= w^2 - 4*w + 4, so
     adding 1 to phigh is safe */
#if 0
  /* Slower? */
  ularith_add_ul_2ul (&plow, &phigh, tlow);
#else
  phigh += (plow != 0UL ? 1UL : 0UL);
#endif

  modul_add (r, &phigh, &thigh, m);

#if defined(MODTRACE)
  printf (" == %lu /* PARI */ \n", r[0]);
#endif
}
                         

/* FIXME: check for overflow if b > m */
MAYBE_UNUSED
static inline void
modul_mulredc_ul (residueul_t r, const residueul_t a, const unsigned long b,
                  const unsigned long invm, const modulusul_t m)
{
  unsigned long plow, phigh, tlow, thigh;
  ASSERT_EXPENSIVE (m[0] % 2 != 0);
#if defined(MODTRACE)
  printf ("(%lu * %lu / 2^%ld) %% %lu", 
          a[0], b, 8 * sizeof(unsigned long), m[0]);
#endif
  
  ularith_mul_ul_ul_2ul (&plow, &phigh, a[0], b);
  tlow = plow * invm;
  ularith_mul_ul_ul_2ul (&tlow, &thigh, tlow, m[0]);
  phigh += (plow != 0UL ? 1UL : 0UL);
  r[0] = (phigh >= m[0] - thigh) ? (phigh - (m[0] - thigh)) : (phigh + thigh);
  
#if defined(MODTRACE)
  printf (" == %lu /* PARI */ \n", r[0]);
#endif
}
                         

/* Computes (a * b + c)/ 2^wordsize % m. Requires that 
   a * b + c < 2^wordsize * m */

MAYBE_UNUSED
static inline void
modul_muladdredc_ul (residueul_t r, const residueul_t a, const unsigned long b,
                     const unsigned long c, const unsigned long invm, 
                     const modulusul_t m)
{
  unsigned long plow, phigh, tlow, thigh;
  ASSERT_EXPENSIVE (m[0] % 2 != 0);
#if defined(MODTRACE)
  printf ("(%lu * %lu / 2^%ld) %% %lu", 
          a[0], b, 8 * sizeof(unsigned long), m[0]);
#endif
  
  ularith_mul_ul_ul_2ul (&plow, &phigh, a[0], b);
  ularith_add_ul_2ul (&plow, &phigh, c);
  tlow = plow * invm;
  ularith_mul_ul_ul_2ul (&tlow, &thigh, tlow, m[0]);
  phigh += (plow != 0UL ? 1UL : 0UL);
  r[0] = (phigh >= m[0] - thigh) ? (phigh - (m[0] - thigh)) : (phigh + thigh);
  
#if defined(MODTRACE)
  printf (" == %lu /* PARI */ \n", r[0]);
#endif
}
                         

MAYBE_UNUSED
static inline void
modul_div2 (residueul_t r, const residueul_t a, const modulusul_t m)
{
#if (defined(__i386__) || defined(__x86_64__)) && defined(__GNUC__)
  unsigned long s = a[0], t = m[0];
  ASSERT_EXPENSIVE (m[0] % 2UL != 0UL);

  __asm__(
	  "add %1, %0\n\t"
	  "rcr $1, %0\n\t"
	  "shr $1, %1\n\t"
	  "cmovnc %1, %0\n"
	  : "+&r" (t), "+&r" (s)
	  : : "cc"
	  );
  r[0] = t;
#else
  ASSERT_EXPENSIVE (m[0] % 2UL != 0UL);
  if (a[0] % 2UL == 0UL)
    r[0] = a[0] / 2UL;
  else
    r[0] = a[0] / 2UL + m[0] / 2UL + 1UL;
#endif
}


MAYBE_UNUSED
static inline int
modul_next (residueul_t r, const modulusul_t m)
{
    return (++r[0] == m[0]);
}


MAYBE_UNUSED
static inline int
modul_finished (const residueul_t r, const modulusul_t m)
{
    return (r[0] == m[0]);
}


/* prototypes of non-inline functions */
int modul_div3 (residueul_t, const residueul_t, const modulusul_t);
int modul_div5 (residueul_t, const residueul_t, const modulusul_t);
int modul_div7 (residueul_t, const residueul_t, const modulusul_t);
int modul_div11 (residueul_t, const residueul_t, const modulusul_t);
int modul_div13 (residueul_t, const residueul_t, const modulusul_t);
void modul_gcd (modintul_t, const residueul_t, const modulusul_t);
void modul_pow_ul (residueul_t, const residueul_t, const unsigned long, 
		   const modulusul_t);
void modul_2pow_ul (residueul_t, const unsigned long, const modulusul_t);
void modul_pow_mp (residueul_t, const residueul_t, const unsigned long *,
                   const int, const modulusul_t);
void modul_2pow_mp (residueul_t, const unsigned long *, const int, 
                    const modulusul_t);
void modul_V_ul (residueul_t, const residueul_t, const unsigned long, 
                 const modulusul_t);
void modul_V_mp (residueul_t, const residueul_t, const unsigned long *,
		 const int, const modulusul_t);
int modul_sprp (const residueul_t, const modulusul_t);
int modul_sprp2 (const modulusul_t);
int modul_isprime (const modulusul_t);
int modul_inv (residueul_t, const residueul_t, const modulusul_t);
int modul_jacobi (const residueul_t, const modulusul_t);

/* Cruft: this belongs in modredc_*.h, not here  */
unsigned long modul_invmodlong (unsigned long n);

#endif  /* MOD_UL_H */
