/* Some functions for modular arithmetic with residues and modulus in
   unsigned long variables. Residues are stored in Montgomery form,
   reduction after multiplication is done with REDC. Due to inlining, 
   this file must be included in the caller's source code with #include */

/* Naming convention: all function start with modredcul, for 
   MODulus REDC Unsigned Long, followed by underscore, functionality of 
   function (add, mul, etc), and possibly underscore and specification of 
   what argument types the function takes (_ul, etc). */

#ifndef __MODREDC_UL_H

#define __MODREDC_UL_H

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
/* Helper macros */

/* A macro for function renaming. All functions here start with modul_ */
#define MODREDCUL_RENAME(x) modredcul_##x

#define MODREDCUL_SIZE 1
#define MODREDCUL_MAXBITS LONG_BIT

typedef unsigned long residueredcul_t[MODREDCUL_SIZE];
typedef unsigned long modintredcul_t[MODREDCUL_SIZE];
typedef struct { 
  unsigned long m;
  unsigned long invm;
  residueredcul_t one;
} __modulusredcul_t;
typedef __modulusredcul_t modulusredcul_t[1];


/* ==================== Functions used internally ==================== */

static inline void
modredcul_add (residueredcul_t, const residueredcul_t, 
               const residueredcul_t, const modulusredcul_t);
/* Computes (a * 2^wordsize) % m */
MAYBE_UNUSED
static inline void
modredcul_tomontgomery (residueredcul_t r, const residueredcul_t a, 
                        const modulusredcul_t m)
{
  ASSERT_EXPENSIVE (a[0] < m[0].m);
  ularith_div_2ul_ul_ul_r (r, 0UL, a[0], m[0].m);
}


/* Computes (a / 2^wordsize) % m. Assumes a < m */
MAYBE_UNUSED
static inline void
modredcul_frommontgomery (residueredcul_t r, const residueredcul_t a, 
                          const modulusredcul_t m)
{
  unsigned long tlow, thigh;
  tlow = a[0] * m[0].invm;
  ularith_mul_ul_ul_2ul (&tlow, &thigh, tlow, m[0].m);
  r[0] = thigh + ((a[0] != 0UL) ? 1UL : 0UL);
}

/* Computes ((phigh*2^wordsize + plow) / 2^wordsize) % m.
   Requires phigh < m */
MAYBE_UNUSED
static inline void
modredcul_redc (residueredcul_t r, const unsigned long plow, 
                const unsigned long phigh, const modulusredcul_t m)
{
  unsigned long tlow, thigh, t = phigh;

  ASSERT_EXPENSIVE (phigh < m[0].m);

  tlow = plow * m[0].invm;
  ularith_mul_ul_ul_2ul (&tlow, &thigh, tlow, m[0].m);
  /* Let w = 2^wordsize. We know (phigh * w + plow) + (thigh * w + tlow) 
     == 0 (mod w) so either plow == tlow == 0, or plow !=0 and tlow != 0. 
     In the former case we want phigh + thigh + 1, in the latter 
     phigh + thigh. Since t = phigh < m, and modredcul_add can handle the 
     case where the second operand is equal to m, adding 1 is safe */

  t += (plow != 0UL) ? 1UL : 0UL; /* Does not depend on the mul */ 

  modredcul_add (r, &t, &thigh, m);
}


/* Requires a < m, then result == a+b (mod m).
   If a + b < m, then r = a + b, otherwise r = a + b - m.
   Implies r < b if b >= m. */
MAYBE_UNUSED
static inline void
modredcul_add_semi (residueredcul_t r, const residueredcul_t a, 
		    const residueredcul_t b, const modulusredcul_t m)
{
  ASSERT_EXPENSIVE (a[0] < m[0].m);

#if (defined(__i386__) || defined(__x86_64__)) && defined(__GNUC__)
  {
    unsigned long t = a[0] - m[0].m, tr = a[0] + b[0];
    
    __asm__ (
      "add %2, %1\n\t"   /* t = t + b ( t == a - m + b (mod w)) */
      "cmovc %1, %0\n\t"  /* if (cy) tr = t */
      : "+r" (tr), "+&r" (t)
      : ULARITH_CONSTRAINT_G (b[0])
      : "cc"
    );
    r[0] = tr;
  }
#else
  r[0] = (b[0] >= m[0].m - a[0]) ? (b[0] - (m[0].m - a[0])) : (a[0] + b[0]);
#endif
}


/* Requires phigh < ULONG_MAX. If phigh:plow = a*b with a,b < ULONG_MAX,
   phigh <= ULONG_MAX - 1, so that works. 
   If phigh < m, then r < m, otherwise r <= phigh. */
MAYBE_UNUSED
static inline void
modredcul_redc_semi (residueredcul_t r, const unsigned long plow, 
		     const unsigned long phigh, const modulusredcul_t m)
{
  unsigned long tlow, thigh, t = phigh;

  ASSERT_EXPENSIVE (phigh < ULONG_MAX);

  tlow = plow * m[0].invm;
  ularith_mul_ul_ul_2ul (&tlow, &thigh, tlow, m[0].m);
  t += (plow != 0UL) ? 1UL : 0UL;
  /* We have thigh < m since thigh = trunc (m * (something % w) / w), 
     so this add always produces the correct residue class (although 
     it may have r >= m if phigh >= m. This add is slightly slower
     than in the modredcul_redc() function, since thigh depends on the 
     previous mul and is used first in the add. */
  modredcul_add_semi (r, &thigh, &t, m);
}


/* ================= Functions that are part of the API ================= */

/* Some functions for integers of the same width as the modulus */

MAYBE_UNUSED
static inline void
modredcul_intset (modintredcul_t r, const modintredcul_t s)
{
  r[0] = s[0];
}


MAYBE_UNUSED
static inline void
modredcul_intset_ul (modintredcul_t r, const unsigned long s)
{
  r[0] = s;
}


MAYBE_UNUSED
static inline int
modredcul_intequal (const modintredcul_t a, const modintredcul_t b)
{
  return (a[0] == b[0]);
}


MAYBE_UNUSED
static inline int
modredcul_intequal_ul (const modintredcul_t a, const unsigned long b)
{
  return (a[0] == b);
}


MAYBE_UNUSED
static inline int
modredcul_intcmp (const modintredcul_t a, const modintredcul_t b)
{
  return (a[0] < b[0]) ? -1 : (a[0] == b[0]) ? 0 : 1;
}


MAYBE_UNUSED
static inline int
modredcul_intcmp_ul (const modintredcul_t a, const unsigned long b)
{
  return (a[0] < b) ? -1 : (a[0] == b) ? 0 : 1;
}

MAYBE_UNUSED
static inline int
modredcul_intfits_ul (const modintredcul_t a MAYBE_UNUSED)
{
  return 1;
}

MAYBE_UNUSED
static inline void
modredcul_intadd (modintredcul_t r, const modintredcul_t a, 
                  const modintredcul_t b)
{
  r[0] = a[0] + b[0];
}

MAYBE_UNUSED
static inline void
modredcul_intsub (modintredcul_t r, const modintredcul_t a, 
                  const modintredcul_t b)
{
  r[0] = a[0] - b[0];
}


/* Returns the number of bits in a, that is, floor(log_2(a))+1. 
   For a == 0 returns 0. */
MAYBE_UNUSED
static inline int
modredcul_intbits (const modintredcul_t a)
{
  if (a[0] == 0)
    return 0;
  return LONG_BIT - ularith_clz (a[0]);
}

MAYBE_UNUSED
static inline void
modredcul_intshr (modintredcul_t r, const modintredcul_t s, const int i)
{
  r[0] = s[0] >> i;
}


MAYBE_UNUSED
static inline void
modredcul_intshl (modintredcul_t r, const modintredcul_t s, const int i)
{
  r[0] = s[0] << i;
}


/* r = n/d. We require d|n */
MAYBE_UNUSED
static inline void
modredcul_intdivexact (modintredcul_t r, const modintredcul_t n, 
                       const modintredcul_t d)
{
  /* ularith_invmod() is faster than a DIV */
  r[0] = n[0] * ularith_invmod(d[0]); 
}


/* r = n%d */
MAYBE_UNUSED
static inline void
modredcul_intmod (modintredcul_t r, const modintredcul_t n, 
                  const modintredcul_t d)
{
  r[0] = n[0] % d[0]; 
}


/* Functions for the modulus */

MAYBE_UNUSED
static inline void
modredcul_initmod_ul (modulusredcul_t m, const unsigned long s)
{
  m[0].m = s;
  m[0].invm = -ularith_invmod (s);
  if (m[0].m == 1UL)
    m[0].one[0] = 0UL;
  else
    {
      m[0].one[0] = 1UL;
      modredcul_tomontgomery (m[0].one, m[0].one, m);
    }
}


MAYBE_UNUSED
static inline void
modredcul_initmod_uls (modulusredcul_t m, const modintredcul_t s)
{
  m[0].m = s[0];
  m[0].invm = -ularith_invmod (s[0]);
  if (m[0].m == 1UL)
    m[0].one[0] = 0UL;
  else
    {
      m[0].one[0] = 1UL;
      modredcul_tomontgomery (m[0].one, m[0].one, m);
    }
}


MAYBE_UNUSED
static inline unsigned long
modredcul_getmod_ul (const modulusredcul_t m)
{
  return m[0].m;
}


MAYBE_UNUSED
static inline void
modredcul_getmod_uls (modintredcul_t r, const modulusredcul_t m)
{
  r[0] = m[0].m;
}


MAYBE_UNUSED
static inline void
modredcul_clearmod (modulusredcul_t m MAYBE_UNUSED)
{
  return;
}


/* Functions for residues */

/* Initialises a residue_t type and sets it to zero */
MAYBE_UNUSED
static inline void
modredcul_init (residueredcul_t r, const modulusredcul_t m MAYBE_UNUSED)
{
  r[0] = 0UL;
}


/* Initialises a residue_t type, but does not set it to zero. For fixed length
   residue_t types, that leaves nothing to do at all. */
MAYBE_UNUSED
static inline void
modredcul_init_noset0 (residueredcul_t r MAYBE_UNUSED, 
                       const modulusredcul_t m MAYBE_UNUSED)
{
  return;
}


MAYBE_UNUSED
static inline void
modredcul_clear (residueredcul_t r MAYBE_UNUSED, 
                 const modulusredcul_t m MAYBE_UNUSED)
{
  return;
}


MAYBE_UNUSED
static inline void
modredcul_set (residueredcul_t r, const residueredcul_t s, 
               const modulusredcul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (s[0] < m[0].m);
  r[0] = s[0];
}


MAYBE_UNUSED
static inline void
modredcul_set_ul (residueredcul_t r, const unsigned long s, 
                  const modulusredcul_t m)
{
  unsigned long plow, phigh;
  ularith_mul_ul_ul_2ul (&plow, &phigh, s, m[0].one[0]);
  modredcul_redc (r, plow, phigh, m);
  modredcul_tomontgomery (r, r, m);
}


/* Sets the residue_t to the class represented by the integer s. Assumes that
   s is reduced (mod m), i.e. 0 <= s < m */

MAYBE_UNUSED
static inline void
modredcul_set_ul_reduced (residueredcul_t r, const unsigned long s, 
                          const modulusredcul_t m MAYBE_UNUSED)
{
  ASSERT (s < m[0].m);
  r[0] = s;
  modredcul_tomontgomery (r, r, m);
}


MAYBE_UNUSED
static inline void
modredcul_set_uls (residueredcul_t r, const modintredcul_t s, 
		   const modulusredcul_t m)
{
  r[0] = s[0] % m[0].m;
  modredcul_tomontgomery (r, r, m);
}


MAYBE_UNUSED
static inline void
modredcul_set_uls_reduced (residueredcul_t r, const modintredcul_t s, 
			   const modulusredcul_t m)
{
  ASSERT (s[0] < m[0].m);
  r[0] = s[0];
  modredcul_tomontgomery (r, r, m);
}


/* This one is so trivial that we don't really require m in the
 * interface. For interface homogeneity we make it take the m parameter 
 * anyway.
 */
MAYBE_UNUSED 
static inline void 
modredcul_set0 (residueredcul_t r, const modulusredcul_t m MAYBE_UNUSED) 
{ 
  r[0] = 0UL;
}


MAYBE_UNUSED 
static inline void 
modredcul_set1 (residueredcul_t r, const modulusredcul_t m) 
{
  r[0] = m[0].one[0];
}


/* Exchanges the values of the two arguments */

MAYBE_UNUSED
static inline void
modredcul_swap (residueredcul_t a, residueredcul_t b, 
                const modulusredcul_t m MAYBE_UNUSED)
{
  unsigned long t;
  ASSERT_EXPENSIVE (a[0] < m[0].m && b[0] < m[0].m);
  t = a[0];
  a[0] = b[0];
  b[0] = t;
}


MAYBE_UNUSED
static inline unsigned long
modredcul_get_ul (const residueredcul_t s, 
	          const modulusredcul_t m MAYBE_UNUSED)
{
  unsigned long t;
  ASSERT_EXPENSIVE (s[0] < m[0].m);
  modredcul_frommontgomery (&t, s, m);
  return t;
}


MAYBE_UNUSED
static inline void
modredcul_get_uls (modintredcul_t r, const residueredcul_t s, 
		   const modulusredcul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (s[0] < m[0].m);
  modredcul_frommontgomery (r, s, m);
}


MAYBE_UNUSED
static inline int
modredcul_equal (const residueredcul_t a, const residueredcul_t b, 
             const modulusredcul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (a[0] < m[0].m && b[0] < m[0].m);
  return (a[0] == b[0]);
}


MAYBE_UNUSED
static inline int
modredcul_is0 (const residueredcul_t a, const modulusredcul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (a[0] < m[0].m);
  return (a[0] == 0UL);
}


MAYBE_UNUSED
static inline int
modredcul_is1 (const residueredcul_t a, const modulusredcul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (a[0] < m[0].m);
  return (a[0] == m[0].one[0]);
}


/* Requires a < m and b <= m, then r == a+b (mod m) and r < m */
MAYBE_UNUSED
static inline void
modredcul_add (residueredcul_t r, const residueredcul_t a, 
               const residueredcul_t b, const modulusredcul_t m)
{
  ASSERT_EXPENSIVE (a[0] < m[0].m && b[0] <= m[0].m);
#ifdef MODTRACE
  printf ("modul_add: a = %lu, b = %lu", a[0], b[0]);
#endif

#if (defined(__i386__) || defined(__x86_64__)) && defined(__GNUC__)
  {
    unsigned long t = a[0] + b[0], tr = a[0] - m[0].m;
  
    __asm__ (
      "add %2, %0\n\t"   /* tr += b */
      "cmovnc %1, %0\n\t"  /* if (!cy) tr = t */
      : "+&r" (tr)
      : "rm" (t), ULARITH_CONSTRAINT_G (b[0])
      : "cc"
    );
    ASSERT_EXPENSIVE (tr == ((a[0] >= m[0].m - b[0]) ? (a[0] - (m[0].m - b[0])) : (a[0] + b[0])));
    r[0] = tr;
  }
#else
  r[0] = (b[0] >= m[0].m - a[0]) ? (b[0] - (m[0].m - a[0])) : (a[0] + b[0]);
#endif

#ifdef MODTRACE
  printf (", r = %lu\n", r[0]);
#endif
  ASSERT_EXPENSIVE (r[0] < m[0].m);
}


MAYBE_UNUSED
static inline void
modredcul_add_ul (residueredcul_t r, const residueredcul_t a,
                  const unsigned long b, const modulusredcul_t m)
{
  residueredcul_t t;

  /* TODO: speed up */
  ASSERT_EXPENSIVE (a[0] < m[0].m);
  modredcul_init_noset0 (t, m);
  modredcul_set_ul (t, b, m);
  modredcul_add (r, a, t, m);
  modredcul_clear (t, m);
}


MAYBE_UNUSED
static inline void
modredcul_sub (residueredcul_t r, const residueredcul_t a, 
               const residueredcul_t b, const modulusredcul_t m)
{
  ASSERT_EXPENSIVE (a[0] < m[0].m && b[0] < m[0].m);
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
      : ULARITH_CONSTRAINT_G (b[0]), "r" (m[0].m)
      : "cc"
    );
    r[0] = tr;
  }
#elif 1
  /* Seems to be faster than the one below */
  {
    unsigned long t = 0UL, tr;
    if ((tr = a[0] - b[0]) > a[0])
      t = m[0].m;
    r[0] = tr + t;
  }
#else
  r[0] = (a[0] < b[0]) ? (a[0] - b[0] + m[0].m) : (a[0] - b[0]);
#endif

#ifdef MODTRACE
  printf (", r = %lu\n", r[0]);
#endif
  ASSERT_EXPENSIVE (r[0] < m[0].m);
}


MAYBE_UNUSED
static inline void
modredcul_sub_ul (residueredcul_t r, const residueredcul_t a,
                  const unsigned long b, const modulusredcul_t m)
{
  residueredcul_t t;

  /* TODO: speed up */
  ASSERT_EXPENSIVE (a[0] < m[0].m);
  modredcul_init_noset0 (t, m);
  modredcul_set_ul (t, b, m);
  modredcul_sub (r, a, t, m);
  modredcul_clear (t, m);
}


MAYBE_UNUSED
static inline void
modredcul_neg (residueredcul_t r, const residueredcul_t a, 
               const modulusredcul_t m)
{
  ASSERT_EXPENSIVE (a[0] < m[0].m);
  if (a[0] == 0UL)
    r[0] = a[0];
  else
    r[0] = m[0].m - a[0];
}


MAYBE_UNUSED
static inline void
modredcul_mul (residueredcul_t r, const residueredcul_t a, 
               const residueredcul_t b, const modulusredcul_t m)
{
  unsigned long plow, phigh;

  ASSERT_EXPENSIVE (m[0].m % 2 != 0);
  ASSERT_EXPENSIVE (a[0] < m[0].m && b[0] < m[0].m);
#if defined(MODTRACE)
  printf ("(%lu * %lu / 2^%ld) %% %lu", a[0], b[0], LONG_BIT, m[0].m);
#endif
  
  ularith_mul_ul_ul_2ul (&plow, &phigh, a[0], b[0]);

#if defined(__x86_64__) && defined(__GNUC__)

  /* TODO: are the register constraints watertight? 
     %rax gets modified but putting tlow as an output constraint with "+"
     will keep r from getting allocated in %rax, which is a shame
     since we often want the result in %rax for the next multiply. */
  
  __asm__ (
    "imulq %[invm], %%rax\n\t"
    "cmpq $1, %%rax \n\t"                /* if plow != 0, increase phigh */
    "sbbq $-1, %[phigh]\n\t"
    "mulq %[m]\n\t"
    "lea (%[phigh],%%rdx,1), %[r]\n\t"  /* compute (rdx + thigh) mod m */
    "subq %[m], %[phigh]\n\t"
    "addq %%rdx, %[phigh]\n\t"
    "cmovcq %[phigh], %[r]\n\t"
    : [phigh] "+&r" (phigh), [r] "=&r" (r[0])
    : [invm] "rm" (m[0].invm), [m] "rm" (m[0].m), "a" (plow)
    : "%rdx", "cc"
  );
#else
  modredcul_redc (r, plow, phigh, m);
#endif

#if defined(MODTRACE)
  printf (" == %lu /* PARI */ \n", r[0]);
#endif
}


MAYBE_UNUSED
static inline void
modredcul_sqr (residueredcul_t r, const residueredcul_t a, 
               const modulusredcul_t m)
{
  unsigned long plow, phigh;

  ASSERT_EXPENSIVE (m[0].m % 2 != 0);
  ASSERT_EXPENSIVE (a[0] < m[0].m);

  ularith_sqr_ul_2ul (&plow, &phigh, a[0]);
#if defined(__x86_64__) && defined(__GNUC__)

  /* TODO: are the register constraints watertight? 
     %rax gets modified but putting tlow as an output constraint with "+"
     will keep r from getting allocated in %rax, which is a shame
     since we often want the result in %rax for the next multiply. */
  
  __asm__ (
    "imulq %[invm], %%rax\n\t"
    "cmpq $1, %%rax \n\t"                /* if plow != 0, increase phigh */
    "sbbq $-1, %[phigh]\n\t"
    "mulq %[m]\n\t"
    "lea (%[phigh],%%rdx,1), %[r]\n\t"  /* compute (rdx + thigh) mod m */
    "subq %[m], %[phigh]\n\t"
    "addq %%rdx, %[phigh]\n\t"
    "cmovcq %[phigh], %[r]\n\t"
    : [phigh] "+&r" (phigh), [r] "=&r" (r[0])
    : [invm] "rm" (m[0].invm), [m] "rm" (m[0].m), "a" (plow)
    : "%rdx", "cc"
  );
#else
  modredcul_redc (r, plow, phigh, m);
#endif
}


MAYBE_UNUSED
static inline void
modredcul_div2 (residueredcul_t r, const residueredcul_t a, 
                const modulusredcul_t m)
{
#if (defined(__i386__) || defined(__x86_64__)) && defined(__GNUC__)
  unsigned long s = a[0], t = m[0].m;
  ASSERT_EXPENSIVE (m[0].m % 2UL != 0UL);

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
  ASSERT_EXPENSIVE (m[0].m % 2UL != 0UL);
  if (a[0] % 2UL == 0UL)
    r[0] = a[0] / 2UL;
  else
    r[0] = a[0] / 2UL + m[0].m / 2UL + 1UL;
#endif
}


MAYBE_UNUSED
static inline int
modredcul_next (residueredcul_t r, const modulusredcul_t m)
{
    return (++r[0] == m[0].m);
}


MAYBE_UNUSED
static inline int
modredcul_finished (const residueredcul_t r, const modulusredcul_t m)
{
    return (r[0] == m[0].m);
}


/* prototypes of non-inline functions */
int modredcul_div3 (residueredcul_t, const residueredcul_t, 
                     const modulusredcul_t);
int modredcul_div5 (residueredcul_t, const residueredcul_t, 
                     const modulusredcul_t);
int modredcul_div7 (residueredcul_t, const residueredcul_t, 
                     const modulusredcul_t);
int modredcul_div11 (residueredcul_t, const residueredcul_t, 
                      const modulusredcul_t);
int modredcul_div13 (residueredcul_t, const residueredcul_t, 
                      const modulusredcul_t);
void modredcul_gcd (modintredcul_t, const residueredcul_t, 
                    const modulusredcul_t);
void modredcul_pow_ul (residueredcul_t, const residueredcul_t, 
                       const unsigned long, const modulusredcul_t);
void modredcul_2pow_ul (residueredcul_t, const unsigned long, 
                        const modulusredcul_t);
void modredcul_pow_mp (residueredcul_t, const residueredcul_t, 
                   const unsigned long *, const int, const modulusredcul_t);
void modredcul_2pow_mp (residueredcul_t, const unsigned long *, const int, 
                        const modulusredcul_t);
void modredcul_V_ul (residueredcul_t, const residueredcul_t, 
		     const unsigned long, const modulusredcul_t);
void modredcul_V_mp (residueredcul_t, const residueredcul_t, 
		     const unsigned long *, const int, const modulusredcul_t);
int modredcul_sprp (const residueredcul_t, const modulusredcul_t);
int modredcul_sprp2 (const modulusredcul_t);
int modredcul_isprime (const modulusredcul_t);
int modredcul_inv (residueredcul_t, const residueredcul_t, 
		   const modulusredcul_t);
int modredcul_jacobi (const residueredcul_t, const modulusredcul_t);

#endif
