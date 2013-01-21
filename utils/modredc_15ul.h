/* Some functions for modular arithmetic with residues and modulus in
   unsigned long variables. The modulus can be up to 1.5 unsigned longs 
   in size (meaning 48 bits if unsigned long has 32 bits, or 96 bits if 
   unsigned long has 64 bits). Residues are stored in Montgomery form,
   reduction after multiplication is done with REDC. Due to inlining, 
   this file must be included in the caller's source code with #include */

/* Naming convention: all function start with modredc15ul, for 
   MODulus REDC 1.5 Unsigned Longs, followed by underscore, functionality of 
   function (add, mul, etc), and possibly underscore and specification of 
   what argument types the function takes (_ul, etc). */

#ifndef MODREDC_15UL_H
#define MODREDC_15UL_H

/**********************************************************************/
#include <assert.h>
#if defined(MODTRACE)
#include <stdio.h>
#endif
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

/* A macro for function renaming. All functions here start with 
   modredc15ul_ */
#define MODREDC15UL_RENAME(x) modredc15ul_##x

#define MODREDC15UL_SIZE 2
#define MODREDC15UL_MINBITS LONG_BIT
#define MODREDC15UL_MAXBITS (LONG_BIT + LONG_BIT/2)

typedef unsigned long residueredc15ul_t[MODREDC15UL_SIZE];
typedef unsigned long modintredc15ul_t[MODREDC15UL_SIZE];
typedef struct { 
  modintredc15ul_t m;
  residueredc15ul_t one;
  unsigned long invm;
} __modulusredc15ul_t;
typedef __modulusredc15ul_t modulusredc15ul_t[1];


/* ==================== Functions used internally ==================== */

static inline void
modredc15ul_add (residueredc15ul_t r, const residueredc15ul_t a, 
		 const residueredc15ul_t b, const modulusredc15ul_t m);
static inline void
modredc15ul_get_uls (modintredc15ul_t r, const residueredc15ul_t s, 
		     const modulusredc15ul_t m MAYBE_UNUSED);

MAYBE_UNUSED
static inline void
modredc15ul_tomontgomery (residueredc15ul_t r, const residueredc15ul_t s,
			  const modulusredc15ul_t m)
{
  int i;

  r[0] = s[0];
  r[1] = s[1];
  /* TODO FIXME: ridiculously slow */
  for (i = 0; i < 2 * LONG_BIT; i++)
    modredc15ul_add (r, r, r, m);
}


MAYBE_UNUSED
static inline void
modredc15ul_frommontgomery (residueredc15ul_t r, const residueredc15ul_t s,
			    const modulusredc15ul_t m)
{
  unsigned long t[5], k;
  
  /* Do two REDC steps */
  k = s[0] * m[0].invm;
  ularith_mul_ul_ul_2ul (&(t[0]), &(t[1]), k, m[0].m[0]);
  if (t[0] != 0UL)
    t[1]++;
  t[2] = 0;
  ularith_add_ul_2ul (&(t[1]), &(t[2]), s[1]); /* t[2] <= 1 */
  ularith_mul_ul_ul_2ul (&(t[0]), &(t[4]), k, m[0].m[1]); /* t[4] < 2^(w/2) */
  ularith_add_2ul_2ul (&(t[1]), &(t[2]), t[0], t[4]);
  /* Here, t2:t1:0 = s + k*m */
  k = t[1] * m[0].invm;
  ularith_mul_ul_ul_2ul (&(t[0]), &(t[4]), k, m[0].m[0]);
  if (t[0] != 0UL)
    t[4]++;
  t[3] = 0UL;
  ularith_add_ul_2ul (&(t[2]), &(t[3]), t[4]); /* t[3] <= 1 */
  ularith_mul_ul_ul_2ul (&(t[0]), &(t[4]), k, m[0].m[1]);
  ularith_add_2ul_2ul (&(t[2]), &(t[3]), t[0], t[4]);
  r[0] = t[2];
  r[1] = t[3];
}

/* Do a one-word REDC, i.e., divide by 2^LONG_BIT */
MAYBE_UNUSED
static inline void
modredc15ul_redc1 (residueredc15ul_t r, const residueredc15ul_t s,
		   const modulusredc15ul_t m)
{
  unsigned long t[4], k;
  
  k = s[0] * m[0].invm;
  ularith_mul_ul_ul_2ul (&(t[0]), &(t[1]), k, m[0].m[0]);
  if (t[0] != 0UL)
    t[1]++;
  t[2] = 0;
  ularith_add_ul_2ul (&(t[1]), &(t[2]), s[1]); /* t[2] <= 1 */
  ularith_mul_ul_ul_2ul (&(t[0]), &(t[3]), k, m[0].m[1]); /* t[3] < 2^(w/2) */
  ularith_add_2ul_2ul (&(t[1]), &(t[2]), t[0], t[3]);

  /* r = (k*m + s) div wb. k <= wb-1, s<m, thus r < m */
  r[0] = t[1];
  r[1] = t[2];
}

/* ================= Functions that are part of the API ================= */

/* Some functions for integers of the same width as the modulus */

MAYBE_UNUSED
static inline void
modredc15ul_intset (modintredc15ul_t r, const modintredc15ul_t s)
{
  r[0] = s[0];
  r[1] = s[1];
}

MAYBE_UNUSED
static inline void
modredc15ul_intset_ul (modintredc15ul_t r, const unsigned long s)
{
  r[0] = s;
  r[1] = 0UL;
}

MAYBE_UNUSED
static inline int
modredc15ul_intequal (const modintredc15ul_t a, const modintredc15ul_t b)
{
  return (a[0] == b[0] && a[1] == b[1]);
}

MAYBE_UNUSED
static inline int
modredc15ul_intequal_ul (const modintredc15ul_t a, const unsigned long b)
{
  return (a[0] == b && a[1] == 0UL);
}

/* Returns 1 if a < b, 0 otherwise */
MAYBE_UNUSED
static inline int
modredc15ul_intlt (const modintredc15ul_t a, const modintredc15ul_t b)
{
    modintredc15ul_t t;

    modredc15ul_intset (t, a);
    return ularith_sub_2ul_2ul_cy (&(t[0]), &(t[1]), b[0], b[1]);
}

MAYBE_UNUSED
static inline int
modredc15ul_intcmp (const modintredc15ul_t a, const modintredc15ul_t b)
{
  if (a[1] < b[1])
    return -1;
  if (a[1] > b[1])
    return 1;
  return (a[0] < b[0]) ? -1 : (a[0] == b[0]) ? 0 : 1;
}

MAYBE_UNUSED
static inline int
modredc15ul_intcmp_ul (const modintredc15ul_t a, const unsigned long b)
{
  if (a[1] > 0UL)
    return 1;
  return (a[0] < b) ? -1 : (a[0] == b) ? 0 : 1;
}

MAYBE_UNUSED
static inline int
modredc15ul_intfits_ul (const modintredc15ul_t a)
{
  return (a[1] == 0UL);
}

MAYBE_UNUSED
static inline void
modredc15ul_intadd (modintredc15ul_t r, const modintredc15ul_t a,
		    const modintredc15ul_t b)
{
  modintredc15ul_t t;
  modredc15ul_intset (t, a);
  ularith_add_2ul_2ul (&(t[0]), &(t[1]), b[0], b[1]);
  modredc15ul_intset (r, t);
}

MAYBE_UNUSED
static inline void
modredc15ul_intsub (modintredc15ul_t r, const modintredc15ul_t a,
		    const modintredc15ul_t b)
{
  modintredc15ul_t t;
  modredc15ul_intset (t, a);
  ularith_sub_2ul_2ul (&(t[0]), &(t[1]), b[0], b[1]);
  modredc15ul_intset (r, t);
}

/* Returns the number of bits in a, that is, floor(log_2(a))+1. 
   For a == 0 returns 0. */
MAYBE_UNUSED
static inline int
modredc15ul_intbits (const modintredc15ul_t a)
{
  if (a[1] > 0UL)
    return 2*LONG_BIT - ularith_clz (a[1]);

  if (a[0] > 0UL)
    return LONG_BIT - ularith_clz (a[0]);
  
  return 0;
}


MAYBE_UNUSED
static inline void
modredc15ul_intshr (modintredc15ul_t r, const modintredc15ul_t s, const int i)
{
#if (__GNUC__ == 4 && __GNUC_MINOR__ <= 2)
  /* gcc 4.2.1 on 32-bit Intel CPUs seems to get confused over register 
     allocation when s == r, when using the shrd instruction via inline
     assembly in ularith_shrd(), and when modredc15ul_intshr() itself gets 
     inlined. This extra variable seems to fix it. */
  unsigned long t = s[1];
  r[0] = s[0];
  ularith_shrd (&(r[0]), t, i);
  r[1] = t >> i;
#else
  r[0] = s[0];
  ularith_shrd (&(r[0]), s[1], i);
  r[1] = s[1] >> i;
#endif
}


MAYBE_UNUSED
static inline void
modredc15ul_intshl (modintredc15ul_t r, const modintredc15ul_t s, const int i)
{
  r[1] = s[1];
  ularith_shld (&(r[1]), s[0], i);
  r[0] = s[0] << i;
}


/* r = n/d. We require d|n */
MAYBE_UNUSED
static inline void
modredc15ul_intdivexact (modintredc15ul_t r, const modintredc15ul_t n,
                         const modintredc15ul_t d)
{
  modintredc15ul_t n1, d1;
  unsigned long invf, r0, k0, k1;
#ifdef WANT_ASSERT_EXPENSIVE
  unsigned long s0 = n[0], s1 = n[1];
#endif
  
  n1[0] = n[0];
  n1[1] = n[1];
  d1[0] = d[0];
  d1[1] = d[1];

  /* Make d odd. TODO: use ctzl */
  while (d1[0] % 2 == 0UL)
    {
      ASSERT_EXPENSIVE (n1[0] % 2 == 0UL);
      ularith_shrd (&(d1[0]), d1[1], 1);
      d1[1] >>= 1;
      ularith_shrd (&(n1[0]), n1[1], 1);
      n1[1] >>= 1;
    }
  
  invf = ularith_invmod (d1[0]);
  r0 = invf * n1[0];
  ularith_mul_ul_ul_2ul (&k0, &k1, r0, d1[0]);
  ularith_sub_2ul_2ul (&(n1[0]), &(n1[1]), k0, k1);
  n1[1] -= r0 * d1[1];
  ASSERT (n1[0] == 0UL);
  r[0] = r0;
  r[1] = invf * n1[1];
#ifdef WANT_ASSERT_EXPENSIVE
  ASSERT_EXPENSIVE (d[1] == 0UL || r[1] == 0UL);
  ularith_mul_ul_ul_2ul (&k0, &k1, r[0], d[0]);
  ularith_sub_2ul_2ul (&s0, &s1, k0, k1);
  ASSERT_EXPENSIVE (s0 == 0UL);
  ularith_mul_ul_ul_2ul (&k0, &k1, r[0], d[1]);
  ASSERT_EXPENSIVE (k1 == 0UL);
  s1 -= k0;
  ularith_mul_ul_ul_2ul (&k0, &k1, r[1], d[0]);
  ASSERT_EXPENSIVE (k1 == 0UL);
  s1 -= k0;
  ASSERT_EXPENSIVE (s1 == 0UL);
#endif
}

/* r = n%d */
MAYBE_UNUSED
static inline void
modredc15ul_intmod (modintredc15ul_t r, const modintredc15ul_t n,
                    const modintredc15ul_t d)
{
  if (d[1] == 0UL)
    {
      if (n[1] < d[0])
        {
          unsigned long dummy;
          ularith_div_2ul_ul_ul (&dummy, &(r[0]), n[0], n[1], d[0]);
          r[1] = 0UL;
        }
      else
        {
          unsigned long t, dummy;
          ularith_div_2ul_ul_ul (&dummy, &t, n[1], 0UL, d[0]);
          ularith_div_2ul_ul_ul (&dummy, &(r[0]), n[0], t, d[0]);
          r[1] = 0UL;
        }
    }
  else    
    {
      modintredc15ul_t t1, t2;
      unsigned long q, dummy;
      int i;

      /* Divide both by 2^k s.t. 2^(LONG_BIT-1) <= d/2^k < 2^LONG_BIT */
      modredc15ul_intshr (t1, n, 1);
      modredc15ul_intshr (t2, d, 1);
      if (t2[1] != 0UL)
        {
          i = LONG_BIT - ularith_clz (t2[1]);
          modredc15ul_intshr (t1, t1, i);
          modredc15ul_intshr (t2, t2, i);
        }
      ASSERT(t2[1] == 0);
      ASSERT((t2[0] & (1UL << (LONG_BIT - 1))) != 0UL);
      ASSERT (t1[1] < t2[0]);
      
      ularith_div_2ul_ul_ul (&q, &dummy, t1[0], t1[1], t2[0]);
      ularith_mul_ul_ul_2ul (&(t1[0]), &(t1[1]), q, d[0]);
      t1[1] += q * d[1];
      /* printf ("n=%lu*2^%d + %lu; d=%lu*2^%d + %lu; q=%lu; "
              "t1=%lu*2^%d + %lu; ", 
              n[1], LONG_BIT, n[0], d[1], LONG_BIT, d[0], q,
              t1[1], LONG_BIT, t1[0]); */
      
      if (modredc15ul_intlt (n, t1))
        modredc15ul_intsub (t1, t1, d);
      ASSERT(!modredc15ul_intlt (n, t1));
      modredc15ul_intsub(r, n, t1);
      /* printf ("r=%lu*2^%d + %lu;\n", r[1], LONG_BIT, r[0]); */
    }
}


/* Functions for the modulus */

/* Init the modulus from a multi-word integer. s must point to an array of
   at least two unsigned longs, where s[0] is the low word of the modulus, 
   and s[1] is the high word. */
MAYBE_UNUSED
static inline void
modredc15ul_initmod_uls (modulusredc15ul_t m, const modintredc15ul_t s)
{
  ASSERT (s[1] > 0UL);
  ASSERT (s[1] < (1UL << (LONG_BIT / 2)));
  modredc15ul_intset (m[0].m, s);
  m[0].invm = -ularith_invmod (s[0]);

  modredc15ul_intset_ul (m[0].one, 0UL);
  modredc15ul_intsub (m[0].one, m[0].one, m[0].m); /* 2^128 - m */
  modredc15ul_intmod (m[0].one, m[0].one, m[0].m);
  
#ifdef WANT_ASSERT_EXPENSIVE
  {
    modintredc15ul_t t;
    modredc15ul_get_uls (t, m[0].one, m);
    ASSERT_EXPENSIVE (modredc15ul_intequal_ul (t, 1UL));
  }
#endif
}


/* Returns the modulus to an array of unsigned longs. */
MAYBE_UNUSED
static inline void
modredc15ul_getmod_uls (modintredc15ul_t r, const modulusredc15ul_t m)
{
  modredc15ul_intset (r, m[0].m);
}


MAYBE_UNUSED
static inline void
modredc15ul_clearmod (modulusredc15ul_t m MAYBE_UNUSED)
{
  return;
}


/* Functions for residues */

/* Initialises a residueredc15ul_t and sets it to zero */
MAYBE_UNUSED
static inline void
modredc15ul_init (residueredc15ul_t r, const modulusredc15ul_t m MAYBE_UNUSED)
{
  modredc15ul_intset_ul (r, 0UL);
}


/* Initialises a residueredc15ul_t, but does not set it to zero. For fixed 
   length residueredc15ul_t, that leaves nothing to do at all. */
MAYBE_UNUSED
static inline void
modredc15ul_init_noset0 (residueredc15ul_t r MAYBE_UNUSED, 
			 const modulusredc15ul_t m MAYBE_UNUSED)
{
  return;
}


MAYBE_UNUSED
static inline void
modredc15ul_clear (residueredc15ul_t r MAYBE_UNUSED, 
		   const modulusredc15ul_t m MAYBE_UNUSED)
{
  return;
}


MAYBE_UNUSED
static inline void
modredc15ul_set (residueredc15ul_t r, const residueredc15ul_t s, 
		 const modulusredc15ul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (modredc15ul_intlt (s, m[0].m));
  modredc15ul_intset (r, s);
}


MAYBE_UNUSED
static inline void
modredc15ul_set_ul (residueredc15ul_t r, const unsigned long s, 
		    const modulusredc15ul_t m)
{
  modredc15ul_intset_ul (r, s);
  modredc15ul_tomontgomery (r, r, m);
#ifdef WANT_ASSERT_EXPENSIVE
  {
    modintredc15ul_t t;
    modredc15ul_get_uls (t, r, m);
    ASSERT_EXPENSIVE (t[0] == s && t[1] == 0UL);
  }
#endif
}


/* Sets the residueredc15ul_t to the class represented by the integer s. 
   Assumes that s is reduced (mod m), i.e. 0 <= s < m */

MAYBE_UNUSED
static inline void
modredc15ul_set_ul_reduced (residueredc15ul_t r, const unsigned long s, 
			    const modulusredc15ul_t m MAYBE_UNUSED)
{
  modredc15ul_set_ul (r, s, m);
}


MAYBE_UNUSED
static inline void
modredc15ul_set_uls (residueredc15ul_t r, const modintredc15ul_t s, 
		     const modulusredc15ul_t m)
{
  if (!modredc15ul_intlt (s, m[0].m))
    modredc15ul_intmod (r, s, m[0].m);
  else
    modredc15ul_intset (r, s);

  modredc15ul_tomontgomery (r, r, m);
}


MAYBE_UNUSED
static inline void
modredc15ul_set_uls_reduced (residueredc15ul_t r, const modintredc15ul_t s, 
			     const modulusredc15ul_t m)
{
  ASSERT (modredc15ul_intlt (s, m[0].m));
  modredc15ul_intset (r, s);
  modredc15ul_tomontgomery (r, r, m);
}


MAYBE_UNUSED 
static inline void 
modredc15ul_set0 (residueredc15ul_t r, const modulusredc15ul_t m MAYBE_UNUSED) 
{ 
  modredc15ul_intset_ul (r, 0UL);
}


MAYBE_UNUSED 
static inline void 
modredc15ul_set1 (residueredc15ul_t r, const modulusredc15ul_t m) 
{ 
  modredc15ul_intset (r, m[0].one);
}


/* Exchanges the values of the two arguments */

MAYBE_UNUSED
static inline void
modredc15ul_swap (residueredc15ul_t a, residueredc15ul_t b, 
		  const modulusredc15ul_t m MAYBE_UNUSED)
{
  modintredc15ul_t t;
  ASSERT_EXPENSIVE (modredc15ul_intlt (a, m[0].m));
  ASSERT_EXPENSIVE (modredc15ul_intlt (b, m[0].m));
  modredc15ul_intset (t, a);
  modredc15ul_intset (a, b);
  modredc15ul_intset (b, t);
}


/* Returns the least significant unsigned long of the residue. How to signal
   if the residue does not fit in one unsigned long? */

MAYBE_UNUSED
static inline unsigned long
modredc15ul_get_ul (const residueredc15ul_t s, 
		    const modulusredc15ul_t m MAYBE_UNUSED)
{
  unsigned long t[2];
  ASSERT_EXPENSIVE (modredc15ul_intlt (s, m[0].m));
  modredc15ul_frommontgomery (t, s, m);
  ASSERT (t[1] == 0UL);
  return t[0];
}


/* Returns the residue into an array of unsigned longs */

MAYBE_UNUSED
static inline void
modredc15ul_get_uls (modintredc15ul_t r, const residueredc15ul_t s, 
		     const modulusredc15ul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (modredc15ul_intlt (s, m[0].m));
  modredc15ul_frommontgomery (r, s, m);
}


MAYBE_UNUSED
static inline int
modredc15ul_equal (const residueredc15ul_t a, const residueredc15ul_t b, 
		   const modulusredc15ul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (modredc15ul_intlt (a, m[0].m));
  ASSERT_EXPENSIVE (modredc15ul_intlt (b, m[0].m));
  return (modredc15ul_intequal(a, b));
}


MAYBE_UNUSED
static inline int
modredc15ul_is0 (const residueredc15ul_t a, const modulusredc15ul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (modredc15ul_intlt (a, m[0].m));
  return (modredc15ul_intequal_ul(a, 0UL));
}


MAYBE_UNUSED
static inline int
modredc15ul_is1 (const residueredc15ul_t a, const modulusredc15ul_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (modredc15ul_intlt (a, m[0].m));
  return (a[0] == m[0].one[0] && a[1] == m[0].one[1]);
}


MAYBE_UNUSED
static inline void
modredc15ul_add (residueredc15ul_t r, const residueredc15ul_t a, 
		 const residueredc15ul_t b, const modulusredc15ul_t m)
{
  const unsigned long t0 = b[0], t1 = b[1]; /* r, a, and/or b may overlap */
  ASSERT_EXPENSIVE (modredc15ul_intlt (a, m[0].m));
  ASSERT_EXPENSIVE (modredc15ul_intlt (b, m[0].m));

  modredc15ul_intset (r, a);
  ularith_add_2ul_2ul (&(r[0]), &(r[1]), t0, t1);
  ularith_sub_2ul_2ul_ge (&(r[0]), &(r[1]), m[0].m[0], m[0].m[1]);
  ASSERT_EXPENSIVE (modredc15ul_intlt (r, m[0].m));
}


MAYBE_UNUSED
static inline void
modredc15ul_sub (residueredc15ul_t r, const residueredc15ul_t a, 
		 const residueredc15ul_t b, const modulusredc15ul_t m)
{
  ASSERT_EXPENSIVE (modredc15ul_intlt (a, m[0].m));
  ASSERT_EXPENSIVE (modredc15ul_intlt (b, m[0].m));

#ifdef HAVE_GCC_STYLE_AMD64_ASM
  {
    unsigned long s1 = m[0].m[0], s2 = m[0].m[1], t1 = a[0], t2 = a[1];
    
    __asm__ (
	     "subq %4, %0\n\t"
	     "sbbq %5, %1\n\t"    /* r -= b */
	     "cmovncq %6, %2\n\t" /* If !carry, s = 0 */
	     "cmovncq %6, %3\n"
	     : "+&r" (t1), "+&r" (t2), "+&r" (s1), "+r" (s2)
	     : ULARITH_CONSTRAINT_G (b[0]), 
	       ULARITH_CONSTRAINT_G (b[1]), "rm" (0UL)
	     : "cc"
	     );
    ularith_add_2ul_2ul (&t1, &t2, s1, s2);
    r[0] = t1;
    r[1] = t2;
  }
#else
  {
    unsigned long t1 = a[0], t2 = a[1];
    ularith_sub_2ul_2ul (&t1, &t2, b[0], b[1]);
    
    if (t2 > a[1] || (t2 == a[1] && t1 > a[0]))
      ularith_add_2ul_2ul (&t1, &t2, m[0].m[0], m[0].m[1]);

    r[0] = t1;
    r[1] = t2;
  }
#endif
}


MAYBE_UNUSED
static inline void
modredc15ul_add1 (residueredc15ul_t r, const residueredc15ul_t a, 
		  const modulusredc15ul_t m)
{
  modredc15ul_add(r, a, m[0].one, m);
}


MAYBE_UNUSED
static inline void
modredc15ul_add_ul (residueredc15ul_t r, const residueredc15ul_t a,
		    const unsigned long b, const modulusredc15ul_t m)
{
  residueredc15ul_t t;

  /* TODO: speed up */
  ASSERT_EXPENSIVE (modredc15ul_intlt (a, m[0].m));
  modredc15ul_init_noset0 (t, m);
  modredc15ul_set_ul (t, b, m);
  modredc15ul_add (r, a, t, m);
  modredc15ul_clear (t, m);
}


MAYBE_UNUSED
static inline void
modredc15ul_sub_ul (residueredc15ul_t r, const residueredc15ul_t a,
		    const unsigned long b, const modulusredc15ul_t m)
{
  residueredc15ul_t t;

  /* TODO: speed up */
  ASSERT_EXPENSIVE (modredc15ul_intlt (a, m[0].m));
  modredc15ul_init_noset0 (t, m);
  modredc15ul_set_ul (t, b, m);
  modredc15ul_sub (r, a, t, m);
  modredc15ul_clear (t, m);
}


MAYBE_UNUSED
static inline void
modredc15ul_neg (residueredc15ul_t r, const residueredc15ul_t a, 
		 const modulusredc15ul_t m)
{
  ASSERT_EXPENSIVE (modredc15ul_intlt (a, m[0].m));
  if (a[0] == 0UL && a[1] == 0UL)
    modredc15ul_set (r, a, m);
  else
    {
      unsigned long t1 = m[0].m[0], t2 = m[0].m[1];
      ularith_sub_2ul_2ul (&t1, &t2, a[0], a[1]);
      r[0] = t1;
      r[1] = t2;
    }
}


MAYBE_UNUSED
static inline void
modredc15ul_div2 (residueredc15ul_t r, const residueredc15ul_t a, 
		  const modulusredc15ul_t m)
{
  ASSERT_EXPENSIVE (m[0].m[0] % 2UL != 0UL);
  ASSERT_EXPENSIVE (modredc15ul_intlt (a, m[0].m));
  modredc15ul_intset (r, a);
  if (r[0] % 2UL == 1UL)
    ularith_add_2ul_2ul (&(r[0]), &(r[1]), m[0].m[0], m[0].m[1]);
  ularith_shrd (&(r[0]), r[1], 1);
  r[1] >>= 1;
}


#ifdef WANT_ASSERT_EXPENSIVE
#if defined(__x86_64__)
#define ABORT_IF_CY "jnc 1f\n\tlea _GLOBAL_OFFSET_TABLE_(%%rip), %%rbx\n\tcall abort@plt\n1:\n\t"
#elif defined(__i386__)
#define ABORT_IF_CY "jnc 1f\n\tcall abort\n1:\n\t"
#endif
#else
#define ABORT_IF_CY
#endif

MAYBE_UNUSED
static inline void
modredc15ul_mul (residueredc15ul_t r, const residueredc15ul_t a, 
                 const residueredc15ul_t b, const modulusredc15ul_t m)
{
#ifdef HAVE_GCC_STYLE_AMD64_ASM

  ASSERT_EXPENSIVE (modredc15ul_intlt (a, m[0].m));
  ASSERT_EXPENSIVE (modredc15ul_intlt (b, m[0].m));
#if defined(MODTRACE)
  printf ("((%lu * 2^%d + %lu) * (%lu * 2^%d + %lu) / 2^%d) %% "
	  "(%lu * 2^%d + %lu)", 
          a[1], LONG_BIT, a[0], b[1], LONG_BIT, b[0], 2 * LONG_BIT, 
	  m[0].m[1], LONG_BIT, m[0].m[0]);
#endif

/* Since m1>0, m*u is maximal for m0=1 and u=2^64-1, so
   u*m is bounded by (2^96 - 2^64 + 1)*(2^64 - 1) = 
   2^160 - 2^128 - 2^96 - 1. Doesn't really save anything, tho */

  unsigned long dummy;
  __asm__ (
    /* Product of low words */
    "movq %[a0], %%rax\n\t"
    "mulq %[b0]\n\t"         /* rdx:rax = a0*b0 */
    "movq %%rdx, %[t0]\n\t"
    /* Compute u0*m, add to t0:rax */
    "imulq %[invm], %%rax\n\t"
    "movq %%rax, %[t2]\n\t"  /* t2 = u0 */
    "xorl %k[t1], %k[t1]\n\t"
    "mulq %[m0]\n\t"         /* rdx:rax = u0*m0 <= (2^64-1)^2 */
    "negq %%rax\n\t"         /* if low word != 0, carry to high word */
    "movq %[t2], %%rax\n\t"  /* independent, goes in pipe 0 */
    "adcq %%rdx, %[t0]\n\t"
    "setc %b[t1]\n\t"        /* t1:t0 = (a0*b0+u0*m0)/2^64 <= (2^64-1)^2/2^64 = 2*2^64-4 */
    "mulq %[m1]\n\t"         /* rdx:rax <= (2^64-1)*(2^32-1) = 2^96-2^64-2^32+1 */
    "addq %%rax, %[t0]\n\t"
    "movq %[a0], %%rax\n\t"  /* independent, goes in pipe 0 */
    "adcq %%rdx, %[t1]\n\t"  /* t0:t1 = (a0*b0+u0*m)/2^64 */
    ABORT_IF_CY              /* <= ((2^64-1)^2 + 2^160 - 2^128 - 2^96 - 1)/2^64
                                <= 2^96 - 2^32 - 2 */
    
    /* 2 products of low and high word */
    "xorl %k[t2], %k[t2]\n\t"
    "mulq %[b1]\n\t"         /* rdx:rax = a0*b1 <= (2^64-1)*(2^32-1) */
    "addq %%rax, %[t0]\n\t"
    "movq %[a1], %%rax\n\t"  /* independent, goes in pipe 0 */
    "adcq %%rdx, %[t1]\n\t"  /* t1:t0 = (a0*b0+u0*m)/2^64 + a0*b1 */
    ABORT_IF_CY              /* <= 2^96 - 2^32 - 2 + 2*(2^64-1)*(2^32-1) 
                               = 2*2^96 - 2^64 - 2*2^32 - 1 */
    /* Free slot here */
    "mulq %[b0]\n\t"         /* rdx:rax = a1*b0 <= (2^64-1)*(2^32-1) */
    "addq %%rax, %[t0]\n\t"
    "movq %[a1], %%rax\n\t"  /* independent, goes in pipe 0 */
    "adcq %%rdx, %[t1]\n\t"  /* t1:t0 = (a0*b0+u0*m)/2^64 + a0*b1 + a1*b0 */
    ABORT_IF_CY              /* <= 2*2^96 - 2^64 - 2*2^32 - 1 + (2^64-1)*(2^32-1)
                                =  3*2^96 - 2*2^64 - 3*2^32 */
    /* Free slot here */
    /* Product of high words */
    "imulq %[b1], %%rax\n\t" /* rax = a1*b1 <= (2^32-1)^2 = 2^64 - 2*2^32 + 1 */
    "addq %%rax, %[t1]\n\t"
    "setc %b[t2]\n\t"        /* t2:t1:rax = (a0*b0+u0*m)/2^64 + a0*b1 + a1*b0 + a1*b1*2^64
                                <= 3*2^96 - 2*2^64 - 3*2^32 + (2^32-1)^2*2^64
                                = 3*2^96 - 2*2^64 - 3*2^32 + 2^128 - 2*2^96 + 2^64
                                = 2^128 + 2^96 - 2^64 - 3*2^32 */
    "movq %[t0], %%rax\n\t"
    /* Compute u1*m, add to t2:t1:t0 */
    "imulq %[invm], %%rax\n\t"
    "movq %%rax, %[t0]\n\t" /* t0 = u1 */
    "mulq %[m0]\n\t"       /* rdx:rax = u1*m0 <= (2^64-1)^2 = 2^128 - 2*2^64 + 1 */
    "negq %%rax\n\t"       /* if low word != 0, carry to high word */
    "movq %[t0], %%rax\n\t"
    "adcq %%rdx, %[t1]\n\t"
    "adcq $0,%[t2]\n\t"    /* t2:t1:0 = (a0*b0+u0*m)/2^64 + a0*b1 + a1*b0 + a1*b1*2^64 + u1*m0 */
    ABORT_IF_CY            /* <= 2^128 + 2^96 - 2^64 - 3*2^32 + 2^128 - 2*2^64 + 1
                              = 2*2^128 + 2^96 - 3*2^64 - 3*2^32 + 1 */
                           /* t2:t1 = ((a0*b0+u0*m)/2^64 + a0*b1 + a1*b0  + u1*m0)/2^64 + a1*b1
                              <= 2*2^64 + 2^32 - 4 */

    "mulq %[m1]\n\t"       /* rdx:rax = u1*m1 <= (2^64-1)*(2^32-1) */
    "addq %%rax, %[t1]\n\t"
    "adcq %%rdx, %[t2]\n\t"/* t2:t1 = ((a0*b0+u0*m)/2^64 + a0*b1 + a1*b0 + u1*m)/2^64 + a1*b1 */
    ABORT_IF_CY            /* <= (2^128 + 2^96 - 2^64 - 3*2^32 + 2^160 - 2^128 - 2^96 - 1)/2^64
                              = (2^160 - 2^64 - 3*2^32 - 1)/2^64
                              <= 2^96 - 2 */

  /* t2:t1:0:0 = a*b + u0*m + u1*m*2^64
     t2:t1 <= (a*b + u0*m + u1*m*2^64) / 2^128
     <= (m^2 + 2^64*m + 2^64*(2^64-1)*m) / 2^128
     =  (m^2 + 2^64*m + (2^128-2^64)*m)/2^128
     =  m + (m^2)/2^128
     <= m + (2^96*m)/2^128
     <= m + m/2^32 */
    "movq %[t1], %%rax\n\t" /* See if result > m */
    "movq %[t2], %%rdx\n\t"
    "subq %[m0], %[t1]\n\t"
    "sbbq %[m1], %[t2]\n\t"
    "cmovc %%rax, %[t1]\n\t" /* No carry -> copy new result */
    "cmovc %%rdx, %[t2]\n\t"
    : [t0] "=&r" (dummy), [t1] "=&r" (r[0]), [t2] "=&r" (r[1])
    : [a0] "rme" (a[0]), [a1] "rme" (a[1]), [b0] "rm" (b[0]), [b1] "rm" (b[1]),
      [m0] "rm" (m[0].m[0]), [m1] "rm" (m[0].m[1]), [invm] "rm" (m[0].invm)
    : "%rax", "%rdx", "cc"
  );
#else
  unsigned long pl, ph, t[4], k;

  ASSERT_EXPENSIVE (modredc15ul_intlt (a, m[0].m));
  ASSERT_EXPENSIVE (modredc15ul_intlt (b, m[0].m));
#if defined(MODTRACE)
  printf ("((%lu * 2^%d + %lu) * (%lu * 2^%d + %lu) / 2^%d) %% "
	  "(%lu * 2^%d + %lu)", 
          a[1], LONG_BIT, a[0], b[1], LONG_BIT, b[0], 2 * LONG_BIT, 
	  m[0].m[1], LONG_BIT, m[0].m[0]);
#endif
  
  /* Product of the two low words */
  ularith_mul_ul_ul_2ul (&(t[0]), &(t[1]), a[0], b[0]); /* t1:t0 = a[0]*b[0] <= W^2 - 2W + 1 */

  /* One REDC step */
  k = t[0] * m[0].invm; /* k <= W-1 */
  ularith_mul_ul_ul_2ul (&pl, &ph, k, m[0].m[0]); /* ph:pl = k*m[0] <= W^2 - 2W + 1 */
  /* t[0] + pl == 0 (mod W) */
  if (pl != 0UL)
    ph++; /* ph <= W-1 */
  t[2] = 0UL;
  ularith_add_ul_2ul (&(t[1]), &(t[2]), ph); /* t2:t1:0 = a[0]*b[0] + k*m[0] <= 2*W^2 - 4W + 2, so
                                                t2:t1 = (a[0]*b[0] + k*m[0]) / W <= 2*W - 4 */
  ularith_mul_ul_ul_2ul (&pl, &ph, k, m[0].m[1]); /* ph:pl <= (W^(1/2)-1)*(W-1) = W^(3/2) - W - W^(1/2) + 1 */
  ularith_add_2ul_2ul (&(t[1]), &(t[2]), pl, ph); /* t2:t1 <= W^(3/2) + W - W^(1/2) - 3 */

  /* Products of one low and one high word  */
  ularith_mul_ul_ul_2ul (&pl, &ph, a[1], b[0]);   /* ph:pl <= W^(3/2) - W - W^(1/2) + 1 */
  ularith_add_2ul_2ul (&(t[1]), &(t[2]), pl, ph); /* t2:t1 <= 2W^(3/2) - 2W^(1/2) - 2 */
  ularith_mul_ul_ul_2ul (&pl, &ph, a[0], b[1]);   /* ph:pl <= W^(3/2) - W - W^(1/2) + 1 */
  ularith_add_2ul_2ul (&(t[1]), &(t[2]), pl, ph); /* t2:t1 <= 3W^(3/2) - 3W^(1/2) - W - 1 */
  t[3] = 0UL;
  pl = a[1] * b[1];                               /* pl <= (W^(1/2)-1)^2 = W - 2W^(1/2) + 1 */
  ularith_add_ul_2ul (&(t[2]), &(t[3]), pl);      /* t3:t2:t1 <= W^2 + W^(3/2) - 3W^(1/2) - 1 */
  k = t[1] * m[0].invm;
  ularith_mul_ul_ul_2ul (&pl, &ph, k, m[0].m[0]); /* ph:pl <= W^2 - 2W + 1 */
  if (pl != 0UL)
    ph++;
  ularith_add_ul_2ul (&(t[2]), &(t[3]), ph);      /* t3:t2:t1 <= 2W^2 + W^(3/2) - 2W - 3W^(1/2), t1 = 0 so
						     t3:t2 <= 2W + W^(1/2) - 3 */
  ularith_mul_ul_ul_2ul (&pl, &ph, k, m[0].m[1]); /* ph:pl <= W^(3/2) - W - W^(1/2) + 1 */
  ularith_add_2ul_2ul (&(t[2]), &(t[3]), pl, ph); /* t3:t2 <= W^(3/2) + W - 2 */

  /* Result may be larger than m, but is < 2*m */

  ularith_sub_2ul_2ul_ge (&(t[2]), &(t[3]), m[0].m[0], m[0].m[1]);

  r[0] = t[2];
  r[1] = t[3];
#endif

#if defined(MODTRACE)
  printf (" == (%lu * 2^%d + %lu) /* PARI */ \n", r[1], LONG_BIT, r[0]);
#endif
  ASSERT_EXPENSIVE (modredc15ul_intlt (r, m[0].m));
}


MAYBE_UNUSED
static inline void
modredc15ul_sqr (residueredc15ul_t r, const residueredc15ul_t a, 
                 const modulusredc15ul_t m)
{
#ifdef HAVE_GCC_STYLE_AMD64_ASM
  ASSERT_EXPENSIVE (modredc15ul_intlt (a, m[0].m));
#if defined(MODTRACE)
  printf ("((%lu * 2^%d + %lu)^2 / 2^%d) %% (%lu * 2^%d + %lu)", 
          a[1], LONG_BIT, a[0], 2 * LONG_BIT, m[0].m[1], LONG_BIT, m[0].m[0]);
#endif
  
  unsigned long dummy;
  __asm__ (
    /* Product of low words */
    "movq %[a0], %%rax\n\t"
    "mulq %%rax\n\t"         /* rdx:rax = a0*a0 */
    "movq %%rdx, %[t0]\n\t"
    /* Compute u0*m, add to t0:rax */
    "imulq %[invm], %%rax\n\t"
    "movq %%rax, %[t2]\n\t"  /* t2 = u0 */
    "xorl %k[t1], %k[t1]\n\t"
    "mulq %[m0]\n\t"         /* rdx:rax = u0*m0 <= (2^64-1)^2 */
    "negq %%rax\n\t"         /* if low word != 0, carry to high word */
    "movq %[t2], %%rax\n\t"  /* independent, goes in pipe 0 */
    "adcq %%rdx, %[t0]\n\t"
    "setc %b[t1]\n\t"        /* t1:t0 = (a0*a0+u0*m0)/2^64 <= (2^64-1)^2/2^64 = 2*2^64-4 */
    "mulq %[m1]\n\t"         /* rdx:rax <= (2^64-1)*(2^32-1) = 2^96-2^64-2^32+1 */
    "addq %%rax, %[t0]\n\t"
    "movq %[a0], %%rax\n\t"  /* independent, goes in pipe 0 */
    "adcq %%rdx, %[t1]\n\t"  /* t0:t1 = (a0*a0+u0*m)/2^64 */
    ABORT_IF_CY              /* <= ((2^64-1)^2 + 2^160 - 2^96 - 2^128 - 1)/2^64
                                <= 2^96 - 2^32 - 2 */
    
    /* Product of low and high word */
    "xorl %k[t2], %k[t2]\n\t"
    "mulq %[a1]\n\t"         /* rdx:rax = a0*a1 <= (2^64-1)*(2^32-1) */
    "shlq $1,%%rax\n\t"
    "rclq $1,%%rdx\n\t"
    ABORT_IF_CY
    "addq %%rax, %[t0]\n\t"
    "adcq %%rdx, %[t1]\n\t"  /* t1:t0 = (a0*a0+u0*m)/2^64 + 2*a0*a1 */
    ABORT_IF_CY              /* <= 2^96 - 2^32 - 2 + 2*(2^64-1)*(2^32-1)
                                =  3*2^96 - 2*2^64 - 3*2^32 */
    "movq %[a1], %%rax\n\t"  /* independent, goes in pipe 0 */
    /* Free slot here */
    /* Product of high words */
    "imulq %%rax, %%rax\n\t" /* rax = a1*a1 <= (2^32-1)^2 = 2^64 - 2*2^32 + 1 */
    "addq %%rax, %[t1]\n\t"
    "setc %b[t2]\n\t"        /* t2:t1:rax = (a0*a0+u0*m)/2^64 + 2*a0*a1 + a1*a1*2^64
                                <= 3*2^96 - 2*2^64 - 3*2^32 + (2^32-1)^2*2^64
                                = 3*2^96 - 2*2^64 - 3*2^32 + 2^128 - 2*2^96 + 2^64
                                = 2^128 + 2^96 - 2^64 - 3*2^32 */
    "movq %[t0], %%rax\n\t"
    /* Compute u1*m, add to t2:t1:t0 */
    "imulq %[invm], %%rax\n\t"
    "movq %%rax, %[t0]\n\t" /* t0 = u1 */
    "mulq %[m0]\n\t"       /* rdx:rax = u1*m0 <= (2^64-1)^2 = 2^128 - 2*2^64 + 1 */
    "negq %%rax\n\t"       /* if low word != 0, carry to high word */
    "movq %[t0], %%rax\n\t"
    "adcq %%rdx, %[t1]\n\t"
    "adcq $0,%[t2]\n\t"    /* t2:t1:0 = (a0*a0+u0*m)/2^64 + 2*a0*a1 + a1*a1*2^64 + u1*m0 */
    ABORT_IF_CY            /* <= 2^128 + 2^96 - 2^64 - 3*2^32 + 2^128 - 2*2^64 + 1
                              = 2*2^128 + 2^96 - 3*2^64 - 3*2^32 + 1 */
                           /* t2:t1 = ((a0*a0+u*m)/2^64 + a0*a1 + a1*a0  + u*m0)/2^64 + a1*a1
                              <= 2*2^64 + 2^32 - 4 */

    "mulq %[m1]\n\t"       /* rdx:rax = u1*m1 <= (2^64-1)*(2^32-1) */
    "addq %%rax, %[t1]\n\t"
    "adcq %%rdx, %[t2]\n\t"/* t2:t1 = ((a0*a0+u0*m)/2^64 + 2*a0*a1 + u1*m)/2^64 + a1*a1 */
    ABORT_IF_CY            /* <= (2^128 + 2^96 - 2^64 - 3*2^32 + 2^160 - 2^128 - 2^96 - 1)/2^64
                              = (2^160 - 2^64 - 3*2^32 - 1)/2^64
                              <= 2^96 - 2 */

  /* t2:t1:0:0 = a*b + u0*m + u1*m*2^64
     t2:t1 <= (a*b + u0*m + u1*m*2^64) / 2^128
     <= (m^2 + 2^64*m + 2^64*(2^64-1)*m)/2^128
     =  (m^2 + 2^64*m + (2^128-2^64)*m)/2^128
     =  m + (m^2)/2^128
     <= m + (2^96*m)/2^128
     <= m + m/2^32 */
    "movq %[t1], %%rax\n\t" /* See if result > m */
    "movq %[t2], %%rdx\n\t"
    "subq %[m0], %[t1]\n\t"
    "sbbq %[m1], %[t2]\n\t"
    "cmovc %%rax, %[t1]\n\t" /* No carry -> copy new result */
    "cmovc %%rdx, %[t2]\n\t"
    : [t0] "=&r" (dummy), [t1] "=&r" (r[0]), [t2] "=&r" (r[1])
    : [a0] "rme" (a[0]), [a1] "rme" (a[1]), [m0] "rm" (m[0].m[0]), [m1] "rm" (m[0].m[1]), 
      [invm] "rm" (m[0].invm)
    : "%rax", "%rdx", "cc"
  );
#else
  ASSERT_EXPENSIVE (modredc15ul_intlt (a, m[0].m));
#if defined(MODTRACE)
  printf ("((%lu * 2^%d + %lu)^2 / 2^%d) %% (%lu * 2^%d + %lu)", 
          a[1], LONG_BIT, a[0], 2 * LONG_BIT, m[0].m[1], LONG_BIT, m[0].m[0]);
#endif
  
  unsigned long pl, ph, t[4], k;
  
  /* Square of low word */
  ularith_mul_ul_ul_2ul (&(t[0]), &(t[1]), a[0], a[0]); /* t1:t0 = a[0]*a[0] <= W^2 - 2W + 1 */

  /* One REDC step */
  k = t[0] * m[0].invm; /* k <= W-1 */
  ularith_mul_ul_ul_2ul (&pl, &ph, k, m[0].m[0]); /* ph:pl = k*m[0] <= W^2 - 2W + 1 */
  /* t[0] + pl == 0 (mod W) */
  if (pl != 0UL)
    ph++; /* ph <= W-1 */
  t[2] = 0UL;
  ularith_add_ul_2ul (&(t[1]), &(t[2]), ph); /* t2:t1:0 = a[0]*a[0] + k*m[0] <= 2*W^2 - 4W + 2, so
                                                t2:t1 = (a[0]*a[0] + k*m[0]) / W <= 2*W - 4 */
  ularith_mul_ul_ul_2ul (&pl, &ph, k, m[0].m[1]); /* ph:pl <= (W^(1/2)-1)*(W-1) = W^(3/2) - W - W^(1/2) + 1 */
  ularith_add_2ul_2ul (&(t[1]), &(t[2]), pl, ph); /* t2:t1 <= W^(3/2) + W - W^(1/2) - 3 */

  /* Product of low and high word  */
  ularith_mul_ul_ul_2ul (&pl, &ph, a[1], a[0]);   /* ph:pl <= W^(3/2) - W - W^(1/2) + 1 */
  ularith_add_2ul_2ul (&(t[1]), &(t[2]), pl, ph); /* t2:t1 <= 2W^(3/2) - 2W^(1/2) - 2 */
  ularith_add_2ul_2ul (&(t[1]), &(t[2]), pl, ph); /* t2:t1 <= 3W^(3/2) - 3W^(1/2) - W - 1 */
  t[3] = 0UL;

  /* Square of high word */
  pl = a[1] * a[1];                               /* pl <= (W^(1/2)-1)^2 = W - 2W^(1/2) + 1 */
  ularith_add_ul_2ul (&(t[2]), &(t[3]), pl);      /* t3:t2:t1 <= W^2 + W^(3/2) - 3W^(1/2) - 1 */
  k = t[1] * m[0].invm;
  ularith_mul_ul_ul_2ul (&pl, &ph, k, m[0].m[0]); /* ph:pl <= W^2 - 2W + 1 */
  if (pl != 0UL)
    ph++;
  ularith_add_ul_2ul (&(t[2]), &(t[3]), ph);      /* t3:t2:t1 <= 2W^2 + W^(3/2) - 2W - 3W^(1/2), t1 = 0 so
						     t3:t2 <= 2W + W^(1/2) - 3 */
  ularith_mul_ul_ul_2ul (&pl, &ph, k, m[0].m[1]); /* ph:pl <= W^(3/2) - W - W^(1/2) + 1 */
  ularith_add_2ul_2ul (&(t[2]), &(t[3]), pl, ph); /* t3:t2 <= W^(3/2) + W - 2 */

  /* Result may be larger than m, but is < 2*m */

  ularith_sub_2ul_2ul_ge (&(t[2]), &(t[3]), m[0].m[0], m[0].m[1]);

  r[0] = t[2];
  r[1] = t[3];
#endif

#if defined(MODTRACE)
  printf (" == (%lu * 2^%d + %lu) /* PARI */ \n", r[1], LONG_BIT, r[0]);
#endif
  ASSERT_EXPENSIVE (modredc15ul_intlt (r, m[0].m));
}


MAYBE_UNUSED
static inline int
modredc15ul_next (residueredc15ul_t r, const modulusredc15ul_t m)
{
  ularith_add_ul_2ul (&(r[0]), &(r[1]), 1UL);
  return (modredc15ul_intequal (r, m[0].m));
}


MAYBE_UNUSED
static inline int
modredc15ul_finished (const residueredc15ul_t r, const modulusredc15ul_t m)
{
  return (modredc15ul_intequal (r, m[0].m));
}


/* Division by small integer n, where (n-1)*m may NOT overflow the most 
   significant word. Returns 1 if n is invertible modulo m, 0 if not. 
   
   w_mod_n is word base (e.g., 2^32 or  2^64) mod n
   inv_n contains -1/i (mod n) if i is coprime to n, or 0 if i is not coprime 
   to n, for 0 <= i < n
   c = n^(-1) (mod word base)
*/

static inline int
modredc15ul_divn (residueredc15ul_t r, const residueredc15ul_t a, 
		  const unsigned long n, const unsigned long w_mod_n, 
		  const unsigned long *inv_n, const unsigned long c,
		  const modulusredc15ul_t m)
{
  const unsigned long an = ((a[1] % n)*w_mod_n + a[0] % n) % n;
  const unsigned long mn = ((m[0].m[1] % n)*w_mod_n + m[0].m[0] % n) % n;
  unsigned long k;
  residueredc15ul_t t;
#ifdef WANT_ASSERT_EXPENSIVE
  residueredc15ul_t a_backup;

  modredc15ul_init_noset0 (a_backup, m);
  modredc15ul_set (a_backup, a, m);
#endif
  
  modredc15ul_init_noset0 (t, m);
  t[1] = a[1];
  t[0] = a[0];
  
  if (inv_n[mn] == 0)
    {
      modredc15ul_clear (t, m);
      return 0;
    }

  /* Make t[1]:t[0] divisible by n */
  /* We want a+km == 0 (mod n), so k = -a*m^{-1} (mod n) */
  k = (inv_n[mn] * an) % n;
  ularith_mul_ul_ul_2ul (&(t[0]), &(t[1]), m[0].m[0], k);
  t[1] += m[0].m[1] * k;
  ularith_add_2ul_2ul (&(t[0]), &(t[1]), a[0], a[1]);
  
  /* Now t[1]:t[0] is divisible by n */
  ASSERT_EXPENSIVE (((t[1] % n)*w_mod_n + t[0] % n) % n == 0UL);
  
  /* May overwrite a */
  r[1] = t[1] / n;
  r[0] = t[0] * c;

#ifdef WANT_ASSERT_EXPENSIVE
  {
    unsigned long i;
    modredc15ul_set (t, r, m);
    for (i = 1; i < n; i++)
      modredc15ul_add (t, t, r, m);
    ASSERT_EXPENSIVE (modredc15ul_equal (t, a_backup, m));
    modredc15ul_clear (a_backup, m);
  }
#endif

  modredc15ul_clear (t, m);
  return 1;
}



/* prototypes of non-inline functions */
int modredc15ul_div3 (residueredc15ul_t, const residueredc15ul_t, 
		      const modulusredc15ul_t);
int modredc15ul_div5 (residueredc15ul_t, const residueredc15ul_t, 
		      const modulusredc15ul_t);
int modredc15ul_div7 (residueredc15ul_t, const residueredc15ul_t, 
		      const modulusredc15ul_t);
int modredc15ul_div11 (residueredc15ul_t, const residueredc15ul_t, 
		       const modulusredc15ul_t);
int modredc15ul_div13 (residueredc15ul_t, const residueredc15ul_t, 
		       const modulusredc15ul_t);
void modredc15ul_gcd (modintredc15ul_t, const residueredc15ul_t, 
		      const modulusredc15ul_t);
void modredc15ul_pow_ul (residueredc15ul_t, const residueredc15ul_t, 
			 const unsigned long, const modulusredc15ul_t);
void modredc15ul_2pow_ul (residueredc15ul_t, const unsigned long, 
                          const modulusredc15ul_t);
void modredc15ul_pow_mp (residueredc15ul_t, const residueredc15ul_t, 
			 const unsigned long *, const int, 
			 const modulusredc15ul_t);
void modredc15ul_2pow_mp (residueredc15ul_t, const unsigned long *, const int, 
			  const modulusredc15ul_t);
void modredc15ul_V_ul (residueredc15ul_t, const residueredc15ul_t, 
		       const unsigned long, const modulusredc15ul_t);
void modredc15ul_V_mp (residueredc15ul_t, const residueredc15ul_t, 
		       const unsigned long *, const int, 
		       const modulusredc15ul_t);
int modredc15ul_sprp (const residueredc15ul_t, const modulusredc15ul_t);
int modredc15ul_sprp2 (const modulusredc15ul_t);
int modredc15ul_isprime (const modulusredc15ul_t);
int modredc15ul_inv (residueredc15ul_t, const residueredc15ul_t, 
		     const modulusredc15ul_t);
int modredc15ul_jacobi (const residueredc15ul_t, const modulusredc15ul_t);
#endif  /* MODREDC_15UL_H */
