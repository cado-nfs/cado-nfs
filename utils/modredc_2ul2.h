/* Some functions for modular arithmetic with residues and modulus in
   unsigned long variables. The modulus can be up to 2 unsigned longs 
   in size with the two most significant bits zero (meaning 62 bits if 
   unsigned long has 32 bits, or 126 bits if unsigned long has 64 bits).
   Moduli must be odd and have the upper word non-zero. Residues are stored 
   in Montgomery form, reduction after multiplication is done with REDC. 
   Due to inlining, this file must be included in the caller's source code 
   with #include */

/* Naming convention: all function start with modredc2ul2, for 
   MODulus REDC 2 Unsigned Longs minus 2 bits, followed by underscore, 
   functionality of function (add, mul, etc), and possibly underscore and 
   specification of what argument types the function takes (_ul, etc). */

#ifndef MODREDC_2UL2_H
#define MODREDC_2UL2_H

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
   modredc2ul2_ */
#define MODREDC2UL2_RENAME(x) modredc2ul2_##x

#define MODREDC2UL2_SIZE 2
#define MODREDC2UL2_MINBITS LONG_BIT
#define MODREDC2UL2_MAXBITS (2 * LONG_BIT - 2)

typedef unsigned long residueredc2ul2_t[MODREDC2UL2_SIZE];
typedef unsigned long modintredc2ul2_t[MODREDC2UL2_SIZE];
typedef struct { 
  modintredc2ul2_t m;
  residueredc2ul2_t one;
  unsigned long invm;
} __modulusredc2ul2_t;
typedef __modulusredc2ul2_t modulusredc2ul2_t[1];


/* ==================== Functions used internally ==================== */

static inline int
modredc2ul2_intlt (const modintredc2ul2_t a, const modintredc2ul2_t b);

static inline void
modredc2ul2_add (residueredc2ul2_t r, const residueredc2ul2_t a, 
		 const residueredc2ul2_t b, const modulusredc2ul2_t m);
static inline void
modredc2ul2_get_uls (modintredc2ul2_t r, const residueredc2ul2_t s, 
		     const modulusredc2ul2_t m MAYBE_UNUSED);
static inline void
modredc2ul2_intset (modintredc2ul2_t r, const modintredc2ul2_t s);

MAYBE_UNUSED
static inline void
modredc2ul2_tomontgomery (residueredc2ul2_t r, const residueredc2ul2_t s,
			  const modulusredc2ul2_t m)
{
  int i;

  modredc2ul2_intset (r, s);
  /* TODO FIXME: ridiculously slow */
  for (i = 0; i < (MODREDC2UL2_SIZE * LONG_BIT); i++)
    modredc2ul2_add (r, r, r, m);
}


/* Do a one-word REDC, i.e., r == s / w (mod m), w = 2^LONG_BIT. 
   If m > w, r < 2m. If s<m, then r<m */
MAYBE_UNUSED
static inline void
modredc2ul2_redc1 (residueredc2ul2_t r, const residueredc2ul2_t s,
		   const modulusredc2ul2_t m)
{
  unsigned long t[4], k;
  
  k = s[0] * m[0].invm;
  ularith_mul_ul_ul_2ul (&(t[0]), &(t[1]), k, m[0].m[0]);
  if (s[0] != 0UL)
    t[1]++;
  t[2] = 0;
  ularith_add_ul_2ul (&(t[1]), &(t[2]), s[1]); /* t[2] <= 1 */
  ularith_mul_ul_ul_2ul (&(t[0]), &(t[3]), k, m[0].m[1]); /* t[3] < 2^w-1 */
  ularith_add_2ul_2ul (&(t[1]), &(t[2]), t[0], t[3]);     /* t[2] < 2^w */

  /* r = (k*m + s) / w, k <= w-1. If s < m, then r < m */
  r[0] = t[1];
  r[1] = t[2];
}

/* Converts s out of Montgomery form by dividing by 2^(2*LONG_BIT).
   Requires s < m. */
MAYBE_UNUSED
static inline void
modredc2ul2_frommontgomery (residueredc2ul2_t r, const residueredc2ul2_t s,
			    const modulusredc2ul2_t m)
{
  residueredc2ul2_t t;
  
  /* Do two REDC steps */
  modredc2ul2_redc1 (t, s, m);
  modredc2ul2_redc1 (r, t, m);
}


/* ================= Functions that are part of the API ================= */

/* Some functions for integers of the same width as the modulus */

MAYBE_UNUSED
static inline void
modredc2ul2_intset (modintredc2ul2_t r, const modintredc2ul2_t s)
{
  r[0] = s[0];
  r[1] = s[1];
}

MAYBE_UNUSED
static inline void
modredc2ul2_intset_ul (modintredc2ul2_t r, const unsigned long s)
{
  r[0] = s;
  r[1] = 0UL;
}

MAYBE_UNUSED
static inline int
modredc2ul2_intequal (const modintredc2ul2_t a, const modintredc2ul2_t b)
{
  return (a[0] == b[0] && a[1] == b[1]);
}

MAYBE_UNUSED
static inline int
modredc2ul2_intequal_ul (const modintredc2ul2_t a, const unsigned long b)
{
  return (a[0] == b && a[1] == 0UL);
}

/* Returns 1 if a < b, 0 otherwise */
MAYBE_UNUSED
static inline int
modredc2ul2_intlt (const modintredc2ul2_t a, const modintredc2ul2_t b)
{
    modintredc2ul2_t t;

    modredc2ul2_intset (t, a);
    return ularith_sub_2ul_2ul_cy (&(t[0]), &(t[1]), b[0], b[1]);
}

MAYBE_UNUSED
static inline int
modredc2ul2_intcmp (const modintredc2ul2_t a, const modintredc2ul2_t b)
{
  if (a[1] < b[1])
    return -1;
  if (a[1] > b[1])
    return 1;
  return (a[0] < b[0]) ? -1 : (a[0] == b[0]) ? 0 : 1;
}

MAYBE_UNUSED
static inline int
modredc2ul2_intcmp_ul (const modintredc2ul2_t a, const unsigned long b)
{
  if (a[1] > 0UL)
    return 1;
  return (a[0] < b) ? -1 : (a[0] == b) ? 0 : 1;
}

MAYBE_UNUSED
static inline int
modredc2ul2_intfits_ul (const modintredc2ul2_t a)
{
  return (a[1] == 0UL);
}

MAYBE_UNUSED
static inline void
modredc2ul2_intadd (modintredc2ul2_t r, const modintredc2ul2_t a,
		    const modintredc2ul2_t b)
{
  modintredc2ul2_t t;
  modredc2ul2_intset (t, a);
  ularith_add_2ul_2ul (&(t[0]), &(t[1]), b[0], b[1]);
  modredc2ul2_intset (r, t);
}

MAYBE_UNUSED
static inline void
modredc2ul2_intsub (modintredc2ul2_t r, const modintredc2ul2_t a,
		    const modintredc2ul2_t b)
{
  modintredc2ul2_t t;
  modredc2ul2_intset (t, a);
  ularith_sub_2ul_2ul (&(t[0]), &(t[1]), b[0], b[1]);
  modredc2ul2_intset (r, t);
}

/* Returns the number of bits in a, that is, floor(log_2(a))+1. 
   For a == 0 returns 0. */
MAYBE_UNUSED
static inline int
modredc2ul2_intbits (const modintredc2ul2_t a)
{
  if (a[1] > 0UL)
    return 2*LONG_BIT - ularith_clz (a[1]);

  if (a[0] > 0UL)
    return LONG_BIT - ularith_clz (a[0]);
  
  return 0;
}


MAYBE_UNUSED
static inline void
modredc2ul2_intshr (modintredc2ul2_t r, const modintredc2ul2_t s, const int i)
{
  r[0] = s[0];
  ularith_shrd (&(r[0]), s[1], i);
  r[1] = s[1] >> i;
}


MAYBE_UNUSED
static inline void
modredc2ul2_intshl (modintredc2ul2_t r, const modintredc2ul2_t s, const int i)
{
  r[1] = s[1];
  ularith_shld (&(r[1]), s[0], i);
  r[0] = s[0] << i;
}


/* r = n/d. We require d|n */
MAYBE_UNUSED
static inline void
modredc2ul2_intdivexact (modintredc2ul2_t r, const modintredc2ul2_t n,
                         const modintredc2ul2_t d)
{
  modintredc2ul2_t n1, d1;
  unsigned long invf, r0, k0, k1;
  int i;
#ifdef WANT_ASSERT_EXPENSIVE
  unsigned long s0 = n[0], s1 = n[1];
#endif
  
  modredc2ul2_intset (n1, n);
  modredc2ul2_intset (d1, d);

  /* Make d odd */
  if (d1[0] == 0UL)
    {
      ASSERT (n1[0] == 0UL);
      d1[0] = d1[1];
      d1[1] = 0UL;
      n1[0] = n1[1];
      n1[1] = 0UL;
    }
  ASSERT(d1[0] != 0UL);
  i = ularith_ctz (d1[0]);
  ularith_shrd (&(d1[0]), d1[1], i);
  d1[1] >>= i;
  ASSERT((n1[0] & ((1UL << i) - 1UL)) == 0UL);
  ularith_shrd (&(n1[0]), n1[1], i);
  n1[1] >>= i;
  
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
modredc2ul2_intmod (modintredc2ul2_t r, const modintredc2ul2_t n,
                    const modintredc2ul2_t d)
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
      modintredc2ul2_t t1, t2;
      unsigned long q, dummy;
      int i;

      /* Divide both by 2^k s.t. 2^(LONG_BIT-1) <= d/2^k < 2^LONG_BIT */
      modredc2ul2_intshr (t1, n, 1);
      modredc2ul2_intshr (t2, d, 1);
      if (t2[1] != 0UL)
        {
          i = LONG_BIT - ularith_clz (t2[1]);
          modredc2ul2_intshr (t1, t1, i);
          modredc2ul2_intshr (t2, t2, i);
        }
      ASSERT(t2[1] == 0);
      ASSERT((t2[0] & (1UL << (LONG_BIT - 1))) != 0UL);
      ASSERT (t1[1] < t2[0]);
      
      ularith_div_2ul_ul_ul (&q, &dummy, t1[0], t1[1], t2[0]);
      ularith_mul_ul_ul_2ul (&(t1[0]), &(t1[1]), q, d[0]);
      t1[1] += q * d[1];
      /* printf ("n=%lu*2^%d + %lu; d=%lu*2^%d + %lu; q=%lu; "
              "t1=%lu*2^%d + %lu;\n", 
              n[1], LONG_BIT, n[0], d[1], LONG_BIT, d[0], q,
              t1[1], LONG_BIT, t1[0]); */
      
      if (modredc2ul2_intlt (n, t1))
        modredc2ul2_intsub (t1, t1, d);
      ASSERT(!modredc2ul2_intlt (n, t1));
      modredc2ul2_intsub(r, n, t1);
    }
}


/* Functions for the modulus */

/* Init the modulus from a multi-word integer. s must point to an array of
   at least two unsigned longs, where s[0] is the low word of the modulus, 
   and s[1] is the high word. */
MAYBE_UNUSED
static inline void
modredc2ul2_initmod_uls (modulusredc2ul2_t m, const modintredc2ul2_t s)
{
  ASSERT (s[1] > 0UL);
  ASSERT (s[1] < (1UL << (LONG_BIT - 2)));
  modredc2ul2_intset (m[0].m, s);
  m[0].invm = -ularith_invmod (s[0]);

  modredc2ul2_intset_ul (m[0].one, 0UL);
  modredc2ul2_intsub (m[0].one, m[0].one, m[0].m); /* 2^128 - m */
  modredc2ul2_intmod (m[0].one, m[0].one, m[0].m);
  
#ifdef WANT_ASSERT_EXPENSIVE
  {
    modintredc2ul2_t t;
    modredc2ul2_get_uls (t, m[0].one, m);
    ASSERT_EXPENSIVE (modredc2ul2_intequal_ul (t, 1UL));
  }
#endif
}


/* Returns the modulus to an array of unsigned longs. */
MAYBE_UNUSED
static inline void
modredc2ul2_getmod_uls (modintredc2ul2_t r, const modulusredc2ul2_t m)
{
  modredc2ul2_intset (r, m[0].m);
}


MAYBE_UNUSED
static inline void
modredc2ul2_clearmod (modulusredc2ul2_t m MAYBE_UNUSED)
{
  return;
}


/* Functions for residues */

/* Initialises a residueredc2ul2_t and sets it to zero */
MAYBE_UNUSED
static inline void
modredc2ul2_init (residueredc2ul2_t r, const modulusredc2ul2_t m MAYBE_UNUSED)
{
  modredc2ul2_intset_ul (r, 0UL);
}


/* Initialises a residueredc2ul2_t, but does not set it to zero. For fixed 
   length residueredc2ul2_t, that leaves nothing to do at all. */
MAYBE_UNUSED
static inline void
modredc2ul2_init_noset0 (residueredc2ul2_t r MAYBE_UNUSED, 
			 const modulusredc2ul2_t m MAYBE_UNUSED)
{
  return;
}


MAYBE_UNUSED
static inline void
modredc2ul2_clear (residueredc2ul2_t r MAYBE_UNUSED, 
		   const modulusredc2ul2_t m MAYBE_UNUSED)
{
  return;
}


MAYBE_UNUSED
static inline void
modredc2ul2_set (residueredc2ul2_t r, const residueredc2ul2_t s, 
		 const modulusredc2ul2_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (modredc2ul2_intlt (s, m[0].m));
  modredc2ul2_intset (r, s);
}


MAYBE_UNUSED
static inline void
modredc2ul2_set_ul (residueredc2ul2_t r, const unsigned long s, 
		    const modulusredc2ul2_t m)
{
  modredc2ul2_intset_ul (r, s);
  modredc2ul2_tomontgomery (r, r, m);
#ifdef WANT_ASSERT_EXPENSIVE
  {
    modintredc2ul2_t t;
    modredc2ul2_get_uls (t, r, m);
    ASSERT_EXPENSIVE (t[0] == s && t[1] == 0UL);
  }
#endif
}


/* Sets the residueredc2ul2_t to the class represented by the integer s. 
   Assumes that s is reduced (mod m), i.e. 0 <= s < m */

MAYBE_UNUSED
static inline void
modredc2ul2_set_ul_reduced (residueredc2ul2_t r, const unsigned long s, 
			    const modulusredc2ul2_t m MAYBE_UNUSED)
{
  modredc2ul2_set_ul (r, s, m);
}


MAYBE_UNUSED
static inline void
modredc2ul2_set_uls (residueredc2ul2_t r, const modintredc2ul2_t s, 
		     const modulusredc2ul2_t m)
{
  if (!modredc2ul2_intlt (s, m[0].m))
    modredc2ul2_intmod (r, s, m[0].m);
  else
    modredc2ul2_intset (r, s);

  modredc2ul2_tomontgomery (r, r, m);
}


MAYBE_UNUSED
static inline void
modredc2ul2_set_uls_reduced (residueredc2ul2_t r, const modintredc2ul2_t s, 
			     const modulusredc2ul2_t m)
{
  ASSERT (modredc2ul2_intlt (s, m[0].m));
  modredc2ul2_intset (r, s);
  modredc2ul2_tomontgomery (r, r, m);
}


MAYBE_UNUSED 
static inline void 
modredc2ul2_set0 (residueredc2ul2_t r, const modulusredc2ul2_t m MAYBE_UNUSED) 
{ 
  modredc2ul2_intset_ul (r, 0UL);
}


MAYBE_UNUSED 
static inline void 
modredc2ul2_set1 (residueredc2ul2_t r, const modulusredc2ul2_t m) 
{ 
  modredc2ul2_intset (r, m[0].one);
}


/* Exchanges the values of the two arguments */

MAYBE_UNUSED
static inline void
modredc2ul2_swap (residueredc2ul2_t a, residueredc2ul2_t b, 
		  const modulusredc2ul2_t m MAYBE_UNUSED)
{
  modintredc2ul2_t t;
  ASSERT_EXPENSIVE (modredc2ul2_intlt (a, m[0].m));
  ASSERT_EXPENSIVE (modredc2ul2_intlt (b, m[0].m));
  modredc2ul2_intset (t, a);
  modredc2ul2_intset (a, b);
  modredc2ul2_intset (b, t);
}


/* Returns the least significant unsigned long of the residue. How to signal
   if the residue does not fit in one unsigned long? */

MAYBE_UNUSED
static inline unsigned long
modredc2ul2_get_ul (const residueredc2ul2_t s, 
		    const modulusredc2ul2_t m MAYBE_UNUSED)
{
  residueredc2ul2_t t;
  ASSERT_EXPENSIVE (modredc2ul2_intlt (s, m[0].m));
  modredc2ul2_frommontgomery (t, s, m);
  ASSERT (t[1] == 0UL);
  return t[0];
}


/* Returns the residue into an array of unsigned longs */

MAYBE_UNUSED
static inline void
modredc2ul2_get_uls (modintredc2ul2_t r, const residueredc2ul2_t s, 
		     const modulusredc2ul2_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (modredc2ul2_intlt (s, m[0].m));
  modredc2ul2_frommontgomery (r, s, m);
}


MAYBE_UNUSED
static inline int
modredc2ul2_equal (const residueredc2ul2_t a, const residueredc2ul2_t b, 
		   const modulusredc2ul2_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (modredc2ul2_intlt (a, m[0].m));
  ASSERT_EXPENSIVE (modredc2ul2_intlt (b, m[0].m));
  return (modredc2ul2_intequal(a, b));
}


MAYBE_UNUSED
static inline int
modredc2ul2_is0 (const residueredc2ul2_t a, const modulusredc2ul2_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (modredc2ul2_intlt (a, m[0].m));
  return (modredc2ul2_intequal_ul(a, 0UL));
}


MAYBE_UNUSED
static inline int
modredc2ul2_is1 (const residueredc2ul2_t a, const modulusredc2ul2_t m MAYBE_UNUSED)
{
  ASSERT_EXPENSIVE (modredc2ul2_intlt (a, m[0].m));
  return (modredc2ul2_intequal(a, m[0].one));
}


MAYBE_UNUSED
static inline void
modredc2ul2_add (residueredc2ul2_t r, const residueredc2ul2_t a, 
		 const residueredc2ul2_t b, const modulusredc2ul2_t m)
{
  const unsigned long t0 = b[0], t1 = b[1]; /* r, a, and/or b may overlap */
  ASSERT_EXPENSIVE (modredc2ul2_intlt (a, m[0].m));
  ASSERT_EXPENSIVE (modredc2ul2_intlt (b, m[0].m));

  modredc2ul2_intset (r, a);
  ularith_add_2ul_2ul (&(r[0]), &(r[1]), t0, t1);
  ularith_sub_2ul_2ul_ge (&(r[0]), &(r[1]), m[0].m[0], m[0].m[1]);
  ASSERT_EXPENSIVE (modredc2ul2_intlt (r, m[0].m));
}


MAYBE_UNUSED
static inline void
modredc2ul2_sub (residueredc2ul2_t r, const residueredc2ul2_t a, 
		 const residueredc2ul2_t b, const modulusredc2ul2_t m)
{
  ASSERT_EXPENSIVE (modredc2ul2_intlt (a, m[0].m));
  ASSERT_EXPENSIVE (modredc2ul2_intlt (b, m[0].m));

#ifdef HAVE_GCC_STYLE_AMD64_ASM
  {
    unsigned long s1 = m[0].m[0], s2 = m[0].m[1], t1 = a[0], t2 = a[1];
    
    __asm__ (
	     "subq %4, %0\n\t"
	     "sbbq %5, %1\n\t"    /* t -= b */
	     "cmovncq %6, %2\n\t" /* If !carry, s = 0 */
	     "cmovncq %6, %3\n"
	     : "+&r" (t1), "+&r" (t2), "+&r" (s1), "+r" (s2)
	     : "g" (b[0]), "g" (b[1]), "rm" (0UL)
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
modredc2ul2_add1 (residueredc2ul2_t r, const residueredc2ul2_t a, 
		  const modulusredc2ul2_t m)
{
  modredc2ul2_add(r, a, m[0].one, m);
}


MAYBE_UNUSED
static inline void
modredc2ul2_add_ul (residueredc2ul2_t r, const residueredc2ul2_t a,
		    const unsigned long b, const modulusredc2ul2_t m)
{
  residueredc2ul2_t t;

  /* TODO: speed up */
  ASSERT_EXPENSIVE (modredc2ul2_intlt (a, m[0].m));
  modredc2ul2_init_noset0 (t, m);
  modredc2ul2_set_ul (t, b, m);
  modredc2ul2_add (r, a, t, m);
  modredc2ul2_clear (t, m);
}


MAYBE_UNUSED
static inline void
modredc2ul2_sub_ul (residueredc2ul2_t r, const residueredc2ul2_t a,
		    const unsigned long b, const modulusredc2ul2_t m)
{
  residueredc2ul2_t t;

  /* TODO: speed up */
  ASSERT_EXPENSIVE (modredc2ul2_intlt (a, m[0].m));
  modredc2ul2_init_noset0 (t, m);
  modredc2ul2_set_ul (t, b, m);
  modredc2ul2_sub (r, a, t, m);
  modredc2ul2_clear (t, m);
}


MAYBE_UNUSED
static inline void
modredc2ul2_neg (residueredc2ul2_t r, const residueredc2ul2_t a, 
		 const modulusredc2ul2_t m)
{
  ASSERT_EXPENSIVE (modredc2ul2_intlt (a, m[0].m));
  if (modredc2ul2_is0 (a, m))
    modredc2ul2_set (r, a, m);
  else
    modredc2ul2_intsub (r, m[0].m, a);
}


MAYBE_UNUSED
static inline void
modredc2ul2_div2 (residueredc2ul2_t r, const residueredc2ul2_t a, 
		  const modulusredc2ul2_t m)
{
  ASSERT_EXPENSIVE (m[0].m[0] % 2UL != 0UL);
  ASSERT_EXPENSIVE (modredc2ul2_intlt (a, m[0].m));
  modredc2ul2_intset (r, a);
  if (r[0] % 2UL == 1UL)
    modredc2ul2_intadd (r, r, m[0].m);
  modredc2ul2_intshr (r, r, 1);
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
modredc2ul2_mul (residueredc2ul2_t r, const residueredc2ul2_t a, 
               const residueredc2ul2_t b, const modulusredc2ul2_t m)
{
#ifdef HAVE_GCC_STYLE_AMD64_ASM

  ASSERT_EXPENSIVE (modredc2ul2_intlt (a, m[0].m));
  ASSERT_EXPENSIVE (modredc2ul2_intlt (b, m[0].m));
#if defined(MODTRACE)
  printf ("((%lu * 2^%d + %lu) * (%lu * 2^%d + %lu) / 2^%d) %% "
	  "(%lu * 2^%d + %lu)", 
          a[1], LONG_BIT, a[0], b[1], LONG_BIT, b[0], 2 * LONG_BIT, 
	  m[0].m[1], LONG_BIT, m[0].m[0]);
#endif

  unsigned long dummy;
  __asm__ (
    /* Product of low words */
    "movq %[a0], %%rax\n\t"
    "mulq %[b0]\n\t"         /* rdx:rax = a0*b0 <= (2^64-1)^2 */
    "movq %%rdx, %[t0]\n\t"
    /* Compute u0*m, add to t0:rax */
    "imulq %[invm], %%rax\n\t"
    "movq %%rax, %[t2]\n\t"  /* t2 = u0 */
    "xorl %k[t1], %k[t1]\n\t"
    "mulq %[m0]\n\t"         /* rdx:rax = u0*m0 <= (2^64-1)^2 */
    "negq %%rax\n\t"         /* if low word != 0, carry to high word */
    "movq %[t2], %%rax\n\t"  /* independent, goes in pipe 0 */
    "adcq %%rdx, %[t0]\n\t"
    "setc %b[t1]\n\t"        /* t1:t0 = (a0*b0 + u0*m0) / 2^64 <= 2*2^64 - 4 */
    "mulq %[m1]\n\t"         
    "addq %%rax, %[t0]\n\t"
    "movq %[a0], %%rax\n\t"  /* independent, goes in pipe 0 */
    "adcq %%rdx, %[t1]\n\t"  /* t1:t0 = (a0*b0+u0*m)/2^64 */
    ABORT_IF_CY              /* <= 2^126 - 2^62 */
                             
    
    /* 2 products of low and high word */
    "xorl %k[t2], %k[t2]\n\t"
    "mulq %[b1]\n\t"         /* rdx:rax = a0*b1 <= (2^64-1)*(2^63-2^32-1) */
    "addq %%rax, %[t0]\n\t"
    "movq %[a1], %%rax\n\t"  /* independent, goes in pipe 0 */
    "adcq %%rdx, %[t1]\n\t"  /* t1:t0 = (a0*b0+u0*m)/2^64 + a0*b1 */
    ABORT_IF_CY              /* <= 2^126 - 2^62 + (2^64-1)*(2^63-2^32-1)
    			        = 2^127 + 2^126 - 2^96 ... */
                             
    /* Free slot here */
    "mulq %[b0]\n\t"         /* rdx:rax = a1*b0 <= (2^63-2^32-1)*(2^64-1) */
    "addq %%rax, %[t0]\n\t"
    "movq %[a1], %%rax\n\t"  /* independent, goes in pipe 0 */
    "adcq %%rdx, %[t1]\n\t"  
    "setc %b[t2]\n\t"        /* t2:t1:t0 = (a0*b0+u0*m)/2^64 + a0*b1 + a1*b0 */
    			     /* <= 2^126 - 2^62 + 2*(2^64-1)*(2^63-2^32-1)
                                = 2^128 + 2^126 - 2*2^96 ... */
    /* Product of high words */
    "mulq %[b1]\n\t"         /* rdx:rax = a1*b1 <= (2^63-2^32-1)^2 */
    "addq %%rax, %[t1]\n\t"
    "movq %[t0], %%rax\n\t"
    "adcq %%rdx, %[t2]\n\t"  /* t2:t1:t0 = (a*b+u0*m)/2^64 */
    ABORT_IF_CY              /* <= ((2^127-2^96-1)^2+(2^64-1)*(2^126-2^64+1))/2^64 
                                = 2^190 - 2^160 ... */
    /* Free slot here */
    /* Compute u1*m, add to t2:t1:t0 */
    "imulq %[invm], %%rax\n\t"
    "movq %%rax, %[t0]\n\t" /* t0 = u1 */
    /* Free slot here */
    "mulq %[m0]\n\t"        /* rdx:rax = m0*u1 <= (2^64-1)^2 */
    "negq %%rax\n\t"        /* if low word != 0, carry to high word */
    "movq %[t0], %%rax\n\t"
    "adcq %%rdx, %[t1]\n\t"
    "adcq $0,%[t2]\n\t"     /* t2:t1:0 = (a*b+u0*m)/2^64 + u1*m0 */
    ABORT_IF_CY             /* <= 2^190 - 2^160 + 2*2^128 + 2^126 ... */
                            
    "mulq %[m1]\n\t"        /* rdx:rax = u1*m1 */
    "addq %%rax, %[t1]\n\t"
    "adcq %%rdx, %[t2]\n\t" /* t2:t1 = ((a*b+u0*m)/2^64 + u1*m)/2^64 */
    ABORT_IF_CY             /* <= 2^127 - 2^96 - 1 */
                           
    "movq %[t1], %%rax\n\t" /* See if result > m */
    "movq %[t2], %%rdx\n\t"
    "subq %[m0], %[t1]\n\t"
    "sbbq %[m1], %[t2]\n\t"
    "cmovc %%rax, %[t1]\n\t" /* No carry -> copy new result */
    "cmovc %%rdx, %[t2]\n\t"
    : [t0] "=&r" (dummy), [t1] "=&r" (r[0]), [t2] "=&r" (r[1])
    : [a0] "g" (a[0]), [a1] "g" (a[1]), [b0] "rm" (b[0]), [b1] "rm" (b[1]),
      [m0] "rm" (m[0].m[0]), [m1] "rm" (m[0].m[1]), [invm] "rm" (m[0].invm)
    : "%rax", "%rdx", "cc"
  );
#else /* HAVE_GCC_STYLE_AMD64_ASM */

  unsigned long pl, ph, t[4], k;
  
  ASSERT_EXPENSIVE (modredc2ul2_intlt (a, m[0].m));
  ASSERT_EXPENSIVE (modredc2ul2_intlt (b, m[0].m));
#if defined(MODTRACE)
  printf ("((%lu * 2^%d + %lu) * (%lu * 2^%d + %lu) / 2^%d) %% "
	  "(%lu * 2^%d + %lu)", 
          a[1], LONG_BIT, a[0], b[1], LONG_BIT, b[0], 2 * LONG_BIT, 
	  m[0].m[1], LONG_BIT, m[0].m[0]);
#endif

  /* m < 1/4 W^2,  a,b < m */
  
  /* Product of the two low words */
  ularith_mul_ul_ul_2ul (&(t[0]), &(t[1]), a[0], b[0]); 

  /* One REDC step */
  modredc2ul2_redc1 (t, t, m); /* t < 2m < 1/2 W^2 */

  /* Products of one low and one high word  */
  ularith_mul_ul_ul_2ul (&pl, &ph, a[1], b[0]);   /* ph:pl < 1/4 W^2 */
  ularith_add_2ul_2ul (&(t[0]), &(t[1]), pl, ph); /* t1:t0 < 3/4 W^2 */
  ularith_mul_ul_ul_2ul (&pl, &ph, a[0], b[1]);   /* ph:pl < 1/4 W^2 */
  ularith_add_2ul_2ul (&(t[0]), &(t[1]), pl, ph); /* t1:t0 < W^2 */

  /* Product of the two high words */
  ularith_mul_ul_ul_2ul (&pl, &(t[2]), a[1], b[1]); /* t2:pl < 1/16 W^2 */
  ularith_add_ul_2ul (&(t[1]), &(t[2]), pl); /* t2:t1:t0 < 1/16 W^3 + W^2 */

  /* Compute t2:t1:t0 := t2:t1:t0 + km, km < Wm < 1/4 W^3 */
  k = t[0] * m[0].invm;
  ularith_mul_ul_ul_2ul (&pl, &ph, k, m[0].m[0]);
  if (t[0] != 0UL)
    ph++; /* t[0] = 0 */
  ularith_add_ul_2ul (&(t[1]), &(t[2]), ph);
  ularith_mul_ul_ul_2ul (&pl, &ph, k, m[0].m[1]); /* ph:pl < 1/4 W^2 */
  ularith_add_2ul_2ul (&(t[1]), &(t[2]), pl, ph);
  /* t2:t1:0 < 1/16 W^3 + W^2 + 1/4 W^3 < 5/16 W^3 + W^2 */

  /* Result may be larger than m, but is < 2*m */

  ularith_sub_2ul_2ul_ge (&(t[1]), &(t[2]), m[0].m[0], m[0].m[1]);

  r[0] = t[1];
  r[1] = t[2];
#endif
#if defined(MODTRACE)
  printf (" == (%lu * 2^%d + %lu) /* PARI */ \n", r[1], LONG_BIT, r[0]);
#endif
  ASSERT_EXPENSIVE (modredc2ul2_intlt (r, m[0].m));
}


MAYBE_UNUSED
static inline void
modredc2ul2_sqr (residueredc2ul2_t r, const residueredc2ul2_t a, 
                 const modulusredc2ul2_t m)
{
#ifdef HAVE_GCC_STYLE_AMD64_ASM

  ASSERT_EXPENSIVE (modredc2ul2_intlt (a, m[0].m));
#if defined(MODTRACE)
  printf ("((%lu * 2^%d + %lu) * (%lu * 2^%d + %lu) / 2^%d) %% "
	  "(%lu * 2^%d + %lu)", 
          a[1], LONG_BIT, a[0], b[1], LONG_BIT, b[0], 2 * LONG_BIT, 
	  m[0].m[1], LONG_BIT, m[0].m[0]);
#endif

/* m <= 2^126-1
   Since m1>0, m*u is maximal for m0=1 and u=2^64-1, so
   u*m is bounded by (2^126 - 2^64 + 1)*(2^64 - 1) = 
   2^190 - 2^128 - 2^126 + 2*2^64 - 1.
   If a,b <= 2^127-2^96-1, then
   ((a*b+u0*m)/2^64 + u1*m)/2^64 <=  2^127-2^96-1
   If we allow non-canonical residues up to 2^127-2^96-1, we can skip
   the final conditional subtraction. These residues are still < 2^127,
   so an addition does not overflow */

  unsigned long dummy;
  __asm__ (
    /* Product of low words */
    "movq %[a0], %%rax\n\t"
    "mulq %[a0]\n\t"         /* rdx:rax = a0^2 <= (2^64-1)^2 */
    "movq %%rdx, %[t0]\n\t"
    /* Compute u0*m, add to t0:rax */
    "imulq %[invm], %%rax\n\t"
    "movq %%rax, %[t2]\n\t"  /* t2 = u0 */
    "xorl %k[t1], %k[t1]\n\t"
    "mulq %[m0]\n\t"         /* rdx:rax = u0*m0 <= (2^64-1)^2 */
    "negq %%rax\n\t"         /* if low word != 0, carry to high word */
    "movq %[t2], %%rax\n\t"  /* independent, goes in pipe 0 */
    "adcq %%rdx, %[t0]\n\t"
    "setc %b[t1]\n\t"        /* t1:t0 = (a0^2 + u0*m0) / 2^64 <= 2*2^64 - 4 */
    "mulq %[m1]\n\t"         
    "addq %%rax, %[t0]\n\t"
    "movq %[a0], %%rax\n\t"  /* independent, goes in pipe 0 */
    "adcq %%rdx, %[t1]\n\t"  /* t1:t0 = (a0^2+u0*m)/2^64 */
    ABORT_IF_CY              /* <= 2^126 - 2^62 */
                             
    /* 2 products of low and high word */
    "xorl %k[t2], %k[t2]\n\t"
    "mulq %[a1]\n\t"         /* rdx:rax = a0*a1 <= (2^64-1)*(2^63-2^32-1) */
    "addq %%rax, %[t0]\n\t"
    "adcq %%rdx, %[t1]\n\t"  /* t1:t0 = (a0^2+u0*m)/2^64 + a0*a1 */
    ABORT_IF_CY              /* <= 2^126 - 2^62 + (2^64-1)*(2^63-2^32-1)
    			        = 2^127 + 2^126 - 2^96 ... */
    "addq %%rax, %[t0]\n\t"
    "adcq %%rdx, %[t1]\n\t"  
    "movq %[a1], %%rax\n\t"  /* independent, goes in pipe 0 */
    "setc %b[t2]\n\t"        /* t2:t1:t0 = (a0*a0+u0*m)/2^64 + 2*a0*a1 */
    			     /* <= 2^126 - 2^62 + 2*(2^64-1)*(2^63-2^32-1)
                                = 2^128 + 2^126 - 2*2^96 ... */
    
    /* Product of high words */
    "mulq %[a1]\n\t"         /* rdx:rax = a1^2 <= (2^63-2^32-1)^2 */
    "addq %%rax, %[t1]\n\t"
    "movq %[t0], %%rax\n\t"
    "adcq %%rdx, %[t2]\n\t"  /* t2:t1:t0 = (a^2+u0*m)/2^64 */
    ABORT_IF_CY              /* <= ((2^127-2^96-1)^2+(2^64-1)*(2^126-2^64+1))/2^64 
                                = 2^190 - 2^160 ... */
    /* Free slot here */
    /* Compute u1*m, add to t2:t1:t0 */
    "imulq %[invm], %%rax\n\t"
    "movq %%rax, %[t0]\n\t" /* t0 = u1 */
    /* Free slot here */
    "mulq %[m0]\n\t"        /* rdx:rax = m0*u1 <= (2^64-1)^2 */
    "negq %%rax\n\t"        /* if low word != 0, carry to high word */
    "movq %[t0], %%rax\n\t"
    "adcq %%rdx, %[t1]\n\t"
    "adcq $0,%[t2]\n\t"     /* t2:t1:0 = (a*a+u0*m)/2^64 + u1*m0 */
    ABORT_IF_CY             /* <= 2^190 - 2^160 + 2*2^128 + 2^126 ... */
                            
    "mulq %[m1]\n\t"        /* rdx:rax = u1*m1 */
    "addq %%rax, %[t1]\n\t"
    "adcq %%rdx, %[t2]\n\t" /* t2:t1 = ((a*a+u0*m)/2^64 + u1*m)/2^64 */
    ABORT_IF_CY             /* <= 2^127 - 2^96 - 1 */
                           
    "movq %[t1], %%rax\n\t" /* See if result > m */
    "movq %[t2], %%rdx\n\t"
    "subq %[m0], %[t1]\n\t"
    "sbbq %[m1], %[t2]\n\t"
    "cmovc %%rax, %[t1]\n\t" /* No carry -> copy new result */
    "cmovc %%rdx, %[t2]\n\t"
    : [t0] "=&r" (dummy), [t1] "=&r" (r[0]), [t2] "=&r" (r[1])
    : [a0] "g" (a[0]), [a1] "g" (a[1]), 
      [m0] "rm" (m[0].m[0]), [m1] "rm" (m[0].m[1]), [invm] "rm" (m[0].invm)
    : "%rax", "%rdx", "cc"
  );
#else /* HAVE_GCC_STYLE_AMD64_ASM */

  unsigned long pl, ph, t[4], k;
  
  ASSERT_EXPENSIVE (modredc2ul2_intlt (a, m[0].m));
#if defined(MODTRACE)
  printf ("((%lu * 2^%d + %lu)^2 %% (%lu * 2^%d + %lu)", 
          a[1], LONG_BIT, a[0], 2 * LONG_BIT, m[0].m[1], LONG_BIT, m[0].m[0]);
#endif

  /* m < 1/4 W^2,  a < m */
  
  /* Square of the low word */
  ularith_mul_ul_ul_2ul (&(t[0]), &(t[1]), a[0], a[0]); 

  /* One REDC step */
  modredc2ul2_redc1 (t, t, m); /* t < 2m < 1/2 W^2 */

  /* Products of low and high word  */
  ularith_mul_ul_ul_2ul (&pl, &ph, a[1], a[0]);   /* ph:pl < 1/4 W^2 */
  ularith_add_2ul_2ul (&(t[0]), &(t[1]), pl, ph); /* t1:t0 < 3/4 W^2 */
  ularith_add_2ul_2ul (&(t[0]), &(t[1]), pl, ph); /* t1:t0 < W^2 */

  /* Square of high word */
  ularith_mul_ul_ul_2ul (&pl, &(t[2]), a[1], a[1]); /* t2:pl < 1/16 W^2 */
  ularith_add_ul_2ul (&(t[1]), &(t[2]), pl); /* t2:t1:t0 < 1/16 W^3 + W^2 */

  /* Compute t2:t1:t0 := t2:t1:t0 + km, km < Wm < 1/4 W^3 */
  k = t[0] * m[0].invm;
  ularith_mul_ul_ul_2ul (&pl, &ph, k, m[0].m[0]);
  if (t[0] != 0UL)
    ph++; /* t[0] = 0 */
  ularith_add_ul_2ul (&(t[1]), &(t[2]), ph);
  ularith_mul_ul_ul_2ul (&pl, &ph, k, m[0].m[1]); /* ph:pl < 1/4 W^2 */
  ularith_add_2ul_2ul (&(t[1]), &(t[2]), pl, ph);
  /* t2:t1:0 < 1/16 W^3 + W^2 + 1/4 W^3 < 5/16 W^3 + W^2 */

  /* Result may be larger than m, but is < 2*m */

  ularith_sub_2ul_2ul_ge (&(t[1]), &(t[2]), m[0].m[0], m[0].m[1]);

  r[0] = t[1];
  r[1] = t[2];
#endif
#if defined(MODTRACE)
  printf (" == (%lu * 2^%d + %lu) /* PARI */ \n", r[1], LONG_BIT, r[0]);
#endif
  ASSERT_EXPENSIVE (modredc2ul2_intlt (r, m[0].m));
}


MAYBE_UNUSED
static inline int
modredc2ul2_next (residueredc2ul2_t r, const modulusredc2ul2_t m)
{
  ularith_add_ul_2ul (&(r[0]), &(r[1]), 1UL);
  return (modredc2ul2_intequal (r, m[0].m));
}


MAYBE_UNUSED
static inline int
modredc2ul2_finished (const residueredc2ul2_t r, const modulusredc2ul2_t m)
{
  return (modredc2ul2_intequal (r, m[0].m));
}

/* Division by small integer n, where (n-1)*m may overflow the most 
   significant word. Returns 1 if n is invertible modulo m, 0 if not. 
   
   w_mod_n is word base (e.g., 2^32 or  2^64) mod n
   inv_n contains -1/i (mod n) if i is coprime to n, or 0 if i is not coprime 
   to n, for 0 <= i < n
   c = n^(-1) (mod word base)
*/

MAYBE_UNUSED
static inline int
modredc2ul2_divn (residueredc2ul2_t r, const residueredc2ul2_t a, 
		  const unsigned long n, const unsigned long w_mod_n, 
		  const unsigned long *inv_n, const unsigned long c,
		  const modulusredc2ul2_t m)
{
  const unsigned long an = ((a[1] % n) * w_mod_n + a[0] % n) % n;
  const unsigned long mn = ((m[0].m[1] % n) * w_mod_n + m[0].m[0] % n) % n;

  residueredc2ul2_t t, t2;
  unsigned long k;
#ifdef WANT_ASSERT_EXPENSIVE
  residueredc2ul2_t a_backup;
#endif
  
  if (inv_n[mn] == 0)
    return 0;

#ifdef WANT_ASSERT_EXPENSIVE
  modredc2ul2_init_noset0 (a_backup, m);
  modredc2ul2_set (a_backup, a, m);
#endif
  modredc2ul2_init_noset0 (t, m);
  modredc2ul2_init_noset0 (t2, m);
  modredc2ul2_set (t, a, m);
  
  /* Make t[1]:t[0] == a+km (mod w^2) with a+km divisible by n */
  /* We want a+km == 0 (mod n), so k = -a*m^{-1} (mod n) */
  k = (inv_n[mn] * an) % n;
  ASSERT_EXPENSIVE ((an + k*mn) % n == 0);
  ularith_mul_ul_ul_2ul (&(t[0]), &(t[1]), m[0].m[0], k);
  t[1] += m[0].m[1] * k;
  ularith_add_2ul_2ul (&(t[0]), &(t[1]), a[0], a[1]);

  /* We want r = (a+km)/n. */

  /* May overwrite a */
  r[0] = t[0] * c;

  /* r0 == (a+km)/n (mod w) 
     (r1*w + r0) * n = (a+km)
     (r1*w + r0) * n == t (mod w^2)
     r1*w*n == t - n*r0 (mod w^2)
                            t - n*r0 == 0 (mod w), thus
     r1*n == (t - n*r0)/w (mod w) */

  ularith_mul_ul_ul_2ul (&(t2[0]), &(t2[1]), r[0], n);
  ularith_sub_2ul_2ul (&(t[0]), &(t[1]), t2[0], t2[1]);
  ASSERT_EXPENSIVE (t[0] == 0UL);
  r[1] = t[1] * c;

#ifdef WANT_ASSERT_EXPENSIVE
  {
    unsigned long i;
    modredc2ul2_set (t, r, m);
    for (i = 1; i < n; i++)
      modredc2ul2_add (t, t, r, m);
    ASSERT_EXPENSIVE (modredc2ul2_equal (t, a_backup, m));
    modredc2ul2_clear (a_backup, m);
  }
#endif

  modredc2ul2_clear (t, m);
  modredc2ul2_clear (t2, m);

  return 1;
}


/* prototypes of non-inline functions */
int modredc2ul2_div3 (residueredc2ul2_t, const residueredc2ul2_t, 
		      const modulusredc2ul2_t);
int modredc2ul2_div5 (residueredc2ul2_t, const residueredc2ul2_t, 
		       const modulusredc2ul2_t);
int modredc2ul2_div7 (residueredc2ul2_t, const residueredc2ul2_t, 
		      const modulusredc2ul2_t);
int modredc2ul2_div11 (residueredc2ul2_t, const residueredc2ul2_t, 
		       const modulusredc2ul2_t);
int modredc2ul2_div13 (residueredc2ul2_t, const residueredc2ul2_t, 
		       const modulusredc2ul2_t);
void modredc2ul2_gcd (modintredc2ul2_t, const residueredc2ul2_t, 
		      const modulusredc2ul2_t);
void modredc2ul2_pow_ul (residueredc2ul2_t, const residueredc2ul2_t, 
			 const unsigned long, const modulusredc2ul2_t);
void modredc2ul2_2pow_ul (residueredc2ul2_t, const unsigned long, 
                          const modulusredc2ul2_t);
void modredc2ul2_pow_mp (residueredc2ul2_t, const residueredc2ul2_t, 
			 const unsigned long *, const int, 
			 const modulusredc2ul2_t);
void modredc2ul2_2pow_mp (residueredc2ul2_t, const unsigned long *, const int, 
			  const modulusredc2ul2_t);
void modredc2ul2_V_ul (residueredc2ul2_t, const residueredc2ul2_t, 
		       const unsigned long, const modulusredc2ul2_t);
void modredc2ul2_V_mp (residueredc2ul2_t, const residueredc2ul2_t, 
		       const unsigned long *, const int, 
		       const modulusredc2ul2_t);
int modredc2ul2_sprp (const residueredc2ul2_t, const modulusredc2ul2_t);
int modredc2ul2_sprp2 (const modulusredc2ul2_t);
int modredc2ul2_isprime (const modulusredc2ul2_t);
int modredc2ul2_inv (residueredc2ul2_t, const residueredc2ul2_t, 
		     const modulusredc2ul2_t);
int modredc2ul2_jacobi (const residueredc2ul2_t, const modulusredc2ul2_t);
#endif  /* MODREDC_2UL2_H */
