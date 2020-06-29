/* Some functions for modular arithmetic with residues and modulus in
   unsigned long variables. Residues are stored in Montgomery form,
   reduction after multiplication is done with REDC. Due to inlining,
   this file must be included in the caller's source code with #include */

/* Naming convention: all function start with modredcul, for
   MODulus REDC Unsigned Long, followed by underscore, functionality of
   function (add, mul, etc), and possibly underscore and specification of
   what argument types the function takes (_ul, etc). */

#ifndef MODREDC_UL_H
#define MODREDC_UL_H

/**********************************************************************/
#include "cado_config.h"  // just because we're a header.
#include <stddef.h>
#include <limits.h>
#include <stdint.h>
#include "macros.h"
#include "ularith.h"

/*********************************************************************/
/* Helper macros */

/* A macro for function renaming. All functions here start with modul_ */
#define MODREDCUL_RENAME(x) modredcul_##x

#define MODREDCUL_SIZE 1
#define MODREDCUL_MAXBITS LONG_BIT

typedef unsigned long residueredcul_t[MODREDCUL_SIZE];
typedef unsigned long modintredcul_t[MODREDCUL_SIZE];
struct __modulusredcul_t {
  unsigned long m;
  unsigned long invm;
  residueredcul_t one;
};
typedef struct __modulusredcul_t modulusredcul_t[1];


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
  ularith_redc(r, plow, phigh, m[0].m, m[0].invm);
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

#if (defined(__i386__) && defined(__GNUC__)) || defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
  {
    unsigned long t = a[0] - m[0].m, tr = a[0] + b[0];

    __asm__ __VOLATILE (
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
modredcul_intinit (modintredcul_t r)
{
  r[0] = 0;
}


MAYBE_UNUSED
static inline void
modredcul_intclear (modintredcul_t r MAYBE_UNUSED)
{
  return;
}


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
static inline void
modredcul_intset_uls (modintredcul_t r, const unsigned long *s, const size_t n)
{
  ASSERT_ALWAYS(n <= MODREDCUL_SIZE);
  if (n == 0)
    r[0] = 0;
  else
    r[0] = s[0];
}


MAYBE_UNUSED
static inline unsigned long
modredcul_intget_ul (const residueredcul_t s)
{
  return s[0];
}


MAYBE_UNUSED
static inline size_t  
modredcul_intget_uls (unsigned long *r, const residueredcul_t s)
{
  r[0] = s[0];
  return 1;
}


MAYBE_UNUSED
static inline double
modredcul_intget_double (const residueredcul_t s)
{
  return (double) s[0];
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
modredcul_intcmp_uint64 (const modintredcul_t a, const uint64_t b)
{
  if (b > ULONG_MAX)
    return -1;
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
static inline size_t 
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
      /* We want to compute 2^b % m, where b is the number of bits in an
         unsigned long. We know m is odd, so the remainder is not zero.
         Thus if we compute (2^b - 1) % m + 1, the +1 will not make the
         result equal to m and thus will produce the correct result while
         using only a single-word division. */
      m[0].one[0] = (-1UL) % m[0].m + 1;
    }
}

/* same as modredcul_initmod_ul, but does not compute m[0].one */
MAYBE_UNUSED
static inline void
modredcul_initmod_ul_raw (modulusredcul_t m, const unsigned long s)
{
  m[0].m = s;
  m[0].invm = -ularith_invmod (s);
}

MAYBE_UNUSED
static inline void
modredcul_initmod_int (modulusredcul_t m, const modintredcul_t s)
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
modredcul_getmod_int (modintredcul_t r, const modulusredcul_t m)
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


/* Puts in r the value of s * beta mod m, where beta is the word base.
   Note: s can be any unsigned long, in particular can be larger than m.
   When 0 <= s < m, use modredcul_set_ul_reduced for better efficiency. */
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
modredcul_set_int (residueredcul_t r, const modintredcul_t s,
		   const modulusredcul_t m)
{
  r[0] = s[0] % m[0].m;
  modredcul_tomontgomery (r, r, m);
}


MAYBE_UNUSED
static inline void
modredcul_set_int_reduced (residueredcul_t r, const modintredcul_t s,
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
modredcul_get_int (modintredcul_t r, const residueredcul_t s,
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
#ifdef MODTRACE
  printf ("modul_add: a = %lu, b = %lu", a[0], b[0]);
#endif

  ularith_addmod_ul_ul(r, a[0], b[0], m[0].m);

#ifdef MODTRACE
  printf (", r = %lu\n", r[0]);
#endif
}


MAYBE_UNUSED
static inline void
modredcul_add1 (residueredcul_t r, const residueredcul_t a,
                const modulusredcul_t m)
{
  modredcul_add(r, a, m[0].one, m);
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
#ifdef MODTRACE
  printf ("submod_ul: a = %lu, b = %lu", a[0], b[0]);
#endif

  ularith_submod_ul_ul(r, a[0], b[0], m[0].m);

#ifdef MODTRACE
  printf (", r = %lu\n", r[0]);
#endif
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
  printf ("(%lu * %lu / 2^%d) %% %lu", a[0], b[0], LONG_BIT, m[0].m);
#endif

  ularith_mul_ul_ul_2ul (&plow, &phigh, a[0], b[0]);
  modredcul_redc (r, plow, phigh, m);

#if defined(MODTRACE)
  printf (" == %lu /* PARI */ \n", r[0]);
#endif
}


/* For a residue class a (mod m) and non-negative integer b, set r to
   the smallest non-negative integer in the residue class a*b (mod m). */

MAYBE_UNUSED
static inline void
modredcul_mul_ul_ul (unsigned long *r, const residueredcul_t a,
                     const unsigned long b, const modulusredcul_t m)
{
  unsigned long plow, phigh;

  ASSERT_EXPENSIVE (m[0].m % 2 != 0);
  ASSERT_EXPENSIVE (a[0] < m[0].m);
#if defined(MODTRACE)
  printf ("(%lu * %lu / 2^%d) %% %lu", a[0], b, LONG_BIT, m[0].m);
#endif

  ularith_mul_ul_ul_2ul (&plow, &phigh, a[0], b);
  /* We have a <= m-1, b <= 2^LONG_BIT - 1. Thus the product
     phigh:plow <= (m-1)*(2^LONG_BIT - 1) = m*2^LONG_BIT - 2^LONG_BIT - m + 1,
     and with m >= 1,
     phigh:plow <= m*2^LONG_BIT - 2^LONG_BIT, so phigh < m. */
  modredcul_redc (r, plow, phigh, m);

#if defined(MODTRACE)
  printf (" == %lu /* PARI */ \n", *r);
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

#if defined(MODTRACE)
  printf ("(%lu^2 / 2^%d) %% %lu", a[0], LONG_BIT, m[0].m);
#endif

  ularith_sqr_ul_2ul (&plow, &phigh, a[0]);
  modredcul_redc (r, plow, phigh, m);

#if defined(MODTRACE)
  printf (" == %lu /* PARI */ \n", r[0]);
#endif
}


MAYBE_UNUSED
static inline void
modredcul_div2 (residueredcul_t r, const residueredcul_t a,
                const modulusredcul_t m)
{
  r[0] = ularith_div2mod(a[0], m[0].m);
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
#ifdef __cplusplus
extern "C" {
#endif

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
int modredcul_intinv (residueredcul_t, const residueredcul_t,
                      const modulusredcul_t);
int modredcul_batchinv (residueredcul_t *, const residueredcul_t *,
                        const size_t, const residueredcul_t,
                        const modulusredcul_t);
int modredcul_batchinv_ul (unsigned long *, const unsigned long *,
                           unsigned long, const size_t,
                           const modulusredcul_t);
int modredcul_batch_Q_to_Fp (unsigned long *, unsigned long,
                             unsigned long, unsigned long,
                             const unsigned long *, size_t);
int modredcul_jacobi (const residueredcul_t, const modulusredcul_t);

#ifdef __cplusplus
}
#endif

#endif  /* MODREDC_UL_H */
