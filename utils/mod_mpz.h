#ifndef MOD_MPZ_H
#define MOD_MPZ_H

#include <assert.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <gmp.h>
#include "macros.h"

#ifndef ASSERT
#define ASSERT(x)	assert(x)
#endif

/* Even simple assertions are relatively expensive in very simple functions.
   If we want them anyway to hunt a bug, define WANT_ASSERT_EXPENSIVE */
#ifdef WANT_ASSERT_EXPENSIVE
#include <valgrind/memcheck.h>
#undef ASSERT_EXPENSIVE
#define SET_MPZ_UNDEF(x) VALGRIND_MAKE_MEM_UNDEFINED(x, sizeof(x[0]));
#define ASSERT_EXPENSIVE(x) ASSERT(x)
#define ASSERT_DEFINED(x) VALGRIND_CHECK_VALUE_IS_DEFINED(x); 
#define ASSERT_MPZ(x) ASSERT_DEFINED(x[0]); VALGRIND_CHECK_MEM_IS_DEFINED(x->_mp_d, llabs(x->_mp_size)); VALGRIND_CHECK_MEM_IS_ADDRESSABLE(x->_mp_d, llabs(x->_mp_alloc));
#define ASSERT_MPZ_RES(x) ASSERT_MPZ(x); ASSERT_EXPENSIVE (mpz_cmp(m, x) > 0);
#define ASSERT_UL_RES(x) ASSERT_DEFINED(x); ASSERT_EXPENSIVE (mpz_cmp_ui(m, x) > 0);
#else
#define SET_MPZ_UNDEF(x)
#define ASSERT_EXPENSIVE(x)
#define ASSERT_DEFINED(x)
#define ASSERT_MPZ(x)
#define ASSERT_MPZ_RES(x)
#define ASSERT_UL_RES(x)
#endif

#define MODMPZ_MAXBITS LONG_MAX

/* A macro for function renaming. All functions here start with modmpz_ */
#define MODMPZ_RENAME(x) modmpz_##x

typedef mpz_t residuempz_t;
typedef mpz_t modintmpz_t;
typedef mpz_t modulusmpz_t;


MAYBE_UNUSED
static inline void
modmpz_intinit (modintmpz_t r)
{
  mpz_init (r);
}


MAYBE_UNUSED
static inline void
modmpz_intclear (modintmpz_t r)
{
  ASSERT_MPZ(r);
  mpz_clear (r);
  SET_MPZ_UNDEF(r);
}


MAYBE_UNUSED
static inline void
modmpz_intset (modintmpz_t r, const modintmpz_t s)
{
  ASSERT_MPZ(s);
  mpz_set (r, s);
}


MAYBE_UNUSED
static inline void
modmpz_intset_ul (modintmpz_t r, const unsigned long s)
{
  ASSERT_DEFINED(s);
  mpz_set_ui (r, s);
}


MAYBE_UNUSED
static inline void
modmpz_intset_uls (modintmpz_t r, const unsigned long *s, const size_t n)
{
  mpz_import (r, n, -1, sizeof(*s), 0, 0, s);
}


/* Get the least significant unsigned long of r */
MAYBE_UNUSED
static inline unsigned long 
modmpz_intget_ul (const modintmpz_t s)
{
  ASSERT_MPZ(s);
  return mpz_get_ui(s);
}


MAYBE_UNUSED
static inline size_t  
modmpz_intget_uls (unsigned long *r, const modintmpz_t s)
{
  size_t count;
  ASSERT_MPZ(s);
  if (mpz_sgn(s) == 0) {
    r[0] = 0;
    return 1;
  }
  mpz_export (r, &count, -1, sizeof(*r), 0, 0, s);
  return count;
}

MAYBE_UNUSED
static inline double
modmpz_intget_double (const modintmpz_t s)
{
  ASSERT_MPZ(s);
  return mpz_get_d (s);
}

MAYBE_UNUSED
static inline int
modmpz_intequal (const modintmpz_t a, const modintmpz_t b)
{
  ASSERT_MPZ(a);
  ASSERT_MPZ(b);
  return (mpz_cmp(a, b) == 0);
}


MAYBE_UNUSED
static inline int
modmpz_intequal_ul (const modintmpz_t a, const unsigned long b)
{
  ASSERT_MPZ(a);
  ASSERT_DEFINED(b);
  return (mpz_cmp_ui (a, b) == 0);
}


MAYBE_UNUSED
static inline int
modmpz_intcmp (const modintmpz_t a, const modintmpz_t b)
{
  ASSERT_MPZ(a);
  ASSERT_MPZ(b);
  return (mpz_cmp(a, b));
}


MAYBE_UNUSED
static inline int
modmpz_intcmp_ul (const modintmpz_t a, const unsigned long b)
{
  ASSERT_MPZ(a);
  ASSERT_DEFINED(b);
  return (mpz_cmp_ui(a, b));
}


MAYBE_UNUSED
static inline int
modmpz_intcmp_uint64 (const modintmpz_t a, const uint64_t b)
{
  ASSERT(ULONG_MAX == UINT32_MAX || ULONG_MAX == UINT64_MAX);
  ASSERT_MPZ(a);
  ASSERT_DEFINED(b);
  if (ULONG_MAX == UINT32_MAX) {
    mpz_t t;
    int r;
    mpz_init(t);
    mpz_import (t, 1, -1, sizeof(uint64_t), 0, 0, &b);
    r = mpz_cmp (a, t);
    mpz_clear(t);
    return r;
  } else {
    return (mpz_cmp_ui(a, b));
  }
}


MAYBE_UNUSED
static inline int
modmpz_intfits_ul (const modintmpz_t a)
{
  ASSERT_MPZ(a);
  return (mpz_fits_ulong_p(a));
}


MAYBE_UNUSED
static inline void
modmpz_intadd (modintmpz_t r, const modintmpz_t a, const modintmpz_t b)
{
  ASSERT_MPZ(a);
  ASSERT_MPZ(b);
  mpz_add (r, a, b);
}


MAYBE_UNUSED
static inline void
modmpz_intsub (modintmpz_t r, const modintmpz_t a, const modintmpz_t b)
{
  ASSERT_MPZ(a);
  ASSERT_MPZ(b);
  mpz_sub (r, a, b);
}


MAYBE_UNUSED
static inline void
modmpz_intshr (modintmpz_t r, const modintmpz_t s, const int i)
{
  ASSERT_MPZ(s);
  ASSERT_DEFINED(i);
  mpz_tdiv_q_2exp (r, s, i);
}


MAYBE_UNUSED
static inline void
modmpz_intshl (modintmpz_t r, const modintmpz_t s, const int i)
{
  ASSERT_MPZ(s);
  ASSERT_DEFINED(i);
  mpz_mul_2exp (r, s, i);
}


/* Returns the number of bits in a, that is, floor(log_2(n))+1. 
   For n==0 returns 0. */
MAYBE_UNUSED
static inline size_t 
modmpz_intbits (const modintmpz_t a)
{
  ASSERT_MPZ(a);
  if (mpz_sgn(a) == 0)
    return 0;
  return (mpz_sizeinbase(a, 2));
}


/* r = n/d. We require d|n */
MAYBE_UNUSED
static inline void
modmpz_intdivexact (modintmpz_t r, const modintmpz_t n, const modintmpz_t d)
{
  ASSERT_MPZ(n);
  ASSERT_MPZ(d);
  ASSERT_EXPENSIVE(mpz_divisible_p(n, d));
  mpz_divexact (r, n, d);
}


/* r = n%d */
MAYBE_UNUSED
static inline void
modmpz_intmod (modintmpz_t r, const modintmpz_t n, 
              const modintmpz_t d)
{
  ASSERT_MPZ(n);
  ASSERT_MPZ(d);
  mpz_mod (r, n, d);
}


/* Functions for the modulus */

MAYBE_UNUSED
static inline void
modmpz_initmod_ul (modulusmpz_t m, const unsigned long s)
{
  ASSERT_DEFINED(s);
  mpz_init_set_ui (m, s);
}


MAYBE_UNUSED
static inline void
modmpz_initmod_int (modulusmpz_t m, const modintmpz_t s)
{
  ASSERT_MPZ(s);
  mpz_init_set (m, s);
}


MAYBE_UNUSED
static inline unsigned long
modmpz_getmod_ul (const modulusmpz_t m)
{
  ASSERT_MPZ(m);
  ASSERT(mpz_fits_ulong_p(m));
  return mpz_get_ui (m);
}


MAYBE_UNUSED
static inline void
modmpz_getmod_int (modintmpz_t r, const modulusmpz_t m)
{
  ASSERT_MPZ(m);
  mpz_set (r, m);
}


MAYBE_UNUSED
static inline void
modmpz_clearmod (modulusmpz_t m MAYBE_UNUSED)
{
  ASSERT_MPZ(m);
  mpz_clear(m);
  SET_MPZ_UNDEF(m);
}


/* Functions for residues */

/* Initialises a residue_t type and sets it to zero */
MAYBE_UNUSED
static inline void
modmpz_init (residuempz_t r, const modulusmpz_t m MAYBE_UNUSED)
{
  ASSERT_MPZ(m);
  mpz_init(r);
}


/* Initialises a residue_t type, but does not set it to zero. For fixed length
   residue_t types, that leaves nothing to do at all. */
MAYBE_UNUSED
static inline void
modmpz_init_noset0 (residuempz_t r MAYBE_UNUSED, 
		   const modulusmpz_t m MAYBE_UNUSED)
{
  ASSERT_MPZ(m);
  mpz_init(r);
}


MAYBE_UNUSED
static inline void
modmpz_clear (residuempz_t r MAYBE_UNUSED, 
             const modulusmpz_t m MAYBE_UNUSED)
{
  ASSERT_MPZ_RES(r);
  ASSERT_MPZ(m);
  mpz_clear (r);
  SET_MPZ_UNDEF(r);
}


MAYBE_UNUSED
static inline void
modmpz_set (residuempz_t r, const residuempz_t s, 
           const modulusmpz_t m MAYBE_UNUSED)
{
  ASSERT_MPZ_RES (s);
  ASSERT_MPZ(m);
  mpz_set (r, s);
}


MAYBE_UNUSED
static inline void
modmpz_set_ul (residuempz_t r, const unsigned long s, const modulusmpz_t m)
{
  ASSERT_MPZ(m);
  mpz_set_ui (r, s);
  mpz_mod (r, r, m);
}


/* Sets the residue_t to the class represented by the integer s. Assumes that
   s is reduced (mod m), i.e. 0 <= s < m */

MAYBE_UNUSED
static inline void
modmpz_set_ul_reduced (residuempz_t r, const unsigned long s, 
                      const modulusmpz_t m MAYBE_UNUSED)
{
  ASSERT_UL_RES (s);
  ASSERT_MPZ(m);
  mpz_set_ui (r, s);
}


MAYBE_UNUSED
static inline void
modmpz_set_int (residuempz_t r, const modintmpz_t s, 
		const modulusmpz_t m)
{
  ASSERT_MPZ(m);
  mpz_set (r, s);
  mpz_mod (r, r, m);
}


MAYBE_UNUSED
static inline void
modmpz_set_int_reduced (residuempz_t r, const modintmpz_t s, 
		       const modulusmpz_t m MAYBE_UNUSED)
{
  ASSERT_MPZ_RES(s);
  ASSERT_MPZ(m);
  mpz_set (r, s);
}


MAYBE_UNUSED 
static inline void 
modmpz_set0 (residuempz_t r, const modulusmpz_t m MAYBE_UNUSED) 
{ 
  ASSERT_MPZ(m);
  mpz_set_ui (r, 0);
}


MAYBE_UNUSED 
static inline void 
modmpz_set1 (residuempz_t r, const modulusmpz_t m MAYBE_UNUSED) 
{
  ASSERT_MPZ(m);
  if (mpz_cmp_ui (m, 1) == 0)
    mpz_set_ui (r, 0);
  else
    mpz_set_ui (r, 1);
}


/* Exchanges the values of the two arguments */
MAYBE_UNUSED
static inline void
modmpz_swap (residuempz_t a, residuempz_t b, 
            const modulusmpz_t m MAYBE_UNUSED)
{
  ASSERT_MPZ(a);
  ASSERT_MPZ(b);
  ASSERT_MPZ(m);
  mpz_swap(a, b);
}


MAYBE_UNUSED
static inline unsigned long
modmpz_get_ul (const residuempz_t s, const modulusmpz_t m MAYBE_UNUSED)
{
  ASSERT_MPZ_RES (s);
  ASSERT_MPZ(m);
  ASSERT_EXPENSIVE (mpz_fits_ulong_p(s));
  return mpz_get_ui (s);
}


MAYBE_UNUSED
static inline void
modmpz_get_int (modintmpz_t r, const residuempz_t s, 
		   const modulusmpz_t m MAYBE_UNUSED)
{
  ASSERT_MPZ_RES (s);
  ASSERT_MPZ(m);
  mpz_set (r, s);
}


MAYBE_UNUSED
static inline int
modmpz_equal (const residuempz_t a, const residuempz_t b, 
	      const modulusmpz_t m MAYBE_UNUSED)
{
  ASSERT_MPZ_RES (a);
  ASSERT_MPZ_RES (b);
  ASSERT_MPZ(m);
  return (mpz_cmp(a, b) == 0);
}


MAYBE_UNUSED
static inline int
modmpz_is0 (const residuempz_t a, const modulusmpz_t m MAYBE_UNUSED)
{
  ASSERT_MPZ_RES (a);
  ASSERT_MPZ(m);
  return (mpz_sgn(a) == 0);
}


MAYBE_UNUSED
static inline int
modmpz_is1 (const residuempz_t a, const modulusmpz_t m MAYBE_UNUSED)
{
  ASSERT_MPZ_RES (a);
  ASSERT_MPZ(m);
  return (mpz_cmp_ui(a, 1) == 0);
}


MAYBE_UNUSED
static inline void
modmpz_add (residuempz_t r, const residuempz_t a, const residuempz_t b, 
	   const modulusmpz_t m)
{
  ASSERT_MPZ_RES (a);
  ASSERT_MPZ_RES (b);
  ASSERT_MPZ(m);
  mpz_add (r, a, b);
  if (mpz_cmp(r, m) >= 0)
    mpz_sub (r, r, m);
}


MAYBE_UNUSED
static inline void
modmpz_add_ul (residuempz_t r, const residuempz_t a, const unsigned long b, 
	      const modulusmpz_t m)
{
  ASSERT_MPZ_RES (a);
  ASSERT_MPZ(m);
  mpz_add_ui (r, a, b);
  mpz_mod (r, r, m);
}


MAYBE_UNUSED
static inline void
modmpz_add1 (residuempz_t r, const residuempz_t a, const modulusmpz_t m)
{
  ASSERT_MPZ_RES (a);
  ASSERT_MPZ(m);
  modmpz_add_ul (r, a, 1, m);
}


MAYBE_UNUSED
static inline void
modmpz_sub (residuempz_t r, const residuempz_t a, const residuempz_t b, 
	    const modulusmpz_t m)
{
  ASSERT_MPZ_RES (a);
  ASSERT_MPZ_RES (b);
  ASSERT_MPZ(m);
  mpz_sub (r, a, b);
  if (mpz_sgn(r) < 0)
    mpz_add (r, r, m);
}


MAYBE_UNUSED
static inline void
modmpz_sub_ul (residuempz_t r, const residuempz_t a, const unsigned long b, 
	       const modulusmpz_t m)
{
  ASSERT_MPZ_RES (a);
  ASSERT_MPZ(m);
  mpz_sub_ui (r, a, b);
  mpz_mod (r, r, m);
}


MAYBE_UNUSED
static inline void
modmpz_neg (residuempz_t r, const residuempz_t a, const modulusmpz_t m)
{
  ASSERT_MPZ_RES (a);
  ASSERT_MPZ(m);
  if (mpz_sgn(a) == 0)
    mpz_set_ui (r, 0);
  else
    mpz_sub (r, m, a);
}


MAYBE_UNUSED
static inline void
modmpz_mul (residuempz_t r, const residuempz_t a, const residuempz_t b, 
           const modulusmpz_t m)
{
  ASSERT_MPZ_RES (a);
  ASSERT_MPZ_RES (b);
  ASSERT_MPZ (r);
  ASSERT_MPZ(m);
  mpz_mul (r, a, b);
  mpz_mod (r, r, m);
}


MAYBE_UNUSED
static inline void
modmpz_sqr (residuempz_t r, const residuempz_t a, const modulusmpz_t m)
{
  ASSERT_MPZ_RES (a);
  ASSERT_MPZ(m);
  mpz_mul (r, a, a);
  mpz_mod (r, r, m);
}


MAYBE_UNUSED
static inline void
modmpz_div2 (residuempz_t r, const residuempz_t a, const modulusmpz_t m)
{
  ASSERT_MPZ_RES (a);
  ASSERT_MPZ(m);
  ASSERT_EXPENSIVE (mpz_odd_p (m));
  if (mpz_even_p(a)) {
    mpz_tdiv_q_2exp (r, a, 1);
  } else {
    mpz_add (r, a, m);
    mpz_tdiv_q_2exp (r, r, 1);
  }
}


MAYBE_UNUSED
static inline int
modmpz_next (residuempz_t r, const modulusmpz_t m)
{
  ASSERT_MPZ_RES (r);
  ASSERT_MPZ(m);
  mpz_add_ui (r, r, 1);
  return (mpz_cmp(r, m) == 0);
}


MAYBE_UNUSED
static inline int
modmpz_finished (const residuempz_t a, const modulusmpz_t m)
{
  ASSERT_MPZ_RES (a);
  ASSERT_MPZ(m);
  return (mpz_cmp(a, m) == 0);
}


static inline int 
modmpz_div_ul (residuempz_t r, const residuempz_t a, const unsigned long b, 
              const modulusmpz_t m)
{
  mpz_t t;
  ASSERT_MPZ_RES (a);
  ASSERT_MPZ(m);
  mpz_init (t);
  mpz_set_ui (t, b);
  if (mpz_invert (t, t, m) == 0) {
    mpz_clear (t);
    return 0;
  }
  modmpz_mul (r, a, t, m);
  mpz_clear (t);
  return 1;
}

static inline int 
modmpz_div3 (residuempz_t r, const residuempz_t a, const modulusmpz_t m)
{
  return modmpz_div_ul(r, a, 3, m);
}

static inline int 
modmpz_div5 (residuempz_t r, const residuempz_t a, const modulusmpz_t m)
{
  return modmpz_div_ul(r, a, 5, m);
}

static inline int 
modmpz_div7 (residuempz_t r, const residuempz_t a, const modulusmpz_t m)
{
  return modmpz_div_ul(r, a, 7, m);
}

static inline int 
modmpz_div11 (residuempz_t r, const residuempz_t a, const modulusmpz_t m)
{
  return modmpz_div_ul(r, a, 11, m);
}

static inline int 
modmpz_div13 (residuempz_t r, const residuempz_t a, const modulusmpz_t m)
{
  return modmpz_div_ul(r, a, 13, m);
}

static inline void 
modmpz_gcd (modintmpz_t r, const residuempz_t a, const modulusmpz_t m)
{
  ASSERT_MPZ_RES (a);
  ASSERT_MPZ(m);
  mpz_gcd (r, a, m);
}

#define MOD_NO_SHARED_MOD_POW_UL 1
static inline void 
modmpz_pow_ul (residuempz_t r, const residuempz_t a, const unsigned long e, 
              const modulusmpz_t m)
{
  ASSERT_MPZ_RES (a);
  ASSERT_MPZ(m);
  mpz_powm_ui (r, a, e, m);
}

#define MOD_NO_SHARED_MOD_POW_MP 1
static inline void 
modmpz_2pow_ul (residuempz_t r, const unsigned long e, const modulusmpz_t m)
{
  mpz_set_ui (r, 2);
  ASSERT_MPZ(m);
  mpz_powm_ui (r, r, e, m);
}

static inline void 
modmpz_pow_mp (residuempz_t r, const residuempz_t a, const unsigned long *e,
              const int l, const modulusmpz_t m)
{
  mpz_t t;
  ASSERT_MPZ_RES (a);
  ASSERT_MPZ(m);
  mpz_init (t);
  mpz_import (t, l, -1, sizeof(unsigned long), 0, 0, e);
  mpz_powm (r, a, t, m);
  mpz_clear (t);
}

static inline void 
modmpz_2pow_mp (residuempz_t r, const unsigned long *e, const int l, 
               const modulusmpz_t m)
{
  ASSERT_MPZ(m);
  mpz_set_ui (r, 2);
  modmpz_pow_mp (r, r, e, l, m);
}

static inline int 
modmpz_sprp (const residuempz_t a MAYBE_UNUSED, const modulusmpz_t m)
{
  ASSERT_MPZ(m);
  /* FIXME */
  return mpz_probab_prime_p(m, 3);
}

static inline int 
modmpz_sprp2 (const modulusmpz_t m)
{
  ASSERT_MPZ(m);
  return mpz_probab_prime_p(m, 1);
}

static inline int 
modmpz_isprime (const modulusmpz_t m)
{
  ASSERT_MPZ(m);
  return mpz_probab_prime_p(m, 5);
}

static inline int 
modmpz_inv (residuempz_t r, const residuempz_t a, const modulusmpz_t m)
{
  ASSERT_MPZ(m);
  return mpz_invert (r, a, m);
}

static inline int 
modmpz_jacobi (const residuempz_t a, const modulusmpz_t m)
{
  ASSERT_MPZ(m);
  return mpz_jacobi(a, m);
}

/* Prototypes */
void 
modmpz_V_mp (residuempz_t, const residuempz_t, const unsigned long *,
	     int, const modulusmpz_t);
void 
modmpz_V_ul (residuempz_t, const residuempz_t, unsigned long, 
	     const modulusmpz_t);

int modmpz_batchinv (residuempz_t *, const residuempz_t *,
                     size_t n, const residuempz_t, const modulusmpz_t);

#endif  /* MOD_MPZ_H */
