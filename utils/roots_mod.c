/* Routines to solve x^k = a (mod p).

   Reference:
   [1] Adleman-Manders-Miller Root Extraction Method Revisited, Zhengjun Cao,
   Qian Sha, Xiao Fan, 2011, http://arxiv.org/abs/1111.4877.
*/

#include "cado.h"
#include <stdlib.h>
#include <stdint.h>
#include <gmp.h>

#include "gmp_aux.h"
#include "rootfinder.h"
#include "modredc_ul.h"
#include "mod_ul.h"
#include "modredc_ul_default.h"

static int 
mod_roots (residue_t *, residue_t, int, modulus_t);

/* For i < 50, isprime_table[i] == 1 iff i is prime */
static unsigned char isprime_table[] = {
0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 
0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 
0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 
0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0};

static const size_t isprime_table_size = 
    sizeof(isprime_table) / sizeof(isprime_table[0]);

/* return 1/a mod b */
static uint64_t
invert_mod (uint64_t a, uint64_t b)
{
  modulusul_t bb;
  residueul_t aa;
  uint64_t r;
  
  /* Inverse 1/a mod b */
  modul_initmod_ul (bb, b);
  modul_init (aa, bb);
  modul_set_ul (aa, a, bb);
  modul_inv (aa, aa, bb);

  r = modul_get_ul (aa, bb);

  modul_clear (aa, bb);
  modul_clearmod (bb);

  return r;
}

static inline int
isprime_ul (unsigned long n)
{
  int r;
  modulus_t nn;
  mod_initmod_ul (nn, n);
  r = mod_isprime (nn);
  mod_clearmod (nn);
  return r;
}


/* Factor n, put prime factors into p, their exponents into e */
unsigned char
factor_ul (unsigned long *p, unsigned char *e, const unsigned long n)
{
  unsigned long c = n;
  unsigned char ndiv = 0;
  
  if (c % 2 == 0) {
    p[ndiv] = 2;
    e[ndiv] = ularith_ctz(c);
    c >>= e[ndiv];
    ndiv++;
  }

  if (c % 3 == 0) {
    p[ndiv] = 3;
    e[ndiv] = 0;
    do {
      e[ndiv]++;
      c /= 3;
    } while (c % 3 == 0);
    ndiv++;
  }

  if (c % 5 == 0) {
    p[ndiv] = 5;
    e[ndiv] = 0;
    do {
      e[ndiv]++;
      c /= 5;
    } while (c % 5 == 0);
    ndiv++;
  }

  if (c >= 49 && !isprime_ul(c)) {
    unsigned long q, maxq = ularith_sqrt (c);
    for (q = 7; q <= maxq; q += 2) {
      if (c % q == 0) {
        p[ndiv] = q;
        e[ndiv] = 0;
        do {
          e[ndiv]++;
          c /= q;
        } while (c % q == 0);
        ndiv++;
        maxq = ularith_sqrt (c);
        if (q > maxq || isprime_ul(c))
          break;
      }
    }
  }

  if (c > 1) {
    p[ndiv] = c;
    e[ndiv] = 1;
    ndiv++;
  }

  return ndiv;
}


/* Initialise enumeration of all proper divisors of n */
void 
enumeratediv_init (enumeratediv_t *r, const unsigned long n)
{
  unsigned int i;

  r->ndiv = factor_ul (r->p, r->e, n);

  for (i = 0; i < r->ndiv; i++) {
    r->c[i] = 0;
  }
}


/* Iterate through the proper divisors of n. *Not* necessarily in increasing order.
   Return code of 0 means no more divisors. */
unsigned long 
enumeratediv (enumeratediv_t *r)
{
  unsigned int i = 0;
  unsigned long P = 1;
  
  while (i < r->ndiv && ++r->c[i] > r->e[i]) {
    r->c[i] = 0;
    i++;
  }
  if (i == r->ndiv)
    return 0;
  
  for (i = 0; i < r->ndiv; i++) {
    int j;
    for (j = 0; j < r->c[i]; j++)
      P *= r->p[i];
  }
  return P;
}


/* Return in o a k-th primitive root of unity modulo p, and in b a residue 
   that is a q-th non-power for each q|k. Assumes that k | p-1. */
void
omega (residue_t o, residue_t b, const unsigned long k, const modulus_t pp) 
{
  residue_t pow;
  const unsigned long p = mod_getmod_ul (pp);
  unsigned long pprime;
  unsigned long kdivq[15] = {}; /* k/q for each prime q | k */
  unsigned char exp[15] = {};
  unsigned long c = k;
  unsigned int i, ndiv;

  ASSERT (k > 0);
  if (k == 1) {
    mod_set1 (o, pp);
    return;
  }
  
  pprime = (p-1) / k;
  ASSERT (pprime * k == p-1);

  if (c % 2 == 0) {
    /* Not added to prime factor list, we use Jacobi symbol for check */
    c >>= ularith_ctz (c);
  }

  ndiv = factor_ul (kdivq, exp, c);
  for (i = 0; i < ndiv; i++)
    kdivq[i] = k / kdivq[i];

  mod_init (pow, pp);
  mod_set1 (b, pp);
  for (i = 2; i < p; i++) {
    unsigned int j;
    mod_add1 (b, b, pp);
    if ((ndiv == 0 || (ndiv == 1 && k % 2 == 1)) && !isprime_table[i])
      continue;
    if (k % 2 == 0 && mod_jacobi(b, pp) == 1)
      continue;
    mod_pow_ul (o, b, pprime, pp);
    for (j = 0; j < ndiv; j++) {
      mod_pow_ul (pow, o, kdivq[j], pp);
      if (mod_is1(pow, pp))
        break;
    }
    if (j == ndiv)
      break;
  }

#if 0
  /* Very slow but simple test */
  for (c = 1; c < k; c++) {
    mod_pow_ul (pow, o, c, pp);
    ASSERT_ALWAYS (!mod_is1(pow, pp));
    mod_pow_ul (pow, o, k, pp);
    ASSERT_ALWAYS (mod_is1(pow, pp));
  }
#endif
  mod_clear (pow, pp);
}


/************************** square roots *************************************/

/* Uses Tonelli-Shanks, more precisely Algorithm from Table 1 in [1] */
static int
tonelli_shanks(residue_t *rr, int d, uint64_t n, const modulus_t pp)
{
  uint64_t q, s, i, j, k = 0, l;
  residue_t aa, hh, delta, dd, bb, zz;
  const uint64_t p =  mod_getmod_ul(pp);

  /* write p-1 = q*2^s with q odd */
  for (q = p-1, s = 0; (q&1) == 0; q/=2, s++);
  
  mod_init (zz, pp);
  mod_set1 (zz, pp);
  i = 1;
  do {
    /* zz is equal to i (mod pp) */
    mod_add1 (zz, zz, pp);
    i++;
  } while (i < isprime_table_size && (!isprime_table[i] || mod_jacobi (zz, pp) != -1));
  ASSERT_ALWAYS(i < isprime_table_size);
  mod_init (hh, pp);
  mod_init (aa, pp);
  mod_init (delta, pp);
  mod_init (dd, pp);
  mod_init (bb, pp);
  for (i = 0; i < n; i++)
    {
      /* For d divisible by 4, we have to check the Legendre symbol again.
         Consider for example x^4 = 10 (mod 13): x^2 = 10 (mod 13) has two
         roots (6 and 7) but none of them has a root mod 13. */
      if ((d & 3) == 0 && mod_jacobi (rr[d/2+i], pp) != 1)
        continue;

      mod_pow_ul (aa, zz, q, pp);
      mod_pow_ul (bb, rr[d/2+i], q, pp);
      mod_set1 (hh, pp);
      for (j = 1; j < s; j++)
        {
          mod_set (dd, bb, pp);
          for (l = 0; l < s - 1 - j; l++)
            mod_sqr (dd, dd, pp);
          if (!mod_is1 (dd, pp))
            {
              mod_mul (hh, hh, aa, pp);
              mod_sqr (aa, aa, pp);
              mod_mul (bb, bb, aa, pp);
            }
          else
            mod_sqr (aa, aa, pp);
        }
      mod_pow_ul (delta, rr[d/2+i], (q + 1) >> 1, pp);
      mod_mul (hh, hh, delta, pp);
      mod_set (rr[2*k], hh, pp);
      mod_neg (rr[2*k+1], rr[2*k], pp);
      k++;
    }

  mod_clear (hh, pp);
  mod_clear (aa, pp);
  mod_clear (delta, pp);
  mod_clear (dd, pp);
  mod_clear (bb, pp);

  return k;
}


static void
one_root2_V (residue_t rr, const residue_t aa, const modulus_t pp)
{
  residue_t xx, two, discr;
  unsigned long c;
  const unsigned long p = mod_getmod_ul (pp);

  mod_init (xx, pp);
  mod_init (discr, pp);
  mod_init (two, pp);
  mod_set1 (two, pp);
  mod_add1 (two, two, pp); /* two = 2 */
  ASSERT_ALWAYS (p % 4 == 1);

  mod_sub (xx, aa, two, pp); /* x = a - 2 */
  mod_sqr (discr, xx, pp);
  mod_sub (discr, discr, two, pp); /* discr = x^2 - 2 */
  mod_sub (discr, discr, two, pp); /* discr = x^2 - 4 */

  if (mod_jacobi(discr, pp) == -1) {
    mod_V_ul (rr, xx, (p+3)/4, pp);
  } else {
    /* Find a suitable c value */
    residue_t cc;
    mod_init (cc, pp);
    mod_set1 (cc, pp);
    for (c = 2; c < p; c++) {
      mod_add1 (cc, cc, pp); /* cc = c */
      mod_sqr (xx, cc, pp);
      mod_mul (xx, xx, aa, pp);
      mod_sub (xx, xx, two, pp);
      mod_sqr (discr, xx, pp);
      mod_sub (discr, discr, two, pp); /* discr = x^2 - 2 */
      mod_sub (discr, discr, two, pp); /* discr = x^2 - 4 */
      if (mod_jacobi(discr, pp) == -1)
        break;
    }
    ASSERT_ALWAYS(c < p);
    mod_V_ul (xx, xx, (p+3)/4, pp);
    /* Divide out c, using hard-coded functions for small cases */
    while (c % 2 == 0) {
      mod_div2 (xx, xx, pp);
      mod_div2 (cc, cc, pp);
      c /= 2;
    }
    while (c % 3 == 0) {
      c /= 3;
      mod_div3 (xx, xx, pp);
      mod_div3 (cc, cc, pp);
    }
    while (c % 5 == 0) {
      mod_div5 (xx, xx, pp);
      mod_div5 (cc, cc, pp);
      c /= 5;
    }
    while (c % 11 == 0) {
      mod_div11 (xx, xx, pp);
      mod_div11 (cc, cc, pp);
      c /= 11;
    }
    if (c > 1) {
      /* About 1/2^11 + 1/2^13 + 1/2^17 + ... ~= 0.06% of cases get here */
      mod_inv (cc, cc, pp);
      mod_mul (xx, xx, cc, pp);
    }
    mod_set (rr, xx, pp);
    mod_clear (cc, pp);
  }
  
  mod_clear (xx, pp);
  mod_clear (discr, pp);
  mod_clear (two, pp);
}

static int
roots2_V (residue_t *rr, int d, uint64_t n, const modulus_t pp)
{
  uint64_t i, k = 0;

  for (i = 0; i < n; i++)
    {
      if ((d & 3) == 0 && mod_jacobi (rr[d/2+i], pp) != 1)
        continue;
      one_root2_V (rr[2*k], rr[d/2+i], pp);
      mod_neg (rr[2*k + 1], rr[2*k], pp);
      k++;
    }
  return k;
}


/* Roots of x^d = a (mod p) for d even and a 64-bit word */
static int
roots2 (residue_t *rr, residue_t aa, int d, modulus_t pp)
{
  uint64_t n, i, k = 0;
  const uint64_t p =  mod_getmod_ul(pp);

  if (mod_jacobi (aa, pp) != 1)
    return 0;

  /* find the roots of x^(d/2) = a (mod p) */
  n = mod_roots (rr + d / 2, aa, d / 2, pp);

  if (n == 0)
    return n;

  if (p % 4 == 3)
    {
      /* solutions are +/-a^((p+1)/4) */
      for (i = 0; i < n; i++)
        {
          /* If d = 2 (mod 4), then every root of x^(d/2) = a (mod p) gives
             two roots of x^d = a (mod p) since a is a quadratic residue.
             However if d = 0 (mod 4), then a root of x^(d/2) = a (mod p) can
             give no root of x^d = a (mod p), thus we have to check the
             Legendre symbol again.
             Consider for example x^4 = 3 (mod 11): x^2 = 3 mod 11 has
             two roots (5 and 6), then x^2-5 mod 11 has two roots (4 and 7)
             but x^2-6 mod 11 has no root. */
          
          if ((d & 3) == 0 && mod_jacobi (rr[d/2+i], pp) != 1)
            continue;

          mod_pow_ul (rr[2*k], rr[d/2+i], (p + 1) >> 2, pp);
          mod_neg (rr[2*k+1], rr[2*k], pp);
          k ++;
        }
      return 2*k;
    }

  /* case p = 1 (mod 4). */ 
  if (0)
    k = tonelli_shanks (rr, d, n, pp);
  else
    k = roots2_V (rr, d, n, pp);

  return 2*k;
}

/************************** cubic roots **************************************/

static void 
one_cubic_root_2mod3 (residue_t rr, residue_t ddelta, modulus_t pp)
{
  /* when p = 2 (mod 3), then 1/3 = (2p-1)/3 mod (p-1), thus a cubic root
     is delta^((2p-1)/3) mod p. We rewrite exponent as (p+1)/3*2-1 to 
     avoid overflow */

  mod_pow_ul (rr, ddelta, (mod_getmod_ul(pp) + 1) / 3 * 2 - 1, pp);
}


/* Put in rr a root of x^3 = ddelta (mod pp), assuming one exists,
   and in zeta a primitive 3-rd root of unity.
   Use the algorithm from Table 3 in reference [1].
   rr and ddelta overlapping is permissible. */
static void 
one_cubic_root_1mod3 (residue_t rr, residue_t zeta, residue_t ddelta, modulus_t pp)
{
  residue_t rho, a, aprime, b, h, d;
  uint64_t i, j, s, t, smod3, l;
  const uint64_t p = mod_getmod_ul (pp);

  s = (p - 1) / 3;
  t = 1;
  while ((s % 3) == 0)
    s /= 3, t++;
  mod_init (rho, pp);
  mod_init (a, pp);
  mod_init (aprime, pp);
  mod_init (b, pp);
  mod_init (h, pp);
  mod_init (d, pp);
  
  mod_set1 (rho, pp);
  smod3 = s % 3;
  for (i = 2; i < p; i++)
    {
      mod_add1 (rho, rho, pp);
      /* rho is equal to i (mod pp) */
      /* no point in testing i = k*l where we know both k and l are cubic 
         residues, so we do only primes */
      ASSERT_ALWAYS(i < isprime_table_size);
      if (!isprime_table[i])
        continue;
      mod_pow_ul (a, rho, s, pp); /* a = rho^s */
      mod_set (aprime, a, pp);
      for (j = 0; j < t - 1; j++) {
        mod_sqr (b, aprime, pp); /* use b as temp */
        mod_mul (aprime, aprime, b, pp);
      }
      /* aprime = rho^(3^(t-1)*s) = rho^((p-1)/3) */
      if (!mod_is1 (aprime, pp))
        break;
    }
  /* aprime  =  rho^((p-1)/3)  !=  1, so it is a 3-rd rood of 1 */
  mod_set (zeta, aprime, pp);
  /* see below to explain why we start from delta^(2s) for s = 1 mod 3 */
  mod_pow_ul (b, ddelta, (smod3 == 1) ? 2*s : s, pp);
  mod_set1 (h, pp);
  for (i = 1; i < t ; i++)
    {
      mod_set (d, b, pp);
      for (j = 0; j < t - 1 - i; j++)
        mod_pow_ul (d, d, 3, pp);
      if (mod_is1 (d, pp))
        {
          mod_pow_ul (a, a, 3, pp);
        }
      else if (mod_equal (d, aprime, pp))
        {
          mod_sqr (d, a, pp);
          mod_mul (h, h, d, pp);
          mod_mul (a, d, a, pp);
          mod_sqr (d, a, pp);
          mod_mul (b, b, d, pp);
        }
      else
        {
          mod_mul (h, h, a, pp);
          mod_pow_ul (a, a, 3, pp);
          mod_mul (b, b, a, pp);
        }
    }
  /* in the case s = 3l+1, instead of computing r = delta^l*h and then
     inverting r mod p as in Table 3, we use the generic algorithm from
     Table 4, which computes r as delta^alpha*h where alpha = 1/3 mod s,
     i.e., alpha = 2l+1. However this needs to start from b = delta^(2s). */
  l = (s + 1) / 3;
  mod_pow_ul (d, ddelta, (smod3 == 1) ? 2*l+1 : l, pp);
  mod_mul (h, h, d, pp);

  mod_set (rr, h, pp);

  mod_clear (rho, pp);
  mod_clear (a, pp);
  mod_clear (aprime, pp);
  mod_clear (b, pp);
  mod_clear (h, pp);
  mod_clear (d, pp);
}

/* return 1 iff a is a cube mod p */
static int
is_cube (residue_t aa, modulus_t pp)
{
  if ((mod_getmod_ul (pp) % 3) == 1)
    {
      residue_t cc;
      int r;
      mod_init_noset0 (cc, pp);
      mod_pow_ul (cc, aa, (mod_getmod_ul(pp) - 1) / 3, pp);
      r = mod_is1 (cc, pp);
      mod_clear (cc, pp);
      return r;
    }
  else /* p = 2 (mod 3) */
    return 1;
}


/* Roots of x^d = a (mod p) for d divisible by 3 and a 64-bit word */
static int
roots3 (residue_t *rr, residue_t aa, int d, modulus_t pp)
{
  uint64_t i, n, k;
  const uint64_t p = mod_getmod_ul(pp);

  if (!is_cube (aa, pp))
    return 0;

  /* find the roots of x^(d/3) = a (mod p) */
  n = mod_roots (rr + 2 * (d / 3), aa, d / 3, pp);

  if (n == 0)
    return n;

  if ((p % 3) == 1)
    {
      residue_t zeta, t;

      mod_init (zeta, pp);
      mod_init (t, pp);
      mod_set1 (t, pp);
      i = 1;

      for (i = k = 0; i < n; i++)
        {
          /* Note: if d is divisible by 9, then we must check again if each
             root of x^(d/3) = a (mod p) is a cube. For example for d = 9,
             a = 9943082 and p = 20000047, we have three roots of x^3 = a,
             namely 10169532, 11660661 and 18169901, but only 11660661 is a
             cube mod p */
          if ((d % 9) == 0 && is_cube (rr[2 * (d / 3) + i], pp) == 0)
            continue;
          one_cubic_root_1mod3 (rr[k], zeta, rr[2 * (d / 3) + i], pp);
          /* zeta is a cubic root of 1 */
          mod_mul (rr[k+1], rr[3*i], zeta, pp);
          mod_mul (rr[k+2], rr[3*i+1], zeta, pp);
          k += 3;
        }
      mod_clear (zeta, pp);
      mod_clear (t, pp);
      return k;
    }
  else /* p = 2 (mod 3): exactly one root each */
    {
      for (i = 0; i < n; i++)
        one_cubic_root_2mod3 (rr[i], rr[2 * (d / 3) + i], pp);
      return n;
    }
}

/************************** r-th roots ***************************************/

/* Put in rop a r-th root of delta (mod p), assuming one exists,
   using the algorithm from Table 4 in reference [1]. */
static void
one_rth_root (residue_t rop, uint64_t r, residue_t delta, modulus_t pp)
{
  residue_t a, b, c, d, h;
  const uint64_t p = mod_getmod_ul (pp);
  uint64_t i, j, s, t, alpha, rt;

  /* when p-1 is not divisible by r, there is exactly one root */
  if (((p - 1) % r) != 0)
    {
      for (i = 1; (i * (p - 1) + 1) % r != 0; i++);
      mod_pow_ul (rop, delta, (i * (p - 1) + 1) / r, pp);
      return;
    }

  /* now p-1 is divisible by r */

  mod_init (a, pp);
  mod_init (b, pp);
  mod_init (c, pp);
  mod_init (d, pp);
  mod_init (h, pp);
  /* Write s * r^t = p-1, t maximal, and rt = r^t */
  for (s = (p - 1) / r, t = 1, rt = r; (s % r) == 0; s /= r, t++, rt *= r);

#if 1
  omega(c, a, rt, pp);
  mod_pow_ul (a, c, rt/r, pp);
#else  
  for (i = 2; i < p; i++)
    {
      /* c = i (mod p) */
      mod_set_ul (c, i, pp);
      mod_pow_ul (c, c, s, pp);
      mod_set (a, c, pp);
      for (j = 0; j < t - 1; j++)
        mod_pow_ul (a, a, r, pp);
      if (mod_is1 (a, pp) != 1) {
        /* Now c is an r^t-th primitive root of unity, and a is an r-th 
           primitive root of unity */
        break;
      }
    }
#endif
  alpha = invert_mod (s, r);
  /* we have alpha = 1/s mod r, thus alpha*s = 1 + beta*r,
     where beta = (alpha*s - 1)/r, and 1/r = -beta mod s */
  alpha = alpha * s - 1;
  ASSERT(alpha % r == 0);
  /* since alpha is a modular inverse, it cannot be 0 */
  alpha = s - (alpha / r);
  mod_pow_ul (b, delta, r * alpha - 1, pp);
  mod_set1 (h, pp);
  for (i = 1; i < t; i++)
    {
      mod_set (d, b, pp);
      for (j = 0; j < t - 1 - i; j++)
        mod_pow_ul (d, d, r, pp);
      if (mod_is1 (d, pp))
        mod_pow_ul (c, c, r, pp);
      else
        {
          for (j = 0; mod_is1 (d, pp) == 0; mod_mul (d, d, a, pp), j++);
          mod_pow_ul (d, c, j, pp);
          mod_mul (h, h, d, pp);
          mod_pow_ul (c, c, r, pp);
          mod_pow_ul (d, c, j, pp);
          mod_mul (b, b, d, pp);
        }
    }
  mod_pow_ul (rop, delta, alpha, pp);
  mod_mul (rop, rop, h, pp);
  mod_clear (a, pp);
  mod_clear (b, pp);
  mod_clear (c, pp);
  mod_clear (d, pp);
  mod_clear (h, pp);
}

static int
is_rth_power (residue_t a, uint64_t r, modulus_t pp)
{
  const uint64_t p = mod_getmod_ul (pp);
  residue_t t;
  int ret = 1;

  if ((p % r) == 1)
    {
      mod_init (t, pp);
      mod_pow_ul (t, a, (p - 1) / r, pp);
      ret = mod_is1 (t, pp);
      mod_clear (t, pp);
    }
  return ret;
}

/* Roots of x^d = a (mod p), assuming d is not divisible by 2 nor 3.
   (This code works in fact for d divisible by 3 too, if one starts by r = 3
   in the loop below.) */
static int
roots (residue_t *rr, residue_t a, int d, modulus_t pp)
{
  uint64_t r, n, i, j, k;
  const uint64_t p = mod_getmod_ul (pp);
  residue_t z, *rr0;

  /* first find the smallest prime r dividing d (r can be d) */
  for (r = 5; d % r; r += 2);

  if (is_rth_power (a, r, pp) == 0)
    return 0;

  /* find the roots of x^(d/r) = a (mod p) */
  n = mod_roots (rr0 = rr + d - d / r, a, d / r, pp);

  /* For the roots of x^r = a (mod p):
     either p = 1 (mod r), then
     (a) if a^((p-1)/r) <> 1, there is no solution;
     (b) if a^((p-1)/r) == 1, there are r solutions. If we have one of them,
         the others are x*zeta^i for 1 <= i < r, with zeta^((p-1)/r) <> 1;
     or p <> 1 (mod r), then there is only one solution, which is given by
     a^(e/r) where e = i*(p-1)+1 is divisible by r. */

  if ((p % r) != 1)
    {
      uint64_t e;

      for (e = p; e % r; e += p - 1);
      for (i = 0; i < n; i++)
        mod_pow_ul (rr[i], rr0[i], e / r, pp);
      return n;
    }

  /* now p = 1 (mod r) */
  k = 0;
  /* get a primitive r-th root of unity z */
  mod_init (z, pp);
  for (i = 2; i < p; i++)
    {
      mod_set_ul (z, i, pp);
      mod_pow_ul (z, z, (p - 1) / r, pp);
      if (mod_is1 (z, pp) == 0)
        break;
    }
  for (i = k = 0; i < n; i++)
    {
      /* check rr0[i] is a r-th power */
      mod_pow_ul (rr[k], rr0[i], (p - 1) / r, pp);
      if (mod_is1 (rr[k], pp) == 0)
        continue;
      /* get one r-th root */
      one_rth_root (rr[k], r, rr0[i], pp);
      for (j = 1, k++; j < r; j++, k++)
        mod_mul (rr[k], rr[k-1], z, pp);
    }
  mod_clear (z, pp);
  return k;
}

/*****************************************************************************/

static int 
mod_roots (residue_t *rr, residue_t aa, int d, modulus_t pp)
{
  if (d == 1)
    {
      mod_set (rr[0], aa, pp);
      return 1;
    }
  else if ((d & 1) == 0) /* d is even */
    {
      return roots2 (rr, aa, d, pp);
    }
  else if (d % 3 == 0)
    {
      return roots3 (rr, aa, d, pp);
    }
  else
    return roots (rr, aa, d, pp);
}


/* sort the roots r[0], ..., r[n-1] in increasing order */
static void
sort_roots (uint64_t *r, int n)
{
  int i, j;
  uint64_t t;

  for (i = 1; i < n; i++)
    {
      t = r[i];
      for (j = i; j > 0 && r[j-1] > t; j--)
        r[j] = r[j-1];
      r[j] = t;
    }
}

#define MAX_DEGREE 10

/* put in r[0], r[1], ... the roots of x^d = a (mod p),
   and return the number of roots.
   Assumes 0 <= a < p.
*/
int
roots_mod_uint64 (uint64_t *r, uint64_t a, int d, uint64_t p)
{
  int n = -1, i;

  ASSERT_ALWAYS(d <= MAX_DEGREE);

  if (d == 1)
    {
      r[0] = a;
      return 1;
    }

  if (sizeof (unsigned long) == 8)
    {
      modulus_t pp;
      residue_t aa, rr[MAX_DEGREE];

      mod_initmod_ul (pp, p);
      mod_init (aa, pp);
      mod_set_ul (aa, a, pp);
      for (i = 0; i < d; i++)
        mod_init (rr[i], pp);

      n = mod_roots (rr, aa, d, pp);

      for (i = 0; i < n; i++) {
        r[i] = mod_get_ul (rr[i], pp);
#ifndef NDEBUG
        /* Check that it's a d-th root of a */
        mod_pow_ul (aa, rr[i], d, pp);
        ASSERT (mod_get_ul (aa, pp) == a);
#endif
      }
      sort_roots (r, n);
#ifndef NDEBUG
      for (i = 1; i < n; i++) {
        /* Check for dupes */
        if (r[i-1] >= r[i]) {
          fprintf (stderr, "%lu^(1/%d) (mod %lu), r[%d]: %lu >= r[%d]: %lu\n", 
                   a, d, p, i-1, r[i-1], i, r[i]);
          ASSERT(r[i-1] < r[i]);
        }
      }
#endif

      for (i = 0; i < d; i++)
        mod_clear (rr[i], pp);
      mod_clear (aa, pp);
      mod_clearmod (pp);
    }
  else
    {
      mpz_t *f;
      f = (mpz_t*) malloc ((d + 1) * sizeof (mpz_t));
      for (i = 0; i <= d; i++)
        mpz_init (f[i]);
      mpz_set_ui (f[d], 1);
      mpz_set_uint64 (f[0], p - a);
      n = poly_roots_uint64 (r, f, d, p);
      for (i = 0; i <= d; i++)
        mpz_clear (f[i]);
      free (f);
      sort_roots (r, n);
    }

  return n;
}
