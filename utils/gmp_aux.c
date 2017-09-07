/* auxiliary routines on GMP data-types */

#include "cado.h"
#include <stdint.h>
#include <string.h>
#include <limits.h>     /* for INT_MAX */
#include "gmp_aux.h"
#include "macros.h"

/* old versions of GMP do not provide mpn_neg (was mpn_neg_n) and mpn_xor_n */
#if !GMP_VERSION_ATLEAST(5,0,0)
mp_limb_t
mpn_neg (mp_limb_t *rp, const mp_limb_t *sp, mp_size_t n)
{
  mp_size_t i;

  for (i = 0; i < n; i++)
    rp[i] = ~sp[i];
  return mpn_add_1 (rp, rp, n, 1);
}

void
mpn_xor_n (mp_limb_t *rp, const mp_limb_t *s1p, const mp_limb_t *s2p,
	   mp_size_t n)
{
  mp_size_t i;

  for (i = 0; i < n; i++)
    rp[i] = s1p[i] ^ s2p[i];
}

void
mpn_zero(mp_limb_t *rp, mp_size_t n)
{
    memset(rp, 0, n * sizeof(mp_limb_t));
}
void
mpn_copyi(mp_limb_t *rp, const mp_limb_t *up, mp_size_t n)
{
    memmove(rp, up, n * sizeof(mp_limb_t));
}
void
mpn_copyd(mp_limb_t *rp, const mp_limb_t *up, mp_size_t n)
{
    memmove(rp, up, n * sizeof(mp_limb_t));
}
#endif

/* Set z to q. Warning: on 32-bit machines, we cannot use mpz_set_ui! */
void
mpz_set_uint64 (mpz_t z, uint64_t q)
{
  if (sizeof (unsigned long) == 8)
    mpz_set_ui (z, (unsigned long) q);
  else
    {
      ASSERT_ALWAYS (sizeof (unsigned long) == 4);
      mpz_set_ui (z, (unsigned long) (q >> 32));
      mpz_mul_2exp (z, z, 32);
      /* The & here should be optimized into a direct cast to a 32-bit
       * register in most cases (TODO: check) */
      mpz_add_ui (z, z, (unsigned long) (q & 4294967295UL));
    }
}

void
mpz_set_int64 (mpz_t z, int64_t q)
{
  if (sizeof (long) == 8)
    mpz_set_si (z, (long) q);
  else if (q >= 0)
    mpz_set_uint64 (z, q);
  else
    {
      mpz_set_uint64 (z, -q);
      mpz_neg (z, z);
    }
}

uint64_t
mpz_get_uint64 (mpz_srcptr z)
{
    uint64_t q;

    if (sizeof (unsigned long) == 8)
        q = mpz_get_ui (z);
    else
    {
        ASSERT_ALWAYS (sizeof (unsigned long) == 4);
        ASSERT_ALWAYS (sizeof (mp_limb_t) == 4);
        ASSERT_ALWAYS (GMP_LIMB_BITS == 32);
        q = mpz_get_ui (z); /* get the low word of z */
        q += ((uint64_t) mpz_getlimbn(z,1)) << 32;
    }
    return q;
}

int64_t
mpz_get_int64 (mpz_srcptr z)
{
    return mpz_get_uint64(z) * (int64_t) mpz_sgn(z);
}

int
mpz_fits_int64_p (mpz_srcptr z)
{
    int l = mpz_sizeinbase(z, 2);
    if (l <= 63) return 1;
    /* Also accept -2^63, which is INT64_MIN */
    if (mpz_sgn(z) < 0 && l == 64  && mpz_scan1(z, 0) == 63) return 1;
    return 0;
}

void
mpz_mul_uint64 (mpz_t a, mpz_srcptr b, uint64_t c)
{
  if (sizeof (unsigned long) >= sizeof (uint64_t))
    mpz_mul_ui (a, b, (unsigned long) c);
  else
    {
      mpz_t d;
      mpz_init (d);
      mpz_set_uint64 (d, c);
      mpz_mul (a, b, d);
      mpz_clear (d);
    }
}

void
mpz_addmul_uint64 (mpz_t a, mpz_srcptr b, uint64_t c)
{
  if (sizeof (unsigned long) >= sizeof (uint64_t))
    mpz_addmul_ui (a, b, (unsigned long) c);
  else
    {
      mpz_t d;
      mpz_init (d);
      mpz_set_uint64 (d, c);
      mpz_addmul (a, b, d);
      mpz_clear (d);
    }
}

void
mpz_submul_uint64 (mpz_t a, mpz_srcptr b, uint64_t c)
{
  if (sizeof (unsigned long) >= sizeof (uint64_t))
    mpz_submul_ui (a, b, (unsigned long) c);
  else
    {
      mpz_t d;
      mpz_init (d);
      mpz_set_uint64 (d, c);
      mpz_submul (a, b, d);
      mpz_clear (d);
    }
}

void
mpz_submul_int64 (mpz_t a, mpz_srcptr b, int64_t c)
{
  if (c >= 0)
    mpz_submul_uint64 (a, b, (uint64_t) c);
  else
    mpz_addmul_uint64 (a, b, (uint64_t) (-c));
}

void
mpz_divexact_uint64 (mpz_t a, mpz_srcptr b, uint64_t c)
{
  if (sizeof (unsigned long) >= sizeof (uint64_t))
    mpz_divexact_ui (a, b, (unsigned long) c);
  else
    {
      mpz_t d;
      mpz_init (d);
      mpz_set_uint64 (d, c);
      mpz_divexact (a, b, d);
      mpz_clear (d);
    }
}

int
mpz_divisible_uint64_p (mpz_t a, uint64_t c)
{
  if (sizeof (unsigned long) >= sizeof (uint64_t))
    return mpz_divisible_ui_p (a, (unsigned long) c);
  else
    {
      mpz_t d;
      int ret;

      mpz_init (d);
      mpz_set_uint64 (d, c);
      ret = mpz_divisible_p (a, d);
      mpz_clear (d);
      return ret;
    }
}

void
mpz_mul_int64 (mpz_t a, mpz_srcptr b, int64_t c)
{
  if (sizeof (long) == 8)
    mpz_mul_si (a, b, (long) c);
  else
    {
      mpz_t d;
      mpz_init (d);
      mpz_set_int64 (d, c);
      mpz_mul (a, b, d);
      mpz_clear (d);
    }
}

/* a <- a + b * c */
void
mpz_addmul_int64 (mpz_t a, mpz_srcptr b, int64_t c)
{
  if (sizeof (long) == 8)
    mpz_addmul_si (a, b, (long) c);
  else
    {
      mpz_t d;
      mpz_init (d);
      mpz_set_int64 (d, c);
      mpz_addmul (a, b, d);
      mpz_clear (d);
    }
}

/* returns the smallest prime p, q < p <= ULONG_MAX, or 0 if no such prime
   exists. guaranteed correct for q < 300M */
unsigned long
ulong_nextprime (unsigned long q)
{
  mpz_t p;
  unsigned long s[] = {16661633, 18790021, 54470491, 73705633, 187546133,
                       300164287, (unsigned long) (-1L)};
  int i;

  mpz_init (p);
  mpz_set_ui (p, q);
  do
    {
      mpz_nextprime (p, p);
      if (mpz_fits_ulong_p (p))
        q = mpz_get_ui (p);
      else {
        q = 0;
        break;
      }
      for (i = 0; q > s[i]; i++);
    }
  while (q == s[i]);
  mpz_clear (p);
  return q;
}

uint64_t
uint64_nextprime (uint64_t q)
{
  mpz_t z;

  mpz_init (z);
  mpz_set_uint64 (z, q);
  mpz_nextprime (z, z);
  q = mpz_get_uint64 (z);
  mpz_clear (z);
  return q;
}

#define REPS 3 /* number of Miller-Rabin tests in isprime */

/* For GMP 6.0.0:
   with REPS=1, the smallest composite reported prime is 1537381
   with REPS=2, it is 1943521
   with REPS=3, it is 465658903
   See also https://en.wikipedia.org/wiki/Miller–Rabin_primality_test
*/
int
ulong_isprime (unsigned long p)
{
  mpz_t P;
  int res;
  
  mpz_init_set_ui (P, p);
  res = mpz_probab_prime_p (P, REPS);
  mpz_clear (P);
  return res;
}

/* returns the smallest composite >= q with factors >= pmin */
unsigned long
ulong_nextcomposite (unsigned long q, unsigned long pmin)
{
  ASSERT_ALWAYS(q >= 2);
  for (;;q++)
    {
      if (ulong_isprime (q))
        continue;
      unsigned long p;
      for (p = 2; p < pmin && q % p; p += 1 + (p > 2));
      if (p >= pmin)
        break;
    }
  return q;
}

/* return the number of bits of p, counting from the least significant end */
int nbits (uintmax_t p)
{
  int k;

  for (k = 0; p != 0; p >>= 1, k ++);
  return k;
}

/* q <- n/d rounded to nearest, assuming d <> 0
   r <- n - q*d
   Output: -d/2 <= r < d/2
*/
void
mpz_ndiv_qr (mpz_t q, mpz_t r, mpz_t n, const mpz_t d)
{
  int s;

  ASSERT (mpz_cmp_ui (d, 0) != 0);
  mpz_fdiv_qr (q, r, n, d); /* round towards -inf, r has same sign as d */
  mpz_mul_2exp (r, r, 1);
  s = mpz_cmpabs (r, d);
  mpz_div_2exp (r, r, 1);
  if (s > 0) /* |r| > |d|/2 */
    {
      mpz_add_ui (q, q, 1);
      mpz_sub (r, r, d);
    }
}

void
mpz_ndiv_qr_ui (mpz_t q, mpz_t r, mpz_t n, unsigned long int d)
{
  int s;

  ASSERT (d != 0);
  mpz_fdiv_qr_ui (q, r, n, d); /* round towards -inf, r has same sign as d */
  mpz_mul_2exp (r, r, 1);
  s = mpz_cmpabs_ui (r, d);
  mpz_div_2exp (r, r, 1);
  if (s > 0) /* |r| > |d|/2 */
    {
      mpz_add_ui (q, q, 1);
      mpz_sub_ui (r, r, d);
    }
}

/* q <- n/d rounded to nearest, assuming d <> 0
   Output satisfies |n-q*d| <= |d|/2
*/
void
mpz_ndiv_q (mpz_t q, mpz_t n, const mpz_t d)
{
  mpz_t r;

  mpz_init (r);
  mpz_ndiv_qr (q, r, n, d);
  mpz_clear (r);
}

void
mpz_ndiv_q_ui (mpz_t q, mpz_t n, unsigned long int d)
{
  mpz_t r;

  mpz_init (r);
  mpz_ndiv_qr_ui (q, r, n, d);
  mpz_clear (r);
}

/* Return non-zero if a and b are coprime, else return 0 if gcd(a, b) != 1 */
int
mpz_coprime_p (mpz_t a, mpz_t b)
{
  mpz_t g;
  mpz_init(g);
  mpz_gcd (g, a, b);
  int ret = mpz_cmp_ui (g, 1);
  mpz_clear(g);
  return (ret == 0) ? 1 : 0;
}

long double
mpz_get_ld (mpz_t z)
{
  long double ld;
  double d;
  mpz_t t;

  d = mpz_get_d (z);
  mpz_init (t);
  mpz_set_d (t, d);
  mpz_sub (t, z, t);
  ld = (long double) d + (long double) mpz_get_d (t);
  mpz_clear (t);
  return ld;
}

/* returns the p-valuation of a, where p is expected to be prime. INT_MAX
 * is returned if a==0 */
int mpz_p_valuation(mpz_srcptr a, mpz_srcptr p)
{
    mpz_t c;
    mpz_init(c);
    int v = 0;
    if (mpz_size(a) == 0) return INT_MAX;
    mpz_set(c, a);
    for( ; mpz_divisible_p(c, p) ; v++)
        mpz_fdiv_q(c, c, p);
    mpz_clear(c);
    return v;
}

/* returns the p-valuation of a, where p is expected to be prime. INT_MAX
 * is returned if a==0 */
int mpz_p_valuation_ui(mpz_srcptr a, unsigned long p)
{
    mpz_t c;
    mpz_init(c);
    int v = 0;
    if (mpz_size(a) == 0) return INT_MAX;
    mpz_set(c, a);
    for( ; mpz_divisible_ui_p(c, p) ; v++)
        mpz_fdiv_q_ui(c, c, p);
    mpz_clear(c);
    return v;
}

