/* auxiliary routines on GMP data-types */

#include "cado.h" // IWYU pragma: keep

#include <stddef.h>
#include <stdbool.h>
#include <stdint.h>
#include <limits.h>     /* for INT_MAX */
#include <math.h>

#include <gmp.h>

#include "gmp_aux.h"
#include "getprime.h"  // for getprime_mt, prime_info_clear, prime_info_init
#include "macros.h"    // for ASSERT_ALWAYS // IWYU pragma: keep

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

#if !GMP_VERSION_ATLEAST(6,1,0)
int mpn_zero_p(const mp_limb_t *rp, mp_size_t n)
{
    for( ; n-- ; ) if (rp[n]) return 0;
    return 1;
}
#endif

#if ULONG_BITS < 64
/* Set z to q. Warning: on 32-bit machines, we cannot use mpz_set_ui! */
void
mpz_set_uint64 (mpz_ptr z, uint64_t q)
{
  if (q <= ULONG_MAX)
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
mpz_set_int64 (mpz_ptr z, int64_t q)
{
  if (LONG_MIN <= q && q <= LONG_MAX)
    mpz_set_si (z, (long) q);
  else if (q >= 0)
    mpz_set_uint64 (z, (uint64_t) q);
  else
    {
      mpz_set_uint64 (z, -(uint64_t)q);
      mpz_neg (z, z);
    }
}

void mpz_init_set_uint64 (mpz_ptr z, uint64_t x)
{
    mpz_init2(z, 64);
    mpz_set_uint64(z, x);
}

void mpz_init_set_int64 (mpz_ptr z, int64_t x)
{
    mpz_init2(z, 64);
    mpz_set_int64(z, x);
}

uint64_t
mpz_get_uint64 (mpz_srcptr z)
{
    return mpz_getlimbn_uint64(z, 0);
}

uint64_t
mpz_getlimbn_uint64 (mpz_srcptr z, unsigned int i)
{
    uint64_t q;

    if (sizeof (unsigned long) == 8)
        q = mpz_getlimbn (z, i);
    else
    {
        ASSERT_ALWAYS (sizeof (unsigned long) == 4);
        ASSERT_ALWAYS (sizeof (mp_limb_t) == 4);
        ASSERT_ALWAYS (GMP_LIMB_BITS == 32);
        q = mpz_getlimbn (z, 2*i); /* get the low word of z */
        q += ((uint64_t) mpz_getlimbn(z, 2 * i + 1)) << 32;
    }
    return q;
}

int
mpz_fits_uint64_p (mpz_srcptr z)
{
    if (mpz_sgn(z) < 0)
        return 0;
    size_t l = mpz_sizeinbase(z, 2);
    return (l <= 64);
}

int64_t
mpz_get_int64 (mpz_srcptr z)
{
    const uint64_t l = mpz_get_uint64(z);
    if (mpz_sgn(z) >= 0)
        return l & (uint64_t) INT64_MAX;
    /* We want 1 -> -1, 2 -> -2, ...,
       2^63-1 -> -2^63+1, 2^63 -> -2^63, 2^63+1 -> -1 */
    return (int64_t) -(((l - 1) & (uint64_t) INT64_MAX) + 1);
}

int
mpz_fits_sint64_p (mpz_srcptr z)
{
    size_t l = mpz_sizeinbase(z, 2);
    if (l <= 63) return 1;
    /* Also accept -2^63, which is INT64_MIN */
    if (mpz_sgn(z) < 0 && l == 64  && mpz_scan1(z, 0) == 63) return 1;
    return 0;
}

void
mpz_mul_uint64 (mpz_ptr a, mpz_srcptr b, uint64_t c)
{
  if (c <= ULONG_MAX)
    mpz_mul_ui (a, b, (unsigned long) c);
  else
    {
      mpz_t d;
      mpz_init_set_uint64 (d, c);
      mpz_mul (a, b, d);
      mpz_clear (d);
    }
}

int
mpz_cmp_uint64 (mpz_srcptr a, uint64_t c)
{
  if (c <= ULONG_MAX)
    return mpz_cmp_ui (a, (unsigned long) c);
  else
    {
      mpz_t d;
      mpz_init_set_uint64 (d, c);
      int r = mpz_cmp (a, d);
      mpz_clear (d);
      return r;
    }
}

int
mpz_cmp_int64 (mpz_srcptr a, int64_t c)
{
  if (LONG_MIN <= c && c <= LONG_MAX)
    return mpz_cmp_si (a, (long) c);
  else
    {
      mpz_t d;
      mpz_init_set_int64 (d, c);
      int r = mpz_cmp (a, d);
      mpz_clear (d);
      return r;
    }
}

void
mpz_add_uint64 (mpz_ptr a, mpz_srcptr b, uint64_t c)
{
  if (c <= ULONG_MAX)
    mpz_add_ui (a, b, (unsigned long) c);
  else
    {
      mpz_t d;
      mpz_init_set_uint64 (d, c);
      mpz_add (a, b, d);
      mpz_clear (d);
    }
}

void
mpz_add_int64 (mpz_ptr a, mpz_srcptr b, int64_t c)
{
  if (c >= 0)
    mpz_add_uint64 (a, b, (uint64_t) c);
  else
    mpz_sub_uint64 (a, b, -(uint64_t) c);
}

void
mpz_sub_uint64 (mpz_ptr a, mpz_srcptr b, uint64_t c)
{
  if (c <= ULONG_MAX)
    mpz_sub_ui (a, b, (unsigned long) c);
  else
    {
      mpz_t d;
      mpz_init_set_uint64 (d, c);
      mpz_sub (a, b, d);
      mpz_clear (d);
    }
}

void
mpz_sub_int64 (mpz_ptr a, mpz_srcptr b, int64_t c)
{
  if (c >= 0)
    mpz_sub_uint64 (a, b, (uint64_t) c);
  else
    mpz_add_uint64 (a, b, -(uint64_t) c);
}

void
mpz_addmul_uint64 (mpz_ptr a, mpz_srcptr b, uint64_t c)
{
  if (c <= ULONG_MAX)
    mpz_addmul_ui (a, b, (unsigned long) c);
  else
    {
      mpz_t d;
      mpz_init_set_uint64 (d, c);
      mpz_addmul (a, b, d);
      mpz_clear (d);
    }
}

void
mpz_submul_uint64 (mpz_ptr a, mpz_srcptr b, uint64_t c)
{
  if (c <= ULONG_MAX)
    mpz_submul_ui (a, b, (unsigned long) c);
  else
    {
      mpz_t d;
      mpz_init_set_uint64 (d, c);
      mpz_submul (a, b, d);
      mpz_clear (d);
    }
}

void
mpz_divexact_uint64 (mpz_ptr a, mpz_srcptr b, uint64_t c)
{
  if (c <= ULONG_MAX)
    mpz_divexact_ui (a, b, (unsigned long) c);
  else
    {
      mpz_t d;
      mpz_init_set_uint64 (d, c);
      mpz_divexact (a, b, d);
      mpz_clear (d);
    }
}

int
mpz_divisible_uint64_p (mpz_ptr a, uint64_t c)
{
  if (c <= ULONG_MAX)
    return mpz_divisible_ui_p (a, (unsigned long) c);
  else
    {
      mpz_t d;
      int ret;

      mpz_init_set_uint64 (d, c);
      ret = mpz_divisible_p (a, d);
      mpz_clear (d);
      return ret;
    }
}

void
mpz_mul_int64 (mpz_ptr a, mpz_srcptr b, int64_t c)
{
  if (LONG_MIN <= c && c <= LONG_MAX)
    mpz_mul_si (a, b, (long) c);
  else
    {
      mpz_t d;
      mpz_init_set_int64 (d, c);
      mpz_mul (a, b, d);
      mpz_clear (d);
    }
}

uint64_t
uint64_nextprime (uint64_t q)
{
  mpz_t z;

  mpz_init_set_uint64 (z, q);
  mpz_nextprime (z, z);
  q = mpz_get_uint64 (z);
  mpz_clear (z);
  return q;
}

uint64_t
mpz_tdiv_qr_uint64 (mpz_ptr Q, mpz_ptr R, mpz_srcptr N, uint64_t d)
{
    mpz_t D;
    mpz_init_set_uint64 (D, d);
    mpz_tdiv_qr (Q, R, N, D);
    uint64_t r = mpz_get_uint64 (R);
    mpz_clear (D);
    return r;
}

uint64_t
mpz_tdiv_q_uint64 (mpz_ptr Q, mpz_srcptr N, uint64_t d)
{
    mpz_t R;
    mpz_init2 (R, 64);
    uint64_t r = mpz_tdiv_qr_uint64 (Q, R, N, d);
    mpz_clear (R);
    return r;
}

uint64_t
mpz_tdiv_r_uint64 (mpz_ptr R, mpz_srcptr N, uint64_t d)
{
    mpz_t D;
    mpz_init_set_uint64 (D, d);
    mpz_tdiv_r (R, N, D);
    uint64_t r = mpz_get_uint64 (R);
    mpz_clear (D);
    return r;
}

uint64_t
mpz_tdiv_uint64 (mpz_srcptr N, uint64_t d)
{
    mpz_t R;
    mpz_init2 (R, 64);
    uint64_t r = mpz_tdiv_r_uint64 (R, N, d);
    mpz_clear (R);
    return r;
}

#endif


void
mpz_submul_int64 (mpz_ptr a, mpz_srcptr b, int64_t c)
{
  if (c >= 0)
    mpz_submul_uint64 (a, b, (uint64_t) c);
  else
    mpz_addmul_uint64 (a, b, -(uint64_t) c);
}

/* a <- a + b * c */
void
mpz_addmul_int64 (mpz_ptr a, mpz_srcptr b, int64_t c)
{
  if (c >= 0)
    mpz_addmul_uint64 (a, b, (uint64_t) c);
  else
    mpz_submul_uint64 (a, b, -(uint64_t) c);
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

#define REPS 3 /* number of Miller-Rabin tests in isprime */

/* For GMP 6.1.2:
   with REPS=1, the smallest composite reported prime is 1537381
                (then 1803601, 1943521, 2237017, 3604201, 5095177, ...);
   with REPS=2, it is 1943521 (then 16661633, 18790021, 54470491, ...);
   with REPS=3, it is 465658903 (then 2242724851, 5969607379, 6635692801, ...);
   with REPS=4, it is 239626837621 (then 277376554153, 537108528133,
                      898547205403, and no other <= 10^12).
   See also https://en.wikipedia.org/wiki/Miller–Rabin_primality_test
   and http://www.trnicely.net/misc/mpzspsp.html.
   GMP 6.2.1:
   This implementation of the Baillie-PSW test was checked up to 31*2^46
   (https://gmplib.org/repo/gmp-6.2/file/tip/mpz/millerrabin.c#l11)
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
      for (p = 2; p < pmin && q % p; p += 1UL + (p > 2));
      if (p >= pmin)
        break;
    }
  return q;
}

/* Put in r the smallest legitimate value that it at least s + diff (note
   that if s+diff is already legitimate, then r = s+diff will result.

   Here, legitimate means prime or squarefree composite, with the constraint
   that all the prime factors must be in [pmin, pmax[ .

   The prime factors of r are put in factor_r, and the number of them is
   returned. The caller must have allocated factor_r with enough space.
   */
int 
next_mpz_with_factor_constraints(mpz_ptr r,
        unsigned long factor_r[],
        mpz_srcptr s,
        unsigned long diff,
        unsigned long pmin,
        unsigned long pmax)
{
    ASSERT_ALWAYS(pmin > 2);
    mpz_add_ui(r, s, diff);
    if (mpz_cmp_ui(r, pmin) < 0)
        mpz_set_ui(r, pmin);
    if (mpz_even_p(r))
        mpz_add_ui(r, r, 1);
    for( ; ; mpz_add_ui(r, r, 2)) { // only odd values are considered.
        /* This function is used for the composite special-q code, and it
         * is really important that no composite special-q is considered
         * prime.*/
        if (mpz_probab_prime_p(r, 10)) {
            if (mpz_cmp_ui(r, pmax) < 0) {
                factor_r[0] = mpz_get_ui(r);
                return 1;
            }
        } else { // r is composite. Check its factorization.
            // Work with unsigned long
            unsigned long rr = mpz_get_ui(r);
            int nf = 0;
            prime_info pi;
            prime_info_init(pi);
            unsigned long p;
            bool ok = true;
            for( ; (p = getprime_mt(pi)) < pmin ; ) {
                if ((rr % p) == 0) {
                    ok = false;
                    break;
                }
            }
            if (!ok) {
                // r is divisible by a prime < pmin.
                // try next candidate.
                prime_info_clear(pi);
                continue;
            }
            // Compute the factorization
            for( ; rr > 1 && p < pmax ; p = getprime_mt(pi)) {
                if (rr % p)
                    continue;

                int c = 0;
                for( ; rr % p == 0 ; rr /= p, c++) ;
                // p divides exactly c times r. Note that 2 or
                // more is too many.
                if (c > 1)
                    break;
                factor_r[nf++] = p;
                // check primality of cofactor.
                bool cofac_prime;
                {
                    mpz_t xxx;
                    mpz_init_set_ui(xxx, rr);
                    cofac_prime = mpz_probab_prime_p(xxx, 10);
                    mpz_clear(xxx);
                }
                
                if (!cofac_prime)
                    continue;

                if (rr >= pmax)
                    break;

                // We have a winner.
                factor_r[nf++] = rr;
                prime_info_clear(pi);
                return nf;
            }
            prime_info_clear(pi);
            // We have still a composite in rr, it must have a prime > pmax.
            // This is a fail. Continue with next candidate.
        }
    }
    ASSERT_ALWAYS(0); // Should never get there.
}

/* return the number of bits of p, counting from the least significant end.
 * 0 has 0 bit
 * 1 has 1 bit
 * 2 and 3 have 2 bits,
 * 4 to 7 have 3 bits, etc
 */
int nbits (uintmax_t p)
{
  int k;

  for (k = 0; p != 0; p >>= 1U, k ++);
  return k;
}

/* q <- n/d rounded to nearest, assuming d <> 0
   r <- n - q*d
   Output: -d/2 <= r < d/2
*/
void
mpz_ndiv_qr (mpz_ptr q, mpz_ptr r, mpz_srcptr n, mpz_srcptr d)
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
mpz_ndiv_qr_ui (mpz_ptr q, mpz_ptr r, mpz_srcptr n, unsigned long int d)
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
mpz_ndiv_q (mpz_ptr q, mpz_srcptr n, mpz_srcptr d)
{
  mpz_t r;

  mpz_init (r);
  mpz_ndiv_qr (q, r, n, d);
  mpz_clear (r);
}

void
mpz_ndiv_q_ui (mpz_ptr q, mpz_srcptr n, unsigned long int d)
{
  mpz_t r;

  mpz_init (r);
  mpz_ndiv_qr_ui (q, r, n, d);
  mpz_clear (r);
}

/* a <- b cmod c with -c/2 <= a < c/2 */
void
mpz_ndiv_r (mpz_ptr a, mpz_srcptr b, mpz_srcptr c)
{
  mpz_mod (a, b, c); /* now 0 <= a < c */

  mp_size_t n = (mp_size_t) mpz_size (c);
  mp_limb_t aj, cj;
  int sub = 0;
  unsigned int sh = GMP_NUMB_BITS - 1;

  if (n == 0) {
      mpz_set_ui (a, 0);
      return;
  }

  if (mpz_getlimbn (a, n-1) >= (mp_limb_t) 1U << sh)
    sub = 1;
  else
    {
      while (n-- > 0)
	{
	  cj = mpz_getlimbn (c, n);
	  aj = mpz_getlimbn (a, n) << 1U;
	  if (n > 0)
	    aj |= mpz_getlimbn (a, n-1) >> sh;
	  if (aj > cj)
	    {
	      sub = 1;
	      break;
	    }
	  else if (aj < cj)
	    break;
	}
    }
  if (sub)
    mpz_sub (a, a, c);
}

int mpz_rdiv_q(mpz_ptr q, mpz_srcptr a, mpz_srcptr b)
{
    /* Return the relative integer q which is closest to a/b. We
     * guarantee -1/2<=a/b-q<1/2.
     * It's a pity that this is not supported directly by gmp... */
    mpz_t r;
    mpz_init(r);

    mpz_fdiv_qr(q,r,a,b);
    /* b>0: We want -b/2 <= a-bq < b/2 */
    /* b<0: We want  b/2 < a-bq <= -b/2 */
    mpz_mul_2exp(r,r,1);
    if (mpz_cmp(r,b) * mpz_sgn(b) >= 0) {
        mpz_add_ui(q, q, 1);
    }
    mpz_clear(r);

    return 1;
}

/* Return non-zero if a and b are coprime, else return 0 if gcd(a, b) != 1 */
int
mpz_coprime_p (mpz_srcptr a, mpz_srcptr b)
{
  mpz_t g;
  mpz_init(g);
  mpz_gcd (g, a, b);
  int ret = mpz_cmp_ui (g, 1);
  mpz_clear(g);
  return (ret == 0) ? 1 : 0;
}

/* returns the p-valuation of a, where p is expected to be prime. INT_MAX
 * is returned if a==0 */
int mpz_p_valuation(mpz_srcptr a, mpz_srcptr p)
{
    mpz_t c;
    int v = 0;
    if (mpz_size(a) == 0) return INT_MAX;
    mpz_init_set(c, a);
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
    int v = 0;
    if (mpz_size(a) == 0) return INT_MAX;
    mpz_init_set(c, a);
    for( ; mpz_divisible_ui_p(c, p) ; v++)
        mpz_fdiv_q_ui(c, c, p);
    mpz_clear(c);
    return v;
}

void memfill_random(void *p, size_t s, gmp_randstate_t rstate)
{
    size_t w = sizeof(mp_limb_t);
    mp_limb_t * pw = (mp_limb_t *)p;
    for( ; s >= w ; pw++, s -= w) *pw = gmp_urandomb_ui(rstate, GMP_LIMB_BITS);
    char * pc = (char*) pw;
    for( ; s ; pc++, s--) *pc = gmp_urandomb_ui(rstate, CHAR_BIT);
}

/* Adapted from mpz_ln2 in utils/usp.c */
double mpz_log(mpz_srcptr a)
{
    size_t l = mpz_sizeinbase (a, 2);
    if (l <= 1024) { /* a fits in a double */
        return log(fabs(mpz_get_d(a)));
    } else {
        mpz_t b;
        mpz_init (b);
        mpz_tdiv_q_2exp (b, a, l - 900U);
        double r = mpz_get_d(b);
        r = ((double) l - 900.0)*log(2.0) + log(fabs(r));
        mpz_clear (b);
        return r;
    }
}
