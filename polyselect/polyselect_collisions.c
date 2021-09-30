#include "cado.h"		// IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>		// malloc ...
#include <stdint.h>		// uint64_t
#include <limits.h>		/* for CHAR_BIT */
#include <math.h>               // log
#include <gmp.h>
#include "auxiliary.h"
#include "cado_poly.h"
#include "getprime.h"		// getprime
#include "gmp_aux.h"		// mpz_set_uint64
#include "macros.h"		// ASSERT
#include "memusage.h"
#include "misc.h"
#include "modredc_ul.h"
#include "mpz_poly.h"
#include "polyselect_arith.h"
#include "polyselect_shash.h"
#include "polyselect_hash.h"
#include "polyselect_locals.h"
#include "polyselect_norms.h"
#include "polyselect_main_queue.h"
#include "portability.h"
#include "roots_mod.h"
#include "size_optimization.h"
#include "polyselect_collisions.h"
#include "timing.h"		// for seconds
#include "usp.h"		// usp_root_data
#include "verbose.h"		// verbose_output_print

static void
check_divexact_ui(mpz_ptr r, mpz_srcptr d, const char *d_name MAYBE_UNUSED,
		  const unsigned long q, const char *q_name MAYBE_UNUSED)
{
#ifdef DEBUG_POLYSELECT
  if (mpz_divisible_ui_p(d, q) == 0)
    {
      gmp_fprintf(stderr, "Error: %s=%Zd not divisible by %s=%lu\n",
		  d_name, d, q_name, q);
      exit(1);
    }
#endif
  mpz_divexact_ui(r, d, q);
}

static void
check_divexact(mpz_ptr r, mpz_srcptr d, const char *d_name MAYBE_UNUSED,
	       const mpz_srcptr q, const char *q_name MAYBE_UNUSED)
{
#ifdef DEBUG_POLYSELECT
  if (mpz_divisible_p(d, q) == 0)
    {
      gmp_fprintf(stderr, "Error: %s=%Zd not divisible by %s=%Zd\n",
		  d_name, d, q_name, q);
      exit(1);
    }
#endif
  mpz_divexact(r, d, q);
}

/* rq is a root of N = (m0 + rq)^d mod (q^2) */
/* This is called by polyselect_hash_add
 *
 * If rq is NULL, then q is 1, and we simply have a (p1,p2) match
 */
void
match(unsigned long p1, unsigned long p2, const int64_t i,
      uint64_t q, 
      mpz_srcptr rq,
      polyselect_thread_locals_ptr loc)
{
  polyselect_poly_header_srcptr header = loc->header;
  mpz_t l, mtilde, m, adm1, t, k;
  mpz_poly f, g, f_raw, g_raw;
  int cmp, did_optimize;
  double skew, logmu;

  /* the expected rotation space is S^5 for degree 6 */
#ifdef DEBUG_POLYSELECT
  gmp_printf("Found match: (%lu,%lld) (%lu,%lld) for "
	     "ad=%Zd, q=%llu, rq=%Zd\n",
	     p1, (long long) i, p2, (long long) i, header->ad,
	     (unsigned long long) q, rq);
  gmp_printf("m0=%Zd\n", header->m0);
#endif

  mpz_init(l);
  mpz_init(m);
  mpz_init(t);
  mpz_init(k);
  mpz_init(adm1);
  mpz_init(mtilde);

  mpz_poly_init(f, header->d);
  mpz_poly_init(g, 1);
  mpz_poly_init(f_raw, header->d);
  mpz_poly_init(g_raw, 1);
  /* we have l = p1*p2*q */
  mpz_set_ui(l, p1);
  mpz_mul_ui(l, l, p2);
  mpz_mul_ui(l, l, q);
  /* mtilde = header->m0 + rq + i*q^2 */
  mpz_set_si(mtilde, i);
  if (rq) {
      mpz_mul_ui(mtilde, mtilde, q);
      mpz_mul_ui(mtilde, mtilde, q);
      mpz_add(mtilde, mtilde, rq);
  }
  mpz_add(mtilde, mtilde, header->m0);
  /* we should have Ntilde - mtilde^d = 0 mod {p1^2,p2^2,q^2} */

  /* Small improvement: we have Ntilde = mtilde^d + l^2*R with R small.
     If p^2 divides R, with p prime to d*ad, then we can accumulate p into l,
     which will give an even smaller R' = R/p^2.
     Note: this might produce duplicate polynomials, since a given p*l
     might be found in different ways. For example with revision b5a1635 and
     polyselect -P 60000 -N 12939597433839929710052817774007139127064894178566832462175875720079522272519444917218095639720802504629187785806903263303 -degree 5 -t 1 -admin 780 -admax 840 -incr 60 -nq 2317
     the polynomial with Y1 = 35641965604484971 is found four times:
     * once with q = 92537 = 37 * 41 * 61
     * then with q = 182573 = 41 * 61 * 73
     * then with q = 110741 = 37 * 41 * 73
     * and finally with q = 164761 = 37 * 61 * 73
     As a workaround, we only allow p > qmax, the largest prime factor of q.
   */

  /* compute the largest prime factor of q */
  unsigned long qmax = 1;
  for (unsigned long j = 0; j < LEN_SPECIAL_Q - 1; j++)
    if ((q % SPECIAL_Q[j]) == 0)
      qmax = SPECIAL_Q[j];

  mpz_mul_ui(m, header->ad, header->d);
  mpz_pow_ui(m, m, header->d);
  mpz_divexact(m, m, header->ad);
  mpz_mul(m, m, header->N);		/* m := Ntilde = d^d*ad^(d-1)*N */
  mpz_pow_ui(t, mtilde, header->d);
  mpz_sub(t, m, t);
  mpz_divexact(t, t, l);
  mpz_divexact(t, t, l);
  unsigned long p;

  prime_info pi;
  prime_info_init(pi);
  /* Note: we could find p^2 dividing t in a much more efficient way, for
     example by precomputing the product of all primes < 2*P, then doing
     a gcd with t, which gives say g, then computing gcd(t, t/g).
     But if P is small, it would gain little with respect to the naive loop
     below, and if P is large, we have only a few hits, thus the global
     overhead will be small too. */
  for (p = 2; p <= loc->main->Primes[loc->main->lenPrimes - 1]; p = getprime_mt(pi))
    {
      if (p <= qmax || polyselect_poly_header_skip(header, p))
	continue;
      while (mpz_divisible_ui_p(t, p * p))
	{
	  mpz_mul_ui(l, l, p);
	  mpz_divexact_ui(t, t, p * p);
	}
    }
  prime_info_clear(pi);
  /* end of small improvement */

  /* we want mtilde = d*ad*m + a_{d-1}*l with -d*ad/2 <= a_{d-1} < d*ad/2.
     We have a_{d-1} = mtilde/l mod (d*ad). */
  mpz_mul_ui(m, header->ad, header->d);
  if (mpz_invert(adm1, l, m) == 0)
    {
      fprintf(stderr, "Error in 1/l mod (d*ad)\n");
      exit(1);
    }
  mpz_mul(adm1, adm1, mtilde);
  mpz_mod(adm1, adm1, m);	/* m is d*ad here */

  /* we make -d*ad/2 <= adm1 < d*ad/2 */
  mpz_mul_2exp(t, adm1, 1);
  if (mpz_cmp(t, m) >= 0)
    mpz_sub(adm1, adm1, m);

  mpz_mul(m, adm1, l);
  mpz_sub(m, mtilde, m);
  check_divexact_ui(m, m, "m-a_{d-1}*l", header->d, "d");
  check_divexact(m, m, "(m-a_{d-1}*l)/d", header->ad, "ad");
  mpz_set(g->coeff[1], l);
  mpz_neg(g->coeff[0], m);
  mpz_set(f->coeff[header->d], header->ad);
  mpz_pow_ui(t, m, header->d);
  mpz_mul(t, t, header->ad);
  mpz_sub(t, header->N, t);
  mpz_set(f->coeff[header->d - 1], adm1);
  check_divexact(t, t, "t", l, "l");
  mpz_pow_ui(mtilde, m, header->d - 1);
  mpz_mul(mtilde, mtilde, adm1);
  mpz_sub(t, t, mtilde);
  for (unsigned long j = header->d - 2; j > 0; j--)
    {
      check_divexact(t, t, "t", l, "l");
      /* t = a_j*m^j + l*R thus a_j = t/m^j mod l */
      mpz_pow_ui(mtilde, m, j);
      /* fdiv rounds toward -infinity: adm1 = floor(t/mtilde) */
      mpz_fdiv_q(adm1, t, mtilde);	/* t -> adm1 * mtilde + t */
      mpz_invert(k, mtilde, l);	/* search adm1 + k such that
				   t = (adm1 + k) * m^j mod l */
      mpz_mul(k, k, t);
      mpz_sub(k, k, adm1);
      mpz_mod(k, k, l);

      mpz_mul_2exp(k, k, 1);
      cmp = mpz_cmp(k, l);
      mpz_div_2exp(k, k, 1);
      if (cmp >= 0)
	mpz_sub(k, k, l);
      mpz_add(adm1, adm1, k);
      mpz_set(f->coeff[j], adm1);
      /* subtract adm1*m^j */
      mpz_submul(t, mtilde, adm1);
    }
  check_divexact(t, t, "t", l, "l");
  mpz_set(f->coeff[0], t);

  /* As noticed by Min Yang, Qingshu Meng, Zhangyi Wang, Lina Wang and
     Huanguo Zhang in "Polynomial Selection for the Number Field Sieve in an
     Elementary Geometric View" (https://eprint.iacr.org/2013/583),
     if the coefficient of degree d-2 is of the same sign as the leading
     coefficient, the size optimization will not work well, thus we simply
     discard those polynomials. */
  if (mpz_sgn(f->coeff[header->d]) * mpz_sgn(f->coeff[header->d - 2]) > 0)
    {
      loc->stats->discarded1++;
      goto end;
    }


  mpz_poly_cleandeg(f, header->d);
  ASSERT_ALWAYS(mpz_poly_degree(f) == (int) header->d);
  mpz_poly_cleandeg(g, 1);
  ASSERT_ALWAYS(mpz_poly_degree(g) == (int) 1);

  mpz_poly_set(g_raw, g);
  mpz_poly_set(f_raw, f);

  /* _raw lognorm */
  skew = L2_skewness(f, SKEWNESS_DEFAULT_PREC);
  logmu = L2_lognorm(f, skew);

  {
    /* information on all polynomials */
    loc->stats->collisions++;
    loc->stats->tot_found++;
    polyselect_data_series_add(loc->stats->raw_lognorm, logmu);
    polyselect_data_series_add(loc->stats->raw_proj_alpha,
				 get_alpha_projective(f, get_alpha_bound()));
  }

  /* if the polynomial has small norm, we optimize it */
  did_optimize = optimize_raw_poly(f, g, loc->main);

  /* print optimized (maybe size- or size-root- optimized) polynomial */
  if (did_optimize && loc->main->verbose >= 0)
      output_polynomials(f_raw, g_raw, header->N, f, g, loc);

end:
  mpz_clear(l);
  mpz_clear(m);
  mpz_clear(t);
  mpz_clear(k);
  mpz_clear(adm1);
  mpz_clear(mtilde);
  mpz_poly_clear(f);
  mpz_poly_clear(g);
  mpz_poly_clear(f_raw);
  mpz_poly_clear(g_raw);
}


/*{{{ polyselect_proots_dispatch_to_shash_flat_ugly */

/* This code dispatches congruence classes of the roots_per_prime[] table,
 * modulo the primes that are listed in Primes, within the range
 * [-umax...umax], into the quick hash table H.
 *
 * This overly complicated code provides some tiny performance gain, but
 * the simpler version below is fine too.
 *
 * the value number_of_roots_per_prime[(number of primes)] must be 0xff
 *
 * This means that the following must typically appear before the call to
 * this function.
 *
 * loc->R->nr[loc->R->size] = 0xff;
 *
 */
void polyselect_proots_dispatch_to_shash_flat_ugly(
        polyselect_shash_ptr H,
        const uint32_t * Primes,
        const unsigned long * roots_per_prime,
        const uint8_t * number_of_roots_per_prime,
        int64_t umax)
{
  uint64_t **cur1, **cur2, *ccur1, *ccur2;
  long *pc, *epc;
  int64_t ppl, neg_umax, v1, v2, nv;
  uint8_t vpnr;

#if polyselect_SHASH_NBUCKETS == 256
#define CURRENT(V) (H->current + (uint8_t) (V))
#else
#define CURRENT(V) (H->current + ((V) & (polyselect_SHASH_NBUCKETS - 1)))
#endif

  polyselect_shash_reset(H);

  pc = (long *) roots_per_prime;
  nv = *pc;
  const uint32_t * pprimes = Primes - 1;
  const uint8_t * pnr = number_of_roots_per_prime;
  neg_umax = -umax;

  /* This define inserts 2 values v1 and v2 with a interlace.
     The goal is to have a little time to read ccurX from L0
     cache before to use it. The best seems a
     three read interlacing in fact, two seems too short. */
#define INSERT_2I(I1,I2)                                                \
  do {                                                                  \
    cur1 = CURRENT(I1); ccur1 = *cur1;					\
    cur2 = CURRENT(I2); ccur2 = *cur2;					\
    *ccur1++ = I1; __builtin_prefetch(ccur1, 1, 3); *cur1 = ccur1;	\
    *ccur2++ = I2; __builtin_prefetch(ccur2, 1, 3); *cur2 = ccur2;	\
  } while (0)
  /* This version is slow because ccur1 is used immediatly after
     it has been read from L0 cache -> 3 ticks of latency on P4 Nehalem. */
#define INSERT_I(I)						\
  do {								\
    cur1 = CURRENT(I); ccur1 = *cur1; *ccur1++ = I;		\
    __builtin_prefetch(ccur1, 1, 3); *cur1 = ccur1;		\
  } while (0)

  int64_t b;
  b = (int64_t) ((double) umax * 0.3333333333333333);
  do
    {
      do
	{
	  vpnr = *pnr++;
	  pprimes++;
	}
      while (!vpnr);
      if (UNLIKELY(vpnr == 0xff))
	goto bend;
      ppl = *pprimes;
      __builtin_prefetch(((void *) pnr) + 0x040, 0, 3);
      __builtin_prefetch(((void *) pprimes) + 0x80, 0, 3);
      __builtin_prefetch(((void *) pc) + 0x100, 0, 3);
      ppl *= ppl;
      epc = pc + vpnr;
      if (UNLIKELY(ppl > b))
	{
	  b = umax >> 1;
	  goto iter2;
	}
      do
	{
	  v1 = nv;
	  cur1 = CURRENT(v1);
	  ccur1 = *cur1;
	  v2 = v1 - ppl;
	  cur2 = CURRENT(v2);
	  ccur2 = *cur2;
	  nv = *++pc;
	  *ccur1++ = v1;
	  __builtin_prefetch(ccur1, 1, 3);
	  *cur1 = ccur1;
	  *ccur2++ = v2;
	  __builtin_prefetch(ccur2, 1, 3);
	  *cur2 = ccur2;
	  v1 += ppl;
	  cur1 = CURRENT(v1);
	  ccur1 = *cur1;
	  v2 -= ppl;
	  cur2 = CURRENT(v2);
	  ccur2 = *cur2;
	  *ccur1++ = v1;
	  __builtin_prefetch(ccur1, 1, 3);
	  *cur1 = ccur1;
	  *ccur2++ = v2;
	  __builtin_prefetch(ccur2, 1, 3);
	  *cur2 = ccur2;
	  v1 += ppl;
	  cur1 = CURRENT(v1);
	  ccur1 = *cur1;
	  v2 -= ppl;
	  cur2 = CURRENT(v2);
	  ccur2 = *cur2;
	  *ccur1++ = v1;
	  __builtin_prefetch(ccur1, 1, 3);
	  *cur1 = ccur1;
	  *ccur2++ = v2;
	  __builtin_prefetch(ccur2, 1, 3);
	  *cur2 = ccur2;
	  v1 += ppl;
	  v2 -= ppl;
	  if (LIKELY(v1 > umax))
	    {
	      if (UNLIKELY(v2 >= neg_umax))
		INSERT_I(v2);
	  } else if (UNLIKELY(v2 >= neg_umax))
	    INSERT_2I(v1, v2);
	  else
	    INSERT_I(v1);
	}
      while (pc != epc);
    }
  while (1);

  do
    {
      do
	{
	  vpnr = *pnr++;
	  pprimes++;
	}
      while (!vpnr);
      if (UNLIKELY(vpnr == 0xff))
	goto bend;
      ppl = *pprimes;
      __builtin_prefetch(((void *) pnr) + 0x040, 0, 3);
      __builtin_prefetch(((void *) pprimes) + 0x100, 0, 3);
      __builtin_prefetch(((void *) pc) + 0x280, 0, 3);
      ppl *= ppl;
      epc = pc + vpnr;
    iter2:
      if (UNLIKELY(ppl > b))
	goto iter1;
      do
	{
	  v1 = nv;
	  cur1 = CURRENT(v1);
	  ccur1 = *cur1;
	  v2 = v1 - ppl;
	  cur2 = CURRENT(v2);
	  ccur2 = *cur2;
	  nv = *++pc;
	  *ccur1++ = v1;
	  __builtin_prefetch(ccur1, 1, 3);
	  *cur1 = ccur1;
	  *ccur2++ = v2;
	  __builtin_prefetch(ccur2, 1, 3);
	  *cur2 = ccur2;
	  v1 += ppl;
	  cur1 = CURRENT(v1);
	  ccur1 = *cur1;
	  v2 -= ppl;
	  cur2 = CURRENT(v2);
	  ccur2 = *cur2;
	  *ccur1++ = v1;
	  __builtin_prefetch(ccur1, 1, 3);
	  *cur1 = ccur1;
	  *ccur2++ = v2;
	  __builtin_prefetch(ccur2, 1, 3);
	  *cur2 = ccur2;
	  v1 += ppl;
	  v2 -= ppl;
	  if (LIKELY(v1 > umax))
	    {
	      if (UNLIKELY(v2 >= neg_umax))
		INSERT_I(v2);
	  } else if (UNLIKELY(v2 >= neg_umax))
	    INSERT_2I(v1, v2);
	  else
	    INSERT_I(v1);
	}
      while (pc != epc);
    }
  while (1);

  do
    {
      do
	{
	  vpnr = *pnr++;
	  pprimes++;
	}
      while (!vpnr);
      if (UNLIKELY(vpnr == 0xff))
	goto bend;
      ppl = *pprimes;
      __builtin_prefetch(((void *) pnr) + 0x040, 0, 3);
      __builtin_prefetch(((void *) pprimes) + 0x100, 0, 3);
      __builtin_prefetch(((void *) pc) + 0x280, 0, 3);
      ppl *= ppl;
      epc = pc + vpnr;
    iter1:
      do
	{
	  v1 = nv;
	  cur1 = CURRENT(v1);
	  ccur1 = *cur1;
	  v2 = v1 - ppl;
	  cur2 = CURRENT(v2);
	  ccur2 = *cur2;
	  nv = *++pc;
	  *ccur1++ = v1;
	  __builtin_prefetch(ccur1, 1, 3);
	  *cur1 = ccur1;
	  *ccur2++ = v2;
	  __builtin_prefetch(ccur2, 1, 3);
	  *cur2 = ccur2;
	  v1 += ppl;
	  v2 -= ppl;
	  if (LIKELY(v1 > umax))
	    {
	      if (UNLIKELY(v2 >= neg_umax))
		INSERT_I(v2);
	  } else if (UNLIKELY(v2 >= neg_umax))
	    INSERT_2I(v1, v2);
	  else
	    INSERT_I(v1);
	}
      while (pc != epc);
    }
  while (1);

bend:
#undef INSERT_2I
#undef INSERT_I

  ;
}

/*}}}*/

/*{{{ polyselect_proots_dispatch_to_shash_flat */
void polyselect_proots_dispatch_to_shash_flat(
        polyselect_shash_ptr H,
        const uint32_t * Primes,
        size_t lenPrimes,
        const unsigned long * roots_per_prime,
        const uint8_t * number_of_roots_per_prime,
        int64_t umax)
{
    polyselect_shash_reset(H);
    unsigned long c = 0;
    for (unsigned long nprimes = 0; nprimes < lenPrimes; nprimes++)
    {
        unsigned long p = Primes[nprimes];
        int64_t ppl = (int64_t) p *(int64_t) p;
        unsigned long nr = number_of_roots_per_prime[nprimes];
        for (unsigned long j = 0; j < nr; j++, c++)
        {
            // int64_t u0 = (((int64_t) roots_per_prime[c] + umax) % ppl) - umax;
            int64_t u0 = roots_per_prime[c];
            for(int64_t u = u0 ; u < umax ; u += ppl)
                polyselect_shash_add(H, u);
            for(int64_t u = u0 - ppl ; u + umax >= 0; u -= ppl)
                polyselect_shash_add(H, u);
        }
    }

    for (int i = 0; i < polyselect_SHASH_NBUCKETS; i++)
        ASSERT(H->current[i] <= H->base[i + 1]);
}
/*}}}*/

/*{{{ polyselect_proots_dispatch_to_shash_notflat */
/* This is a slight variation around the previous implementation. Here,
 * we expect the roots to be organized into several different tables per
 * prime.
 */
void polyselect_proots_dispatch_to_shash_notflat(
        polyselect_shash_ptr H,
        const uint32_t * Primes,
        size_t lenPrimes,
        uint64_t * const * roots_per_prime,
        const uint8_t * number_of_roots_per_prime,
        int64_t umax)
{
      for (unsigned long nprimes = 0; nprimes < lenPrimes; nprimes++)
	{
	  unsigned long p = Primes[nprimes];
          int64_t ppl = (int64_t) p *(int64_t) p;
	  unsigned long nr = number_of_roots_per_prime[nprimes];
	  for (unsigned long j = 0; j < nr; j++)
	    {
                // int64_t u0 = (((int64_t) roots_per_prime[nprimes][j] + umax) % ppl) - umax;
                int64_t u0 = roots_per_prime[nprimes][j];
                for(int64_t u = u0 ; u < umax ; u += ppl)
                    polyselect_shash_add(H, u);
                for(int64_t u = u0 - ppl ; u + umax >= 0 ; u -= ppl)
                    polyselect_shash_add(H, u);
            }
        }

      for (int i = 0; i < polyselect_SHASH_NBUCKETS; i++)
          ASSERT(H->current[i] <= H->base[i + 1]);
}/*}}}*/

/*{{{ polyselect_proots_dispatch_to_hash_notflat */
/* same as above, but for a hash (not shash) table */
void polyselect_proots_dispatch_to_hash_notflat(
        polyselect_hash_ptr H,
        const uint32_t * Primes,
        size_t lenPrimes,
        uint64_t * const * roots_per_prime,
        const uint8_t * number_of_roots_per_prime,
        int64_t umax,
        unsigned long q,
        mpz_srcptr rq,
        polyselect_thread_locals_ptr loc
        )
{
      for (unsigned long nprimes = 0; nprimes < lenPrimes; nprimes++)
	{
	  unsigned long p = Primes[nprimes];
          int64_t ppl = (int64_t) p *(int64_t) p;
	  unsigned long nr = number_of_roots_per_prime[nprimes];
	  for (unsigned long j = 0; j < nr; j++)
	    {
                // int64_t u0 = (((int64_t) roots_per_prime[nprimes][j] + umax) % ppl) - umax;
                int64_t u0 = roots_per_prime[nprimes][j];
                for(int64_t u = u0 ; u < umax ; u += ppl)
                    polyselect_hash_add(H, p, u, q, rq, loc);
                for(int64_t u = u0 - ppl ; u + umax >= 0 ; u -= ppl)
                    polyselect_hash_add(H, p, u, q, rq, loc);
            }
        }
}/*}}}*/

/*{{{ polyselect_proots_dispatch_to_hash_flat */
void polyselect_proots_dispatch_to_hash_flat(
        polyselect_hash_ptr H,
        const uint32_t * Primes,
        size_t lenPrimes,
        const unsigned long * roots_per_prime,
        const uint8_t * number_of_roots_per_prime,
        int64_t umax,
        unsigned long q,
        mpz_srcptr rq,
        polyselect_thread_locals_ptr loc
        )
{
    unsigned long c = 0;
    for (unsigned long nprimes = 0; nprimes < lenPrimes; nprimes++)
    {
        unsigned long p = Primes[nprimes];
        int64_t ppl = (int64_t) p *(int64_t) p;
        unsigned long nr = number_of_roots_per_prime[nprimes];
        for (unsigned long j = 0; j < nr; j++, c++)
        {
            // int64_t u0 = (((int64_t) roots_per_prime[nprimes][j] + umax) % ppl) - umax;
            int64_t u0 = roots_per_prime[c];
            for(int64_t u = u0 ; u < umax ; u += ppl)
                polyselect_hash_add(H, p, u, q, rq, loc);
            for(int64_t u = u0 - ppl ; u + umax >= 0 ; u -= ppl)
                polyselect_hash_add(H, p, u, q, rq, loc);
        }
    }
}/*}}}*/

static unsigned long compute_and_lift_proots(polyselect_thread_locals_ptr loc)/*{{{*/
{
    polyselect_proots_ptr R = loc->R;
    polyselect_poly_header_srcptr header = loc->header;
    gmp_randstate_ptr rstate = loc->rstate;

    unsigned long tot_roots = 0;
    uint64_t *
        rp = (uint64_t *) malloc(header->d * sizeof(uint64_t));
    if (rp == NULL)
    {
        fprintf(stderr, "Error, cannot allocate memory in collision_on_p\n");
        exit(1);
    }

    for (unsigned long nprimes = 0; nprimes < loc->main->lenPrimes; nprimes++)
    {
        unsigned long p = loc->main->Primes[nprimes];

        /* add fake roots to keep indices */
        if (polyselect_poly_header_skip(header, p))
        {
            R->nr[nprimes] = 0;	// nr = 0.
            R->roots[nprimes] = NULL;
            continue;
        }

        unsigned long nrp = roots_mod_uint64(rp,
                mpz_fdiv_ui(header->Ntilde, p),
                header->d,
                p, rstate);
        tot_roots += nrp;
        nrp = roots_lift(rp, header->Ntilde, header->d, header->m0, p, nrp);
        polyselect_proots_add(R, nrp, rp, nprimes);
    }
    free(rp);
    return tot_roots;
}/*}}}*/

/* find collisions between "P" primes, return number of loops */
unsigned long
collision_on_p(polyselect_shash_ptr H, polyselect_thread_locals_ptr loc)
{
  int found = 0;
  int st = milliseconds();

  int64_t umax = polyselect_main_data_get_M(loc->main);

  /* first compute and lift all roots modulo the primes in
   * loc->main->Primes ; we store that in loc->R
   */
  unsigned long tot_roots = compute_and_lift_proots(loc);

  /* We first store only i's (and not p's), and look for collisions,
   * which occur very rarely.
   * if we find out that there is a collision on i, we run the search
   * again for the p's.
   */

  polyselect_shash_reset(H);
  polyselect_proots_dispatch_to_shash_notflat(H,
          loc->main->Primes,
          loc->main->lenPrimes,
          loc->R->roots,
          loc->R->nr,
          umax);

  st = milliseconds() - st;

  if (loc->main->verbose > 2)
    {
      fprintf(stderr, "# computing %lu p-roots took %dms\n", tot_roots, st);
    }

  st = milliseconds();
  found = polyselect_shash_find_collision(H);


  if (found)
    {				/* do the real work */
      polyselect_hash_t H;

      polyselect_hash_init(H, INIT_FACTOR * loc->main->lenPrimes, match);

      polyselect_proots_dispatch_to_hash_notflat(H,
              loc->main->Primes,
              loc->main->lenPrimes,
              loc->R->roots,
              loc->R->nr,
              umax,
              1, NULL, loc);

#ifdef DEBUG_POLYSELECT
      fprintf(stderr, "# collision_on_p took %lums\n", milliseconds() - st);
      gmp_fprintf(stderr, "# p polyselect_hash_size: %u for ad = %Zd\n",
		  H->size, header->ad);
#endif

#ifdef DEBUG_HASH_TABLE
      fprintf(stderr,
	      "# p polyselect_hash_size: %u, polyselect_hash_alloc: %u\n",
	      H->size, H->alloc);
      fprintf(stderr, "# hash table coll: %lu, all_coll: %lu\n", H->coll,
	      H->coll_all);
#endif
      polyselect_hash_clear(H);
    }

  loc->stats->potential_collisions++;

  return tot_roots;
}

static inline void
collision_on_each_sq(unsigned long q,
		     mpz_srcptr rqqz,
		     unsigned long *inv_qq,
                     polyselect_shash_ptr H,
                     polyselect_thread_locals_ptr loc)
{
  int found;

#ifdef DEBUG_POLYSELECT2
  int st = milliseconds();
#endif

  int64_t umax = polyselect_main_data_get_M(loc->main);

#if 0
  loc->R->nr[loc->R->size] = 0xff;	/* use guard to end */
  polyselect_proots_dispatch_to_shash_flat_ugly(H,
          loc->main->Primes,
          inv_qq,
          loc->R->nr,
          umax);
#else
  /* inv_qq is created by the caller as a flat list.  */
  polyselect_proots_dispatch_to_shash_flat(H,
          loc->main->Primes,
          loc->main->lenPrimes,
          inv_qq,
          loc->R->nr,
          umax);
#endif

  found = polyselect_shash_find_collision(H);

  if (found)
    {				/* do the real work */
      polyselect_hash_t H;

      polyselect_hash_init(H, INIT_FACTOR * loc->main->lenPrimes, match);

      polyselect_proots_dispatch_to_hash_flat(H,
              loc->main->Primes,
              loc->main->lenPrimes,
              inv_qq,
              loc->R->nr,
              umax,
              q, rqqz, loc);

      polyselect_hash_clear(H);
    }

  /* use DEBUG_POLYSELECT2 since this is too verbose */
#ifdef DEBUG_POLYSELECT2
  fprintf(stderr, "# inner collision_on_each_sq took %lums\n",
	  milliseconds() - st);
  fprintf(stderr, "# - q polyselect_hash_alloc (q=%lu): %u\n", q, H->alloc);
#endif

#ifdef DEBUG_HASH_TABLE
  fprintf(stderr,
	  "# p polyselect_hash_size: %u, polyselect_hash_alloc: %u\n",
	  H->size, H->alloc);
  fprintf(stderr, "# hash table coll: %lu, all_coll: %lu\n", H->coll,
	  H->coll_all);
#endif

  loc->stats->potential_collisions++;
}


/* Given p, rp, q, invqq[], for each rq of q, compute (rp - rq) / q^2 */
static inline void
collision_on_each_sq_r(unsigned long q,
		       const mpz_t * rqqz,
		       unsigned long *inv_qq,
		       unsigned long number_pr,
		       int count,
                       polyselect_shash_ptr H,
                       polyselect_thread_locals_ptr loc)
{
  if (count == 0)
    return;

  unsigned long c = 0;
  unsigned long **tinv_qq = malloc(count * sizeof(unsigned long *));

  if (!tinv_qq)
    {
      fprintf(stderr, "Error, cannot allocate memory in %s\n", __FUNCTION__);
      exit(1);
    }
  tinv_qq[0] = malloc((number_pr + 1) * count * sizeof(unsigned long));
  for (int k = 0; k < count; k++)
    {
      /* number_pr + 1 for guard for pre-load in collision_on_each_sq (nv) */
      if (k)
          tinv_qq[k] = tinv_qq[k-1] + number_pr + 1;
      tinv_qq[k][number_pr] = 0;
    }

  int st = milliseconds();

  /* for each rp, compute (rp-rq)*1/q^2 (mod p^2) */
  for (unsigned long nprimes = 0; nprimes < loc->main->lenPrimes; nprimes++)
    {
      if (!loc->R->nr[nprimes])
	continue;
      uint8_t nr = loc->R->nr[nprimes];
      unsigned long p = loc->main->Primes[nprimes];
      uint64_t pp = ((uint64_t) p) * ((uint64_t) p);

      modulusredcul_t modpp;
      residueredcul_t res_rqi, res_rp, res_tmp;
      modredcul_initmod_ul_raw(modpp, pp);
      modredcul_init(res_rqi, modpp);
      modredcul_init(res_rp, modpp);
      modredcul_init(res_tmp, modpp);

      for (int k = 0; k < count; k++)
	{
	  unsigned long rqi = mpz_fdiv_ui(rqqz[k], pp);
	  modredcul_intset_ul(res_rqi, rqi);
	  modredcul_intset_ul(res_tmp, inv_qq[nprimes]);
	  for (uint8_t i = 0; i < nr; i++, c++)
	    {
	      unsigned long rp = loc->R->roots[nprimes][i];
	      modredcul_intset_ul(res_rp, rp);
	      /* rp - rq */
	      modredcul_sub(res_rp, res_rp, res_rqi, modpp);
	      /* res_rp = (rp - rq) / q[i]^2 */
	      modredcul_mul(res_rp, res_rp, res_tmp, modpp);
	      tinv_qq[k][c] = modredcul_intget_ul(res_rp);
	    }
	  c -= nr;
	}
      c += nr;

      modredcul_clear(res_rp, modpp);
      modredcul_clear(res_rqi, modpp);
      modredcul_clear(res_tmp, modpp);
      modredcul_clearmod(modpp);
    }

  if (loc->main->verbose > 2)
    {
      fprintf(stderr,
	      "#  substage: batch %d many (rp-rq)*1/q^2 took %lums\n",
	      count, milliseconds() - st);
      st = milliseconds();
    }

  /* core function to find collisions */
  for (int k = 0; k < count; k++)
    collision_on_each_sq(q, rqqz[k], tinv_qq[k], H, loc);

  if (loc->main->verbose > 2)
    fprintf(stderr,
	    "#  substage: collision-detection %d many rq took %lums\n",
	    count, milliseconds() - st);

  free(tinv_qq[0]);
  free(tinv_qq);
}

/* Next combination */
static inline unsigned int
aux_nextcomb(unsigned int *ind, unsigned int len_q, unsigned int *len_nr)
{
  unsigned int i;

  /* bottom change first */
  for (i = len_q - 1;; i--)
    {
      if (ind[i] < (len_nr[i] - 1))
	{
	  ind[i]++;
	  return 1;
      } else
	{
	  if (i == 0)
	    break;
	  ind[i] = 0;
	}
    }
  return 0;
}


/* Compute crt-ed rq (qqz,rqqz) = (q_1 * ... * q_k,
                                   CRT([r_1, ..., r_k], [q_1, ..., q_k]))
   Well, almost. I think that it's all modulo the q_i^2 (squared!)
 */
static inline void
aux_return_rq(polyselect_qroots_srcptr SQ_R,
	      const unsigned long *idx_q,
	      const unsigned int *idx_nr,
	      unsigned long k, mpz_ptr qqz, mpz_ptr rqqz)
{
  unsigned long i, q[k], rq[k];

  /* q and roots */
  for (i = 0; i < k; i++)
    {
      q[i] = SQ_R->q[idx_q[i]];
      rq[i] = SQ_R->roots[idx_q[i]][idx_nr[i]];
    }

  /* crt roots */
  crt_sq(qqz, rqqz, q, rq, k);

  return;
}

/* Consider each rq which is the product of k pairs (q,r).
 * In this routine the q[i] are fixed, only the roots mod q[i] change.
 *
 * The name is misleading. There is nothing about reentrancy here, this
 * function is just a component of collision_on_sq
 */
static inline void
collision_on_batch_sq_r(polyselect_qroots_srcptr SQ_R,
			unsigned long q,
			const unsigned long *idx_q,
			unsigned long *inv_qq,
			unsigned long number_pr,
			unsigned long *curr_nq,
			unsigned long k,
                        polyselect_shash_ptr H,
                        polyselect_thread_locals_ptr loc)
{
  int count;
  unsigned int ind_qr[k];	/* indices of roots for each small q */
  unsigned int len_qnr[k];	/* for each small q, number of roots */
  unsigned long i;
  mpz_t qqz, rqqz[BATCH_SIZE];

  mpz_init(qqz);
  for (i = 0; i < BATCH_SIZE; i++)
    mpz_init(rqqz[i]);

#if 0
  fprintf(stderr, "q: %lu, ", q);
  for (i = 0; i < k; i++)
    fprintf(stderr, "%u ", SQ_R->q[idx_q[i]]);
  fprintf(stderr, ", ");
  for (i = 0; i < k; i++)
    fprintf(stderr, "%u ", SQ_R->nr[idx_q[i]]);
  fprintf(stderr, "\n");
#endif

  /* we proceed with BATCH_SIZE many rq for each time */
  for (i = 0; i < k; i++)
    {
      ind_qr[i] = 0;
      len_qnr[i] = SQ_R->nr[idx_q[i]];
    }
  for(int re = 1 ; re ; )
    {
      /* compute BATCH_SIZE such many rqqz[] */
      int num_rq = 0;
      for (count = 0; re && count < BATCH_SIZE; count++)
	{
          /* We do multiple inversions over and over again (in crt_sq),
           * even though the set of q's doesn't change... This has no
           * noticeable impact on the running time, though.
           */
	  aux_return_rq(SQ_R, idx_q, ind_qr, k, qqz, rqqz[count]);
	  re = aux_nextcomb(ind_qr, k, len_qnr);
	  (*curr_nq)++;
	  num_rq++;
	  if ((*curr_nq) >= loc->main->nq)
	    re = 0;
	}

      /* core function for a fixed qq and several rqqz[] */
      collision_on_each_sq_r(q, (const mpz_t *) rqqz, inv_qq,
			     number_pr, num_rq, H, loc);
    }

  mpz_clear(qqz);
  for (i = 0; i < BATCH_SIZE; i++)
    mpz_clear(rqqz[i]);
}

static inline void invert_q2_mod_all_p2(/*{{{*/
        unsigned long *invqq,
        unsigned long q,
        const uint32_t * Primes,
        unsigned long lenPrimes,
        polyselect_poly_header_srcptr header,
        const uint8_t * number_of_roots_per_prime)
{
  /* Step 1: inversion; compute 1/q^2 (mod p_i^2) to invqq[i] */
  for (unsigned long nprimes = 0; nprimes < lenPrimes; nprimes++)
    {
      unsigned long p = Primes[nprimes];
      if (polyselect_poly_header_skip(header, p))
	continue;
      unsigned int nr = number_of_roots_per_prime[nprimes];
      if (nr == 0)
	continue;
      uint64_t pp = ((uint64_t) p) * (uint64_t) p;

      modulusredcul_t modpp;
      residueredcul_t qq, tmp;
      modredcul_initmod_ul(modpp, pp);
      modredcul_init(qq, modpp);
      modredcul_init(tmp, modpp);

      /* q^2/B (mod pp). Warning: for large nq, we might have q > p^2, therefore
         we must first reduce q mod p^2 before calling modredcul_intset_ul. */
      modredcul_intset_ul(tmp, q % pp);
      modredcul_sqr(qq, tmp, modpp);
      /* B/q^2 (mod pp) */
      modredcul_intinv(tmp, qq, modpp);
      invqq[nprimes] = modredcul_intget_ul(tmp);

      modredcul_clear(tmp, modpp);
      modredcul_clear(qq, modpp);
      modredcul_clearmod(modpp);
    }

}/*}}}*/

/* collision on special-q, call collision_on_batch_sq */
void
collision_on_sq(unsigned long c, polyselect_shash_ptr H, polyselect_thread_locals_ptr loc)
{
  unsigned long *invqq = malloc(loc->main->lenPrimes * sizeof(unsigned long));
  if (!invqq)
    {
      fprintf(stderr, "Error, cannot allocate memory in %s\n", __FUNCTION__);
      exit(1);
    }

  /* init special-q roots */
  polyselect_qroots_t SQ_R;

  polyselect_qroots_init(SQ_R);
  comp_sq_roots(loc->header, SQ_R, loc->rstate);
  //polyselect_qroots_print (SQ_R);


  /* find a suitable lq */
  unsigned long k;
  unsigned long lq = find_suitable_lq(loc->header, SQ_R, &k, loc->main);

  unsigned long q, idx_q[lq], curr_nq = 0;

  first_comb(k, idx_q);

  for ( ; curr_nq < loc->main->nq ; )
    {
      q = return_q_norq(SQ_R, idx_q, k);

      /* collision batch */
      {
          int st = milliseconds();

          /* Step 1: inversion; compute 1/q^2 (mod p_i^2) to invqq[i] */
          invert_q2_mod_all_p2(invqq, q,
                  loc->main->Primes, loc->main->lenPrimes,
                  loc->header, loc->R->nr);

          if (loc->main->verbose > 2)
              fprintf(stderr,
                      "# stage (1/q^2 inversion) for %lu primes took %lums\n",
                      loc->main->lenPrimes, milliseconds() - st);
      }

      {
          /* Step 2: find collisions on q. */
          int st2 = milliseconds();

          collision_on_batch_sq_r(SQ_R, q, idx_q, invqq, c,
                  &curr_nq, k, H, loc);
          if (loc->main->verbose > 2)
              fprintf(stderr,
                      "#  stage (special-q) for %lu special-q's took %lums\n",
                      curr_nq, milliseconds() - st2);
      }

      unsigned long ret = next_comb(lq, k, idx_q);
      if (ret == k)		/* in case binomial(lq, k) < nq */
	break;
    }

  free(invqq);

  /* clean */
  polyselect_qroots_clear(SQ_R);
}

