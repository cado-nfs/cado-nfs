#include "cado.h"
#include <stdint.h>
#include <stdio.h>
#include <gmp.h>
#include "gmp_aux.h"
#include "mpz_poly.h"
#include "polyselect_alpha.h"
#include "polyselect_arith.h"
#include "polyselect_collisions_gmp.h"
#include "polyselect_hash.h"
#include "polyselect_locals.h"
#include "polyselect_main_queue.h"
#include "polyselect_norms.h"
#include "roots_mod.h"

/* XXX
 * This file will not stay.
 */

/* two functions that are duplicated in polyselect_collisions as well.
 */									\
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
/* this routine is called from polyselect_str.c */
void
gmp_match(unsigned long p1, unsigned long p2, int64_t i, mpz_srcptr m0,
	  mpz_srcptr ad, unsigned long d, mpz_srcptr N, uint64_t q,
	  mpz_srcptr rq,
          polyselect_thread_locals_ptr loc)
{
  mpz_t l, mtilde, m, adm1, t, k, qq, tmp;
  mpz_poly f, g, f_old, g_old;
  int cmp, did_optimize;
  double skew, logmu;

#ifdef DEBUG_POLYSELECT
  gmp_printf("Found match: (%" PRIu32 ",%lld) (%" PRIu32 ",%lld) for "
	     "ad=%Zd, q=%llu, rq=%Zd\n",
	     p1, (long long) i, p2, (long long) i, ad,
	     (unsigned long long) q, rq);
  gmp_printf("m0=%Zd\n", m0);
#endif
  mpz_init(tmp);
  mpz_init(l);
  mpz_init(m);
  mpz_init(t);
  mpz_init(k);
  mpz_init(qq);
  mpz_init(adm1);
  mpz_init(mtilde);

  mpz_poly_init(f, d);
  mpz_poly_init(g, 1);
  mpz_poly_init(f_old, d);
  mpz_poly_init(g_old, 1);

  /* we have l = p1*p2*q */
  mpz_set_ui(l, p1);
  mpz_mul_ui(l, l, p2);
  mpz_set_uint64(tmp, q);
  mpz_mul(l, l, tmp);
  /* mtilde = m0 + rq + i*q^2 */
  mpz_set(qq, tmp);		// qq = q
  mpz_mul(qq, qq, tmp);		// qq = q^2
  if (i >= 0)
    mpz_set_uint64(tmp, (uint64_t) i);
  else
    {
      mpz_set_uint64(tmp, (uint64_t) (-i));
      mpz_neg(tmp, tmp);
    }
  mpz_set(mtilde, tmp);
  mpz_mul(mtilde, mtilde, qq);
  mpz_add(mtilde, mtilde, rq);
  mpz_add(mtilde, mtilde, m0);
  /* we want mtilde = d*ad*m + a_{d-1}*l with 0 <= a_{d-1} < d*ad.
     We have a_{d-1} = mtilde/l mod (d*ad). */
  mpz_mul_ui(m, ad, d);
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
  check_divexact_ui(m, m, "m-a_{d-1}*l", d, "d");
  check_divexact(m, m, "(m-a_{d-1}*l)/d", ad, "ad");
  mpz_set(g->coeff[1], l);
  mpz_neg(g->coeff[0], m);
  mpz_set(f->coeff[d], ad);
  mpz_pow_ui(t, m, d);
  mpz_mul(t, t, ad);
  mpz_sub(t, N, t);
  mpz_set(f->coeff[d - 1], adm1);
  check_divexact(t, t, "t", l, "l");
  mpz_pow_ui(mtilde, m, d - 1);
  mpz_mul(mtilde, mtilde, adm1);
  mpz_sub(t, t, mtilde);
  for (unsigned long j = d - 2; j > 0; j--)
    {
      check_divexact(t, t, "t", l, "l");
      /* t = a_j*m^j + l*R thus a_j = t/m^j mod l */
      mpz_pow_ui(mtilde, m, j);
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

  mpz_poly_cleandeg(f, d);
  ASSERT_ALWAYS(mpz_poly_degree(f) == (int) d);
  mpz_poly_cleandeg(g, 1);
  ASSERT_ALWAYS(mpz_poly_degree(g) == (int) 1);

  mpz_poly_set(g_old, g);
  mpz_poly_set(f_old, f);

  /* _old lognorm */
  skew = L2_skewness(f, SKEWNESS_DEFAULT_PREC);
  logmu = L2_lognorm(f, skew);

#ifdef HAVE_OPENMP
#pragma	omp critical
#endif
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
    {
      output_polynomials(f_old, g_old, N, f, g, loc);
    }

  mpz_clear(tmp);
  mpz_clear(l);
  mpz_clear(m);
  mpz_clear(t);
  mpz_clear(k);
  mpz_clear(qq);
  mpz_clear(adm1);
  mpz_clear(mtilde);
  mpz_poly_clear(f);
  mpz_poly_clear(g);
  mpz_poly_clear(f_old);
  mpz_poly_clear(g_old);
}

/* compared to the integer code (not gmp), this code is practically
 * identical, except that it doesn't do the shash trick.
 *
 * Note the relevance of the _gmp layer is dubious at best, since root
 * finding is done with the uint64_t routines here...
 */

/* find collisions between "P" primes, return number of loops */
unsigned long
gmp_collision_on_p(polyselect_thread_locals_ptr loc)
{
  unsigned long j, p, nrp, c = 0;
  uint64_t *rp;
  int64_t ppl = 0, u, umax;
  mpz_t zero;

  /* init zero */
  mpz_init_set_ui(zero, 0);

  rp = (uint64_t *) malloc(loc->header->d * sizeof(uint64_t));
  if (rp == NULL)
    {
      fprintf(stderr, "Error, cannot allocate memory in collision_on_p\n");
      exit(1);
    }

  polyselect_hash_t H;

  polyselect_hash_init(H, INIT_FACTOR * loc->main->lenPrimes, gmp_match);

#ifdef DEBUG_POLYSELECT
  int st = milliseconds();
#endif

  umax = polyselect_main_data_get_M(loc->main);

  for (unsigned long nprimes = 0; nprimes < loc->main->lenPrimes; nprimes++)
    {
      p = loc->main->Primes[nprimes];
      ppl = (int64_t) p *(int64_t) p;

      /* XXX This is special to the _gmp code. Why do we need this?
       */
      /* add fake roots to keep indices */
      if (polyselect_poly_header_skip(loc->header, p))
	{
	  loc->R->nr[nprimes] = 0;	// nr = 0.
	  loc->R->roots[nprimes] = NULL;
	  continue;
	}

      /* we want p^2 | N - (m0 + i)^d, thus
         (m0 + i)^d = N (mod p^2) or m0 + i = N^(1/d) mod p^2 */
      nrp =
	  roots_mod_uint64(rp,
                  mpz_fdiv_ui(loc->header->Ntilde, p), loc->header->d, p,
                  loc->rstate);
      roots_lift(rp, loc->header->Ntilde, loc->header->d, loc->header->m0, p, nrp);
      polyselect_proots_add(loc->R, nrp, rp, nprimes);
      for (j = 0; j < nrp; j++, c++)
	{
	  for (u = (int64_t) rp[j]; u < umax; u += ppl)
	    polyselect_hash_add(H, p, u,
                                    loc->header->m0, loc->header->ad,
				    loc->header->d, loc->header->N, 1, zero,
                                    loc);
	  for (u = ppl - (int64_t) rp[j]; u < umax; u += ppl)
	    polyselect_hash_add(H, p, -u,
                                    loc->header->m0, loc->header->ad,
				    loc->header->d, loc->header->N, 1, zero,
                                    loc);
	}
    }

#ifdef DEBUG_POLYSELECT
  fprintf(stderr, "# collision_on_p took %lums\n", milliseconds() - st);
  gmp_fprintf(stderr, "# p polyselect_hash_size: %u for ad = %Zd\n",
	      H->size, header->ad);
#endif

  polyselect_hash_clear(H);

  free(rp);
  mpz_clear(zero);

  loc->stats->potential_collisions++;
  return c;
}


/* as before, this _gmp variant does not do the shash trick, which is
 * particularly ugly and disgusting in the integer version.
 * 
 * Note that inv_qq is an array of uint64_t's anyway. It's a set of
 * inverses mod p^2. If that is constrained to uint64_t, why not
 * constrain q to it as well??
 */
/* collision on each special-q, call collision_on_batch_p() */
static inline void
gmp_collision_on_each_sq(uint64_t q, mpz_srcptr rqqz,
                         uint64_t * inv_qq,
                         polyselect_thread_locals_ptr loc)
{
  unsigned int nr, j;
  unsigned long nprimes, p, c = 0;
  uint64_t pp;
  int64_t ppl, u, v, umax;

#ifdef DEBUG_POLYSELECT
  int st = milliseconds();
#endif

  polyselect_hash_t H;

  polyselect_hash_init(H, INIT_FACTOR * loc->main->lenPrimes, gmp_match);

  umax = polyselect_main_data_get_M(loc->main);

  for (nprimes = 0; nprimes < loc->main->lenPrimes; nprimes++)
    {
      p = loc->main->Primes[nprimes];
      if (polyselect_poly_header_skip(loc->header, p))
	continue;

      /* set p, p^2, ppl */
      pp = (uint64_t) p;
      pp *= (uint64_t) p;
      ppl = (int64_t) pp;
      nr = loc->R->nr[nprimes];

      for (j = 0; j < nr; j++, c++)
	{
	  u = (int64_t) inv_qq[c];

	  for (v = u; v < umax; v += ppl)
	    polyselect_hash_add(H, p, v, loc->header->m0, loc->header->ad,
				    loc->header->d, loc->header->N, q, rqqz,
                                    loc);
	  for (v = ppl - u; v < umax; v += ppl)
	    polyselect_hash_add(H, p, -v, loc->header->m0, loc->header->ad,
				    loc->header->d, loc->header->N, q, rqqz,
                                    loc);

	}			// next rp
    }				// next p

#ifdef DEBUG_POLYSELECT
  fprintf(stderr, "# inner collision_on_each_sq took %lums\n",
	  milliseconds() - st);
  fprintf(stderr, "# - q polyselect_hash_size (q=%lu): %u\n", q, H->size);
#endif

  polyselect_hash_clear(H);

  loc->stats->potential_collisions++;
}


/* Batch SQ mode */
static inline void
gmp_collision_on_batch_sq(uint64_t * q,
			  const mpz_t * qqz,
			  const mpz_t * rqqz,
			  unsigned long size, unsigned long number_pr,
                          polyselect_thread_locals_ptr loc)
{
  if (size == 0)
    return;

  unsigned int i, j, nr;
  unsigned long p, c = 0;
  uint64_t pp, **invqq, rp;
  mpz_t qprod[size], modpp, tmp, tmp1, tmp2, rpmp;

  mpz_init(tmp);
  mpz_init(tmp1);
  mpz_init(tmp2);
  mpz_init(rpmp);
  mpz_init(modpp);
  for (i = 0; i < size; i++)
    mpz_init(qprod[i]);

  invqq = (uint64_t **) malloc(size * sizeof(uint64_t *));
  if (invqq)
    {
      for (i = 0; i < size; i++)
	invqq[i] = malloc(number_pr * sizeof(uint64_t));
  } else
    {
      fprintf(stderr, "Error, cannot allocate memory in %s\n", __FUNCTION__);
      exit(1);
    }

  mpz_set(qprod[0], qqz[0]);
  for (i = 1; i < size; i++)
    mpz_mul(qprod[i], qqz[i], qprod[i - 1]);

  /* Step 1: batch inversion */
  for (unsigned long nprimes = 0; nprimes < loc->main->lenPrimes; nprimes++)
    {

      p = loc->main->Primes[nprimes];
      pp = (uint64_t) p;
      pp *= (uint64_t) p;
      if (polyselect_poly_header_skip(loc->header, p))
	continue;
      if (loc->R->nr[nprimes] == 0)
	continue;
      nr = loc->R->nr[nprimes];

      /* inversion */
      mpz_set_uint64(modpp, pp);
      mpz_invert(tmp1, qprod[size - 1], modpp);

      /* for each (q, r) \in a batch */
      for (i = size - 1; i > 0; i--)
	{
	  mpz_mul(tmp, qprod[i - 1], tmp1);
	  /* for each rp, compute (rp-rq)*1/q^2 (mod p^2) */
	  for (j = 0; j < nr; j++, c++)
	    {
	      rp = loc->R->roots[nprimes][j];
	      mpz_set_uint64(rpmp, rp);
	      mpz_sub(rpmp, rpmp, rqqz[i]);
	      mpz_mul(rpmp, rpmp, tmp);
	      mpz_mod(rpmp, rpmp, modpp);
	      invqq[i][c] = mpz_get_uint64(rpmp);
	    }
	  mpz_mul(tmp, qqz[i], tmp1);
	  mpz_mod(tmp1, tmp, modpp);
	  c -= nr;
	}
      /* last q in the batch is in tmp_modul */
      for (j = 0; j < nr; j++, c++)
	{
	  rp = loc->R->roots[nprimes][j];
	  mpz_set_uint64(rpmp, rp);
	  mpz_sub(rpmp, rpmp, rqqz[0]);
	  mpz_fdiv_r(rpmp, rpmp, modpp);
	  mpz_mul(tmp2, rpmp, tmp1);
	  mpz_mod(tmp2, tmp2, modpp);
	  invqq[0][c] = mpz_get_uint64(tmp2);
	}
    }				// next prime p

  /* Step 2: find collisions on q. */
  for (i = 0; i < size; i++)
    {
      //int st2 = milliseconds();
      gmp_collision_on_each_sq(q[i], rqqz[i], invqq[i], loc);
      //printf ("# outer collision_on_each_sq took %dms\n", milliseconds () - st2);
    }

  for (i = 0; i < size; i++)
    free(invqq[i]);
  free(invqq);
  mpz_clear(tmp);
  mpz_clear(tmp1);
  mpz_clear(tmp2);
  mpz_clear(rpmp);
  mpz_clear(modpp);
  for (i = 0; i < size; i++)
    mpz_clear(qprod[i]);
}


/* collision on special-q, call gmp_collision_on_batch_sq */
void
gmp_collision_on_sq(unsigned long c,
                    polyselect_thread_locals_ptr loc)
{
  // init special-q roots
  unsigned long K, lq;
  polyselect_qroots_t SQ_R;
  polyselect_qroots_init(SQ_R);
  comp_sq_roots(loc->header, SQ_R, loc->rstate);
  // polyselect_qroots_print (SQ_R);

  /* find a suitable lq */
  lq = find_suitable_lq(loc->header, SQ_R, &K, loc->main);

  unsigned long N = lq, tot, i, l, idx_q[K];
  uint64_t q[BATCH_SIZE];
  mpz_t *qqz, *rqqz;

  qqz = (mpz_t *) malloc(BATCH_SIZE * sizeof(mpz_t));
  rqqz = (mpz_t *) malloc(BATCH_SIZE * sizeof(mpz_t));
  if (!qqz || !rqqz)
    {
      fprintf(stderr, "Error, cannot allocate memory "
	      "in gmp_collision_on_sq \n");
      exit(1);
    }
  for (l = 0; l < BATCH_SIZE; l++)
    {
      mpz_init(qqz[l]);
      mpz_init(rqqz[l]);
    }

  // less than lq special primes having roots for this ad
  if (N == 0 || N < K)
    {
      gmp_fprintf(stderr, "# Info: binomial(%lu, %lu) error in "
		  "collision_on_sq(). ad=%Zd.\n", N, K, loc->header->ad);
      return;
    }

  tot = binom(N, K);

  if (tot > loc->main->nq)
    tot = loc->main->nq;

  if (tot < BATCH_SIZE)
    tot = BATCH_SIZE;

#ifdef DEBUG_POLYSELECT
  fprintf(stderr, "# Info: n=%lu, k=%lu, (n,k)=%lu"
	  ", maxnq=%lu, nq=%lu\n", N, K, binom(N, K), nq, tot);
#endif

  i = 0;
  while (i <= (tot - BATCH_SIZE))
    {

      l = i;			// why do I use an extra l here?
      if (l == 0)
	{

	  // enumerate first combination
	  first_comb(K, idx_q);
	  //print_comb (K, idx_q);
	  q[l] = return_q_rq(SQ_R, idx_q, K, qqz[l], rqqz[l]);

	  for (l = 1; l < BATCH_SIZE; l++)
	    {
	      next_comb(N, K, idx_q);
	      q[l] = return_q_rq(SQ_R, idx_q, K, qqz[l], rqqz[l]);
	    }
      } else
	{
	  for (l = 0; l < BATCH_SIZE; l++)
	    {
	      next_comb(N, K, idx_q);
	      q[l] = return_q_rq(SQ_R, idx_q, K, qqz[l], rqqz[l]);
	    }
	}

#ifdef DEBUG_POLYSELECT
      unsigned long j;
      for (j = 0; j < BATCH_SIZE; j++)
	gmp_fprintf(stderr, "q: %lu, qq: %Zd, rqq: %Zd\n",
		    q[j], qqz[j], rqqz[j]);
#endif

      // collision batch
      gmp_collision_on_batch_sq(q, (const mpz_t *) qqz,
				(const mpz_t *) rqqz, BATCH_SIZE, c, loc);
      i += BATCH_SIZE;
    }

  // tail batch
  for (l = 0; l < (tot % BATCH_SIZE); l++)
    {
      next_comb(N, K, idx_q);
      q[l] = return_q_rq(SQ_R, idx_q, K, qqz[l], rqqz[l]);

#ifdef DEBUG_POLYSELECT
      gmp_fprintf(stderr, "q: %lu, qq: %Zd, rqq: %Zd\n",
		  q[l], qqz[l], rqqz[l]);
#endif

    }

  gmp_collision_on_batch_sq(q, (const mpz_t *) qqz,
			    (const mpz_t *) rqqz, tot % BATCH_SIZE, c, loc);

  for (l = 0; l < BATCH_SIZE; l++)
    {
      mpz_clear(qqz[l]);
      mpz_clear(rqqz[l]);
    }
  free(qqz);
  free(rqqz);
  polyselect_qroots_clear(SQ_R);
}

