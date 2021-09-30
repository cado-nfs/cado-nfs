#include "cado.h"		// IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>		// malloc ...
#include <stdint.h>		// uint64_t
#include <limits.h>		/* for CHAR_BIT */
#include <math.h>               // log
#include <gmp.h>
#include "auxiliary.h"
#include "cado_poly.h"
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
#include "polyselect_main_queue.h"
#include "portability.h"
#include "roots_mod.h"
#include "polyselect_collisions.h"
#include "timing.h"		// for seconds
#include "usp.h"		// usp_root_data
#include "verbose.h"		// verbose_output_print

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
static inline void polyselect_proots_dispatch_to_hash_notflat(
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
static inline void polyselect_proots_dispatch_to_hash_flat(
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

unsigned long compute_and_lift_proots(polyselect_thread_locals_ptr loc)/*{{{*/
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
            R->nr[nprimes] = 0;        // nr = 0.
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
collision_on_p(polyselect_shash_ptr H, polyselect_hash_match_t match, polyselect_thread_locals_ptr loc)
{
  unsigned long j, p, nrp, tot_roots = 0;
  uint64_t *rp;
  int64_t ppl = 0, u, umax;
  int found = 0;
  int st = milliseconds();

  rp = (uint64_t *) malloc(loc->header->d * sizeof(uint64_t));
  if (rp == NULL)
    {
      fprintf(stderr, "Error, cannot allocate memory in collision_on_p\n");
      exit(1);
    }

  /* We first store only i's (and not p's), and look for collisions,
   * which occur very rarely.
   * if we find out that there is a collision on i, we run the search
   * again for the p's.
   */
  polyselect_shash_reset(H);
  umax = polyselect_main_data_get_M(loc->main);
  for (unsigned long nprimes = 0; nprimes < loc->main->lenPrimes; nprimes++)
    {
      p = loc->main->Primes[nprimes];
      ppl = (int64_t) p *(int64_t) p;

      /* add fake roots to keep indices */
      if (polyselect_poly_header_skip(loc->header, p))
	{
	  loc->R->nr[nprimes] = 0;	// nr = 0.
	  loc->R->roots[nprimes] = NULL;
	  continue;
	}

      nrp = roots_mod_uint64(rp, mpz_fdiv_ui(loc->header->Ntilde, p), loc->header->d,
			     p, loc->rstate);
      tot_roots += nrp;
      nrp = roots_lift(rp, loc->header->Ntilde, loc->header->d, loc->header->m0, p, nrp);
      polyselect_proots_add(loc->R, nrp, rp, nprimes);
      for (j = 0; j < nrp; j++)
	{
	  for (u = (int64_t) rp[j]; u < umax; u += ppl)
	    polyselect_shash_add(H, u);
	  for (u = ppl - (int64_t) rp[j]; u < umax; u += ppl)
	    polyselect_shash_add(H, -u);
	}
    }
  st = milliseconds() - st;

  if (loc->main->verbose > 2)
    {
      fprintf(stderr, "# computing %lu p-roots took %dms\n", tot_roots, st);
      fprintf(stderr,
	      "# polyselect_shash_size (umax = %" PRId64 ", P = %lu"
	      "): %zu\n", umax, loc->main->P, polyselect_shash_size(H));
      fprintf(stderr, "# expected number of pairs: %zu\n",
	      polyselect_main_data_expected_number_of_pairs(loc->main));
    }

  st = milliseconds();
  found = polyselect_shash_find_collision(H);
  free(rp);

  if (loc->main->verbose > 2)
    {
      fprintf(stderr,
	      "# collision found in shash: %d (probability = 1 / %.1f) [took %dms]\n",
	      found, 4 * log(loc->main->P) * log(loc->main->P), st);
    }


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


/* collision on each special-q, call collision_on_batch_p() */
static inline void
collision_on_each_sq(unsigned long q,
		     mpz_srcptr rqqz,
		     unsigned long *inv_qq,
                     polyselect_shash_ptr H,
                     polyselect_hash_match_t match,
                     polyselect_thread_locals_ptr loc)
{
  int found;

#ifdef DEBUG_POLYSELECT2
  int st = milliseconds();
#endif

  int64_t umax = polyselect_main_data_get_M(loc->main);

#if 1
  loc->R->nr[loc->R->size] = 0xff;     /* use guard to end */
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
                       polyselect_hash_match_t match,
                       polyselect_thread_locals_ptr loc)
{
  if (count == 0)
    return;

  uint8_t i, nr, *pnr;
  unsigned long p, c = 0, rp, rqi;
  int k;
  uint64_t pp;
  unsigned long **tinv_qq = malloc(count * sizeof(unsigned long *));

  if (!tinv_qq)
    {
      fprintf(stderr, "Error, cannot allocate memory in %s\n", __FUNCTION__);
      exit(1);
    }
  for (k = 0; k < count; k++)
    {
      /* number_pr + 1 for guard for pre-load in collision_on_each_sq (nv) */
      tinv_qq[k] = malloc((number_pr + 1) * sizeof(unsigned long));
      tinv_qq[k][number_pr] = 0;
    }

  int st = milliseconds();
  pnr = loc->R->nr;

  /* for each rp, compute (rp-rq)*1/q^2 (mod p^2) */
  for (unsigned long nprimes = 0; nprimes < loc->main->lenPrimes; nprimes++)
    {
      if (!pnr[nprimes])
	continue;
      nr = pnr[nprimes];
      p = loc->main->Primes[nprimes];
      pp = p * p;

      modulusredcul_t modpp;
      residueredcul_t res_rqi, res_rp, res_tmp;
      modredcul_initmod_ul_raw(modpp, pp);
      modredcul_init(res_rqi, modpp);
      modredcul_init(res_rp, modpp);
      modredcul_init(res_tmp, modpp);

      for (k = 0; k < count; k++)
	{
	  rqi = mpz_fdiv_ui(rqqz[k], pp);
	  modredcul_intset_ul(res_rqi, rqi);
	  modredcul_intset_ul(res_tmp, inv_qq[nprimes]);
	  for (i = 0; i < nr; i++, c++)
	    {
	      rp = loc->R->roots[nprimes][i];
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
  for (k = 0; k < count; k++)
    collision_on_each_sq(q, rqqz[k], tinv_qq[k], H, match, loc);

  if (loc->main->verbose > 2)
    fprintf(stderr,
	    "#  substage: collision-detection %d many rq took %lums\n",
	    count, milliseconds() - st);

  for (k = 0; k < count; k++)
    free(tinv_qq[k]);
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
                        polyselect_hash_match_t match,
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
			     number_pr, num_rq, H, match, loc);
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
collision_on_sq(unsigned long c, polyselect_shash_ptr H, polyselect_hash_match_t match, polyselect_thread_locals_ptr loc)
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
                  &curr_nq, k, H, match, loc);
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

