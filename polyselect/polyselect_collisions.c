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
#include "polyselect_main_queue.h"
#include "portability.h"
#include "roots_mod.h"
#include "polyselect_collisions.h"
#include "polyselect_thread.h"
#include "timing.h"		// for seconds
#include "usp.h"		// usp_root_data
#include "verbose.h"		// verbose_output_print

/* CCS = computational collision search
 * DCS = decisional collision search
 */

/*{{{ polyselect_proots_dispatch_to_shash_flat_ugly */

/* This code dispatches congruence classes of the roots_per_prime[] table,
 * modulo the primes that are listed in Primes, within the range
 * [-umax...umax], into the quick hash table H.
 *
 * This overly complicated code provides some tiny performance gain, but
 * the simpler version below is fine too.
 *
 */
void polyselect_proots_dispatch_to_shash_flat_ugly(
        polyselect_shash_ptr H,
        const uint32_t * Primes,
        size_t lenPrimes,
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

  // reset has been done by the caller with polyselect_shash_reset_multi
  // polyselect_shash_reset(H);

  pc = (long *) roots_per_prime;
  nv = *pc;
  const uint32_t * pprimes = Primes - 1;
  const uint8_t * pnr = number_of_roots_per_prime;
  const uint8_t * pnr_end = pnr + lenPrimes;
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
      for( ; ; ) {
          if (UNLIKELY(pnr == pnr_end))
              goto bend;
          pprimes++;
          if ((vpnr = *pnr++) != 0)
              break;
      }
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
      for( ; ; ) {
          if (UNLIKELY(pnr == pnr_end))
              goto bend;
          pprimes++;
          if ((vpnr = *pnr++) != 0)
              break;
      }
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
      for( ; ; ) {
          if (UNLIKELY(pnr == pnr_end))
              goto bend;
          pprimes++;
          if ((vpnr = *pnr++) != 0)
              break;
      }
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
    // reset has been done by the caller with polyselect_shash_reset_multi
    // polyselect_shash_reset(H);
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

/*{{{ polyselect_proots_dispatch_to_shash_flat */
void polyselect_proots_dispatch_to_shash2_flat(
        polyselect_shash_ptr H,
        const uint32_t * Primes,
        size_t lenPrimes,
        const unsigned long * roots_per_prime,
        const uint8_t * number_of_roots_per_prime,
        int64_t umax)
{
    // reset has been done by the caller with polyselect_shash_reset_multi
    // polyselect_shash_reset(H);
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
                polyselect_shash2_add(H, u, p);
            for(int64_t u = u0 - ppl ; u + umax >= 0; u -= ppl)
                polyselect_shash2_add(H, u, p);
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
    // reset has been done by the caller with polyselect_shash_reset_multi
    // polyselect_shash_reset(H);
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

/*{{{ polyselect_proots_dispatch_to_shash2_notflat
 *
 * touch H->pmem as well, and store the prime info.
 * */
void polyselect_proots_dispatch_to_shash2_notflat(
        polyselect_shash_ptr H,
        const uint32_t * Primes,
        size_t lenPrimes,
        uint64_t * const * roots_per_prime,
        const uint8_t * number_of_roots_per_prime,
        int64_t umax)
{
    // reset has been done by the caller with polyselect_shash_reset_multi
    // polyselect_shash_reset(H);
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
                polyselect_shash2_add(H, u, p);
            for(int64_t u = u0 - ppl ; u + umax >= 0 ; u -= ppl)
                polyselect_shash2_add(H, u, p);
        }
    }

    for (int i = 0; i < polyselect_SHASH_NBUCKETS; i++)
        ASSERT(H->current[i] <= H->base[i + 1]);
}/*}}}*/

struct polyselect_CCS_subtask_data {
    unsigned long q;
    mpz_srcptr rq;
    /* only for the flat version (sq) */
    unsigned long * invq_roots_per_prime;
};

/*{{{ polyselect_CCS_notflat_subtask
 *
 * This is more costly than the decisional version, and needs only be
 * called rarely.
 *
 * Acts upon all collisions in the table.
 */
void polyselect_CCS_notflat_subtask(polyselect_thread_ptr thread)
{
    struct polyselect_CCS_subtask_data * arg = thread->team->task->arg;
    unsigned long q = arg->q;
    mpz_srcptr rq = arg->rq;

    unsigned int nt = thread->team->task->expected;
    unsigned int it = thread->index_in_sync_zone;
    polyselect_primes_table_srcptr pt = thread->team->league->pt;

    polyselect_shash_ptr SH = thread->team->SH[it];

    polyselect_thread_team_enter_roaming(thread->team, thread);
    pthread_mutex_unlock(&thread->team->lock);
    /********* BEGIN UNLOCKED SECTION **************/
    polyselect_thread_chronogram_chat(thread, "enter dispatch_shash2_nf");

    // reset has been done by the caller with polyselect_shash_reset_multi
    // polyselect_shash_reset(H);
    {
        size_t qt = pt->lenPrimes / nt;
        size_t rt = pt->lenPrimes % nt;
        unsigned long i0 = qt * it + MIN(it, rt);
        unsigned long i1 = i0 + qt + (it < rt);
        polyselect_proots_dispatch_to_shash2_notflat(SH,
                pt->Primes + i0,
                i1 - i0,
                thread->team->R->roots + i0,
                thread->team->R->nr + i0,
                polyselect_main_data_get_M(thread->team->league->main));
    }

    polyselect_thread_chronogram_chat(thread, "leave dispatch_shash2_nf");

    polyselect_thread_team_roaming_barrier(thread->team, thread);
    // barrier_wait(&thread->team->sync_task->barrier, NULL, NULL, NULL);

    polyselect_thread_chronogram_chat(thread, "enter transverse_shash2");
    {
        /* which of the buckets do we have to scan for collisions ? */
        unsigned int wt = polyselect_SHASH_NBUCKETS;
        unsigned int rt = wt % nt;
        unsigned int k0 = it * (wt / nt) + MIN(it, rt);
        unsigned int k1 = k0 + (wt / nt) + (it < rt);

        polyselect_shash2_find_collision_multi(thread->team->SH, nt, k0, k1,
                q, rq, thread);
    }

    polyselect_thread_chronogram_chat(thread, "leave transverse_shash2");
    /********** END UNLOCKED SECTION ***************/
    pthread_mutex_lock(&thread->team->lock);
    polyselect_thread_team_leave_roaming(thread->team, thread);
}/*}}}*/

/*{{{ polyselect_DCS_notflat_subtask
 *
 * quickly tell whether there is a collision or not.
 *
 * boolean result is stored in int * thread->team->sync_task->arg
 */
void polyselect_DCS_notflat_subtask(polyselect_thread_ptr thread)
{
    unsigned int nt = thread->team->task->expected;
    unsigned int it = thread->index_in_sync_zone;
    polyselect_primes_table_srcptr pt = thread->team->league->pt;

    /* thread->team->SH[it] might be different from one call to the next,
     * and this is perhaps a source for bad performance if the team is
     * large, since we're killing the L2 cache. This could be fixed, but
     * my reasoning is that it's not necessarily worth the trouble since
     * the transversal pass thrashes these caches in a similar way
     * anyway.
     */
    polyselect_shash_ptr SH = thread->team->SH[it];

    polyselect_thread_team_enter_roaming(thread->team, thread);
    pthread_mutex_unlock(&thread->team->lock);
    /********* BEGIN UNLOCKED SECTION **************/
    polyselect_thread_chronogram_chat(thread, "enter dispatch_shash_nf");

    // reset has been done by the caller with polyselect_shash_reset_multi
    // polyselect_shash_reset(H);
    {
        size_t qt = pt->lenPrimes / nt;
        size_t rt = pt->lenPrimes % nt;
        unsigned long i0 = qt * it + MIN(it, rt);
        unsigned long i1 = i0 + qt + (it < rt);
        polyselect_proots_dispatch_to_shash_notflat(SH,
                pt->Primes + i0,
                i1 - i0,
                thread->team->R->roots + i0,
                thread->team->R->nr + i0,
                polyselect_main_data_get_M(thread->team->league->main));
    }

    polyselect_thread_chronogram_chat(thread, "leave dispatch_shash_nf");

    polyselect_thread_team_roaming_barrier(thread->team, thread);
    // barrier_wait(&thread->team->sync_task->barrier, NULL, NULL, NULL);

    polyselect_thread_chronogram_chat(thread, "enter transverse_shash");
    int found;

    {
        /* which of the buckets do we have to scan for collisions ? */
        unsigned int wt = polyselect_SHASH_NBUCKETS;
        unsigned int rt = wt % nt;
        unsigned int k0 = it * (wt / nt) + MIN(it, rt);
        unsigned int k1 = k0 + (wt / nt) + (it < rt);

        found = polyselect_shash_find_collision_multi(thread->team->SH, nt, k0, k1);
    }

    polyselect_thread_chronogram_chat(thread, "leave transverse_shash");
    /********** END UNLOCKED SECTION ***************/
    pthread_mutex_lock(&thread->team->lock);
    polyselect_thread_team_leave_roaming(thread->team, thread);

    if (found) {
        *((int*)thread->team->task->arg) = 1;
    }
}
/*}}}*/

/* {{{ collision_on_p_conductor
 *
 * find collisions between "P" primes, return number of loops */
unsigned long
collision_on_p_conductor(polyselect_thread_ptr thread)
{
  /* first compute and lift all roots modulo the primes in
   * thread->main->Primes ; we store that in thread->team->R
   */
  unsigned long tot_roots = polyselect_proots_compute_conductor(thread);

  /* We first store only i's (and not p's), and look for collisions,
   * which occur very rarely. Each thread does that for its range of
   * primes. In a second pass, the buckets of all threads are examined
   * collectively for collisions.
   *
   * if we find out that there is a collision on i, we run the search
   * again for the p's.
   */

  int found = 0;

  /* This is called with the team lock held ! */

  polyselect_shash_reset_multi(thread->team->SH, thread->team->count->sync);
  polyselect_thread_team_post_work(thread->team, thread, polyselect_DCS_notflat_subtask, &found);


  if (found) {/* do the real work */
      struct polyselect_CCS_subtask_data arg[1] = {{ .q = 1, .rq = NULL }};
      polyselect_shash_reset_multi(thread->team->SH, thread->team->count->sync);
      polyselect_thread_team_post_work(thread->team, thread, polyselect_CCS_notflat_subtask, arg);
  }

  thread->stats->potential_collisions++;

  return tot_roots;
}  /* }}} */

struct polyselect_DCS_flat_subtask_data {
    unsigned long * invq_roots_per_prime;
    int found;
};

/*{{{ polyselect_DCS_flat_subtask
 *
 * quickly tell whether there is a collision or not.
 *
 * boolean result is stored in int * thread->team->sync_task->arg
 */
void polyselect_DCS_flat_subtask(polyselect_thread_ptr thread)
{
    struct polyselect_DCS_flat_subtask_data * arg = thread->team->task->arg;
    unsigned long * invq_roots_per_prime = arg->invq_roots_per_prime;

    unsigned int nt = thread->team->task->expected;
    unsigned int it = thread->index_in_sync_zone;
    polyselect_primes_table_srcptr pt = thread->team->league->pt;

    /* thread->team->SH[it] might be different from one call to the next,
     * and this is perhaps a source for bad performance if the team is
     * large, since we're killing the L2 cache. This could be fixed, but
     * my reasoning is that it's not necessarily worth the trouble since
     * the transversal pass thrashes these caches in a similar way
     * anyway.
     */
    polyselect_shash_ptr SH = thread->team->SH[it];

    polyselect_thread_team_enter_roaming(thread->team, thread);
    pthread_mutex_unlock(&thread->team->lock);
    /********* BEGIN UNLOCKED SECTION **************/
    polyselect_thread_chronogram_chat(thread, "enter dispatch_shash_f");

    // reset has been done by the caller with polyselect_shash_reset_multi
    // polyselect_shash_reset(SH);
    {
        size_t qt = pt->lenPrimes / nt;
        size_t rt = pt->lenPrimes % nt;
        unsigned long i0 = qt * it + MIN(it, rt);
        unsigned long i1 = i0 + qt + (it < rt);
        /* This is unfortunate. It's caused by the "flat" layout. */
        unsigned int z = 0;
        for(unsigned int i = 0 ; i < i0 ; i++) {
            z += thread->team->R->nr[i];
        }
        polyselect_proots_dispatch_to_shash_flat(SH,
                pt->Primes + i0,
                i1 - i0,
                invq_roots_per_prime + z,
                thread->team->R->nr + i0,
                polyselect_main_data_get_M(thread->team->league->main));
    }

    polyselect_thread_chronogram_chat(thread, "leave dispatch_shash_f");

    polyselect_thread_team_roaming_barrier(thread->team, thread);
    //// barrier_wait(&thread->team->sync_task->barrier, NULL, NULL, NULL);

    polyselect_thread_chronogram_chat(thread, "enter transverse_shash");
    int found;

    {
        /* which of the buckets do we have to scan for collisions ? */
        unsigned int wt = polyselect_SHASH_NBUCKETS;
        unsigned int rt = wt % nt;
        unsigned int k0 = it * (wt / nt) + MIN(it, rt);
        unsigned int k1 = k0 + (wt / nt) + (it < rt);

        found = polyselect_shash_find_collision_multi(thread->team->SH, nt, k0, k1);
    }

    polyselect_thread_chronogram_chat(thread, "leave transverse_shash");
    /********** END UNLOCKED SECTION ***************/
    pthread_mutex_lock(&thread->team->lock);
    polyselect_thread_team_leave_roaming(thread->team, thread);

    if (found)
        arg->found = 1;
}
/*}}}*/

/*{{{ polyselect_CCS_flat_subtask
 *
 * This is more costly than the decisional version, and needs only be
 * called rarely.
 *
 * Acts upon all collisions in the table.
 */
void polyselect_CCS_flat_subtask(polyselect_thread_ptr thread)
{
    struct polyselect_CCS_subtask_data * arg = thread->team->task->arg;
    unsigned long q = arg->q;
    mpz_srcptr rq = arg->rq;
    unsigned long * invq_roots_per_prime = arg->invq_roots_per_prime;

    unsigned int nt = thread->team->task->expected;
    unsigned int it = thread->index_in_sync_zone;
    polyselect_primes_table_srcptr pt = thread->team->league->pt;

    polyselect_shash_ptr SH = thread->team->SH[it];

    polyselect_thread_team_enter_roaming(thread->team, thread);
    pthread_mutex_unlock(&thread->team->lock);
    /********* BEGIN UNLOCKED SECTION **************/
    polyselect_thread_chronogram_chat(thread, "enter dispatch_shash2_nf");

    // reset has been done by the caller with polyselect_shash_reset_multi
    // polyselect_shash_reset(SH);
    {
        size_t qt = pt->lenPrimes / nt;
        size_t rt = pt->lenPrimes % nt;
        unsigned long i0 = qt * it + MIN(it, rt);
        unsigned long i1 = i0 + qt + (it < rt);
        /* This is unfortunate. It's caused by the "flat" layout. */
        unsigned int z = 0;
        for(unsigned int i = 0 ; i < i0 ; i++) {
            z += thread->team->R->nr[i];
        }
        polyselect_proots_dispatch_to_shash2_flat(SH,
                pt->Primes + i0,
                i1 - i0,
                invq_roots_per_prime + z,
                thread->team->R->nr + i0,
                polyselect_main_data_get_M(thread->team->league->main));
    }

    polyselect_thread_chronogram_chat(thread, "leave dispatch_shash2_nf");

    polyselect_thread_team_roaming_barrier(thread->team, thread);
    // barrier_wait(&thread->team->sync_task->barrier, NULL, NULL, NULL);

    polyselect_thread_chronogram_chat(thread, "enter transverse_shash2");
    {
        /* which of the buckets do we have to scan for collisions ? */
        unsigned int wt = polyselect_SHASH_NBUCKETS;
        unsigned int rt = wt % nt;
        unsigned int k0 = it * (wt / nt) + MIN(it, rt);
        unsigned int k1 = k0 + (wt / nt) + (it < rt);

        polyselect_shash2_find_collision_multi(thread->team->SH, nt, k0, k1,
                q, rq, thread);
    }

    polyselect_thread_chronogram_chat(thread, "leave transverse_shash2");
    /********** END UNLOCKED SECTION ***************/
    pthread_mutex_lock(&thread->team->lock);
    polyselect_thread_team_leave_roaming(thread->team, thread);
}/*}}}*/
/* collision on each special-q, call collision_on_batch_p() */
static inline void
collision_on_each_sq(unsigned long q,
		     mpz_srcptr rq,
		     unsigned long * invq_roots_per_prime,
                     polyselect_thread_ptr thread)
{
    struct polyselect_DCS_flat_subtask_data arg[1] = {{
        .invq_roots_per_prime = invq_roots_per_prime,
        .found = 0
    }};

    polyselect_shash_reset_multi(thread->team->SH, thread->team->count->sync);
    polyselect_thread_team_post_work(thread->team, thread, polyselect_DCS_flat_subtask, arg);

    if (arg->found) {/* do the real work */
      struct polyselect_CCS_subtask_data arg[1] = {{
          .q = q, .rq = rq, .invq_roots_per_prime = invq_roots_per_prime }};
      polyselect_shash_reset_multi(thread->team->SH, thread->team->count->sync);
      polyselect_thread_team_post_work(thread->team, thread, polyselect_CCS_flat_subtask, arg);
  }

    /*
    fprintf(stderr, "thread %d exits scope before STOP for %d sync thread in team %d\n",
            thread->thread_index, thread->team->count->sync, thread->team->team_index);
            */
  thread->stats->potential_collisions++;
}

struct modcalc_arg {
    const mpz_t * rqqz;
    int count;
    const unsigned long *inv_qq;
    unsigned long ** tinv_qq;
};

unsigned long modcalc_nroots_interval(polyselect_proots_srcptr R, unsigned long i0, unsigned long i1)
{
    unsigned long c = 0;
    for (unsigned long i = i0; i < i1; i++)
        c += R->nr[i];
    return c;
}

void modcalc_subtask(polyselect_thread_ptr thread)/*{{{*/
{
    struct modcalc_arg * arg = thread->team->task->arg;
    int count = arg->count;
    const mpz_t * rqqz = arg->rqqz;
    const unsigned long *inv_qq = arg->inv_qq;
    unsigned long **tinv_qq = arg->tinv_qq;
    unsigned long c = 0;

    unsigned int nt = thread->team->task->expected;
    unsigned int it = thread->index_in_sync_zone;
    polyselect_primes_table_ptr pt = thread->team->league->pt;

    polyselect_thread_team_enter_roaming(thread->team, thread);
    pthread_mutex_unlock(&thread->team->lock);
    /********* BEGIN UNLOCKED SECTION **************/
    polyselect_thread_chronogram_chat(thread, "enter modcalc");

#ifdef DEBUG_POLYSELECT_THREADS
    fprintf(stderr, "enter modcalc with %d threads\n", nt);
#endif

    size_t qt = pt->lenPrimes / nt;
    size_t rt = pt->lenPrimes % nt;
    unsigned long i0 = qt * it + MIN(it, rt);
    unsigned long i1 = i0 + qt + (it < rt);

    /* This is unfortunate. It's caused by the "flat" layout. */
    for (unsigned long i = 0; i < i0; i++) {
        uint8_t nr = thread->team->R->nr[i];
        c += nr;
    }
    for (unsigned long i = i0; i < i1; i++) {
        uint8_t nr = thread->team->R->nr[i];
        if (!nr)
            continue;
        uint32_t p = pt->Primes[i];
        uint64_t pp = (int64_t) p *(int64_t) p;

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
            modredcul_intset_ul(res_tmp, inv_qq[i]);
            for (uint8_t j = 0; j < nr; j++, c++)
            {
                unsigned long rp = thread->team->R->roots[i][j];
                modredcul_intset_ul(res_rp, rp);
                /* rp - rq */
                modredcul_sub(res_rp, res_rp, res_rqi, modpp);
                /* res_rp = (rp - rq) / q[j]^2 */
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
    polyselect_thread_chronogram_chat(thread, "leave modcalc");
    /********** END UNLOCKED SECTION ***************/
    pthread_mutex_lock(&thread->team->lock);
    polyselect_thread_team_leave_roaming(thread->team, thread);
}/*}}}*/

/* Given p, rp, q, invqq[], for each rq of q, compute (rp - rq) / q^2 */
static inline void
collision_on_each_sq_r(unsigned long q,
		       const mpz_t * rqqz,
		       const unsigned long *inv_qq,
		       unsigned long number_pr,
		       int count,
                       polyselect_thread_ptr thread)
{
  if (count == 0)
    return;

  polyselect_thread_chronogram_chat(thread, "enter malloc");

  unsigned long **tinv_qq = malloc(count * sizeof(unsigned long *));

  if (!tinv_qq)
    {
      fprintf(stderr, "Error, cannot allocate memory in %s\n", __func__);
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

  struct modcalc_arg arg[1] = {{
          .count = count,
          .rqqz = rqqz,
          .inv_qq = inv_qq,
          .tinv_qq = tinv_qq
  }};
  polyselect_thread_chronogram_chat(thread, "leave malloc");
  polyselect_thread_team_post_work(thread->team, thread, modcalc_subtask, arg);

  /* core function to find collisions */
  for (int k = 0; k < count; k++)
    collision_on_each_sq(q, rqqz[k], tinv_qq[k], thread);

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
                        polyselect_thread_ptr thread)
{
  unsigned int ind_qr[k];	/* indices of roots for each small q */
  unsigned int len_qnr[k];	/* for each small q, number of roots */
  mpz_t qqz, rqqz[BATCH_SIZE];

  mpz_init(qqz);
  for (unsigned long i = 0; i < BATCH_SIZE; i++)
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

  /* we proceed with BATCH_SIZE many rq each time */
  for (unsigned long i = 0; i < k; i++)
    {
      ind_qr[i] = 0;
      len_qnr[i] = SQ_R->nr[idx_q[i]];
    }
  for(int re = 1 ; re ; )
    {
      /* compute BATCH_SIZE such many rqqz[] */
      int num_rq = 0;
      for (int count = 0; re && count < BATCH_SIZE; count++)
	{
          /* We do multiple inversions over and over again (in crt_sq),
           * even though the set of q's doesn't change... This has no
           * noticeable impact on the running time, though.
           */
	  aux_return_rq(SQ_R, idx_q, ind_qr, k, qqz, rqqz[count]);
	  re = aux_nextcomb(ind_qr, k, len_qnr);
	  (*curr_nq)++;
	  num_rq++;
	  if ((*curr_nq) >= thread->team->league->main->nq)
	    re = 0;
	}

      /* core function for a fixed qq and several rqqz[] */
      collision_on_each_sq_r(q, (const mpz_t *) rqqz, inv_qq,
			     number_pr, num_rq, thread);
    }

  mpz_clear(qqz);
  for (unsigned long i = 0; i < BATCH_SIZE; i++)
    mpz_clear(rqqz[i]);
}

struct invert_q2_mod_all_p2_data {
        unsigned long *invqq;
        unsigned long q;
};

static inline void invert_q2_mod_all_p2_subtask(polyselect_thread_ptr thread) /*{{{*/
{
    struct invert_q2_mod_all_p2_data * arg = thread->team->task->arg;
    unsigned long q = arg->q;
    unsigned long *invqq = arg->invqq;

    polyselect_thread_team_enter_roaming(thread->team, thread);
    pthread_mutex_unlock(&thread->team->lock);
    /********* BEGIN UNLOCKED SECTION **************/
    polyselect_thread_chronogram_chat(thread, "enter invert_q2");
    polyselect_primes_table_srcptr pt = thread->team->league->pt;
    unsigned int nt = thread->team->task->expected;
    unsigned int it = thread->index_in_sync_zone;
    size_t qt = pt->lenPrimes / nt;
    size_t rt = pt->lenPrimes % nt;
    unsigned long i0 = qt * it + MIN(it, rt);
    unsigned long i1 = i0 + qt + (it < rt);

    const uint32_t * Primes = pt->Primes;
    polyselect_poly_header_srcptr header = thread->team->header;
    const uint8_t * number_of_roots_per_prime = thread->team->R->nr;
    for (unsigned long i = i0; i < i1; i++)
    {
        unsigned long p = Primes[i];
        if (polyselect_poly_header_skip(header, p))
            continue;
        unsigned int nr = number_of_roots_per_prime[i];
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
        invqq[i] = modredcul_intget_ul(tmp);

        modredcul_clear(tmp, modpp);
        modredcul_clear(qq, modpp);
        modredcul_clearmod(modpp);
    }
    polyselect_thread_chronogram_chat(thread, "leave invert_q2");
    /********** END UNLOCKED SECTION ***************/
    pthread_mutex_lock(&thread->team->lock);
    polyselect_thread_team_leave_roaming(thread->team, thread);
}/*}}}*/

/* collision on special-q, call collision_on_batch_sq */
void
collision_on_sq_conductor(unsigned long c, polyselect_thread_ptr thread)
{
  unsigned long *invqq = malloc(thread->team->league->pt->lenPrimes * sizeof(unsigned long));
  if (!invqq)
    {
      fprintf(stderr, "Error, cannot allocate memory in %s\n", __func__);
      exit(1);
    }

  /* init special-q roots */
  /* This is rather trivial (we only have a few dozen special-q's). At
   * this point, we don't feel a compelling need to parallelize it
   */
  polyselect_qroots_t SQ_R;
  polyselect_qroots_init(SQ_R);
  polyselect_thread_chronogram_chat(thread, "enter qroots");
  comp_sq_roots(thread->team->header, SQ_R, thread->team->rstate);
  polyselect_thread_chronogram_chat(thread, "leave qroots");


  //polyselect_qroots_print (SQ_R);


  /* find a suitable lq */
  unsigned long k;
  unsigned long lq = find_suitable_lq(thread->team->header, SQ_R, &k, thread->team->league->main);

  unsigned long q, idx_q[lq], curr_nq = 0;

  first_comb(k, idx_q);

  for ( ; curr_nq < thread->team->league->main->nq ; )
    {
      q = return_q_norq(SQ_R, idx_q, k);

      /* collision batch */

      /* Step 1: inversion; compute 1/q^2 (mod p_i^2) to invqq[i] */
      struct invert_q2_mod_all_p2_data arg[1] = {{ .q = q, .invqq = invqq }};
      polyselect_thread_team_post_work(thread->team, thread, invert_q2_mod_all_p2_subtask, arg);
 
      {
          /* Step 2: find collisions on q. */
          int st2 = milliseconds();

          collision_on_batch_sq_r(SQ_R, q, idx_q, invqq, c,
                  &curr_nq, k, thread);
          if (thread->team->league->main->verbose > 2)
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
