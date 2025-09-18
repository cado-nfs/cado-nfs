/*
  Polynomial selection using Kleinjung's algorithm (cf slides presented
  at the CADO Workshop in October 2008, Nancy, France).

  [1. Run and parameters]

  To see parameters and their meaning, run without any argument.

  Please report bugs to the Bug Tracking System on:
  https://gitlab.inria.fr/cado-nfs/cado-nfs/-/issues

  [read this first](https://sympa.inria.fr/sympa/arc/cado-nfs/2020-10/msg00006.html).
*/

#include "cado.h" // IWYU pragma: keep

#define EMIT_ADDRESSABLE_shash_add

/* The following avoids to put #ifdef HAVE_OPENMP ... #endif around each
 * OpenMP pragma. It should come after cado.h, which sets -Werror=all.
 *
#ifdef  __GNUC__
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif
 *
 * unfortunately, while it looks like a reasonable thing to do in theory,
 * it's gcc specific. We can't expect such a thing to work with other
 * compilers.
 */
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include <pthread.h>

#include <gmp.h>

#include "macros.h"
#include "mpz_poly.h"
#include "params.h"
#include "polyselect_data_series.h"
#include "polyselect_main_data.h"
#include "polyselect_main_queue.h"
#include "polyselect_poly_header.h"
#include "polyselect_priority_queue.h"
#include "polyselect_thread_league.h"
#include "polyselect_thread_team.h"
#include "polyselect_collisions.h"
#include "polyselect_shash.h"
#include "polyselect_stats.h"
#include "polyselect_match.h"
#include "polyselect_norms.h"
#include "polyselect_thread.h"
#include "polyselect_alpha.h"
#include "polyselect_special_q.h"
#include "portability.h"
#include "size_optimization.h"
#include "timing.h"		// for seconds
#include "verbose.h"		// verbose_output_print
#include "getprime.h"
#include "dllist.h"
#include "auxiliary.h"

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


void
polyselect_process_match_async(polyselect_thread_league_srcptr league, polyselect_stats_ptr stats, polyselect_match_info_ptr job)
{
  polyselect_poly_header_srcptr header = job->header;
  unsigned long p1 = job->p1;
  unsigned long p2 = job->p2;
  const int64_t i = job->i;
  uint64_t q = job->q;
  mpz_srcptr rq = job->rq;


  mpz_t l, mtilde, m, adm1, t, k;
  mpz_poly f, g, f_raw, g_raw;
  int cmp;

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
#ifndef NDEBUG
  {
      mpz_t r,s;
      mpz_init(r);
      mpz_init(s);
      mpz_pow_ui(s, mtilde, header->d);
      mpz_sub(s, header->Ntilde, s);
      ASSERT_ALWAYS(mpz_tdiv_r_uint64(r, s, p1 * p1) == 0);
      ASSERT_ALWAYS(mpz_tdiv_r_uint64(r, s, p2 * p2) == 0);
      mpz_set_uint64(r, q);
      mpz_mul(r, r, r);
      mpz_tdiv_r(r, s, r);
      ASSERT_ALWAYS(mpz_cmp_ui(r, 0) == 0);
      mpz_clear(s);
      mpz_clear(r);
  }
#endif

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
  for (p = 2; p <= league->pt->Primes[league->pt->lenPrimes - 1]; p = getprime_mt(pi))
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
  mpz_set(mpz_poly_coeff(g, 1), l);
  mpz_neg(mpz_poly_coeff(g, 0), m);
  mpz_set(mpz_poly_coeff(f, header->d), header->ad);
  mpz_pow_ui(t, m, header->d);
  mpz_mul(t, t, header->ad);
  mpz_sub(t, header->N, t);
  mpz_set(mpz_poly_coeff(f, header->d - 1), adm1);
  check_divexact(t, t, "t", l, "l");
  mpz_pow_ui(mtilde, m, header->d - 1);
  mpz_mul(mtilde, mtilde, adm1);
  mpz_sub(t, t, mtilde);
  for (int j = header->d - 2; j > 0; j--)
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
      mpz_set(mpz_poly_coeff(f, j), adm1);
      /* subtract adm1*m^j */
      mpz_submul(t, mtilde, adm1);
    }
  check_divexact(t, t, "t", l, "l");
  mpz_set(mpz_poly_coeff(f, 0), t);

  /* As noticed by Min Yang, Qingshu Meng, Zhangyi Wang, Lina Wang and
     Huanguo Zhang in "Polynomial Selection for the Number Field Sieve in an
     Elementary Geometric View" (https://eprint.iacr.org/2013/583),
     if the coefficient of degree d-2 is of the same sign as the leading
     coefficient, the size optimization will not work well, thus we simply
     discard those polynomials. */
  if (mpz_sgn(mpz_poly_coeff_const(f, header->d)) * mpz_sgn(mpz_poly_coeff_const(f, header->d - 2)) > 0)
    {
      stats->discarded1++;
      goto end;
    }


  mpz_poly_cleandeg(f, header->d);
  ASSERT_ALWAYS(mpz_poly_degree(f) == (int) header->d);
  mpz_poly_cleandeg(g, 1);
  ASSERT_ALWAYS(mpz_poly_degree(g) == (int) 1);

  {
    /* information on all polynomials */
    stats->collisions++;
    stats->tot_found++;

    /* _raw lognorm */
    double skew = L2_skewness(f);
    double logmu = L2_lognorm(f, skew);

    polyselect_data_series_add(stats->raw_lognorm, logmu);
    polyselect_data_series_add(stats->raw_proj_alpha,
				 get_alpha_projective(f, get_alpha_bound()));
  }

  /* check that the algebraic polynomial has content 1, otherwise skip it */

  if (!mpz_poly_has_trivial_content(f)) 
      goto end;

  /* enter size optimization */

  {
      double st = seconds_thread();
      mpz_poly_set(g_raw, g);
      mpz_poly_set(f_raw, f);
      size_optimization(f, g, f_raw, g_raw, league->main->sopt_effort, league->main->verbose);
      stats->optimize_time += seconds_thread() - st;
      stats->opt_found++;
  }

  /* polynomials with f[d-1] * f[d-3] > 0 *after* size-optimization
     give worse exp_E values */
  if (mpz_sgn(mpz_poly_coeff_const(f, f->deg - 1)) * mpz_sgn(mpz_poly_coeff_const(f, f->deg - 3)) > 0) {
      stats->discarded2++;
      goto end;
  }

  /* register all stat to the stats object. This is a local
   * object, so no lock needed !
   */
  {
      stats->collisions_good++;

      double skew = L2_skewness(f);
      double logmu = L2_lognorm(f, skew);
      /* expected_rotation_gain() takes into account the
       * projective alpha */
      double exp_E = logmu + expected_rotation_gain(f, g);

      polyselect_priority_queue_push(stats->best_opt_logmu, logmu);
      polyselect_priority_queue_push(stats->best_exp_E, exp_E);
      polyselect_data_series_add(stats->opt_lognorm, logmu);
      polyselect_data_series_add(stats->exp_E, exp_E);
      polyselect_data_series_add(stats->opt_proj_alpha,
              get_alpha_projective(f, get_alpha_bound()));
  }

  /* print optimized (maybe size- or size-root- optimized)
   * polynomial */

  if (league->main->verbose >= 0) {
      {
          static pthread_mutex_t iolock = PTHREAD_MUTEX_INITIALIZER;
          pthread_mutex_lock(&iolock);
          polyselect_fprintf_poly_pair(stdout, header->N, f_raw, g_raw, 1);
          puts("#");
          polyselect_fprintf_poly_pair(stdout, header->N, f, g, 0);
          /* There's a significant carriage return to print. */
          puts("");
          pthread_mutex_unlock(&iolock);
      }
  }


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

#if 0
static void display_expected_memory_usage(polyselect_main_data_srcptr main, int nthreads)
{
    char buf[16];

    /* XXX document exactly how we reach this estimate. It seems quite bogus.
     */
    size_t exp_size =
            (BATCH_SIZE * 2 + INIT_FACTOR) * main->lenPrimes
                * (sizeof(uint32_t) + sizeof(uint64_t));
    exp_size *= nthreads;

    printf("# Info: estimated peak memory=%s (%d thread(s),"
            " batch %d inversions on SQ)\n",
            size_disp(exp_size, buf),
            nthreads, BATCH_SIZE);
}
#endif

static void declare_usage(param_list_ptr pl)
{
  param_list_decl_usage(pl, "degree",
			"(required, alias d) polynomial degree");
  param_list_decl_usage(pl, "n", "(required, alias N) input number");
  param_list_decl_usage(pl, "P",
			"(required) deg-1 coeff of g(x) has two prime factors in [P,2P]\n");

  param_list_decl_usage(pl, "admax", "maximal value for ad (+ 1)");
  param_list_decl_usage(pl, "admin", "minimal value for ad (default 0)");
  param_list_decl_usage(pl, "incr", "increment of ad (default 60)");
  param_list_decl_usage(pl, "maxtime",
			"stop the search after maxtime seconds");

  char str[200];
  snprintf(str, 200, "maximum number of special-q's considered\n"
	   "               for each ad (default %d)", DEFAULT_NQ);
  param_list_decl_usage(pl, "nq", str);
  snprintf(str, 200, "number of polynomials kept (default %d)", DEFAULT_POLYSELECT_KEEP);
  param_list_decl_usage(pl, "keep", str);
  snprintf(str, 200, "size-optimization effort (default %d)",
	   SOPT_DEFAULT_EFFORT);
  param_list_decl_usage(pl, "sopteffort", str);
  param_list_decl_usage(pl, "s", str);
  param_list_decl_usage(pl, "t", "number of threads to use (default 1)");
  param_list_decl_usage(pl, "F", "number of finer-grain threads to use (default 1)");
  param_list_decl_usage(pl, "v", "verbose mode");
  param_list_decl_usage(pl, "q", "quiet mode");
  param_list_decl_usage(pl, "target_E", "target E-value\n");
  param_list_decl_usage(pl, "chronogram", "store chronogram raw data to this file\n");
  verbose_decl_usage(pl);
}

static void usage(const char *argv, const char *missing, param_list_ptr pl)
{
  if (missing)
    {
      fprintf(stderr, "\nError: missing or invalid parameter \"-%s\"\n",
	      missing);
    }
  param_list_print_usage(pl, argv, stderr);
  param_list_clear(pl);
  exit(EXIT_FAILURE);
}

/* This thread loop does not (should not) depend on the way the matches
 * are acted upon. In a sense, we could have different match functions
 * use this same loop.
 */
void * thread_loop(polyselect_thread_ptr thread)
{
    polyselect_thread_team_ptr team = thread->team;
    polyselect_thread_league_ptr league = team->league;
    polyselect_main_data_srcptr main_data = league->main;

    polyselect_thread_bind(thread);

    if (thread->thread_index % main_data->finer_grain_threads == 0) {
        polyselect_thread_team_late_init(team);
    }
    
    polyselect_thread_late_init(thread);

    /* 
     * Asynchronous processing is (currently) single-threaded. Processing
     * an ad coefficient can be done collectively.
     *
     * All threads try to find an asynchronous job. When this fails, only
     * one thread in each team will compete to find a new ad coefficient
     * to process.
     */
    pthread_mutex_t * main_lock = thread->main_lock;
    ASSERT_ALWAYS(main_lock == &main_data->lock);
    pthread_mutex_t * team_lock = &team->lock;
    pthread_mutex_t * league_lock = &league->lock;


    /* a thread is either processing an async job, participating in the
     * processing of a sync job (or attempting to), or holds the team
     * lock
     *
     */
    unsigned int idx_max = polyselect_main_data_number_of_ad_tasks(main_data);

    pthread_mutex_lock(team_lock);

    polyselect_thread_team_i_am_ready(team, thread);

#ifdef DEBUG_POLYSELECT_THREADS
    fprintf(stderr, "thread %d (in team %d of size %d) wants to work\n",
            thread->thread_index,
            thread->team->team_index,
            thread->team->size);
#endif

    for( ; ; ) {
        pthread_mutex_lock(league_lock);
        /* is there an async job ready ? */
        struct dllist_head * ptr = dllist_get_first_node(&league->async_jobs);
        if (ptr) {
            dllist_pop(ptr);
            pthread_mutex_unlock(league_lock);
            polyselect_match_info_ptr job = dllist_entry(ptr, struct polyselect_match_info_s, queue);

            polyselect_thread_team_enter_async(team, thread);
            pthread_mutex_unlock(team_lock);
            /********* BEGIN UNLOCKED SECTION **************/
            polyselect_thread_chronogram_chat(thread, "enter match");
            polyselect_process_match_async(league, thread->stats, job);
            dllist_push_back(&thread->empty_job_slots, &job->queue);
            polyselect_thread_chronogram_chat(thread, "leave match");
            /********** END UNLOCKED SECTION ***************/
            pthread_mutex_lock(team_lock);
            polyselect_thread_team_leave_async(team, thread);
        } else {
            pthread_mutex_unlock(league_lock);
            /* Then we want to contribute to synchronous work. */
            /* note that we have the team lock, at this point.
             * Furthermore, this call is a barrier.
             */
            polyselect_thread_team_enter_sync_zone(team, thread);

#ifdef DEBUG_POLYSELECT_THREADS
            fprintf(stderr, "thread %d is %d-th sync thread in team %d\n",
                    thread->thread_index,
                    thread->index_in_sync_zone,
                    thread->team->team_index);
#endif
            /*
             * This is important because when async jobs come back, we
             * want them to notice that the main crowd has called it a
             * day.
             */

            if (team->done)
                break;

            /* important change. schedule work only when we're __last__
             * to enter this state */
            if (thread->team->count->ready == 0 && thread->team->count->sync2 == 0) {
                pthread_mutex_lock(main_lock);
                unsigned long i = team->main_nonconst->idx;
                team->main_nonconst->idx += (i < idx_max);
                pthread_mutex_unlock(main_lock);

                if (i == idx_max) {
#ifdef DEBUG_POLYSELECT_THREADS
                    fprintf(stderr, "thread %d wants to call it off (sync=%d sync2=%d ready=%d async=%d\n",
                            thread->thread_index,
                            thread->team->count->sync,
                            thread->team->count->sync2,
                            thread->team->count->ready,
                            thread->team->count->async);
#endif
                    /*
                    for( ; team->count->async || !dllist_is_empty(&league->async_jobs) ; ) {
                        pthread_cond_wait(&team->count->w_async_empty, &team->lock);
                    }
                    */

                    team->done = 1;

                    polyselect_thread_team_post_work_stop(team, thread);
                    polyselect_thread_team_leave_sync_zone(team, thread);
                    break;
                }


                /* we're effectively the only one in the team, here, so
                 * we can safely touch these unlocked.
                 */
                polyselect_thread_team_set_idx(team, i);

                polyselect_thread_chronogram_chat(thread, "enter ad, %.0f", mpz_get_d(team->ad));


                unsigned long c = 0;

                thread->stats->number_of_ad_values++;

                /* These calls will temporarily release the team lock */
                c = collision_on_p_conductor(thread);

                if (main_data->nq > 0) {
                    collision_on_sq_conductor(c, thread);
                }

                polyselect_thread_team_post_work_stop(team, thread);

                polyselect_thread_chronogram_chat(thread, "leave ad, %.0f", mpz_get_d(thread->team->ad));

                {
                    /* Not absolutely certain that this printing makes
                     * sense. First, it's unlocked, which sounds quite
                     * dangerous. Also, the semantics of wct0 look quite
                     * awkward. And finally, if we're really that
                     * interested, why not go for the chronogram data
                     * instead?
                     */
                    printf("# thread %u completed ad=%.0f at time=%.2fs ; ad: %.2fs\n",
                            thread->thread_index,
                            mpz_get_d(team->ad),
                            wct_seconds() - main_data->stats->wct0,
                            wct_seconds() - thread->stats->wct0);
                    fflush(stdout);
                }
            } else {
                for( ; ; ) {
                    /* wait for sync tasks to be posted. */
                    pthread_cond_wait(&team->count->w_job, &team->lock);

                    /*
                    fprintf(stderr, "thread %d (%d-th sync thread in team %d of size %d among %d sync threads) wakes up\n",
                            thread->thread_index,
                            thread->index_in_sync_zone,
                            thread->team->team_index,
                            thread->team->task->expected,
                            thread->team->count->ready);
                            */

                    if (team->task->expected == 0) {
                        /*
                        fprintf(stderr, "thread %d (%d-th sync thread in team %d) leaves sync group\n",
                                thread->thread_index,
                                thread->index_in_sync_zone,
                                thread->team->team_index);
                                */
                        break;
                    }

                    if (team->task->expected <= thread->index_in_sync_zone) {
                        /* there is a structural race condition here. We
                         * arrived in the sync section after the leader
                         * thread posted this task. Not much to be said,
                         * we missed the train... All we have to do is
                         * wait for the next one.
                         */
                        continue;
                    }

                    /* This call has the team lock held, but it is
                     * expected that the lock be released during the
                     * call.
                     */
                    polyselect_thread_team_enter_sync_task(team, thread);
                    (*team->task->f)(thread);
                    polyselect_thread_team_leave_sync_task(team, thread);
                }
            }

            pthread_mutex_lock(league_lock);
            /*
            fprintf(stderr, "thread %d moves %zu jobs to async queue (current size: %zu)\n", thread->thread_index,
                    dllist_length(&thread->async_jobs),
                    dllist_length(&league->async_jobs));
                    */
            dllist_bulk_move_back(&league->async_jobs, &thread->async_jobs);
            pthread_mutex_unlock(league_lock);
            polyselect_thread_team_leave_sync_zone(team, thread);

            /*
            fprintf(stderr, "thread %d (%d-th sync thread in team %d) moves on\n",
                    thread->thread_index,
                    thread->index_in_sync_zone,
                    thread->team->team_index);
                    */
        }
        pthread_mutex_lock(main_lock);
        polyselect_main_data_commit_stats_unlocked(team->main_nonconst, thread->stats, team->header->ad);
        pthread_mutex_unlock(main_lock);
    }
    {
        /* XXX TODO acquire a lock, probably main_lock */
        printf("# thread %u exits at time=%.2fs\n",
                thread->thread_index,
                wct_seconds() - main_data->stats->wct0);
        fflush(stdout);
    }
    pthread_mutex_unlock(team_lock);

    return NULL;
}


int main(int argc, char const * argv[])
{
  char const ** argv0 = argv;
  /* nthreads = 0 means: do something automatic */
  int quiet = 0;
  const char * chronogram_file = NULL;



  polyselect_main_data main_data;
  polyselect_main_data_init_defaults(main_data);

  /* read params */
  param_list pl;
  param_list_init(pl);

  declare_usage(pl);

  param_list_configure_switch(pl, "-v", &main_data->verbose);
  param_list_configure_switch(pl, "-q", &quiet);
  param_list_configure_alias(pl, "degree", "-d");
  param_list_configure_alias(pl, "incr", "-i");
  param_list_configure_alias(pl, "n", "-N");

  if (argc == 1)
    usage(argv0[0], NULL, pl);

  argv++, argc--;
  for (; argc;)
    {
      if (param_list_update_cmdline(pl, &argc, &argv))
	continue;
      fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
      usage(argv0[0], NULL, pl);
    }


  polyselect_main_data_parse_Nd(main_data, pl);
  polyselect_main_data_parse_ad_range(main_data, pl);
  polyselect_main_data_parse_P(main_data, pl);
  param_list_parse_ulong(pl, "nq", &main_data->nq);
  chronogram_file = param_list_lookup_string(pl, "chronogram");

  const char * tmp;
  if ((tmp = param_list_lookup_string(pl, "t")) && strcmp(tmp, "auto") == 0) {
      main_data->nthreads = 0;
  } else {
      param_list_parse_uint(pl, "t", &main_data->nthreads);
  }
  param_list_parse_uint(pl, "F", &main_data->finer_grain_threads);
#ifndef HAVE_HWLOC
  if (main_data->nthreads == 0) {
      fprintf(stderr, "Warning: -t auto requires hwloc\n");
      main_data->nthreads = 1;
  }
#endif

  /* size optimization effort that passed to size_optimization */
  param_list_parse_uint(pl, "sopteffort", &main_data->sopt_effort);

  {
      param_list_parse_int(pl, "keep", &main_data->keep);
      polyselect_stats_update_keep(main_data->stats, main_data->keep);
  }

  polyselect_main_data_parse_maxtime_or_target(main_data, pl);

  if (param_list_warn_unused(pl))
    usage(argv0[0], NULL, pl);

  /* print command line */
  verbose_interpret_parameters(pl);

  param_list_print_command_line(stdout, pl);

  /* quiet mode */
  if (quiet == 1)
    main_data->verbose = -1;

  // display_expected_memory_usage(main_data, nthreads);


  polyselect_thread_chronogram_init(chronogram_file);

  /* Try to see if the number of threads that we've been passed makes any
   * sort of sense */
  polyselect_main_data_check_topology(main_data);


  /* Start one league per NUMA node */
  polyselect_main_data_prepare_leagues(main_data);
  polyselect_main_data_prepare_teams(main_data);
  polyselect_main_data_prepare_threads(main_data);

  polyselect_thread_chronogram_init(chronogram_file);

  polyselect_main_data_go_parallel(main_data, thread_loop);

  polyselect_thread_chronogram_clear();
  polyselect_main_data_dispose_threads(main_data);
  polyselect_main_data_dispose_teams(main_data);
  polyselect_main_data_dispose_leagues(main_data);

#ifndef HAVE_RUSAGE_THREAD	/* optimize_time is correct only if
                                   RUSAGE_THREAD works or in mono-thread
                                   mode */
  if (main_data->nthreads != 1)
      main_data->stats->optimize_time = -1;
#endif


  if (main_data->verbose >= 0) {
      /* This is very weird. What's the point? */
      main_data->stats->potential_collisions *= polyselect_main_data_expected_collisions(main_data);
      polyselect_stats_display_final(main_data->stats, main_data->verbose);
  }


  polyselect_main_data_clear(main_data);

  param_list_clear(pl);

  return 0;
}
