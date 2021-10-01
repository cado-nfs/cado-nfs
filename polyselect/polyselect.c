/*
  Polynomial selection using Kleinjung's algorithm (cf slides presented
  at the CADO Workshop in October 2008, Nancy, France).

  [1. Run and parameters]

  The parameters are similar to those in polyselect2.c, except the following,

  "-nq xxx" denotes the number of special-q's trials for each ad;

  Please report bugs to the Bug Tracking System on:
  https://gitlab.inria.fr/cado-nfs/cado-nfs/-/issues
*/

#define EMIT_ADDRESSABLE_shash_add

#include "cado.h"		// IWYU pragma: keep
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
#include <stdbool.h>		// bool
#include <stdio.h>
#include <stdlib.h>		// malloc ...
#include <stdint.h>		// uint64_t
#include <gmp.h>
#include <sys/time.h>
#include "macros.h"		// ASSERT
#include "misc.h"
#include "omp_proxy.h"
#include "params.h"
#include "polyselect_main_queue.h"
#include "polyselect_collisions.h"
#include "polyselect_shash.h"
#include "polyselect_match.h"
#include "polyselect_norms.h"
#include "polyselect_alpha.h"
#include "portability.h"
#include "roots_mod.h"
#include "size_optimization.h"
#include "timing.h"		// for seconds
#include "verbose.h"		// verbose_output_print
#include "getprime.h"
#include "dllist.h"

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


FILE * chronogram = NULL;
pthread_mutex_t chronogram_lock;

static inline void chat_chronogram(const char * fmt, ...)
{
    if (!chronogram) return;

    va_list ap;
    va_start(ap, fmt);
    char * msg;
    int rc = vasprintf(&msg, fmt, ap);
    ASSERT_ALWAYS(rc >= 0);
    struct timeval tv[1];
    gettimeofday(tv, NULL);
    pthread_mutex_lock(&chronogram_lock);
    fprintf(chronogram, "%lu.%06lu %d %s\n", tv->tv_sec, tv->tv_usec, omp_get_thread_num(), msg);
    pthread_mutex_unlock(&chronogram_lock);
    free(msg);
    va_end(ap);
}


void
polyselect_process_match_async(polyselect_main_data_srcptr main, polyselect_stats_ptr stats, polyselect_match_info_ptr job)
{
  polyselect_poly_header_srcptr header = job->header;
  unsigned long p1 = job->p1;
  unsigned long p2 = job->p2;
  const int64_t i = job->i;
  uint64_t q = job->q;
  mpz_srcptr rq = job->rq;


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
  for (p = 2; p <= main->Primes[main->lenPrimes - 1]; p = getprime_mt(pi))
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
      stats->discarded1++;
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
    stats->collisions++;
    stats->tot_found++;
    polyselect_data_series_add(stats->raw_lognorm, logmu);
    polyselect_data_series_add(stats->raw_proj_alpha,
				 get_alpha_projective(f, get_alpha_bound()));
  }

  /* if the polynomial has small norm, we optimize it */
  did_optimize = optimize_raw_poly(f, g, main, stats);

  /* print optimized (maybe size- or size-root- optimized) polynomial */
  if (did_optimize && main->verbose >= 0) {
#ifdef HAVE_OPENMP
#pragma omp critical
#endif
      {
          polyselect_fprintf_poly_pair(stdout, header->N, f_raw, g_raw, 1);
          puts("#");
          polyselect_fprintf_poly_pair(stdout, header->N, f, g, 0);
          /* There's a significant carriage return to print. */
          puts("");
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


int main(int argc, char *argv[])
{
  char **argv0 = argv;
  int quiet = 0, nthreads = 1;
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

  param_list_parse_int(pl, "t", &nthreads);
#ifdef HAVE_OPENMP
  omp_set_num_threads(nthreads);
#else
  if (nthreads > 1)
    {
      fprintf(stderr,
	      "Warning, -t %d ignored because openmp support is missing\n",
	      nthreads);
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


  /* initialize primes in [P,2*P] */
  polyselect_main_data_prepare_primes(main_data);

  display_expected_memory_usage(main_data, nthreads);

  unsigned long idx_max = polyselect_main_data_number_of_ad_tasks(main_data);


  if (idx_max < (unsigned long) nthreads)
    {
      fprintf(stderr,
	      "# Warning: the current admin, admax, incr settings only make it possible to run %lu jobs in parallel, so that we won't be able to do %d-thread parallelism as requested\n",
	      idx_max, nthreads);
    }


  unsigned long idx = 0;

  pthread_mutex_init(&chronogram_lock, NULL);
  if (chronogram_file) {
      chronogram = fopen(chronogram_file, "w");
      setbuf(chronogram, NULL);
  }

#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
  {
      /* The loop starts with the mutex locked */
      pthread_mutex_lock(&main_data->lock);
      for(;;) {
          if (!dllist_is_empty(&main_data->async_jobs)) {
              struct dllist_head * ptr = main_data->async_jobs.next;
              polyselect_match_info_ptr job = dllist_entry(ptr, struct polyselect_match_info_s, queue);
              dllist_pop(ptr);
              pthread_mutex_unlock(&main_data->lock);
              /********* BEGIN UNLOCKED SECTION **************/
              chat_chronogram("enter match");
              polyselect_stats stats;
              polyselect_stats_init(stats, main_data->keep);
              polyselect_process_match_async(main_data, stats, job);
              polyselect_match_info_clear(job);
              free(job);
              chat_chronogram("leave match");
              /********** END UNLOCKED SECTION ***************/
              pthread_mutex_lock(&main_data->lock);
              polyselect_main_data_commit_stats_unlocked(main_data, stats, NULL);
              polyselect_stats_clear(stats);
          } else if (idx < idx_max) {
              unsigned long i = idx++;
              pthread_mutex_unlock(&main_data->lock);
              /********* BEGIN UNLOCKED SECTION **************/
              chat_chronogram("enter ad");

              polyselect_thread_locals loc;
              polyselect_thread_locals_init(loc, main_data, i);

              unsigned long c = 0;

              loc->stats->number_of_ad_values++;

              polyselect_shash_t H;
              polyselect_shash_init(H, 4 * loc->main->lenPrimes);
              c = collision_on_p(H, NULL, loc);
              if (loc->main->nq > 0)
                  collision_on_sq(c, H, NULL, loc);
              polyselect_shash_clear(H);

              {
                  printf("# thread %d completed ad=%.0f at time=%.2fs ; ad: %.2fs\n",
                          omp_get_thread_num(),
                          mpz_get_d(loc->ad),
                          wct_seconds() - main_data->stats->wct0,
                          wct_seconds() - loc->stats->wct0);
                  fflush(stdout);
              }

              chat_chronogram("leave ad, %zu", dllist_length(&loc->async_jobs));
              /********** END UNLOCKED SECTION ***************/
              pthread_mutex_lock(&main_data->lock);

              polyselect_main_data_commit_stats_unlocked(main_data, loc->stats, loc->ad);

              dllist_bulk_move_back(&main_data->async_jobs, &loc->async_jobs);
          } else {
              /* we're done! */
              break;
          }
      }

      {
          printf("# thread %d exits at time=%.2fs\n",
                  omp_get_thread_num(),
                  wct_seconds() - main_data->stats->wct0);
          fflush(stdout);
      }

      pthread_mutex_unlock(&main_data->lock);
  }

  if (chronogram_file)
      fclose(chronogram);

  pthread_mutex_destroy(&chronogram_lock);

#ifndef HAVE_RUSAGE_THREAD	/* optimize_time is correct only if
                                   RUSAGE_THREAD works or in mono-thread
                                   mode */
  if (nthreads == 1)
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
