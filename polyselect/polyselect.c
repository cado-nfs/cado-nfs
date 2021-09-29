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
#include "macros.h"		// ASSERT
#include "misc.h"
#include "omp_proxy.h"
#include "params.h"
#include "polyselect_main_queue.h"
#include "polyselect_collisions.h"
#include "polyselect_collisions_gmp.h"
#include "polyselect_shash.h"
#include "portability.h"
#include "roots_mod.h"
#include "size_optimization.h"
#include "timing.h"		// for seconds
#include "verbose.h"		// verbose_output_print

static void newAlgo(polyselect_thread_locals_ptr loc)
{
  unsigned long c = 0;

  loc->stats->number_of_ad_values++;

  if (sizeof(unsigned long int) == 8)
    {
      polyselect_shash_t H;
      polyselect_shash_init(H, 4 * loc->main->lenPrimes);
      c = collision_on_p(H, loc);
      if (loc->main->nq > 0)
	collision_on_sq(c, H, loc);
      polyselect_shash_clear(H);
  } else
    {
      /* This code is used on 32-bit machines. Do we _really_ have to go
       * through this trouble? Why not use uint64_t's all over the place?
       *
       * (OTOH, it's perhaps good practice to have some gmp code around).
       */
      c = gmp_collision_on_p(loc);
      if (loc->main->nq > 0)
	gmp_collision_on_sq(c, loc);
    }

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
      int keep = DEFAULT_POLYSELECT_KEEP;
      param_list_parse_int(pl, "keep", &keep);
      polyselect_stats_setup_keep_best(main_data->stats, keep);
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


#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
  for (unsigned long idx = 0; idx < idx_max; idx++)
  {
      polyselect_thread_locals loc;
      polyselect_thread_locals_init(loc, main_data, idx);

      newAlgo(loc);

      if (main_data->verbose > 0)
#ifdef HAVE_OPENMP
#pragma omp critical
#endif
      {
          printf("# thread %d completed ad=%.0f at time=%.2fs\n",
                  omp_get_thread_num(),
                  mpz_get_d(loc->ad),
                  wct_seconds() - loc->stats->wct0);
          fflush(stdout);
      }

      polyselect_main_data_commit_stats(main_data, loc->stats);

      polyselect_thread_locals_clear(loc);
  }

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
