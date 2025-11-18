#include "cado.h" // IWYU pragma: keep

/*
 * TODO: describe the algorithm being used.
 *
 *
 * 20210930: this code has been edited to refactor some of the machinery
 * form polyselect_collisions.c ; it seems that for the most part, the
 * code here uses the same general framework, only the work that is done
 * upon a match differs.
 */

#define EMIT_ADDRESSABLE_shash_add
#define NEW_ROOTSIEVE
#define INIT_FACTOR 8UL
//#define DEBUG_POLYSELECT

#define BATCH_SIZE 20

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <pthread.h>
#include <float.h>      // DBL_MAX
#include <math.h> // pow
#include <gmp.h>
#include "mpz_poly.h"
#include "area.h"
#include "omp_proxy.h"
#include "polyselect_arith.h"
#include "polyselect_collisions.h"
#include "polyselect_norms.h"
#include "polyselect_alpha.h"
#include "polyselect_shash.h"
#include "polyselect_poly_header.h"
#include "modredc_ul.h"
#include "mpz_vector.h"
#include "roots_mod.h"  // roots_mod_uint64
#include "timing.h"     // milliseconds
#include "verbose.h"    // verbose_decl_usage
#include "auxiliary.h"  // DEFAULT_INCR


#ifdef NEW_ROOTSIEVE
#include "ropt.h"
#include "macros.h"
#include "params.h"
#endif

double best_E = DBL_MAX; /* combined score E (the smaller the better) */
double aver_E = 0.0;
unsigned long found = 0; /* number of polynomials found so far */

mpz_t maxS; /* maximun skewness. O for default max */

/* check that l/2 <= d*m0/P^2, where l = p1 * p2 * q with P <= p1, p2 <= 2P
   q is the product of special-q primes. It suffices to check that
   q <= d*m0/(2P^4).

   Is it really different here from what we have in polyselect ?

   XXX Presently, we're using polyselect_main_data_check_parameters
   instead (it is called from find_suitable_lq). We're not sure if this
   is correct.

*/

#if 0
int
twocubics_check_parameters (polyselect_main_data_ptr main, mpz_t m0, unsigned long d, unsigned long lq)
{
  double maxq = 1.0, maxP;
  int k = lq;

  while (k > 0)
    maxq *= (double) SPECIAL_Q[LEN_SPECIAL_Q - 1 - (k--)];

  maxP = (double) main->Primes[main->lenPrimes - 1];
  if (2.0 * pow (maxP, 4.0) * maxq >= (double) d * mpz_get_d (m0))
    return 0;

  if (maxq > pow (maxP, 2.0))
    return 0;

  return 1;
}
#endif

/* XXX as a workaround, we do this check on top of twocubics_match. Also,
 * it seems that the code above is overly restrictive, it mixes up the
 * min P and the max P. The version below tries to get this right.
 */
int
twocubics_check_parameters_withq (polyselect_main_data_srcptr main, mpz_t m0, unsigned long d, unsigned long q)
{
  if (2.0 * pow (main->P, 4.0) * q >= (double) d * mpz_get_d (m0))
    return 0;

  if (q > pow (main->P, 2.0))
    return 0;

  return 1;
}


/* Compute maximun skewness, which in floor(N^(1/d^2)) */
void
compute_default_max_skew (mpz_t skew, mpz_t N, int d)
{
  mpz_root (skew, N, (unsigned long) d*d);
}


/* what do we do if we encounter a match in the sense of Kleinjung's 2008
 * algorithm. I guess that this is where the algorithm in this file
 * differs from polyselect.c
 */

/* rq is a root of N = (m0 + rq)^d mod (q^2) */
void
twocubics_match (unsigned long p1, unsigned long p2, const int64_t i,
        uint64_t q, mpz_srcptr rq,
        polyselect_thread_locals_ptr loc)
{
#if 0
  unsigned long j;
  int cmp;
  double skew, logmu, E;
  mpz_poly F;
#endif

  polyselect_poly_header_srcptr header = loc->header;

  if (!twocubics_check_parameters_withq (loc->main, loc->header->m0, loc->header->d, q))
      return;


  mpz_t l, r, k, mprime, Nprime, C, l2, tmp, r1, r0, t, adm1, m, skew, root;
  mpz_vector_t a, b, reduced_a, reduced_b;
  mpz_poly f, g;

  mpz_init (root);
  mpz_init (m);
  mpz_init (adm1);
  mpz_init (l);
  mpz_init (l2);
  mpz_init (r);
  mpz_init (k);
  mpz_init (mprime);
  mpz_init (Nprime);
  mpz_init (C);
  mpz_init (tmp);
  mpz_init (r1);
  mpz_init (t);
  mpz_init (r0);
  mpz_init (skew);

  mpz_vector_init (a, header->d+1);
  mpz_vector_init (b, header->d+1);
  mpz_vector_init (reduced_a, header->d+1);
  mpz_vector_init (reduced_b, header->d+1);

  mpz_poly_init (f, header->d);
  mpz_poly_init (g, header->d);

#ifdef DEBUG_POLYSELECT
  gmp_printf ("#### MATCH ######\nN = %Zd\nd = %d\nad = %Zd\n"
              "p1 = %lu\np2 = %lu\nq = %lu\ni = %" PRId64 "\n"
              "rq = %Zd\n", header->N, header->d, header->ad, p1, p2, q, i, rq);
#endif

  /* l = p1*p2*q */
  mpz_set_ui (l, p1);
  mpz_mul_ui (l, l, p2);
  mpz_mul_ui (l, l, q);
#ifdef DEBUG_POLYSELECT
  gmp_printf ("l = p1 * p2 * q\nl == %Zd # has %lu bits\n", l,
              mpz_sizeinbase(l, 2));
#endif
  mpz_mul (l2, l, l); /* l2 = l^2 */

  /* r = rq + i*q^2 */
  mpz_set_si (r, i);
  if (rq) {
      mpz_mul_ui (r, r, q);
      mpz_mul_ui (r, r, q);
      mpz_add (r, r, rq);
  }
#ifdef DEBUG_POLYSELECT
  gmp_printf ("r = rq + i * q^2 \nr == %Zd # has %lu bits\n", r,
              mpz_sizeinbase(r, 2));
#endif

  /* k = header->d^header->d * header->ad^(header->d-1) */
  mpz_mul_ui (k, header->ad, header->d);
  mpz_pow_ui (k, k, header->d-1);
  mpz_mul_ui (k, k, header->d);
#ifdef DEBUG_POLYSELECT
  gmp_printf ("k = d^d * ad^(d-1)\nk == %Zd\n", k);
#endif

  /* Nprime = k * header->N */
  mpz_mul (Nprime, k, header->N);
#ifdef DEBUG_POLYSELECT
  gmp_printf ("Nprime = k * N\nNprime == %Zd\n", Nprime);
#endif

  /* mprime = header->m0 + r */
  mpz_add (mprime, header->m0, r);
#ifdef DEBUG_POLYSELECT
  gmp_printf ("m0 = %Zd\nmprime = m0 + r\nmprime == %Zd\n", header->m0, mprime);
#endif

  /* C = mprime^header->d - Nprime */
  mpz_pow_ui (C, mprime, header->d);
  mpz_sub (C, C, Nprime);
  ASSERT_ALWAYS (mpz_divisible_p (C, l2));
#ifdef DEBUG_POLYSELECT
  gmp_printf ("(mprime^d - Nprime) %% l^2 == 0\n");
#endif

  /* adm1 is such that mprime = header->d*header->ad*m + adm1*l and -header->d*header->ad/2 <= adm1 < header->d*header->ad/2
     We have adm1 = mprime/l mod (header->d*header->ad). */
  mpz_mul_ui (tmp, header->ad, header->d); /* tmp = header->d*header->ad */
  if (mpz_invert (adm1, l, tmp) == 0)
  {
    fprintf (stderr, "Error in 1/l mod (d*ad)\n");
    abort();
  }
  mpz_mul (adm1, adm1, mprime);
  mpz_mod (adm1, adm1, tmp);
  mpz_mul_2exp (t, adm1, 1);
  if (mpz_cmp (t, tmp) >= 0)
    mpz_sub (adm1, adm1, tmp);
#ifdef DEBUG_POLYSELECT
  gmp_printf ("adm1 = %Zd\n", adm1);
#endif

  /* m = (mprime - adm1 * l)/ (header->d * header->ad) */
  mpz_mul (m, adm1, l);
  mpz_sub (m, mprime, m);
  ASSERT_ALWAYS (mpz_divisible_ui_p (m, header->d));
  mpz_divexact_ui (m, m, header->d);
  ASSERT_ALWAYS (mpz_divisible_p (m, header->ad));
  mpz_divexact (m, m, header->ad);
#ifdef DEBUG_POLYSELECT
  gmp_printf ("m = (mprime - adm1*l) / (d*ad)\nm == %Zd\n", m);
#endif

  /* Set vector a = (-m, l, 0, ..., 0) */
  mpz_neg (tmp, m);
  mpz_vector_setcoordinate (a, 0, tmp); /* a[0] = -m */
  mpz_vector_setcoordinate (a, 1, l); /* a[1] = -l */
  for (unsigned int j = 2; j <= header->d; j++)
    mpz_vector_setcoordinate_ui (a, j, 0); /* a[j] = 0 */

  /* Set vector b = (a0, a1, ..., ai, ..., adm1, header->ad) */
  mpz_vector_setcoordinate (b, header->d, header->ad);  /* b[header->d] = header->ad */
  mpz_vector_setcoordinate (b, header->d-1, adm1); /* b[header->d-1] = adm1 */




  mpz_pow_ui (t, m, header->d);
  mpz_mul (t, t, header->ad);
  mpz_sub (t, header->N, t);
  ASSERT_ALWAYS (mpz_divisible_p (t, l));

  mpz_divexact (t, t, l);
  mpz_pow_ui (tmp, m, header->d-1);
  mpz_mul (tmp, tmp, adm1);
  mpz_sub (t, t, tmp);
  for (int j = header->d - 2; j > 0; j--)
  {
    ASSERT_ALWAYS (mpz_divisible_p (t, l));
    mpz_divexact (t, t, l);
    /* t = a_j*m^j + l*R thus a_j = t/m^j mod l */
    mpz_pow_ui (tmp, m, j);
    /* fdiv rounds toward -infinity: r1 = floor(t/tmp) */
    mpz_fdiv_q (r1, t, tmp); /* t -> r1 * tmp + t */
    mpz_invert (k, tmp, l); /* search r1 + k such that */

    mpz_mul (k, k, t);
    mpz_sub (k, k, r1);
    mpz_mod (k, k, l);

    mpz_mul_2exp (k, k, 1);
    int cmp = mpz_cmp (k, l);
    mpz_div_2exp (k, k, 1);
    if (cmp >= 0)
      mpz_sub (k, k, l);
    mpz_add (r1, r1, k);
    mpz_vector_setcoordinate (b, j, r1);
    /* subtract r1*m^j */
    mpz_submul (t, tmp, r1);
  }
  ASSERT_ALWAYS (mpz_divisible_p (t, l));
  mpz_divexact (t, t, l);
  mpz_vector_setcoordinate (b, 0, t);


  mpz_vector_get_mpz_poly(f, a);
  mpz_vector_get_mpz_poly(g, b);
#ifdef DEBUG_POLYSELECT
  printf ("a = ");
  mpz_poly_fprintf (stdout, f);
  printf ("b = ");
  mpz_poly_fprintf (stdout, g);
#endif

  mpz_vector_reduce_with_max_skew (reduced_a, reduced_b, skew, a, b, maxS, header->d);

  mpz_vector_get_mpz_poly(f, reduced_a);
  mpz_vector_get_mpz_poly(g, reduced_b);

#ifdef DEBUG_POLYSELECT
  gmp_printf ("skew = %Zd\nf = ", skew);
  mpz_poly_fprintf (stdout, f);
  printf ("g = ");
  mpz_poly_fprintf (stdout, g);
#endif

  mpz_invert (root, l, header->N);
  mpz_mul (root, root, m);
  mpz_mod (root, root, header->N);
#ifdef DEBUG_POLYSELECT
  gmp_printf ("root = (m / l) %% N\nroot == %Zd\n", root);

  gmp_printf ("## Begin poly file for ad = %Zd and l = %Zd\n", header->ad, l);
#endif

  double skewness, logmu[2], alpha[2], E;

  skewness = L2_combined_skewness2 (f, g);
  logmu[0] = L2_lognorm (g, skewness);
  logmu[1] = L2_lognorm (f, skewness);
  alpha[0] = get_alpha (g, ALPHA_BOUND_SMALL);
  alpha[1] = get_alpha (f, ALPHA_BOUND_SMALL);
  E = logmu[1] + alpha[1] + logmu[0] + alpha[0];

    {
      static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
      pthread_mutex_lock (&lock);
  found ++;
  aver_E += E;

  if (E < best_E)
    {
      best_E = E;
      gmp_printf("n: %Zd\n", header->N);
      for (int i = 0; i <= f->deg; i++)
        gmp_printf ("c%d: %Zd\n", i, mpz_poly_coeff_const(f, i));
      for (int i = 0; i <= g->deg; i++)
        gmp_printf ("Y%d: %Zd\n", i, mpz_poly_coeff_const(g, i));
      printf ("skew: %1.2f\n", skewness);
      printf ("# f lognorm %1.2f, alpha %1.2f, score %1.2f\n",
              logmu[1], alpha[1], logmu[1] + alpha[1]);
      printf ("# g lognorm %1.2f, alpha %1.2f, score %1.2f\n",
              logmu[0], alpha[0], logmu[0] + alpha[0]);
      printf ("# f+g score %1.2f\n", E);
      printf ("# found %lu polynomial(s) so far, aver. E = %1.2f\n",
              found, aver_E / (double) found);
#ifdef DEBUG_POLYSELECT
      printf ("## End poly file\n");
#else
      printf ("\n");
#endif
    }
    pthread_mutex_unlock (&lock);
  }
  
  mpz_clear (root);
  mpz_clear (skew);
  mpz_clear (adm1);
  mpz_clear (l);
  mpz_clear (l2);
  mpz_clear (C);
  mpz_clear (r);
  mpz_clear (k);
  mpz_clear (mprime);
  mpz_clear (m);
  mpz_clear (Nprime);
  mpz_clear (tmp);
  mpz_clear (r1);
  mpz_clear (r0);
  mpz_clear (t);

  mpz_poly_clear (f);
  mpz_poly_clear (g);

  mpz_vector_clear (a);
  mpz_vector_clear (b);
  mpz_vector_clear (reduced_a);
  mpz_vector_clear (reduced_b);
}

static void newAlgo(polyselect_thread_locals_ptr loc)
{
  unsigned long c = 0;

  loc->stats->number_of_ad_values++;

  polyselect_shash_t H;
  polyselect_shash_init(H, 4 * loc->main->lenPrimes);
  c = collision_on_p(H, twocubics_match, loc);
  if (loc->main->nq > 0)
      collision_on_sq(c, H, twocubics_match, loc);
  polyselect_shash_clear(H);
}


static void
declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "degree", "polynomial degree (2 or 3, "
                                      "default is 3)");
  param_list_decl_usage(pl, "n", "(required, alias N) input number");
  param_list_decl_usage(pl, "P", "(required) deg-1 coeff of g(x) has two prime factors in [P,2P]\n");

  param_list_decl_usage(pl, "admax", "max value for ad (+1)");
  param_list_decl_usage(pl, "admin", "min value for ad (default 0)");
  param_list_decl_usage(pl, "incr", "increment of ad (default 60)");
  param_list_decl_usage(pl, "skewness", "maximun skewness possible "
                                        "(default N^(1/9))");
  param_list_decl_usage(pl, "maxtime", "stop the search after maxtime seconds");

  char str[200];
  snprintf (str, 200, "maximum number of special-q's considered\n"
            "               for each ad (default %d)", INT_MAX);
  param_list_decl_usage(pl, "nq", str);
  param_list_decl_usage(pl, "keep", "number of polynomials kept (default 10)");
  param_list_decl_usage(pl, "t", "number of threads to use (default 1)");
  param_list_decl_usage(pl, "v", "verbose mode");
  param_list_decl_usage(pl, "q", "quiet mode");
  snprintf (str, 200, "sieving area (default %.2e)", AREA);
  param_list_decl_usage(pl, "area", str);
  snprintf (str, 200, "algebraic smoothness bound (default %.2e)", BOUND_F);
  param_list_decl_usage(pl, "Bf", str);
  snprintf (str, 200, "rational smoothness bound (default %.2e)", BOUND_G);
  param_list_decl_usage(pl, "Bg", str);
  verbose_decl_usage(pl);
}

static void
usage (const char *argv, const char * missing, param_list pl)
{
  if (missing) {
    fprintf(stderr, "\nError: missing or invalid parameter \"-%s\"\n",
        missing);
  }
  param_list_print_usage(pl, argv, stderr);
  exit (EXIT_FAILURE);
}

int main(int argc, char const * argv[])
{
  char const ** argv0 = argv;
  int quiet = 0, nthreads = 1;

  polyselect_main_data main_data;

  polyselect_main_data_init_defaults(main_data);

  mpz_init (maxS);

  /* read params */
  param_list pl;
  param_list_init (pl);

  declare_usage(pl);

  param_list_configure_switch (pl, "-v", &main_data->verbose);
  param_list_configure_switch (pl, "-q", &quiet);
  param_list_configure_alias(pl, "degree", "-d");
  param_list_configure_alias(pl, "incr", "-i");
  param_list_configure_alias(pl, "n", "-N");

  if (argc == 1)
    usage (argv0[0], NULL, pl);

  /* Add this as a default */
  param_list_add_key(pl, "degree", "3", PARAMETER_FROM_FILE);

  argv++, argc--;
  for ( ; argc; ) {
    if (param_list_update_cmdline (pl, &argc, &argv)) continue;
    fprintf (stderr, "Unhandled parameter %s\n", argv[0]);
    usage (argv0[0], NULL, pl);
  }

  polyselect_main_data_parse_Nd(main_data, pl);
  ASSERT_ALWAYS (2 <= main_data->d && main_data->d <= 3);
  polyselect_main_data_parse_ad_range(main_data, pl);
  polyselect_main_data_parse_P(main_data, pl);
  param_list_parse_ulong(pl, "nq", &main_data->nq);

  param_list_parse_int (pl, "t", &nthreads);
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

  /* XXX we have no sopt parameter here */

  /* used for computation of E */
  if (param_list_parse_double (pl, "area", &area) == 0) /* no -area */
    area = AREA;
  if (param_list_parse_double (pl, "Bf", &bound_f) == 0) /* no -Bf */
    bound_f = BOUND_F;
  if (param_list_parse_double (pl, "Bg", &bound_g) == 0) /* no -Bg */
    bound_g = BOUND_G;



  if (!param_list_parse_mpz(pl, "skewness", maxS))
    mpz_set_ui (maxS, 0);
  else if (mpz_cmp_ui (maxS, 1) < 0)
  {
    gmp_fprintf(stderr, "Error, skewness (%Zd) should be greater or equal "
                        "to 1\n", maxS);
    abort();
  }

  if (param_list_warn_unused(pl))
    usage (argv0[0], NULL, pl);

  /* print command line */
  verbose_interpret_parameters(pl);
  param_list_print_command_line (stdout, pl);

  /* quiet mode */
  if (quiet == 1)
    main_data->verbose = -1;

  /* Compute maxS: use maxS if maxS argument is greater than 0 and lesser than
     default value */
  {
      mpz_t tmp;
      mpz_init (tmp);
      compute_default_max_skew (tmp, main_data->N, 2);
      if (mpz_cmp_ui(maxS, 0) == 0 || mpz_cmp(maxS, tmp) > 0)
          mpz_set (maxS, tmp);
      mpz_clear(tmp);
  }

  polyselect_main_data_prepare_primes(main_data);

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

      polyselect_main_data_commit_stats(main_data, loc->stats, loc->ad);

      polyselect_thread_locals_clear(loc);
  }


  printf ("# Stat: tried %lu ad-value(s), found %lu polynomial(s)\n",
          main_data->stats->number_of_ad_values,
          found);

  polyselect_stats_display_final(main_data->stats, main_data->verbose);

  if (best_E == DBL_MAX)
    /* This line is required by the script: */
    printf ("# No polynomial found, please increase the ad range or decrease P\n");

  mpz_clear (maxS);
  param_list_clear (pl);

  polyselect_main_data_clear(main_data);

  return 0;
}
