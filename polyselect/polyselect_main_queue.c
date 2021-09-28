#include "cado.h"
#include <stdbool.h>		// bool
#include <float.h>		// DBL_MAX
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include "size_optimization.h"
#include "auxiliary.h"
#include "gcd.h"		// for gcd_ul
#include "timing.h"		// for seconds
#include "polyselect_main_queue.h"
#include "polyselect_locals.h"
#include "polyselect_norms.h"


/* Read-Only */

#if 0
const double exp_rot[] = { 0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 0 };

unsigned int sopt_effort = SOPT_DEFAULT_EFFORT;	/* size optimization effort */
double best_E = 0.0;		/* Murphy's E (the larger the better) */

/* read-write global variables */
double potential_collisions = 0.0;
#endif

/* print poly info
 *
 * XXX this does more than printing, as data is inserted in the stats
 * table as well. This is ugly.
 * */
void
sprintf_poly_info(char *buf,
		  size_t size,
		  mpz_poly_srcptr f,
		  mpz_poly_srcptr g,
		  mpz_srcptr n,
		  const int raw, polyselect_thread_locals_ptr loc)
{
  unsigned int i, nroots;
  double skew, logmu, exp_E;

  size_t np = 0;

  const char *prefix = "";

  if (raw)
    prefix = "# ";

  if (raw)
    {
      np += snprintf(buf + np, size - np, "# Raw polynomial:\n");
  } else
    {
      snprintf(buf + np, size - np, "# Size-optimized polynomial:\n");
    }

  np += gmp_snprintf(buf + np, size - np, "%sn: %Zd\n", prefix, n);
  for (i = mpz_poly_degree(g) + 1; i-- != 0;)
    np +=
	gmp_snprintf(buf + np, size - np, "%sY%u: %Zd\n", prefix, i,
		     g->coeff[i]);
  for (i = mpz_poly_degree(f) + 1; i-- != 0;)
    np +=
	gmp_snprintf(buf + np, size - np, "%sc%u: %Zd\n", prefix, i,
		     f->coeff[i]);
  skew = L2_skewness(f, SKEWNESS_DEFAULT_PREC);
  nroots = mpz_poly_number_of_real_roots(f);
  logmu = L2_lognorm(f, skew);
  exp_E = logmu + expected_rotation_gain(f, g);
  if (raw == 1)
    np += snprintf(buf + np, size - np, "# raw exp_E");
  else
    np += snprintf(buf + np, size - np, "# exp_E");

  np +=
      snprintf(buf + np, size - np,
	       " %1.2f, lognorm %1.2f, skew %1.2f, %u rroots\n", exp_E,
	       logmu, skew, nroots);

  ASSERT_ALWAYS(np < size);

  if (raw)
    return;

  if (loc->main->target_E != 0.0 || loc->main->maxtime < DBL_MAX)
    {
      double beta, eta, prob;
      polyselect_data_series_srcptr exp_E = loc->main->stats->exp_E;

      /* estimate the parameters of a Weibull distribution for E */
      /* XXX we're reading a global state, so we need the lock! */
      pthread_mutex_lock(&loc->main->stats_lock);
      polyselect_data_series_estimate_weibull_moments2(&beta, &eta, exp_E);
      unsigned long collisions_good = loc->main->stats->collisions_good;
      pthread_mutex_unlock(&loc->main->stats_lock);

      if (loc->main->target_E != 0.0)
	{
	  prob = 1.0 - exp(-pow(loc->main->target_E / eta, beta));
	  if (prob == 0)	/* for x small, exp(x) ~ 1+x */
	    prob = pow(loc->main->target_E / eta, beta);
	  np += snprintf(buf + np, size - np,
			 "# E: %lu, min %.2f, avg %.2f, max %.2f, stddev %.2f\n",
			 exp_E->size,
			 exp_E->min,
			 polyselect_data_series_mean(exp_E),
			 exp_E->max, polyselect_data_series_std_dev(exp_E));
	  np +=
	      snprintf(buf + np, size - np,
		       "# loc->main->target_E=%.2f: collisions=%.2e, time=%.2e"
		       " (beta %.2f,eta %.2f)\n", loc->main->target_E, 1.0 / prob,
		       seconds() / (prob * collisions_good), beta, eta);
      } else
	{			/* loc->main->maxtime < DBL_MAX */
	  /* time = seconds () / (prob * collisions_good)
	     where  prob = 1 - exp (-(E/eta)^beta) */
	  unsigned long n = collisions_good;	/* #polynomials found so far */
	  double time_so_far = seconds();
	  double time_per_poly = time_so_far / n;	/* average time per poly */
	  double admin_d = mpz_get_d(loc->main->admin);
	  double ad_d = mpz_get_d(loc->ad);
	  double adrange = (ad_d - admin_d) * (loc->main->maxtime / time_so_far);
	  /* WTF ??? I'm commenting that line, but what does it mean ???
	   * Is it a leftover from something?
	   */
	  // adrange = 2.00e+15 - 99900000000000.0;
	  prob = time_so_far / (loc->main->maxtime * n);
	  double E = eta * pow(-log(1 - prob), 1.0 / beta);

	  /* XXX we're modifying a global state, so we need the lock! */
	  pthread_mutex_lock(&loc->main->stats_lock);
          polyselect_data_series_ptr best_exp_E_Weibull = loc->main->stats->best_exp_E_Weibull;
	  polyselect_data_series_add(best_exp_E_Weibull, E);
	  /* since the values of (eta,beta) fluctuate a lot, because
	     they depend on the random samples in estimate_weibull_moments2,
	     we take the average value for best_exp_E */
	  E = polyselect_data_series_mean(best_exp_E_Weibull);
	  pthread_mutex_lock(&loc->main->stats_lock);

	  np += snprintf(buf + np, size - np,
			 "# %.2fs/poly, eta %.2f, beta %.3f, admax %.2e, best exp_E %.2f\n",
			 time_per_poly, eta, beta, admin_d + adrange, E);
	}
    }

  np += snprintf(buf + np, size - np, "\n");
  ASSERT_ALWAYS(np < size);
}


void
output_polynomials(mpz_poly_srcptr f_old, mpz_poly_srcptr g_old,
		   const mpz_t N,
		   mpz_poly_srcptr f, mpz_poly_srcptr g,
		   polyselect_thread_locals_ptr loc)
{
  size_t sz = mpz_sizeinbase(N, 10);
  int length = sz * 12;
  char *str_old = malloc(length);
  ASSERT_ALWAYS(str_old);
  char *str = malloc(length);
  ASSERT_ALWAYS(str);
  if (f_old && g_old)
    {
      sprintf_poly_info(str_old, length, f_old, g_old, N, 1, loc);
    }
  sprintf_poly_info(str, length, f, g, N, 0, loc);

#ifdef HAVE_OPENMP
#pragma omp critical
#endif
  {
    if (f_old != NULL && g_old != NULL)
      if (str_old != NULL)
	printf("%s", str_old);
    if (str != NULL)
      printf("%s", str);
    fflush(stdout);
  }

  if (str_old != NULL)
    free(str_old);
  if (str != NULL)
    free(str);
}

/* Insert a value into a sorted array of length len.
   Returns 1 if element was inserted, 0 if it was too big */
static int
sorted_insert_double(double *array, const size_t len, const double value)
{
  size_t k;
  int result = 0;
  if (len == 0)
    return 0;
#ifdef HAVE_OPENMP
#pragma omp critical
#endif
  if (value < array[len - 1])
    {
      for (k = len - 1; k > 0 && value < array[k - 1]; k--)
	array[k] = array[k - 1];
      array[k] = value;
      result = 1;
    }
  return result;
}

/* return 1 if the polynomial is ok and among the best ones,
   otherwise return 0

   This modifies both F and g
*/
int optimize_raw_poly(mpz_poly_ptr f, mpz_poly_ptr g,
		      polyselect_main_data_ptr main)
{
  double skew;
  mpz_t t;
  double st, logmu, exp_E;

  polyselect_stats_ptr stats = main->stats;

  /* check that the algebraic polynomial has content 1, otherwise skip it */
  mpz_init(t);
  mpz_poly_content(t, f);
  if (mpz_cmp_ui(t, 1) != 0)
    {
      mpz_clear(t);
      return 0;
    }
  mpz_clear(t);

  /* optimize size */

  st = seconds_thread();
  size_optimization(f, g, f, g, main->sopt_effort, main->verbose);
  st = seconds_thread() - st;
  stats->optimize_time += st;
  stats->opt_found++;

  /* polynomials with f[d-1] * f[d-3] > 0 *after* size-optimization
     give worse exp_E values */
  int d = f->deg;
  if (mpz_sgn(f->coeff[d - 1]) * mpz_sgn(f->coeff[d - 3]) > 0)
    {
      stats->discarded2++;
      return 0;
    }

  skew = L2_skewness(f, SKEWNESS_DEFAULT_PREC);
  logmu = L2_lognorm(f, skew);
  /* expected_rotation_gain() takes into account the projective alpha */
  exp_E = logmu + expected_rotation_gain(f, g);

  sorted_insert_double(stats->best_opt_logmu, stats->keep, logmu);
  sorted_insert_double(stats->best_exp_E, stats->keep, exp_E);

  {
    stats->collisions_good++;
    polyselect_data_series_add(stats->opt_lognorm, logmu);
    polyselect_data_series_add(stats->exp_E, exp_E);
    polyselect_data_series_add(stats->opt_proj_alpha,
            get_alpha_projective(f, get_alpha_bound()));
  }

  return 1;
}
