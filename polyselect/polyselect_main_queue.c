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
#include "portability.h"


/* Read-Only */

#if 0
const double exp_rot[] = { 0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 0 };

unsigned int sopt_effort = SOPT_DEFAULT_EFFORT;	/* size optimization effort */
double best_E = 0.0;		/* Murphy's E (the larger the better) */

/* read-write global variables */
double potential_collisions = 0.0;
#endif

/* print poly info */
static size_t
snprintf_poly_info(char *buf,
		  size_t size,
		  mpz_poly_srcptr f,
		  mpz_poly_srcptr g,
		  mpz_srcptr n,
		  const int raw)
{
  size_t np = 0;

  np += snprintf(buf + np, size - np, "# %s polynomial:\n", raw ? "Raw" : "Size-optimized");
  ASSERT_ALWAYS(np <= size);

  char * str;

  np += gmp_snprintf(buf + np, size - np, "%sn: %Zd\n", raw ? "# " : "", n);
  ASSERT_ALWAYS(np <= size);

  mpz_poly_asprintf_cado_format(&str, g, 'Y', raw ? "# " : "");
  np += strlcat(buf + np, str, size - np);
  ASSERT_ALWAYS(np <= size);
  free(str);

  mpz_poly_asprintf_cado_format(&str, f, 'c', raw ? "# " : "");
  np += strlcat(buf + np, str, size - np);
  ASSERT_ALWAYS(np <= size);
  free(str);

  double skew = L2_skewness(f, SKEWNESS_DEFAULT_PREC);
  unsigned int nroots = mpz_poly_number_of_real_roots(f);
  double logmu = L2_lognorm(f, skew);
  double exp_E = logmu + expected_rotation_gain(f, g);

  np += snprintf(buf + np, size - np, "# %sexp_E", raw ? "raw " : "");
  ASSERT_ALWAYS(np <= size);

  np += snprintf(buf + np, size - np,
	       " %1.2f, lognorm %1.2f, skew %1.2f, %u rroots\n", exp_E,
	       logmu, skew, nroots);
  ASSERT_ALWAYS(np <= size);

  return np;
}

void polyselect_fprintf_poly_pair(FILE * fp, mpz_srcptr N,                    mpz_poly_srcptr f, mpz_poly_srcptr g, int raw)
{
    if (!f && !g) return;
    size_t sz = mpz_sizeinbase(N, 10);
    int length = sz * 12;
    char *str = malloc(length);
    ASSERT_ALWAYS(str);
    snprintf_poly_info(str, length, f, g, N, raw);
    fprintf(fp, "%s", str);
    free(str);
}


/* return 1 if the polynomial is ok and among the best ones,
   otherwise return 0

   This modifies both F and g
*/
int optimize_raw_poly(mpz_poly_ptr f, mpz_poly_ptr g,
                        polyselect_main_data_srcptr main,
                        polyselect_stats_ptr stats)
{
  double skew;
  mpz_t t;
  double st, logmu, exp_E;

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

  /* register all stat to the stats object. This is a local object, so no
   * lock needed !
   */
  {
    stats->collisions_good++;
    polyselect_priority_queue_push(stats->best_opt_logmu, logmu);
    polyselect_priority_queue_push(stats->best_exp_E, exp_E);
    polyselect_data_series_add(stats->opt_lognorm, logmu);
    polyselect_data_series_add(stats->exp_E, exp_E);
    polyselect_data_series_add(stats->opt_proj_alpha,
            get_alpha_projective(f, get_alpha_bound()));
  }

  return 1;
}
