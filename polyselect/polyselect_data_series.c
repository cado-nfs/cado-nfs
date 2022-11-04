#include "cado.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include "polyselect_data_series.h"
#include "macros.h"
#include "gcd.h"

void
polyselect_data_series_init (polyselect_data_series_ptr s)
{
  s->size = s->alloc = 0;
  s->x = NULL;
  s->sum = s->sum2 = 0.0;
  s->min = DBL_MAX;
  s->max = -DBL_MAX;
  s->rstate = NULL;
}

void
polyselect_data_series_reset (polyselect_data_series_ptr s)
{
  s->size = 0;
  s->sum = s->sum2 = 0.0;
  s->min = DBL_MAX;
  s->max = -DBL_MAX;
}

void
polyselect_data_series_clear (polyselect_data_series_ptr s)
{
  free (s->x);
}

void
polyselect_data_series_add (polyselect_data_series_ptr s, double x)
{
  if (s->size == s->alloc)
    {
      s->alloc += 1 + s->alloc / 2;
      s->x = realloc (s->x, s->alloc * sizeof (double));
    }
  s->x[s->size++] = x;
  s->sum += x;
  s->sum2 += x * x;
  if (x < s->min)
    s->min = x;
  if (x > s->max)
    s->max = x;
}

double
polyselect_data_series_mean (polyselect_data_series_srcptr s)
{
  return s->sum / (double) s->size;
}

double
polyselect_data_series_variance (polyselect_data_series_srcptr s)
{
  double m = polyselect_data_series_mean (s);
  return s->sum2 / (double) s->size - m * m;
}

double
polyselect_data_series_std_dev (polyselect_data_series_srcptr s)
{
    return sqrt(polyselect_data_series_variance(s));
}

/* given a distribution with mean m and variance v, estimate the parameters
   beta and eta from a matching Weibull distribution, using the method of
   moments:
   m = eta * gamma (1 + 1/beta)
   v = eta^2 * [gamma (1 + 2/beta) - gamma (1 + 1/beta)^2]
 */
void
polyselect_data_series_estimate_weibull_moments(double *beta, double *eta, polyselect_data_series_srcptr s)
{
  double m = polyselect_data_series_mean(s);
  double v = polyselect_data_series_variance(s);
  double y = sqrt(v) / m;

  y = y * (0.7796968012336761 + y * (0.61970313728462 + 0.0562963108244 * y));
  *beta = 1.0 / y;
  *eta = m * (1.0 + y * (0.57721566490153 - 0.655878071520 * y));
}

int polyselect_data_series_snprintf_summary(char * tmp, size_t size, polyselect_data_series_srcptr s)
{
    return snprintf(tmp, size, 
	      " (nr/min/av/max/std): %lu/%1.2f/%1.2f/%1.2f/%1.2f",
	       s->size,
               s->min,
               polyselect_data_series_mean(s),
               s->max,
               polyselect_data_series_std_dev(s));
}


/* Estimation via extreme values: we cut the total n values into samples of k
   values, and for each sample we keep only the minimum. If the series of
   minimum values satisfies a Weibull distribution with parameters beta and eta,
   then the original one has parameters beta (identical) and eta*k^(1/beta).
   Here we choose k near sqrt(n). */
void
polyselect_data_series_estimate_weibull_moments2(double *beta, double *eta, polyselect_data_series_srcptr s)
{
  unsigned long n = s->size;
  unsigned long i, j, k, p, u;
  polyselect_data_series_t smin;
  double min, eta_min;

  ASSERT_ALWAYS(n > 0);

  /* if s->rstate is NULL, then this series must be adjusted by
   * polyselect_stats_init to receive the pointer to the parent random
   * state */
  ASSERT_ALWAYS(s->rstate != NULL);

  polyselect_data_series_init(smin);

  k = (unsigned long) sqrt((double) n);	/* sample size */
  /* We consider full samples only. Since we call this function several times
     with the same sequence, we perform a random permutation of the sequence
     at each call to avoid side effects due to the particular order of
     elements. In practice instead of considering s[j] we consider
     s[(p*j) % n] where p is random with gcd(p,n)=1. */
  {
    do
      p = gmp_urandomm_ui(s->rstate, n);
    while (gcd_uint64(p, n) != 1);
  }

  for (i = 0; i + k <= n; i += k)
    {
      for (j = i, min = DBL_MAX; j < i + k; j++)
	{
	  u = (p * j) % n;
	  if (s->x[u] < min)
	    min = s->x[u];
	}
      polyselect_data_series_add(smin, min);
    }
  polyselect_data_series_estimate_weibull_moments(beta, &eta_min, smin);
  polyselect_data_series_clear(smin);
  *eta = eta_min * pow((double) k, 1.0 / *beta);
}

/* Add the contents of the data series "from" to the data series "to".
 * The input data series is not changed.
 */
void polyselect_data_series_merge(polyselect_data_series_ptr to, polyselect_data_series_srcptr from)
{
  if (to->size + from->size >= to->alloc)
    {
      to->alloc = (to->size + from->size) + to->alloc / 2;
      to->x = realloc (to->x, to->alloc * sizeof (double));
    }
  for(size_t d = 0 ; d < from->size ; d++)
      to->x[to->size++] = from->x[d];
  to->sum += from->sum;
  to->sum2 += from->sum2;
  to->min = MIN(to->min, from->min);
  to->max = MAX(to->max, from->max);
}

