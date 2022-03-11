#ifndef POLYSELECT_DATA_SERIES_H_
#define POLYSELECT_DATA_SERIES_H_

#include "gmp_aux.h"

#ifdef __cplusplus
extern "C" {
#endif

/* structure for linear series of double-precision real number data. This
 * is meant to do keep track for 1- and 2- order moments so that mean and
 * variance can be accessed in O(1).
 */
struct polyselect_data_series_s {
  unsigned long size;  /* number of values */
  unsigned long alloc; /* allocated size */
  double *x;           /* values */
  double sum;          /* sum = x[0] + ... + x[size-1] */
  double sum2;          /* var = x[0]^2 + ... + x[size-1]^2 */
  double min, max;     /* minimum and maximum values */

  /* This may be initialized to a random state that is safe for use for
   * sampling means for this series. Presently the only api call that
   * uses this is polyselect_data_series_estimate_weibull_moments2, so
   * that only the data series on which this function is actually called
   * needs to have this set up.
   *
   * Because it's such a specific use case, this initialization is not
   * done by default in the ctor. *IF* the code using this interface is
   * going to call polyselect_data_series_estimate_weibull_moments2, then
   * we require that this code made provisions to have this pointer point
   * to a valid random state after calling the ctor.  (polyselect_stats
   * does just that).
   */
  gmp_randstate_ptr rstate;
};
typedef struct polyselect_data_series_s polyselect_data_series_t[1];
typedef struct polyselect_data_series_s * polyselect_data_series_ptr;
typedef const struct polyselect_data_series_s * polyselect_data_series_srcptr;

extern void polyselect_data_series_init (polyselect_data_series_ptr);
extern void polyselect_data_series_clear (polyselect_data_series_ptr);
extern void polyselect_data_series_reset (polyselect_data_series_ptr s);
extern void polyselect_data_series_add (polyselect_data_series_ptr, double);
extern double polyselect_data_series_mean (polyselect_data_series_srcptr);
extern double polyselect_data_series_variance (polyselect_data_series_srcptr);
extern double polyselect_data_series_std_dev (polyselect_data_series_srcptr);

extern void polyselect_data_series_estimate_weibull_moments(double *, double *, polyselect_data_series_srcptr);
extern void polyselect_data_series_estimate_weibull_moments2(double *, double *, polyselect_data_series_srcptr);
extern int polyselect_data_series_snprintf_summary(char * tmp, size_t size, polyselect_data_series_srcptr s);
extern void polyselect_data_series_merge(polyselect_data_series_ptr to, polyselect_data_series_srcptr from);


#ifdef __cplusplus
}
#endif

#endif	/* POLYSELECT_DATA_SERIES_H_ */
