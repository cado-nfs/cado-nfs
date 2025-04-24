#ifndef CADO_POLYSELECT_NORMS_H
#define CADO_POLYSELECT_NORMS_H

#include "mpz_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

#define SKEWNESS_DEFAULT_PREC 10

extern double L2_lognorm (mpz_poly_srcptr, double);
extern double L2_skewness (mpz_poly_srcptr, int);
extern double L2_combined_skewness2 (mpz_poly_srcptr f, mpz_poly_srcptr g, int prec);
extern double L2_skew_lognorm (mpz_poly_srcptr, int);


#ifdef __cplusplus
}
#endif

#endif	/* CADO_POLYSELECT_NORMS_H */
