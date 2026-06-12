#ifndef CADO_POLYSELECT_NORMS_H
#define CADO_POLYSELECT_NORMS_H

#include "mpz_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

extern double L2_lognorm (mpz_poly_srcptr, double);
extern double L2_skewness (mpz_poly_srcptr);
extern double L2_combined_skewness2 (mpz_poly_srcptr f, mpz_poly_srcptr g);
extern double L2_skew_lognorm (mpz_poly_srcptr);


#ifdef __cplusplus
}
#endif

#endif	/* CADO_POLYSELECT_NORMS_H */
