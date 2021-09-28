#ifndef POLYSELECT_ALPHA_H_
#define POLYSELECT_ALPHA_H_

#include "mpz_poly.h"
#ifdef __cplusplus
extern "C" {
#endif

/* default prime bounds for the computation of alpha */
#define ALPHA_BOUND_SMALL  100
#ifndef ALPHA_BOUND /* allows to define ALPHA_BOUND in local.sh */
#define ALPHA_BOUND       2000
#endif

extern double get_alpha (mpz_poly_srcptr f, unsigned long B);

extern double get_alpha_projective (mpz_poly_srcptr f, unsigned long B);

extern void set_alpha_bound (unsigned long bound);

extern unsigned long get_alpha_bound (void);

extern double expected_alpha (double logK);

#ifdef __cplusplus
}
#endif

#endif	/* POLYSELECT_ALPHA_H_ */
