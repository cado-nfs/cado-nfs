#ifndef ALPHA3D_H
#define ALPHA3D_H 

#include "mpz_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

double alpha3d(mpz_poly_srcptr f, unsigned long p, gmp_randstate_t rstate, unsigned int N);

#ifdef __cplusplus
}
#endif


#endif /* ALPHA3D_H */
