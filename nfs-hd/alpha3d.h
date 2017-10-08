#ifndef ALPHA3D_H
#define ALPHA3D_H 

#include "cado.h"
#include "utils.h"
#include "mpz_poly.h"

double alpha3d(mpz_poly_srcptr f, unsigned long p, gmp_randstate_t rstate, unsigned int N);

#endif /* ALPHA3D_H */
