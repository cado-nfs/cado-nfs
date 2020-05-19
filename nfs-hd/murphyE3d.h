#ifndef MURPHYE3D_H
#define MURPHYE3D_H 

#include "cado_poly.h"
#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

double murphyE3d(cado_poly_srcptr f, double * lpb, double volume,
    unsigned int N, int q_side, double Q, double s,
    unsigned long p, gmp_randstate_t rstate, unsigned int N_alpha);

#ifdef __cplusplus
}
#endif


#endif /* MURPHYE3D_H */
