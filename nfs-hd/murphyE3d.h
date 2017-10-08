#ifndef MURPHYE3D_H
#define MURPHYE3D_H 

#include <gmp.h>

double murphyE3d(cado_poly_srcptr f, double * lpb, double volume,
    unsigned int N, int q_side, double Q, double s,
    unsigned long p, gmp_randstate_t rstate, unsigned int N_alpha);

#endif /* MURPHYE3D_H */
