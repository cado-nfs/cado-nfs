#include "cado.h"
#include "bblas.hpp"

void m4ri_plu_tests(int n __attribute__((unused)))
{
#ifdef  HAVE_M4RI
    mzd_t * M;
    mzd_t * LU;
    mzp_t * P, *Q;
    M = mzd_init(n, n);
#if 0
    mzd_set_mem(M, m, n);
    uint64_t * m = new uint64_t[n*n/64];
    memfill_random(m, (64) * sizeof(uint64_t), rstate);
    delete[] m;
#else
    my_mzd_randomize(M);
#endif
    LU = mzd_init(n,n);
    P = mzp_init(n);
    Q = mzp_init(n);
    TIME1N(2, mzd_mypluq, LU, M, P, Q, 0);
    TIME1N(2, mzd_myechelonize_m4ri, LU, M, 0, 0);
    TIME1N(2, mzd_myechelonize_pluq, LU, M, 0);
    mzd_free(M);
    mzd_free(LU);
    mzp_free(P);
    mzp_free(Q);
#endif
}


