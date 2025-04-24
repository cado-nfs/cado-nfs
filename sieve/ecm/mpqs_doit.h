#ifndef CADO_MPQS_DOIT_H
#define CADO_MPQS_DOIT_H

#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

void mpqs_doit (mpz_t f, const mpz_t N0, int verbose);
void smooth_stat (int);

#ifdef __cplusplus
}
#endif

#endif	/* CADO_MPQS_DOIT_H */
