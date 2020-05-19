#ifndef MPQS_DOIT_H_
#define MPQS_DOIT_H_

#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

void mpqs_doit (mpz_t f, const mpz_t N0, int verbose);
void smooth_stat (int);

#ifdef __cplusplus
}
#endif

#endif	/* MPQS_DOIT_H_ */
