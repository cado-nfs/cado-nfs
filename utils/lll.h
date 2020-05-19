#ifndef CADO_UTILS_LLL_H_
#define CADO_UTILS_LLL_H_

#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  mpz_t **coeff;
  int NumRows, NumCols;
} mat_Z;

void LLL_init (mat_Z *, int, int);
void LLL_clear (mat_Z *);
long LLL (mpz_t det, mat_Z B, mat_Z* U, mpz_srcptr a, mpz_srcptr b);

#ifdef __cplusplus
}
#endif

#endif	/* CADO_UTILS_LLL_H_ */

