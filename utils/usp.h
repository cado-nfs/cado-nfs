#ifndef USP_H_
#define USP_H_

#include <gmp.h> /* for mpz_t */
#include "mpz_poly.h"

/* this structure represents the interval [a/2^ka, b/2^kb] */
typedef struct {
  mpz_t a;
  int ka;
  mpz_t b;
  int kb;
} usp_root_data;

#ifdef __cplusplus
extern "C" {
#endif

int numberOfRealRoots (mpz_t *p, int n, double T, int verbose, usp_root_data *R);
double rootRefine (usp_root_data *r, mpz_t *p, int n, double precision);
void usp_root_data_init (usp_root_data *R);
void usp_root_data_clear (usp_root_data *R);

#ifdef __cplusplus
}
#endif

#endif  /* USP_H_ */
