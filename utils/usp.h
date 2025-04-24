#ifndef CADO_USP_H
#define CADO_USP_H

#include "mpz_poly.h"
#include <gmp.h> /* for mpz_t */

/* this structure represents the interval [a/2^ka, b/2^kb] */
struct usp_root_interval_s {
  mpz_t a;
  int ka;
  mpz_t b;
  int kb;
};
typedef struct usp_root_interval_s usp_root_interval[1];
typedef struct usp_root_interval_s * usp_root_interval_ptr;
typedef const struct usp_root_interval_s * usp_root_interval_srcptr;

#ifdef __cplusplus
extern "C" {
#endif

/* This is the same as mpz_poly_number_of_real_roots, but we also fill
 * R[0] to R[nroots-1] with disjoint intervals that locate the real
 * roots. Those intervals can be refined lated on.
 *
 * T (if not zero) is a bound on the absolute value of the real roots.
 */
int mpz_poly_number_of_real_roots_extra(mpz_poly_srcptr f, double T, usp_root_interval * R);

/* Refine the root of f given by the interval r */
double usp_root_interval_refine (usp_root_interval_ptr r, mpz_poly_srcptr P, double precision);

void usp_root_interval_init (usp_root_interval_ptr R);
void usp_root_interval_clear (usp_root_interval_ptr R);

void usp_root_interval_set (usp_root_interval_ptr R, mpz_srcptr a, int ka, mpz_srcptr b, int kb);
void usp_root_interval_set_ui (usp_root_interval_ptr R, unsigned long a, int ka, unsigned long b, int kb);

#ifdef __cplusplus
}
#endif

#endif  /* CADO_USP_H */
