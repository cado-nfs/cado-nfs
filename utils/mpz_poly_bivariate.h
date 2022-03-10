#ifndef MPZ_POLY_BIVARIATE_H
#define MPZ_POLY_BIVARIATE_H 

#include <stdio.h>
#include <gmp.h>
#include "mpz_poly.h"

/*
 * Bivariate polynomial are written as polynomial in y with coefficients in
 * mpz_poly in x.
 */

typedef struct {
  unsigned int alloc;
  int deg_y;
  int deg_x;
  mpz_poly * coeff;
} mpz_poly_bivariate_struct_t;

typedef mpz_poly_bivariate_struct_t mpz_poly_bivariate_t[1];
typedef mpz_poly_bivariate_struct_t * mpz_poly_bivariate_ptr;
typedef const mpz_poly_bivariate_struct_t * mpz_poly_bivariate_srcptr;

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Initialize a mpz_poly_bivariate with alloc = d + 1.
 */
void mpz_poly_bivariate_init(mpz_poly_bivariate_ptr f, int d);
/*
 * Initialize a mpz_poly_bivariate with alloc = d + 1 and alloc of all mpz_poly
 * equals to dx + 1.
 */
void mpz_poly_bivariate_init_y_x(mpz_poly_bivariate_ptr f, int dy, int dx);
/*
 * Increase alloc.
 */
void mpz_poly_bivariate_realloc(mpz_poly_bivariate_ptr f, unsigned int nc);
/*
 * Increase alloc and set alloc of all mpz_poly equals dx + 1.
 */
void mpz_poly_bivariate_realloc_x(mpz_poly_bivariate_ptr f, unsigned int nc, int dx);
/*
 * Clear a mpz_poly_bivariate.
 */
void mpz_poly_bivariate_clear(mpz_poly_bivariate_ptr f);
/*
 * Find deg_y of f.
 */
void mpz_poly_bivariate_cleandeg(mpz_poly_bivariate_ptr f, int deg_y);
/*
 * Set the ith coefficient of f.
 */
void mpz_poly_bivariate_setcoeff(mpz_poly_bivariate_ptr f, int i,
    mpz_poly_srcptr z);
/*
 * Print f in a file.
 */
void mpz_poly_bivariate_fprintf(FILE * fp, mpz_poly_bivariate_srcptr f);
/*
 * Compute res(x) = f(x, y).
 */
void mpz_poly_bivariate_eval_y(mpz_poly_ptr res,
    mpz_poly_bivariate_srcptr f, mpz_srcptr y);
/*
 * Compute res(y) = f(x, y).
 */
void mpz_poly_bivariate_eval_x(mpz_poly_ptr res,
    mpz_poly_bivariate_srcptr f, mpz_srcptr x);

/*
 * Compute resultant(x) = res(f, g).
 */
void mpz_poly_bivariate_resultant_y(mpz_poly_ptr resultant,
    mpz_poly_bivariate_srcptr f, mpz_poly_bivariate_srcptr g);

/*
 * Compute resultant(y) = res(f, g).
 */
void mpz_poly_bivariate_resultant_x(mpz_poly_ptr resultant,
    mpz_poly_bivariate_srcptr f, mpz_poly_bivariate_srcptr g);


#ifdef __cplusplus
}
#endif

#endif /* MPZ_POLY_BIVARIATE_H */
