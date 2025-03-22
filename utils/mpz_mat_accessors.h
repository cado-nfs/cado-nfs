#ifndef CADO_UTILS_MPZ_MAT_ACCESSORS_H
#define CADO_UTILS_MPZ_MAT_ACCESSORS_H

#include <gmp.h>
#include "mpz_mat.h"
#include "mpz_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

/*{{{ conversion of rows and columns to polynomials*/
// static void mpz_mat_row_to_poly(mpz_poly_ptr f, mpz_mat_srcptr M, const unsigned int i);
void mpz_mat_row_to_poly_rev(mpz_poly_ptr f, mpz_mat_srcptr M, const unsigned int i);
// static void mpz_mat_column_to_poly(mpz_poly_ptr f, mpz_mat_srcptr M, const unsigned int j);
void mpq_mat_row_to_poly(mpz_poly_ptr f, mpz_ptr lcm, mpq_mat_srcptr M, const unsigned int i);
void mpq_poly_to_mat_row(mpq_mat_ptr M, const unsigned int i, mpz_poly_srcptr f, mpz_srcptr denom);
// static void mpq_mat_column_to_poly(mpz_poly_ptr f, mpz_ptr lcm, mpq_mat_srcptr M, const unsigned int j);
/* }}} */


#ifdef __cplusplus
}
#endif

#endif	/* CADO_UTILS_MPZ_MAT_ACCESSORS_H */
