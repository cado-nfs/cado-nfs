#ifndef ANT_CONVERSIONS_H_
#define ANT_CONVERSIONS_H_

#include "mpz_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

/*{{{ conversion of rows and columns to polynomials*/
/* in ant.cpp */
void mpz_mat_row_to_poly(mpz_poly_ptr f, mpz_mat_srcptr M, const unsigned int i);
void mpz_mat_row_to_poly_rev(mpz_poly_ptr f, mpz_mat_srcptr M, const unsigned int i);
void mpz_mat_column_to_poly(mpz_poly_ptr f, mpz_mat_srcptr M, const unsigned int j);
void mpq_mat_row_to_poly(mpz_poly_ptr f, mpz_ptr lcm, mpq_mat_srcptr M, const unsigned int i);
void mpq_poly_to_mat_row(mpq_mat_ptr M, const unsigned int i, mpz_poly_srcptr f, mpz_srcptr denom);
void mpq_mat_column_to_poly(mpz_poly_ptr f, mpz_ptr lcm, mpq_mat_srcptr M, const unsigned int j);
/* }}} */

#ifdef __cplusplus
}
#endif

#endif	/* ANT_CONVERSIONS_H_ */
