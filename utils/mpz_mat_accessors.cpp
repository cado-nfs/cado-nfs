#include "cado.h" // IWYU pragma: keep

#include <gmp.h>

#include "macros.h"
#include "mpz_mat.h"
#include "mpz_mat_accessors.h"
#include "mpz_poly.h"
#include "runtime_numeric_cast.hpp"

/* conversion of rows and columns to polynomials*/

#if 0
void mpz_mat_row_to_poly(mpz_poly_ptr f, mpz_mat_srcptr M, const unsigned int i)
{
    mpz_poly_realloc(f,M->n);
    unsigned int j;
    for (j = 0 ; j < M->n; j++){
        mpz_poly_setcoeff(f, j, mpz_mat_entry_const(M,i,j));
    }
    mpz_poly_cleandeg(f, M->n - 1);
}
#endif

void mpz_mat_row_to_poly_rev(mpz_poly_ptr f, mpz_mat_srcptr M, const unsigned int i)
{
    mpz_poly_realloc(f,M->n);
    for (int j = 0 ; j < runtime_numeric_cast<int>(M->n); j++){
        mpz_poly_setcoeff(f, j, mpz_mat_entry_const(M, i, M->n-1-j));
    }
    mpz_poly_cleandeg(f, runtime_numeric_cast<int>(M->n - 1));
}

#if 0
void mpz_mat_column_to_poly(mpz_poly_ptr f, mpz_mat_srcptr M, const unsigned int j)
{
    mpz_poly_realloc(f,M->m);
    unsigned int i;
    for (i = 0 ; i < M->m; i++){
        mpz_poly_setcoeff(f,i,mpz_mat_entry_const(M,i,j));
    }
    mpz_poly_cleandeg(f, M->m - 1);
}
#endif

void mpq_mat_row_to_poly(mpz_poly_ptr f, mpz_ptr lcm, mpq_mat_srcptr M, const unsigned int i)
{
    /* read element w[i] as a polynomial. beware, the element might
     * have rational coordinates, that's life ! */
    unsigned int const n = M->n;
    ASSERT_ALWAYS(i < M->m);
    mpz_set_ui(lcm, 1);
    for (unsigned int j = 0; j < n; j++) {
        mpz_lcm(lcm, lcm, mpq_denref(mpq_mat_entry_const(M, i, j)));
    }
    mpz_poly_realloc(f, n);
    for (int j = 0; j < runtime_numeric_cast<int>(n); j++) {
        mpq_srcptr mij = mpq_mat_entry_const(M, i, j);
        mpz_divexact(mpz_poly_coeff(f, j), lcm, mpq_denref(mij));
        mpz_mul(mpz_poly_coeff(f, j), mpz_poly_coeff_const(f, j), mpq_numref(mij));
    }
    mpz_poly_cleandeg(f, runtime_numeric_cast<int>(n)-1);
}

void mpq_poly_to_mat_row(mpq_mat_ptr M, const unsigned int i, mpz_poly_srcptr f, mpz_srcptr denom)
{
    ASSERT_ALWAYS(f->deg < (int) M->n);
    for (int j = 0 ; j < runtime_numeric_cast<int>(M->n); j++){
        mpq_ptr mij = mpq_mat_entry(M,i,j);
        mpz_set(mpq_numref(mij), mpz_poly_coeff_const(f, j));
        mpz_set(mpq_denref(mij), denom);
        mpq_canonicalize(mij);
    }
}

#if 0
void mpq_mat_column_to_poly(mpz_poly_ptr f, mpz_ptr lcm, mpq_mat_srcptr M, const unsigned int j)
{
    mpz_poly_realloc(f,M->m);
    mpz_set_si(lcm,1);
    for (unsigned int i = 0 ; i < M->m ; i++) {
        mpz_lcm(lcm,lcm,mpq_denref(mpq_mat_entry_const(M,i,j)));
    }
    for (unsigned int i = 0 ; i < M->m; i++){
        mpq_srcptr mij = mpq_mat_entry_const(M, i, j);
        mpz_divexact(mpz_poly_coeff(f, i), lcm, mpq_denref(mij));
        mpz_mul(mpz_poly_coeff(f, i), mpz_poly_coeff_const(f, i), mpq_numref(mij));
    }
    mpz_poly_cleandeg(f, M->m - 1);
}
#endif
/**/


