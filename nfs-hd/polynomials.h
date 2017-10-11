#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <stdint.h>
#include "cado.h"
#include "utils.h"

unsigned int mpz_poly_is_reciprocal(mpz_poly_srcptr f);

/*
  f0: the first polynomial.
  f1: the second polynomial.
  g: a polynomial to build f0 and f1, must have the degree n (extension of the
   field).
  a: the a to build f1. (MNFS)
  b: the b to build f1.
  p: the prime number that define the field.
  h: a poylnomial to build f0 and f1.
  coeff0: lower bound of the coefficients of g.
  coeff1: upper bound of the coefficients of g.
  q: a special-q.
  t: dimension of the lattice.
*/
double function_special_q(mpz_poly_ptr f0, mpz_poly_ptr f1, mpz_poly_ptr g,
    mpz_ptr a, mpz_ptr b, mpz_ptr c, mpz_srcptr p, mpz_poly_ptr h,
    int * coeff, unsigned int q, unsigned int t,
    unsigned int nb_times, double * weight, int c_tol, gmp_randstate_t state,
    unsigned int h_set, unsigned int gal);

void function_classical(mpz_poly_ptr f0, mpz_poly_ptr f1, mpz_srcptr p,
    unsigned int n, int * coeff, unsigned int nb_times,
    double * weight, gmp_randstate_t state);

#endif // FUNCTIONS_H
