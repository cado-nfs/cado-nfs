/**
 * \file cado-nfs/polyselect/gfpkdlpolyselect.h
 * 
 * \date 21/08/2014
 * \author Aurore Guillevic
 * \email guillevic@lix.polytechnique.fr
 * \brief header file of gfpkdlpolyselect.c
 *        Compute two polynomials f, g suitable for discrete logarithm in 
 *        extension fields, with the conjugation method.
 *
 * \test TODO
 */

#ifndef CADO_GFPKDLPOLYSELECT_H
#define CADO_GFPKDLPOLYSELECT_H

#define ERROR_UNSPECIFIED 0
#define POLY_OK 1
#define ERROR_DEGREE 2
#define ERROR_NOT_IRREDUCIBLE 3
#define ERROR_NO_DEGREE_FACTOR 4
#define ERROR_SIGNATURE 5
#define ERROR_SUBFIELD 6
#define ERROR_SMEXP 7
#define ERROR_EASYBADIDEALS 8
#define ERROR_VERYBADIDEALS 9
#define ERROR_MAX 10
#define DEG_PY 2

#if DEG_PY > 2
#error "the code works only for Py of degree <= 2, sorry."
#endif

#include <stdio.h>      // FILE
#include <stdbool.h>    // for bool (in C)
#include <gmp.h>
#include "cado_poly.h"   // for MAX_DEGREE
#include "mpz_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

void gfpk_print_params(int n, mpz_srcptr p, mpz_srcptr ell);

// the function to call for generating a .poly file.
// works only for n=2 at the moment.
// , mpz_t ell, unsigned int mnfs
int gfpkdlpolyselect(int n, mpz_srcptr p, mpz_srcptr ell, int mnfs, const char* out_filename);

#ifdef __cplusplus
}
#endif

#endif // CADO_GFPKDLPOLYSELECT_H
