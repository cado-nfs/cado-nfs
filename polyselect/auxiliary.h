#ifndef POLYSELECT_AUXILIARY_H_
#define POLYSELECT_AUXILIARY_H_

/* header file for auxiliary routines for polyselect

Copyright 2008, 2009, 2010, 2013 Emmanuel Thome, Paul Zimmermann

This file is part of CADO-NFS.

CADO-NFS is free software; you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 2.1 of the License, or (at your option)
any later version.

CADO-NFS is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
details.

You should have received a copy of the GNU Lesser General Public License
along with CADO-NFS; see the file COPYING.  If not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include <stdio.h>      // FILE
#include <gmp.h>
#include "cado_poly.h"
#include "mpz_poly.h"
#include "gmp_aux.h"

/* The polynomial selection algorithms that use a linear polynomial will
 * put it on the side given by the following. */
// FIXME: atm, changing these does not work. It should...
#define RAT_SIDE 0
#define ALG_SIDE 1

extern double bound_f, bound_g, area;

/* The maximum degree supported is MAX_DEGREE, as defined in cado_poly.h */

#define NORM_MARGIN 0.2

typedef struct
{
  double kmin, kmax;
} rotation_space;

#ifdef __cplusplus
extern "C" {
#endif

/* alpha */
double special_valuation (mpz_poly_srcptr f, unsigned long p, mpz_srcptr disc, gmp_randstate_ptr rstate);
double special_valuation_affine (mpz_poly_srcptr f, unsigned long p, mpz_srcptr disc, gmp_randstate_ptr rstate);
double get_alpha (mpz_poly_srcptr, unsigned long);
double get_alpha_projective (mpz_poly_srcptr f, unsigned long B);
double get_alpha_affine (mpz_poly_srcptr f, unsigned long B);
double get_alpha_affine_p (mpz_poly_srcptr f, unsigned long p, gmp_randstate_ptr rstate);

/* poly info, being called in order */
void print_cadopoly_fg (FILE*, mpz_t*, int, mpz_t*, int, mpz_srcptr);
double print_cadopoly (FILE*, cado_poly_srcptr);
void print_cadopoly_extra (FILE*, cado_poly, int, char**, double);
double print_poly_fg (mpz_poly_srcptr, mpz_t*, mpz_t, int);
long rotate_aux (mpz_t *f, mpz_t b, mpz_t m, long k0, long k, unsigned int t);
void rotate_auxg_z (mpz_t*, const mpz_t, const mpz_t, const mpz_t, unsigned int);
void do_translate_z (mpz_poly_ptr f, mpz_t *g, const mpz_t k);


double cado_poly_fprintf_with_info (FILE *, cado_poly_ptr, const char *, int);
double cado_poly_fprintf_with_info_and_MurphyE (FILE *fp, cado_poly_ptr,
                                                double, double, double, double,
                                                const char *);
double expected_rotation_gain (mpz_poly_srcptr f, mpz_poly_srcptr g);
void expected_growth (rotation_space *r, mpz_poly_srcptr f, mpz_poly_srcptr g,
                      int i, double maxlognorm, double skew);

void set_alpha_bound (unsigned long bound);
unsigned long get_alpha_bound (void);

#ifdef __cplusplus
}
#endif

#endif	/* POLYSELECT_AUXILIARY_H_ */

