#ifndef CADO_POLYSELECT_AUXILIARY_HPP
#define CADO_POLYSELECT_AUXILIARY_HPP

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

#include <cstdio>      // FILE
#include <gmp.h>
#include "cado_poly.hpp"
#include "mpz_poly.h"
#include "gmp_aux.h"

/* The polynomial selection algorithms that use a linear polynomial will
 * put it on the side given by the following. */
// FIXME: atm, changing these does not work. It should...
#define RAT_SIDE 0
#define ALG_SIDE 1

extern double bound_f, bound_g, area;

/* This goes with one cado_poly object. Each polynomial has its
 * stats, except that we don't compute them for the rational
 * polynomials.
 */
struct cado_poly_stats_one_poly {
    unsigned int nrroots;
    double lognorm;
    double alpha;
    double alpha_proj;
    /* This is the "expected" value of E after root sieving. It takes
     * into account the potential gain from rotation, given the
     * degrees of freedom that are given by the linear polynomial.
     *
     * (in contrast, the "true" E is just lognorm+alpha)
     *
     * This field is 0 if we've computed the stats with
     * cado_poly_compute_stats (a priori final stats), and >0 only if
     * the computation was done with
     * cado_poly_compute_expected_stats.
     */
    double exp_E;
};

struct cado_poly_stats_s {
    struct cado_poly_stats_one_poly (*pols)[1];
    int nb_polys;
};
typedef struct cado_poly_stats_s cado_poly_stats[1];
typedef struct cado_poly_stats_s * cado_poly_stats_ptr;
typedef const struct cado_poly_stats_s * cado_poly_stats_srcptr;


extern void cado_poly_stats_init(cado_poly_stats_ptr spoly, int nb_polys);
extern void cado_poly_stats_clear(cado_poly_stats_ptr spoly);

/* returns the average exp_E for nonlinear polynomials */
extern double cado_poly_compute_expected_stats(cado_poly_stats_ptr, cxx_cado_poly const & cpoly);

/* returns the average of E=lognorm+alpha for nonlinear polynomials */
extern double cado_poly_compute_stats(cado_poly_stats_ptr, cxx_cado_poly const & cpoly);

extern void cado_poly_fprintf_stats(FILE * fp, const char * prefix, cxx_cado_poly const & cpoly, cado_poly_stats_srcptr spoly);

/* The maximum degree supported is MAX_DEGREE, as defined in cado_poly.hpp */

#define NORM_MARGIN 0.2

typedef struct
{
  double kmin, kmax;
} rotation_space;

/* alpha */
double special_valuation (mpz_poly_srcptr f, unsigned long p, mpz_srcptr disc, gmp_randstate_ptr rstate);
double special_valuation_affine (mpz_poly_srcptr f, unsigned long p, mpz_srcptr disc, gmp_randstate_ptr rstate);

/*
 * Those are in polyselect_alpha.h now
double get_alpha (mpz_poly_srcptr, unsigned long);
double get_alpha_projective (mpz_poly_srcptr f, unsigned long B);
double get_alpha_affine (mpz_poly_srcptr f, unsigned long B);
double get_alpha_affine_p (mpz_poly_srcptr f, unsigned long p, gmp_randstate_ptr rstate);
void set_alpha_bound (unsigned long bound);
unsigned long get_alpha_bound ();
*/

/* poly info, being called in order */
void print_cadopoly_fg (FILE*, mpz_poly_srcptr, mpz_poly_srcptr, mpz_srcptr);
double print_cadopoly (FILE*, cxx_cado_poly const &);
void print_cadopoly_extra (FILE*, cxx_cado_poly const &, int, char const **, double);
double print_poly_fg (mpz_poly_srcptr, mpz_poly_srcptr, mpz_srcptr, int);

void cado_poly_set_skewness_if_undefined(cxx_cado_poly & cpoly);
long rotate_aux (mpz_poly_ptr f, mpz_poly_srcptr, long k0, long k, unsigned int t);
void rotate_auxg_z (mpz_poly_ptr f, mpz_poly_srcptr g, mpz_srcptr k, unsigned int t) __attribute__((deprecated));

double expected_rotation_gain (mpz_poly_srcptr f, mpz_poly_srcptr g);

void expected_growth (rotation_space *r, mpz_poly_srcptr f, mpz_poly_srcptr g,
                      int i, double maxlognorm, double skew);


#endif	/* CADO_POLYSELECT_AUXILIARY_HPP */

