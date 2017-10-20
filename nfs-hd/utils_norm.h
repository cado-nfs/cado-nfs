#ifndef UTILS_NORM_H
#define UTILS_NORM_H

#include "utils.h"
#include "mpz_poly.h"
#include "array.h"
#include "sieving_bound.h"
#include "mat_Z.h"
#include "ideal.h"

#ifdef __cplusplus
extern "C" {
#endif


/*
 * Compute the norm of a in the number field defined by f.
 *
 * res: the result.
 * f: a polynomial.
 * a: the polynomial for which we want to compute the norm.
 */
void norm_poly(mpz_ptr res, mpz_poly_srcptr f, mpz_poly_srcptr a);

#ifdef ASSERT_NORM
/*
 * To print if a norm seems suspicious.
 */
void assert_norm(array_srcptr array, sieving_bound_srcptr H, mpz_poly_srcptr f,
    mat_Z_srcptr matrix, int special_q, double log2_base,double spq_log);
#endif // ASSERT_NORM

#ifndef OLD_NORM
/*
 * Init the norm for a special-q (q, g).
 *
 * array: the array in which the norms are initialized.
 * H: the sieving bound which gives the sieving region.
 * matrix: the MqLLL matrix.
 * f: the f which defines the side.
 * spq: the special-q.
 * special_q: 0 if there is no special-q in this side, else 1.
 */
void init_norm(array_ptr array, unsigned int * norm_max,
    FILE * file, sieving_bound_srcptr H, mat_Z_srcptr matrix, mpz_poly_srcptr f,
    double spq_log, int special_q, double log2_base);

#else // OLD_NORM

/*
 * Precompute some part of the computation of the bound on the resultant
 *  of f and ha. We precompute (deg(f) + 1) ^ (i/2) * (i + 1) ^ (deg(f) / 2) *
 *  infinity_norm(f), with i the degree of ha.
 *
 * pre_compute: an array in which we store the precompute elements. Need to be
 *  alloc for #(dimension of the lattice) coefficients.
 * f: a polynomial that defines the side.
 * d: degree of ha.
 */
void pre_computation(double * pre_compute, mpz_poly_srcptr f, unsigned int d);

/*
 * Init the norm for a special-q (q, g).
 *
 * array: the array in which the norms are initialized.
 * pre_compute: the array of precomputed value obtained with pre_computation.
 * H: the sieving bound which gives the sieving region.
 * matrix: the MqLLL matrix.
 * f: the f which defines the side.
 * spq: the special-q.
 * special_q: 0 if there is no special-q in this side, else 1.
 */
void init_norm(array_ptr array, double * pre_compute,
    sieving_bound_srcptr H, mat_Z_srcptr matrix, mpz_poly_srcptr f,
    ideal_spq_srcptr spq, int special_q);
#endif // OLD_NORM

#ifdef __cplusplus
}
#endif
#endif /* UTILS_NORM_H */
