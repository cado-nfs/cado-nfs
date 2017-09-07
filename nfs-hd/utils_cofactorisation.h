#ifndef UTILS_COFACTORISATION_H
#define UTILS_COFACTORISATION_H

#include "cado.h"
#include "utils.h"
#include "uint64_array.h"
#include "mat_Z.h"
#include "sieving_bound.h"
#include "ideal.h"

//Define an array with all the factors.
typedef struct
{
  mpz_t * factorization;
  unsigned int number;
  unsigned int alloc;
} s_factor_t;

typedef s_factor_t factor_t[1];
typedef s_factor_t * factor_ptr;
typedef const s_factor_t * factor_srcptr;

/*
 * Initialize an array of factors with number, the maximum number of elements.
 *
 * factor: the array of factors.
 * number: the maximum number of elements.
 */
void factor_init(factor_ptr factor, unsigned int alloc);

/*
 * Delete an array of factors.
 *
 * factor: the array of factors.
 */
void factor_clear(factor_ptr factor);

void factor_append(factor_ptr factor, mpz_srcptr z);


unsigned int factor_remove(factor_ptr factor, mpz_srcptr z);

/*
 * TODO: useless.
 * Test if the maximum factor of the array factor is less or equal to B. Return
 *  1 if true, 0 otherwise. If factor is sorted, set sort to 1, 0 otherwise.
 *
 * factor: an array of factors.
 * B: the smoothness bound.
 * sort: 1 if factor is sorted, 0 otherwise.
 */
unsigned int factor_is_smooth(factor_srcptr factor, mpz_t B, unsigned int sort);

/*
 * Return 1 if the factorisation is good, 0 otherwise.
 */
unsigned int factor_assert(factor_srcptr factor, mpz_srcptr z);

/*
 * Print an array of factors.
 *
 * factor: an array of factors.
 */
void factor_fprintf(FILE * file, factor_srcptr factor);

/*
 * Remove all the small factors under a certain bound, and store z_root /
 *  (factors) in z. Return 1 if z_root is entirely factorize.
 */
int brute_force_factorize_ul(factor_ptr factor, mpz_ptr z,
    mpz_srcptr z_root, unsigned long bound);

/*
 * To find the relations.
 *
 * indices: for each number field, positions where the norm is less than the
 *  threshold.
 * number_element: number of elements in the sieving region.
 * lpb: large prime bounds.
 * matrix: MqLLL.
 * f: polynomials that define the number fields.
 * H: sieving bounds.
 * V: number of number fields.
 */
unsigned int find_relations(uint64_array_t * indices, uint64_t number_element,
    unsigned int * lpb, mat_Z_srcptr matrix, const mpz_poly * f,
    sieving_bound_srcptr H, unsigned int V, ideal_spq_srcptr special_q,
    unsigned int q_side, int main, FILE * outstd, unsigned int gal,
    unsigned int gal_version, factor_t * gal_norm_denom, int * nb_curves);

#endif /* UTILS_COFACTORISATION_H */
