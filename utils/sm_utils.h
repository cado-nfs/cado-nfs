#ifndef SM_UTILS_H_
#define SM_UTILS_H_

/* Which SMs are used ?
 * There have been several choices during time. We keep track of the old
 * choices, in case there are still parts of the code that make reference
 * to them.
 *
 * Legacy choice:
 *   Compute modulo f(x), with an exponent corresponding to the LCM of
 *   all factors modulo ell. Then take the coefficients of the resulting
 *   polynomials, starting with the one with largest degree, in order to
 *   avoid jokes with the constant coefficient and Galois effects.
 *
 * Transitional:
 *   Do this computation using CRT modulo all factors of f(x), but
 *   keeping the compatibility with the legacy choice.
 *
 * Current choice:
 *   Select a subset of factors of f(x) mod ell so that their degrees sum
 *   up to (just) at least the number of required SMs. See
 *   sm_side_info_init() for the basic algorithm that chooses these
 *   factors.
 *   Then compute modulo each factor f_i(x) mod ell, with the exponent
 *   attached to the degree of this factor. Each coefficient of the
 *   resulting polynomial contributes directly to the output SMs (no
 *   CRT).
 *
 */


#define MAX_LEN_RELSET 1024

#include "mpz_poly.h"

struct sm_side_info_s {
    int unit_rank;
    int nsm; /* number of SMs that are going to be computed. By default, is
                equal to unitrank but can be modified by the user. */
    mpz_t ell;
    mpz_t ell2;
    mpz_poly_srcptr f0;
    mpz_poly f;       /* monic */
    mpz_poly_factor_list fac;
    int is_factor_used[MAX_DEGREE];
    mpz_t exponent;
    mpz_t * exponents;
};
typedef struct sm_side_info_s sm_side_info[1];
typedef struct sm_side_info_s * sm_side_info_ptr;
typedef const struct sm_side_info_s * sm_side_info_srcptr;

typedef struct {
  mpz_poly num[NB_POLYS_MAX];
  mpz_poly denom[NB_POLYS_MAX];
  int nb_polys;
} sm_relset_struct_t;

typedef sm_relset_struct_t sm_relset_t[1];
typedef sm_relset_struct_t * sm_relset_ptr;
typedef const sm_relset_struct_t * sm_relset_srcptr;

#ifdef __cplusplus
extern "C" {
#endif

void sm_side_info_init(sm_side_info_ptr sm, mpz_poly_srcptr f0, mpz_srcptr ell);
void sm_side_info_clear(sm_side_info_ptr sm);
void sm_side_info_print(FILE * out, sm_side_info_srcptr sm);

void sm_relset_init (sm_relset_t r, int *d, int nb_polys);
void sm_relset_clear (sm_relset_t r, int nb_polys);
void sm_relset_copy (sm_relset_t r, sm_relset_srcptr s);

// (a,b) -> a - b*x
void mpz_poly_init_set_ab (mpz_poly_ptr, int64_t, uint64_t);

// array of rows and exponents -> rational fractions corresponding to the
// combination of these rows, stored in a sm_relset structure.
// If F[0] or F[1] is NULL, then no computation is done on the
// corresponding side.
void sm_build_one_relset (sm_relset_ptr rel, uint64_t *r, int64_t *e, int len,
			  mpz_poly * abpolys, mpz_poly_ptr *F, int nb_polys, 
			  const mpz_t ell2);

// Taking a polynomial modulo F as input, compute the corresponding SM
// as a polynomial.
void compute_sm_straightforward(mpz_poly_ptr dst, mpz_poly_srcptr u, sm_side_info_srcptr sm);

// Print coeffs of the SM polynomial
void print_sm (FILE *, mpz_poly, int, int);
// same, with a delimiter
void print_sm2 (FILE *f, mpz_poly SM, int nSM, int d, const char * delim);

/* This does the same as compute_sm_straightforward, except that it works piecewise on
 * the different components. It is thus noticeably faster. Results are
 * compatible, as the change of basis is precomputed within the
 * sm_side_info structure.
 */
void compute_sm_piecewise(mpz_poly_ptr dst, mpz_poly_srcptr u, sm_side_info_srcptr sm);

#ifdef __cplusplus
}
#endif


#endif /* SM_UTILS_H_ */
