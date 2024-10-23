#ifndef SM_UTILS_HPP_
#define SM_UTILS_HPP_

/* Which SMs are used ?
 * There have been several choices during time. We keep track of the old
 * choices, in case there are still parts of the code that make reference
 * to them.
 *
 * Normally, for compatibility with previous data files, we strive to
 * provide access to the old computation conventions. Be aware, however,
 * that this is inherently fragile, and not always tested. The chosen SM
 * mode must be consistent throughout the computation.
 *
 * Legacy choice (internal name SM_MODE_LEGACY_PRE2018, argument
 * --sm-mode legacy)
 *
 *   Compute modulo f(x), with an exponent corresponding to the LCM of
 *   all factors modulo ell. Then take the coefficients of the resulting
 *   polynomials, starting with the one with largest degree, in order to
 *   avoid jokes with the constant coefficient and Galois effects.
 *
 * Transitional: (same)
 *   Do this computation using CRT modulo all factors of f(x), but
 *   keeping the compatibility with the legacy choice.
 *
 * Choice activated in ec940b55a (summer 2018); (internal name
 * SM_MODE_2018, argument --sm-mode 2018. This is also the current
 * default, so that --sm-mode default is an alias to this)
 *   Select a subset of factors of f(x) mod ell so that their degrees sum
 *   up to (just) at least the number of required SMs. See
 *   sm_side_info_init() for the basic algorithm that chooses these
 *   factors.
 *   Then compute modulo each factor f_i(x) mod ell, with the exponent
 *   attached to the degree of this factor. Each coefficient of the
 *   resulting polynomial, taken in increasing degree order, contributes
 *   to a list of deg(f_1)+...+deg(f_k) coordinates. The nsm furst
 *   coordinates of this list are elected as the SM vector.
 *
 * Alternate mode, maybe better (internal name SM_MODE_2019REV, argument
 * --sm-mode 2019rev)
 *   This is almost the same as the previous one, except that the list of
 *   deg(f_1)+...+deg(f_k) coordinates is built by taking the
 *   coefficients of each of the components in reverse degree order.
 *
 * Potential future extension (internal name
 * SM_MODE_PSEUDORANDOM_COMBINATION, argument --sm-mode pseudorandom).
 *   We might consider implementing some sort of pseudo-random
 *   combination of the coefficients. The exact details are not yet
 *   defined. (not implemented).
 */


#define MAX_LEN_RELSET 1024

#include <stdint.h>
#include <stdio.h>      // FILE
#include <gmp.h>        // mpz_t
#include <vector>
#include "mpz_poly.h"
#include "mpz_mat.h"
#include "cxx_mpz.hpp"
#include "cado_poly.h"  // MAX_DEGREE, NB_POLYS_MAX

enum sm_mode {
    SM_MODE_LEGACY_PRE2018 = 1,
    SM_MODE_2018 = 2,
    SM_MODE_2019REV = 3,
    SM_MODE_PSEUDORANDOM_COMBINATION = 4, /* not implemented yet */
};

struct sm_side_info {
    int unit_rank;
    int nsm; /* number of SMs that are going to be computed. By default, is
                equal to unitrank but can be modified by the user. */
    cxx_mpz ell, ell2, ell3;
    cxx_mpz_poly f0;
    cxx_mpz_poly f;       /* monic */
    cxx_mpz exponent;

    struct piece {
        cxx_mpz_poly g;
        int m;
        bool is_used;
        cxx_mpz exponent;
        piece(cxx_mpz_poly const & g, const int & m, bool is_used, int exponent)
            : g(g)
            , m(m)
            , is_used(is_used)
            , exponent(exponent)
        {}
    };

    std::vector<piece> pieces;

    cxx_mpz_mat matrix; // only if sm->mode == SM_MODE_LEGACY_PRE2018

    enum sm_mode mode;

    sm_side_info(mpz_poly_srcptr f0, mpz_srcptr ell, bool handle_small_ell = false);
    void set_mode(const char * mode_string);
    void print(FILE * out) const;

    /* This does the same as compute_sm_straightforward, except that it works piecewise on
     * the different components. It is thus noticeably faster. Results are
     * compatible, as the change of basis is precomputed within the
     * sm_side_info structure.
     */
    void compute_piecewise(mpz_poly_ptr dst, mpz_poly_srcptr u) const;
};

typedef struct {
  mpz_poly * num;
  mpz_poly * denom;
  int nb_polys;
} sm_relset_struct_t;

typedef sm_relset_struct_t sm_relset_t[1];
typedef sm_relset_struct_t * sm_relset_ptr;
typedef const sm_relset_struct_t * sm_relset_srcptr;

/* For MNFS, just "a,b" is not enough. We want to know the sides that are
 * related to the pair.
 * This is only used by thread_sm
 */
struct pair_and_sides {
    /* mpz_poly because we might think of using higher degree */
    cxx_mpz_poly ab;
    unsigned int active_sides[2];
    pair_and_sides(mpz_srcptr a, mpz_srcptr b, int s0, int s1) {
        mpz_poly_set_mpz_ab(ab, a, b);
        active_sides[0] = s0;
        active_sides[1] = s1;
    }
    pair_and_sides(int64_t a, uint64_t b, int s0, int s1) {
        mpz_poly_set_ab(ab, a, b);
        active_sides[0] = s0;
        active_sides[1] = s1;
    }
    pair_and_sides() = default;
    pair_and_sides(pair_and_sides const &) = default;
    pair_and_sides(pair_and_sides &&) = default;
    pair_and_sides& operator=(pair_and_sides const &) = default;
    pair_and_sides& operator=(pair_and_sides &&) = default;
};

#ifdef __cplusplus
extern "C" {
#endif


void sm_relset_init (sm_relset_t r, const mpz_poly_srcptr * F, int nb_polys);
void sm_relset_clear (sm_relset_t r);
void sm_relset_copy (sm_relset_t r, sm_relset_srcptr s);

// (a,b) -> a - b*x
void mpz_poly_init_set_ab (mpz_poly_ptr, int64_t, uint64_t);

// array of rows and exponents -> rational fractions corresponding to the
// combination of these rows, stored in a sm_relset structure.
// If F[0] or F[1] is NULL, then no computation is done on the
// corresponding side.
void sm_build_one_relset (sm_relset_ptr rel,
                    const uint64_t *r, const int64_t *e, int len,
                    std::vector<pair_and_sides> const & ps,
                    const mpz_poly_srcptr * F, int nb_polys,
		    mpz_srcptr ell2);

// Print coeffs of the SM polynomial
void print_sm (FILE *f, sm_side_info const & S, mpz_poly_srcptr SM);
// same, with a delimiter
void print_sm2 (FILE *f, sm_side_info const & S, mpz_poly_srcptr SM, const char * delim);


#ifdef __cplusplus
}
#endif


#endif /* SM_UTILS_HPP_ */
