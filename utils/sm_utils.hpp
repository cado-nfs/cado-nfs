#ifndef CADO_SM_UTILS_HPP
#define CADO_SM_UTILS_HPP

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

#include <cstdint>
#include <cstdio>      // FILE

#include <vector>
#include <array>

#include <gmp.h>

#include "mpz_poly.h"
#include "mpz_mat.h"
#include "cxx_mpz.hpp"
#include "runtime_numeric_cast.hpp"
#include "numbertheory.hpp"

enum sm_mode {
    SM_MODE_LEGACY_PRE2018 = 1,
    SM_MODE_2018 = 2,
    SM_MODE_2019REV = 3,
    SM_MODE_PSEUDORANDOM_COMBINATION = 4, /* not implemented yet */
};

struct sm_side_info {
    cxx_mpz ell, ell2, ell3;
    cxx_mpz_poly f0;
    cxx_mpz_poly f;       /* monic */
    cxx_mpz exponent;

#if 0
    all_valuations_above_p nt;

    struct one_ideal {
        cxx_mpz_mat H;  /* HNF basis with respect to the order */
        cxx_mpz_mat a;  /* valuation helper with respect to the order */
        cxx_mpq_mat uniformizer;        /* as a polynomial in alpha */
        one_ideal(cxx_mpz_mat const & H, cxx_mpz_mat const & a, cxx_mpq_mat const & u)
            : H(H), a(a), uniformizer(u)
        {}
        one_ideal(one_ideal const &) = default;
        one_ideal(one_ideal &&) = default;
        one_ideal& operator=(one_ideal const &) = default;
        one_ideal& operator=(one_ideal &&) = default;
    };
    std::vector<one_ideal> ideals;
#endif

    number_field K;

    /* ell-maximal order O, and O-ideals above ell. M is the
     * multiplication table of O. We don't accept the case where ell
     * ramifies, so we don't need to store the ramification index. On the
     * other hand, a uniformizing element is stored for each ideal, with
     * respect to the polynomial basis.
     */

    number_field_order O;

    int unit_rank;
    int nsm; /* number of SMs that are going to be computed. By default, is
                equal to unitrank but can be modified by the user. */
    struct piece {
        cxx_mpz_poly g;
        int m;
        bool is_used;
        cxx_mpz exponent;
        number_field_prime_ideal fkp;   /* prime ideal (stores a valuation helper) */
        piece(cxx_mpz_poly g,
                int m,
                bool is_used,
                cxx_mpz exponent,
                number_field_prime_ideal fkp)
            : g(std::move(g))
            , m(m)
            , is_used(is_used)
            , exponent(std::move(exponent))
            , fkp(std::move(fkp))
        {}

        private:
        cxx_mpz_poly sm0(number_field_element const &, sm_side_info const &) const;
        mutable std::unique_ptr<cxx_mpz_poly> cached_sm_of_valuation_helper;   // SM(a)
        cxx_mpz_poly const & sm_of_valuation_helper(sm_side_info const &) const;
        public:
        cxx_mpz_poly sm(number_field_element, sm_side_info const &) const;
    };

    std::vector<piece> pieces;

    cxx_mpz_mat matrix; // only if sm->mode == SM_MODE_LEGACY_PRE2018

    enum sm_mode mode = SM_MODE_2018;

    sm_side_info(mpz_poly_srcptr f0, mpz_srcptr ell, bool handle_small_ell = false);
    void set_mode(const char * mode_string);
    void print(FILE * out) const;

    /* This does the same as compute_sm_straightforward, except that it works piecewise on
     * the different components. It is thus noticeably faster. Results are
     * compatible, as the change of basis is precomputed within the
     * sm_side_info structure.
     */
    void compute_piecewise(cxx_mpz_poly & dst, cxx_mpz_poly const & u) const;
};

struct sm_relset {
    std::vector<cxx_mpz_poly> num, denom;
    explicit sm_relset(size_t nb_polys)
        : num(nb_polys)
        , denom(nb_polys)
    {
    }
    sm_relset() = default;
    int nb_polys() const { return runtime_numeric_cast<int>(num.size()); }
};

/* For MNFS, just "a,b" is not enough. We want to know the sides that are
 * related to the pair.
 * This is only used by thread_sm
 */
struct pair_and_sides {
    /* mpz_poly because we might think of using higher degree */
    cxx_mpz_poly ab;
    std::array<int, 2> active_sides;
    pair_and_sides(mpz_srcptr a, mpz_srcptr b, int s0, int s1)
        : active_sides { s0, s1 }
    {
        mpz_poly_set_mpz_ab(ab, a, b);
    }
    pair_and_sides(int64_t a, uint64_t b, int s0, int s1)
        : active_sides { s0, s1 }
    {
        mpz_poly_set_ab(ab, a, b);
    }
    pair_and_sides() = default;
};

// array of rows and exponents -> rational fractions corresponding to the
// combination of these rows, stored in a sm_relset structure.
// If F[0] or F[1] is NULL, then no computation is done on the
// corresponding side.
sm_relset sm_build_one_relset (
                    const uint64_t *r, const int64_t *e, size_t len,
                    std::vector<pair_and_sides> const & ps,
                    std::vector<mpz_poly_srcptr> const & F,
		    mpz_srcptr ell2);

// Print coeffs of the SM polynomial
void print_sm (FILE *f, sm_side_info const & S, mpz_poly_srcptr SM);
// same, with a delimiter
void print_sm2 (FILE *f, sm_side_info const & S, mpz_poly_srcptr SM, const char * delim);


#endif /* CADO_SM_UTILS_HPP */
