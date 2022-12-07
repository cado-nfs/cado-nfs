#ifndef BADIDEALS_HPP_
#define BADIDEALS_HPP_
#include <string>
#include <vector>
#include <iosfwd>
#include <gmp.h>        // for gmp_randstate_t
#include "cxx_mpz.hpp"
#include "mpz_poly.h"

/* Representation of an element of P^1(Z/p^kZ)
 * -------------------------------------------
 *
 * given u with 0 <= u < p^k,  the integer u represents (u:1)
 * given u with p^k <= u < 2*p^k and p|u  the integer p^k+u represents (1:u)
 *
 */

struct badideal {/*{{{*/
    cxx_mpz p;
    /* The r below is a congruence class which is common to all branches
     * specified in the vector below */
    cxx_mpz r;  /* we have 0<=r<p+1 to account for projective ideals */
    int nbad;   /* number of ideals above this (p,r). Always > 1 by
                   definition */
    std::string comments;    /* Used to debug FM */
    struct branch {
        int k;
        cxx_mpz r;  /* we have 0<=r<2*p^k to account for projective ideals */
        std::vector<int> v;
    };
    std::vector<branch> branches;

    std::vector<std::string> sagemath_string; // nbad strings.

    badideal(cxx_mpz const& p, cxx_mpz const& r, int nbad) : p(p), r(r), nbad(nbad) {}
    badideal(std::istream &);

    std::ostream& operator<<(std::ostream& o) const;
    std::istream& operator>>(std::istream& i);

    std::ostream& print_dot_badideals_file(std::ostream & o, int side) const;

    std::ostream& print_dot_badidealinfo_file(std::ostream& o, int side) const;

    static cxx_mpz r_from_rk(cxx_mpz const & p, int k, cxx_mpz const & rk);
};/*}}}*/

std::vector<badideal> badideals_for_polynomial(cxx_mpz_poly const& f, int side);
std::vector<badideal> badideals_for_polynomial(cxx_mpz_poly const& f, int side, gmp_randstate_t state);
std::vector<badideal> badideals_above_p(cxx_mpz_poly const& f, int side, cxx_mpz const& p);
std::vector<badideal> badideals_above_p(cxx_mpz_poly const& f, int side, cxx_mpz const& p, gmp_randstate_t state);

std::string generic_sagemath_string(cxx_mpz_poly const & f, int side, cxx_mpz const & p, cxx_mpz const & r, gmp_randstate_t state);

inline std::istream& operator>>(std::istream& i, badideal & b)
{
    return b.operator>>(i);
}
inline std::ostream& operator<<(std::ostream& o, badideal const & b)
{
    return b.operator<<(o);
}


#endif	/* BADIDEALS_HPP_ */
