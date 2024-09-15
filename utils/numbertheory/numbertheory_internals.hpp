#ifndef UTILS_NUMBERTHEORY_INTERNALS_HPP_
#define UTILS_NUMBERTHEORY_INTERNALS_HPP_

#include <vector>
#include <utility>
#include <string>
#include <ostream>
#include <gmp.h>
#include "fmt/format.h"
#include "mpz_poly.h"
#include "mpz_mat.h"
#include "cxx_mpz.hpp"

namespace numbertheory_internals {
/* Return a basis for a p-maximal order of the number field defined by
 * the polynomial f.
 *
 * The basis of the order is written in lower triangular positive row hnf
 * (magma returns lower triangular centered row hnf, with sometimes negative
 * coefficients).
 */
cxx_mpq_mat p_maximal_order(cxx_mpz_poly const& f, cxx_mpz const& p);

cxx_mpz_mat p_radical_of_order(cxx_mpz_mat const& M, cxx_mpz const& p);

std::vector<std::pair<cxx_mpz_mat, int> > factorization_of_prime(
        cxx_mpq_mat const & B, cxx_mpz_poly const& g,
        cxx_mpz const& p,
        gmp_randstate_ptr state);

cxx_mpz_mat multiplication_table_of_order(cxx_mpq_mat const& O, cxx_mpz_poly const& g);

std::pair<cxx_mpz_mat, cxx_mpz> generate_ideal(cxx_mpq_mat const& O, cxx_mpz_mat const& M, cxx_mpq_mat const& gens);

cxx_mpz_mat valuation_helper_for_ideal(cxx_mpz_mat const& M, cxx_mpz_mat const& I, cxx_mpz const& p);

int valuation_of_ideal_at_prime_ideal(cxx_mpz_mat const& M, cxx_mpz_mat const& I, cxx_mpz_mat const& a, cxx_mpz const& p);
int valuation_of_ideal_at_prime_ideal(cxx_mpz_mat const& M, std::pair<cxx_mpz_mat,cxx_mpz> const& Id, cxx_mpz_mat const& a, int e, cxx_mpz const& p);

int prime_ideal_inertia_degree(cxx_mpz_mat const& I);

std::pair<cxx_mpz, cxx_mpz_mat> prime_ideal_two_element(cxx_mpq_mat const& O, cxx_mpz_poly const& f, cxx_mpz_mat const& M, cxx_mpz_mat const& I);

std::string write_order_element_as_vector(cxx_mpz_mat const& z);

std::string write_element_as_polynomial(cxx_mpq_mat const& theta_q, std::string const& var);

std::vector<cxx_mpz> write_element_as_list_of_integers(cxx_mpq_mat const& theta_q);

cxx_mpz_mat multiply_elements_in_order(cxx_mpz_mat const& M, cxx_mpz_mat const& E, cxx_mpz_mat const& F);

struct ideal_comparator {
    typedef std::pair<cxx_mpz_mat,int> Im_t;
    bool operator()(Im_t const& a, Im_t const& b) const {
        int r = mpz_mat_cmp(a.first, b.first);
        if (r) return r < 0;
        return a.second < b.second;     /* should never happen */
    }
};

struct all_valuations_above_p {/**/
    cxx_mpz_poly f;
    cxx_mpz p;
private:
    cxx_mpq_mat O;
    cxx_mpz_mat M;
    std::vector<std::pair<cxx_mpz_mat, int> > F;
    std::vector<int> inertia;
    std::vector<int> ramification; // valuation of the ideal in the
                              // factorization of the underlying prime
                              // ideal.
    std::pair<cxx_mpz_mat, cxx_mpz> jjinv;
public:
    std::vector<cxx_mpz_mat> helpers;
private:
    std::vector<int> val_base;

public:
    all_valuations_above_p(cxx_mpz_poly const& f, cxx_mpz const& p, gmp_randstate_t state);
    std::vector<int> operator()(std::pair<cxx_mpz_mat, cxx_mpz> const& Id) const;
    std::string sagemath_string(int k, int side);
    std::vector<cxx_mpz> machine_description(int k);
    void print_info(std::ostream& o, int k, cxx_mpz const& r MAYBE_UNUSED, int side) const;
    std::pair<cxx_mpz_mat, cxx_mpz> generate_ideal(cxx_mpq_mat const& gens) const;
    std::pair<cxx_mpz_mat, cxx_mpz> generate_ideal(cxx_mpz_mat const& gens) const;
    std::vector<int> operator()(int k, cxx_mpz const& r) const;
    std::vector<int> multiply_inertia(std::vector<int> const& v) const;



    // getters for e and f
    inline int get_ramification_index(int i) const { return ramification[i]; }
    inline int get_inertia_degree(int i) const { return inertia[i]; }
};/**/

}


#endif	/* UTILS_NUMBERTHEORY_INTERNALS_HPP_ */
