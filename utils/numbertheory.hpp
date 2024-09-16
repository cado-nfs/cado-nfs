#ifndef NUMBERTHEORY_HPP_
#define NUMBERTHEORY_HPP_

#include "numbertheory/number_field.hpp"
#include "numbertheory/number_field_element.hpp"
#include "numbertheory/number_field_fractional_ideal.hpp"
#include "numbertheory/number_field_order.hpp"
#include "numbertheory/number_field_order_element.hpp"
#include "numbertheory/number_field_prime_ideal.hpp"

/* TODO: do away with it! */
#include "numbertheory/numbertheory_internals.hpp"

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


#endif	/* NUMBERTHEORY_HPP_ */
