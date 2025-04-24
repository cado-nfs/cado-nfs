#ifndef CADO_UTILS_NUMBERTHEORY_ALL_VALUATIONS_ABOVE_P_HPP
#define CADO_UTILS_NUMBERTHEORY_ALL_VALUATIONS_ABOVE_P_HPP

#include "gmp_aux.h"
#include "numbertheory/number_field.hpp"
#include "numbertheory/number_field_order.hpp"
#include "numbertheory/number_field_fractional_ideal.hpp"
#include "numbertheory/number_field_prime_ideal.hpp"

namespace numbertheory_internals {

struct all_valuations_above_p {
    cxx_mpz_poly f;
    cxx_mpz p;
private:
    number_field K;
    number_field_order O;
    
    /* factorization of p */
    std::vector<number_field_prime_ideal> F;

    /* inertia degrees. overkill ? */
    std::vector<int> inertia;

    number_field_fractional_ideal jjinv;

private:
    std::vector<int> val_base;

public:
    all_valuations_above_p(cxx_mpz_poly const&, cxx_mpz const& p, cxx_gmp_randstate &);
    void bless_side(int side);
    std::vector<int> operator()(number_field_fractional_ideal const &I) const;
    std::string sagemath_string(int k, int side);
    std::string machine_description(int k);
    void print_info(std::ostream& o, int k, cxx_mpz const& r MAYBE_UNUSED, int side) const;
    // std::pair<cxx_mpz_mat, cxx_mpz> generate_ideal(cxx_mpq_mat const& gens) const;
    // std::pair<cxx_mpz_mat, cxx_mpz> generate_ideal(cxx_mpz_mat const& gens) const;
    std::vector<int> operator()(int k, cxx_mpz const& r) const;
    std::vector<int> multiply_inertia(std::vector<int> const& v) const;

    // getters for e and f
    inline int get_ramification_index(int i) const { return F[i].ramification_index(); }
    inline int get_inertia_degree(int i) const { return F[i].inertia_degree(); }
};

}

#endif	/* UTILS_NUMBERTHEORY_ALL_VALUATIONS_ABOVE_P_HPP_ */
