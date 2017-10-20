#ifndef LAS_COFACTOR_HPP_
#define LAS_COFACTOR_HPP_

#include "gmp.h"
#include "las-types.hpp"
#include "cxx_mpz.hpp"
#include <array>
#include <vector>

int check_leftover_norm (cxx_mpz const & n, sieve_info const & si, int side);

int factor_both_leftover_norms(
        std::array<cxx_mpz, 2> & norms,
        std::array<std::vector<cxx_mpz>, 2> &,
        sieve_info const &);
#endif
