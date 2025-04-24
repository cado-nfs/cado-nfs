#ifndef CADO_FACUL_HPP
#define CADO_FACUL_HPP

#include <vector>

#include "facul_strategies.hpp"
#include "cxx_mpz.hpp"

/* we should have FACUL_NOT_SMOOTH < 0, FACUL_MAYBE = 0,
   and FACUL_SMOOTH, FACUL_AUX >= 1 */

enum facul_status {
    FACUL_NOT_SMOOTH = -1,
    FACUL_MAYBE      = 0,
    FACUL_SMOOTH     = 1,
};

#define STATS_LEN 128


struct facul_result {
    facul_status status = FACUL_MAYBE;
    std::vector<cxx_mpz> primes;
    facul_result() = default;
    // NOLINTNEXTLINE(hicpp-explicit-conversions)
    facul_result(facul_status const s) : status(s) {}
};

facul_result facul (cxx_mpz const &, facul_strategy_oneside const &);

std::vector<facul_result>
facul_both (std::vector<cxx_mpz> const & , facul_strategies const &);

#endif /* FACUL_HPP_ */
