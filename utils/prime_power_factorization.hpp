#ifndef UTILS_PRIME_POWER_FACTORIZATION_HPP_
#define UTILS_PRIME_POWER_FACTORIZATION_HPP_

#include <map>

#include "cxx_mpz.hpp"

namespace cado {

struct prime_power : public std::pair<cxx_mpz, int> {};
struct prime_power_factorization : public std::map<cxx_mpz, int> {};

} /* namespace cado */

#endif	/* UTILS_PRIME_POWER_FACTORIZATION_HPP_ */
