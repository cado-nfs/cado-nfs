#ifndef FACUL_HPP_
#define FACUL_HPP_

#include <array>
#include <vector>
#include <map>
#include "facul_strategies.hpp"

struct cxx_mpz;

/* we should have FACUL_NOT_SMOOTH < 0, FACUL_MAYBE = 0,
   and FACUL_SMOOTH, FACUL_AUX >= 1 */
#define FACUL_NOT_SMOOTH (-1)
#define FACUL_MAYBE (0)
#define FACUL_SMOOTH (1)
#define FACUL_AUX (2)

#define STATS_LEN 128

std::array<int,2>
facul_both (std::array<std::vector<cxx_mpz>, 2>&,
            std::array<cxx_mpz, 2> & ,
	    facul_strategies const &,
            int *);

#endif /* FACUL_HPP_ */
