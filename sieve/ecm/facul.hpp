#ifndef FACUL_HPP_
#define FACUL_HPP_

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

std::vector<int>
facul_both (std::vector<std::vector<cxx_mpz>>&,
            std::vector<cxx_mpz> & ,
	    facul_strategies const &,
            int *);

#endif /* FACUL_HPP_ */
