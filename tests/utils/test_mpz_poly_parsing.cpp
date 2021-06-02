#include "cado.h" // IWYU pragma: keep
// IWYU pragma: no_include <ext/alloc_traits.h>
#include <cstdlib>
#include <iostream>
#include <sstream> // istringstream // IWYU pragma: keep
#include <vector>
#include <utility>
#include <memory>              // for allocator_traits<>::value_type
#include <string>              // for string, basic_string
#include <gmp.h>               // for mpz_cmp, mpz_set_str, mpz_t
#include "cxx_mpz.hpp"
#include "mpz_poly.h"
#include "macros.h"

/* a shorthand so that we can use user-defined literals */
cxx_mpz operator "" _mpz (const char* str, size_t)
{
    cxx_mpz res;
    mpz_set_str(res, str, 0);
    return res;
}

std::vector<std::pair<std::string, std::vector<cxx_mpz>>> examples {
    { "X^128+(X+1)^2+(X^3+X^2+1)*(X+1)-X^2^7", {2UL, 3UL, 2UL, 2UL, 1UL}},
    { "X+1-(X+1)", {} },
    { "z+1", {1UL,1UL} },
    /* see whether we correctly parse the c85 polynomial */
    { "960*x^4+85585660*x^3+405578084588*x^2+4213006218645262637*x-2287975327041639106629845",
     {"-2287975327041639106629845"_mpz, "4213006218645262637"_mpz, "405578084588"_mpz, 85585660UL, 960UL}},
};

std::vector<std::string> expected_failures {
    "x+y",
    "x+y-y",
};

// coverity[root_function]
int main()
{
    for(auto const & example : examples) {
        cxx_mpz_poly f;
        if (!(std::istringstream(example.first) >> f)) {
            std::cerr << "cannot parse polynomial\n";
            exit(EXIT_FAILURE);
        }
        std::cout << f << std::endl;
        ASSERT_ALWAYS((size_t)(f->deg+1) == example.second.size());
        for(size_t i = 0 ; i < example.second.size() ; ++i)
            ASSERT_ALWAYS(mpz_cmp(f->coeff[i], example.second[i]) == 0);
    }
    for(auto const & example : expected_failures) {
        cxx_mpz_poly f;
        if ((std::istringstream(example) >> f)) {
            std::cerr << "unexpected success while parsing bad polynomial\n";
            exit(EXIT_FAILURE);
        }
    }
}
