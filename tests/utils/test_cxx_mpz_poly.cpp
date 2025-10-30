#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <algorithm>

#include <gmp.h>
#include "fmt/base.h"

#include "gmp_aux.h"
#include "mpz_poly.h"
#include "macros.h"

// coverity[root_function]
int main(int argc, char const * argv[])
{
    cxx_gmp_randstate state;

    if (argc > 1) {
        char * p;
        const unsigned long seed = strtoul(argv[1], &p, 0);
        ASSERT_ALWAYS(*p == '\0');
        gmp_randseed_ui(state, seed);
    }

    std::vector<cxx_mpz_poly> v;

    for(int i = 0 ; i < 10 ; i++) {
        cxx_mpz_poly x;
        int const jmax = gmp_urandomm_ui(state, 16);
        mpz_poly_set_randomb(x, jmax, state, 16,
                MPZ_POLY_SIGNED_COEFFICIENTS |
                MPZ_POLY_RRANDOM |
                MPZ_POLY_DEGREE_EXACT);

        v.push_back(x);
    }
    sort(v.begin(), v.end());

    for(auto const & p : v)
        fmt::print("{}\n", p);
}

