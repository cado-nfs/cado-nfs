#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>

#include <utility>
#include <map>

#include <gmp.h>
#include "fmt/base.h"

#include "gmp_aux.h"
#include "cxx_mpfr.hpp"

// coverity[root_function]
int main(int argc, char const * argv[])
{
    cxx_gmp_randstate state;

    if (argc > 1) {
        unsigned long seed;
        seed = strtoul(argv[1], nullptr, 0);
        gmp_randseed_ui(state, seed);
    }

    std::map<unsigned long, cxx_mpfr> v;

    mpfr_prec_t const prec = 256;

    cxx_mpfr m;
    mpfr_set_prec(m, prec);
    mpfr_const_pi(m, MPFR_RNDN);
    mpfr_frac(m, m, MPFR_RNDN);

    /* create 2000 floating point numbers, all between 1/2 and 1 (thus on
     * average 0.75), place them in a map, and sum all of them
     */
    for(int i = 0 ; i < 2000 ; i++) {
        mpfr_mul_ui(m, m, 1 + gmp_urandomm_ui(state, 255), MPFR_RNDN);
        mpfr_div_ui(m, m, 1 + gmp_urandomm_ui(state, 255), MPFR_RNDN);
        mpfr_set_exp(m, 0);
        // fmt::print("{}\n", m);
        unsigned int t = gmp_urandomm_ui(state, 100);
        if (v.find(t) == v.end()) {
            mpfr_set_prec(v[t], prec);
            mpfr_set_zero(v[t], 1);
        }
        mpfr_add(v[t], v[t], m, MPFR_RNDN);
    }

    mpfr_set_zero(m, 1);
    for(auto const & x : v) {
        mpfr_add(m, m, x.second, MPFR_RNDN);
        fmt::print("[{}] {} {}\n", x.first, x.second, m);
    }
}
