#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>

#include <utility>
#include <map>

#include <gmp.h>
#include <mpfr.h>
#include <mpc.h>
#include "fmt/base.h"

#include "gmp_aux.h"
#include "cxx_mpc.hpp"

// coverity[root_function]
int main(int argc, char const * argv[])
{
    cxx_gmp_randstate state;

    if (argc > 1) {
        unsigned long seed;
        seed = strtoul(argv[1], nullptr, 0);
        gmp_randseed_ui(state, seed);
    }

    std::map<unsigned long, cxx_mpc> v;

    mpfr_prec_t const prec = 256;

    cxx_mpc m;

    mpc_set_prec(m, prec);
    mpc_rootofunity(m, 6, 1, MPC_RNDNN);
    mpfr_frac(mpc_realref(m), mpc_realref(m), MPFR_RNDN);
    mpfr_frac(mpc_imagref(m), mpc_imagref(m), MPFR_RNDN);

    /* create 2000 floating point numbers, all between 1/2 and 1 (thus on
     * average 0.75), place them in a map, and sum all of them
     */
    for(int i = 0 ; i < 2000 ; i++) {
        m *= 1 + gmp_urandomm_ui(state, 255);
        m /= 1 + gmp_urandomm_ui(state, 255);
        cxx_mpc u;
        mpc_set_prec(u, prec);
        mpc_rootofunity(u, 2 + gmp_urandomm_ui(state, 14), gmp_urandomm_ui(state, 16), MPC_RNDNN);
        m *= u;
        mpfr_set_exp(mpc_realref(m), 0);
        mpfr_set_exp(mpc_imagref(m), 0);
        // fmt::print("{}\n", m);
        unsigned int t = gmp_urandomm_ui(state, 100);
        if (v.find(t) == v.end()) {
            mpc_set_prec(v[t], prec);
            mpc_set_ui(v[t], 0, MPC_RNDNN);
        }
        v[t] += m;
    }

    mpc_set_ui(m, 0, MPC_RNDNN);
    for(auto const & x : v) {
        m += x.second;
        fmt::print("[{}] {} {}\n", x.first, x.second, m);
    }
}
