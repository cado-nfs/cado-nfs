#include "cado.h" // IWYU pragma: keep

#include <vector>

#include "fmt/base.h"
#include <gmp.h>

#include "cxx_mpz.hpp"
#include "ecm/facul.hpp"
#include "ecm/facul_strategies.hpp"
#include "macros.h"

int main()
{

    const std::vector<unsigned long> lim { 1 << 23, 1 << 21 };
    const std::vector<unsigned int> lpb { 33, 35 };
    const std::vector<unsigned int> mfb { 66, 105 };

    const facul_strategies strat(lim, lpb, mfb,
            {-1, -1}, true, false);

    std::vector<cxx_mpz> n {
        "266630219"_mpz, "8435925448634188287568341721097"_mpz
    };

    // int res = factor_both_leftover_norms(n, factors, lim, strat);
    auto ret = facul_both(n, strat);

    /* This replicates what we have in factor_both_leftover_norms */
    for(int i = 0 ; i < 2 ; i++) {
        fmt::print("n={} --> {}\n", n[i], int(ret[i].status));
        cxx_mpz z = 1;
        cxx_mpz N = n[i];
        for(auto const & p : ret[i].primes) {
            fmt::print("\t{}\n", p);
            z *= p;
            mpz_divexact(N, N, p);
            /* factor_both_leftover_norms claims that repeated factors
             * shouldn't be a problem, but very clearly they _are_ going
             * to be a problem if p only divides once here and yet
             * appears twice in the list! */
        }
        if (N > 1) {
            fmt::print("\t{}\n", N);
            z *= N;
        }
        fmt::print("\tproduct={}\n", z);

        ASSERT_ALWAYS(n[i] == z);
    }
}
