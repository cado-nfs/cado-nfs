#include "cado.h" // IWYU pragma: keep

#include <tuple>
#include <vector>

#include "fmt/base.h"
#include "fmt/ranges.h"

#include "getprime.h"
#include "macros.h"
#include "misc.h"
#include "tests_common.h"

static bool
test_one_with_ref(cxx_mpz const & N, unsigned long B,
                  std::vector<std::pair<cxx_mpz, int>> const & factors_ref,
                  cxx_mpz const & cf_ref)
{
    cxx_mpz cf;
    auto factors = trial_division(N, B, cf);
    bool bf = factors == factors_ref;
    bool bc = cf == cf_ref;

    if (!bf) {
        fmt::print(stderr, "ERROR: for trial_division({}, {}), expected "
                           "factors {}, got {}\n", N, B, factors_ref, factors);
    }
    if (!bc) {
        fmt::print(stderr, "ERROR: for trial_division({}, {}), expected "
                           "cofactor {}, got {}\n", N, B, cf_ref, cf);
    }

    return bf & bc;
}

static bool
tests_static()
{
    cxx_mpz N, cf;
    bool ok = true;

    ok &= test_one_with_ref(42U, 10, {{2, 1}, {3, 1}, {7, 1}}, 1U);
    ok &= test_one_with_ref(-42, 10, {{2, 1}, {3, 1}, {7, 1}}, -1);

    ok &= test_one_with_ref(294U, 10, {{2, 1}, {3, 1}, {7, 2}}, 1U);
    ok &= test_one_with_ref(-294, 10, {{2, 1}, {3, 1}, {7, 2}}, -1);

    return ok;
}

static bool
test_one_without_ref(cxx_mpz const & N, unsigned long B, cxx_mpz const & prod)
{
    if (!tests_common_get_quiet()) {
        fmt::print("# trial_division({}, {})\n", N, B);
    }

    cxx_mpz cf, n, g;
    auto factors = trial_division(N, B, cf);

    n = cf;
    for (auto const & f: factors) {
        for (int i = 0; i < f.second; ++i) {
            mpz_mul(n, n, f.first);
        }
    }

    mpz_gcd(g, cf, prod);

    bool be = N == n;
    bool bs = mpz_sgn(N) == mpz_sgn(cf);
    bool bg = g == 1U;

    if (!be) {
        fmt::print(stderr, "ERROR: for trial_division({}, {}), product of "
                           "factors and cofactor is {}\n", N, B, n);
    }
    if (!bs) {
        fmt::print(stderr, "ERROR: for trial_division({}, {}), cofactor is "
                           "not the same as sign of input integer ({} instead "
                           "of {})\n", N, B, mpz_sgn(cf), mpz_sgn(N));
    }
    if (!be) {
        fmt::print(stderr, "ERROR: for trial_division({}, {}), cofactor is "
                           "still divisible by prime(s) below the bound: "
                           "gcd(cf, prod(pi, pi prime < B)={}\n", N, B, g);
    }

    bool ba = true;
    for (auto const & f: factors) {
        if (mpz_sgn(f.first) <= 0) {
            fmt::print(stderr, "ERROR: for trial_division({}, {}), factor {} "
                               "is <= 0\n", N, B, f.first);
            ba = false;
        } else if (f.first >= B) {
            fmt::print(stderr, "ERROR: for trial_division({}, {}), factor {} "
                               "is >= B\n", N, B, f.first);
            ba = false;
        }
    }

    return be & bs & bg & ba;
}

static bool
tests_dynamic()
{
    bool ok = true;

    std::vector<unsigned int> Bvalues {10, 1000, 1U << 16U };
    std::vector<unsigned int> Nsizes {10, 100, 500, 1000 };

    unsigned long niter = 50; /* default number of iterations */
    tests_common_get_iter(&niter);
    fmt::print("tests_dynamic: will perform {} iteration(s)\n", niter);

    cxx_mpz N, cf, prod;

    for (unsigned int B: Bvalues) {
        /* compute the product of all primes below B */
        prime_info pinf;
        prime_info_init (pinf);
        prod = 1U;
        for (unsigned long p = 2; p < B; p = getprime_mt(pinf)) {
            mpz_mul_ui(prod, prod, p);
        }
        prime_info_clear (pinf); /* free the tables */

        for (unsigned long nbits: Nsizes) {
            for(unsigned long i = 0; i < niter; ++i) {
              tests_common_urandomb(N, nbits);
              test_one_without_ref(N, B, prod);
              mpz_neg(N, N);
              test_one_without_ref(N, B, prod);
            }
        }
    }

    return ok;
}

int main(int argc, char const * argv[])
{
    tests_common_cmdline(&argc, &argv, PARSE_SEED | PARSE_ITER | PARSE_QUIET);
    bool ok = tests_static();
    ok &= tests_dynamic();

    tests_common_clear();
    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
