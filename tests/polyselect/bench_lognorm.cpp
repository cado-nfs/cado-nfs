#include "cado.h"       // IWYU pragma: keep

#include <cstddef>

#include <algorithm>
#include <vector>
#include <string>
#include <stdexcept>

#include <cstdint>
#include <cmath>

#include "fmt/base.h"

#include "mpz_poly.h"
#include "polyselect_norms.hpp"
#include "gmp_aux.h"
#include "timing.h"
#include "utils_cxx.hpp"
#include "params.hpp"

struct lognorm_tester {
    parameter_with_default<std::vector<std::string>,
        "tests",
        "which tests to run",
        "l2_lognorm,l2_skewness,l2_combined_skewness2,l2_skew_lognorm"> tests;
    parameter_with_default<size_t,
        "unit",
        "number of tests per batch",
        "1000"> unit;
    parameter_with_default<size_t,
        "iter",
        "number of batches to run",
        "1000"> loops;
    parameter<unsigned long,
        "seed",
        "random_seed"> seed;
    parameter_with_default<int,
        "bits",
        "average coeff bits", "60"> bits;
    parameter<int,
        "skewbits",
        "typical skewness bits"> skewbits;
    parameter_with_default<int,
        "deg",
        "degree", "6"> deg;

    static void configure(cxx_param_list & pl)
    {
        decltype(tests)::configure(pl);
        decltype(unit)::configure(pl);
        decltype(loops)::configure(pl);
        decltype(seed)::configure(pl);
        decltype(bits)::configure(pl);
        decltype(skewbits)::configure(pl);
        decltype(deg)::configure(pl);
    }

    explicit lognorm_tester(cxx_param_list & pl)
        : tests(pl)
        , unit(pl)
        , loops(pl)
        , seed(pl)
        , bits(pl)
        , skewbits(pl)
        , deg(pl)
    {
        if (skewbits())
            pl.fail("skewbits != 0 is not yet implemented");
    }

    cxx_gmp_randstate state;
    cxx_mpz_poly f;
    cxx_mpz_poly g;

    template<typename T>
        void do_test(std::string const & name, T f)
        {
            std::vector<uint64_t> res;
            res.reserve(loops);

            double s = 0;
            for(size_t i = 0 ; i < loops ; i++) {
                const uint64_t t0 = microseconds();
                for(size_t j = 0 ; j < unit ; j++) {
                    s += f();
                }
                res.push_back(microseconds() - t0);
            }
            std::ranges::sort(res);

            uint64_t m0 = 0, m1 = 0, m2 = 0;
            for(size_t i = res.size() / 4 ; i < (3 * res.size() / 4) ; i++) {
                const uint64_t c = res[i];
                m0++;
                m1+=c;
                m2+=c*c;
            }
            auto mean = double_ratio(m1, m0);
            auto variance = double_ratio(m2, m0) - mean * mean;
            auto sdev = std::sqrt(variance);

            fmt::print("{} {:g} {} {:g} {:g}\n",
                    name, s, res[res.size() / 2], mean, sdev);
        }


    void run() {
        constexpr auto flags = MPZ_POLY_SIGNED_COEFFICIENTS;
        for(auto const & s : tests()) {
            if (s == "l2_lognorm") {
                do_test(s, [&]() {
                    mpz_poly_set_randomb(f, deg(), state, bits(), flags);
                    return L2_lognorm(f, 1.2345);
                    });
            } else if (s == "l2_skew_lognorm") {
                do_test(s, [&]() {
                    mpz_poly_set_randomb(f, deg(), state, bits(), flags);
                    return L2_skew_lognorm(f);
                });
            } else if (s == "l2_skewness") {
                do_test(s, [&]() {
                    mpz_poly_set_randomb(f, deg(), state, bits(), flags);
                    return L2_skewness(f);
                });
            } else if (s == "l2_combined_skewness2") {
                do_test(s, [&]() {
                    mpz_poly_set_randomb(f, deg(), state, bits(), flags);
                    mpz_poly_set_randomb(g, 1, state, 60, flags);
                    return L2_combined_skewness2(f, g);
                });
            } else {
                throw std::invalid_argument(s);
            }
        }
    }
};


int main(int argc, char const * argv[])
{
    cxx_param_list pl;

    pl.declare_usage_header("Bench for the lognorm computation routines");
    lognorm_tester::configure(pl);
    pl.process_command_line(argc, argv);

    lognorm_tester L(pl);

    if (pl.warn_unused())
        pl.fail("unexpected argument(s)");

    L.run();
};
