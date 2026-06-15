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
#include "polynomial.hpp"
#include "double_poly.h"
#include "gmp_aux.h"
#include "timing.h"
#include "utils_cxx.hpp"
#include "params.hpp"

struct double_poly_roots_tester {
    parameter_with_default<std::vector<std::string>,
        "tests",
        "which tests to run",
        "eval,legacy_eval,positive_roots,legacy_positive_roots"> tests;
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

    explicit double_poly_roots_tester(cxx_param_list & pl)
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
            gmp_randseed_ui(state, seed());
            if (s == "positive_roots") {
                do_test(s, [&]() {
                    mpz_poly_set_randomb(f, deg(), state, bits(), flags);
                    auto c = polynomial<double>(f).positive_roots();
                    double s = 0;
                    for(auto x : c) s += x;
                    return s;
                    });
            } else if (s == "legacy_positive_roots") {
                do_test(s, [&]() {
                    double_poly p;
                    double_poly_init(p, deg());
                    std::vector<double> roots(deg(), 0);
                    mpz_poly_set_randomb(f, deg(), state, bits(), flags);
                    double_poly_set_mpz_poly(p, f);
                    double B = double_poly_bound_roots(p);
                    unsigned int nroots = double_poly_compute_roots(roots.data(), p, B);
                    double s = 0;
                    for(unsigned int i = 0 ; i < nroots ; i++) 
                        s += roots[i];
                    double_poly_clear(p);
                    return s;
                    });
            } else if (s == "eval") {
                do_test(s, [&]() {
                    mpz_poly_set_randomb(f, deg(), state, bits(), flags);
                    return polynomial<double>(f)(1.234);
                    });
            } else if (s == "legacy_eval") {
                do_test(s, [&]() {
                    double_poly p;
                    double_poly_init(p, deg());
                    mpz_poly_set_randomb(f, deg(), state, bits(), flags);
                    double_poly_set_mpz_poly(p, f);
                    double s = double_poly_eval(p, 1.234);
                    double_poly_clear(p);
                    return s;
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
    double_poly_roots_tester::configure(pl);
    pl.process_command_line(argc, argv);

    double_poly_roots_tester D(pl);

    if (pl.warn_unused())
        pl.fail("unexpected argument(s)");

    D.run();
};
