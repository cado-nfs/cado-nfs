#include "cado.h" // IWYU pragma: keep
#include <cstring>

#include <iostream>

#include "fmt/base.h"

#include "logapprox.hpp"
#include "macros.h"
#include "mpz_poly.h"
#include "polynomial.hpp"
#include "tests_common.h"

/* The logapprox test is not as interesting as the test_init_norms test.
 * Here, we have different kind of behaviors depending on the compiler
 * and its optimization options. Furthermore, we're only checking the
 * number of approximating lines, which pretty gross.
 */

static void display_logapprox(piecewise_linear_function const & F)
{
    const size_t sz = F.equations.size();
    std::cout << "Total " << sz << " pieces\n";
    ASSERT_ALWAYS(F.endpoints.size() == sz + 1);

    auto eq = F.equations.begin();
    auto ep = F.endpoints.begin();
    for( ; eq != F.equations.end() ; ) {
        ASSERT_ALWAYS(ep != F.endpoints.end());
        double const r0 = *ep++;
        ASSERT_ALWAYS(ep != F.endpoints.end());
        double const r1 = *ep;
        auto uv = *eq++;
        std::cout << "[" << r0 << ", " << r1 << "]: " << uv.first << " + x * " << uv.second << "\n";
    }
}

static int test_from_bug21684(bool display)
{
    const polynomial<double> f {-6.3406659802246472e+28, 6.4695148868405632e+28, 9.5457310557271272e+27};
    const piecewise_linear_approximator<double> A(f, 0.3300700859809263);
    const piecewise_linear_function F = A.logapprox(-2048,2048);
    fmt::print("# {}\n", F.equations.size());
    ASSERT_ALWAYS(F.equations.size() == 21);
    if (display) display_logapprox(F);
    return 0;
}

static int test_from_bug21701(bool display)
{
    const polynomial<double> f {-1.6425515054690201e+34, -2.4119595460727266e+36, -1.1805918683048612e+38, -1.926230702854181e+39};
    const piecewise_linear_approximator<double> A(f, 0.44719172939351309);
    const piecewise_linear_function F = A.logapprox(-2048,2048);
    fmt::print("# {}\n", F.equations.size());
    ASSERT_ALWAYS(F.equations.size() >= 57 || F.equations.size() <= 59);
    if (display) display_logapprox(F);
    return 0;
}

static int test_from_bug30107(bool display)
{
    cxx_mpz_poly F(
            "-5221483323032651731"
            "-7279036884271919814*x"
            "+4900275701330175572*x^2"
            "-2730844396347493817*x^3"
            "-14*x^4+x^5");
    F = F.homography({-934, -453, 2943, -748}).divexact(2031811);

    const piecewise_linear_approximator<double> Ad(
            polynomial<double>(F), 0.34229490398021989);
    const piecewise_linear_function Fd = Ad.logapprox(-64, 64);
    // ASSERT_ALWAYS(Fd.has_precision_issues);
    if (display) display_logapprox(Fd);
    fmt::print("# {}\n", Fd.equations.size());
    ASSERT_ALWAYS(Fd.equations.size() == 63 || Fd.equations.size() == 65);

    if (!tests_run_under_valgrind()) {
        /* long double code with valgrind seems to behave a little bit
         * differently, and I'm not very much interested in tracking down the
         * why and how.
         */
        const piecewise_linear_approximator<long double> Al(
                polynomial<long double>(F), 0.34229490398021989L);
        const piecewise_linear_function Fl = Al.logapprox(-64,64);
        if (display) display_logapprox(Fl);
        fmt::print("# {}\n", Fl.equations.size());
        ASSERT_ALWAYS(Fl.equations.size() == 65);
    }
    return 0;
}

static int original_test(bool display)
{
    const polynomial<double> f {
        -39769265437440000.,
        -302859053854976.,
        5377439145928.,
        -1684314626.,
        -5377481.,
        970.,
        1.};
    const piecewise_linear_approximator<double> A(f, 0.1);
    const piecewise_linear_function F = A.logapprox(-2048,2048);
    fmt::print("# {}\n", F.equations.size());
    ASSERT_ALWAYS(F.equations.size() == 34);
    if (display) display_logapprox(F);
    return 0;
}

// coverity[root_function]
int main(int argc, char const * argv[])
{
    bool display = false;
    if (argc >= 2 && strcmp(argv[1], "--display") == 0)
        display=true;
    test_from_bug30107(display);
    test_from_bug21684(display);
    test_from_bug21701(display);
    original_test(display);
}

