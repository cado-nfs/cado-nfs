#include "cado.h" // IWYU pragma: keep
#include <cstring>
#include <initializer_list>  // for initializer_list
#include <list>              // for list
#include <utility>           // for pair
#include <iostream>
#include "double_poly.h"     // for cxx_double_poly, double_poly_s
#include "logapprox.hpp"
#include "macros.h"

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
    cxx_double_poly f(3);
    for(double const x : {-6.3406659802246472e+28, 6.4695148868405632e+28, 9.5457310557271272e+27}) f->coeff[++f->deg]=x;
    piecewise_linear_approximator A(f, 0.3300700859809263);
    const piecewise_linear_function F = A.logapprox(-2048,2048);
    ASSERT_ALWAYS(F.equations.size() == 21);
    if (display) display_logapprox(F);
    return 0;
}

static int test_from_bug21701(bool display)
{
    cxx_double_poly f(4);
    for(double const x : {-1.6425515054690201e+34, -2.4119595460727266e+36, -1.1805918683048612e+38, -1.926230702854181e+39}) f->coeff[++f->deg]=x;
    piecewise_linear_approximator A(f, 0.44719172939351309);
    const piecewise_linear_function F = A.logapprox(-2048,2048);
    ASSERT_ALWAYS(F.equations.size() == 57);
    if (display) display_logapprox(F);
    return 0;
}

static int test_from_bug30107(bool display)
{
    cxx_double_poly f(6);
    for(double const x: {
            -3.060951772284165e+28,
            -2.0427520002332558e+28,
            6.3458797760215409e+27,
            3.5641203631093495e+27,
            -2.2704812001735614e+26,
            2.1907358313465321e+26})
        f->coeff[++f->deg]=x;
    piecewise_linear_approximator A(f, 0.34229490398021989);
    const piecewise_linear_function F = A.logapprox(-64,64);
    if (display) display_logapprox(F);
    ASSERT_ALWAYS(F.equations.size() == 63);
    return 0;
}

static int original_test(bool display)
{
    cxx_double_poly f(6);
    for(double const x : {-39769265437440000., -302859053854976., 5377439145928., -1684314626., -5377481., 970., 1.}) f->coeff[++f->deg]=x;
    piecewise_linear_approximator A(f, 0.1);
    const piecewise_linear_function F = A.logapprox(-2048,2048);
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
    original_test(display);
    test_from_bug21684(display);
    test_from_bug21701(display);
    test_from_bug30107(display);
}

