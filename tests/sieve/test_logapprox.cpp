#include "cado.h" // IWYU pragma: keep
#include <cstring>
#include <initializer_list>  // for initializer_list
#include <list>              // for list
#include <utility>           // for pair
#include <iostream>
#include "double_poly.h"     // for cxx_double_poly, double_poly_s
#include "logapprox.hpp"
#include "macros.h"

void display_logapprox(piecewise_linear_function & F)
{
    std::cout << "Total " << F.equations.size() << " pieces\n";
    for( ; !F.equations.empty() ; ) {
        double r0 = F.endpoints.front();
        F.endpoints.pop_front();
        double r1 = F.endpoints.front();
        std::pair<double,double> uv = F.equations.front();
        F.equations.pop_front();
        std::cout << "[" << r0 << ", " << r1 << "]: " << uv.first << " + x * " << uv.second << "\n";
    }
}

int test_from_bug21684(bool display)
{
    cxx_double_poly f(3);
    for(double x : {-6.3406659802246472e+28, 6.4695148868405632e+28, 9.5457310557271272e+27}) f->coeff[++f->deg]=x;
    piecewise_linear_approximator A(f, 0.3300700859809263);
    piecewise_linear_function F = A.logapprox(-2048,2048);
    ASSERT_ALWAYS(F.equations.size() == 21);
    if (display) display_logapprox(F);
    return 0;
}

int test_from_bug21701(bool display)
{
    cxx_double_poly f(4);
    for(double x : {-1.6425515054690201e+34, -2.4119595460727266e+36, -1.1805918683048612e+38, -1.926230702854181e+39}) f->coeff[++f->deg]=x;
    piecewise_linear_approximator A(f, 0.44719172939351309);
    piecewise_linear_function F = A.logapprox(-2048,2048);
    ASSERT_ALWAYS(F.equations.size() == 57);
    if (display) display_logapprox(F);
    return 0;
}

int original_test(bool display)
{
    cxx_double_poly f(6);
    for(double x : {-39769265437440000., -302859053854976., 5377439145928., -1684314626., -5377481., 970., 1.}) f->coeff[++f->deg]=x;
    piecewise_linear_approximator A(f, 0.1);
    piecewise_linear_function F = A.logapprox(-2048,2048);
    ASSERT_ALWAYS(F.equations.size() == 34);
    if (display) display_logapprox(F);
    return 0;
}

// coverity[root_function]
int main(int argc, char * argv[])
{
    bool display = false;
    if (argc >= 2 && strcmp(argv[1], "--display") == 0)
        display=true;
    original_test(display);
    test_from_bug21684(display);
    test_from_bug21701(display);
}

