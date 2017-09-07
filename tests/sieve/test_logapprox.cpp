#include "cado.h"
#include "logapprox.hpp"
#include <iostream>
int main()
{
    cxx_double_poly f(6);
    for(double x : {-39769265437440000., -302859053854976., 5377439145928., -1684314626., -5377481., 970., 1.}) f->coeff[++f->deg]=x;
    piecewise_linear_approximator A(f, 0.1);
    piecewise_linear_function F = A.logapprox(-2048,2048);
    ASSERT_ALWAYS(F.equations.size() == 34);
    /*
    std::cout << "Total " << F.equations.size() << " pieces\n";
    for( ; !F.equations.empty() ; ) {
        double r0 = F.endpoints.front();
        F.endpoints.pop_front();
        double r1 = F.endpoints.front();
        std::pair<double,double> uv = F.equations.front();
        F.equations.pop_front();
        std::cout << "[" << r0 << ", " << r1 << "]: " << uv.first << " + x * " << uv.second << "\n";
    }
    */
}
