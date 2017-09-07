#include <vector>
#include <math.h>
#include <list>
#include "double_poly.h"

struct piecewise_linear_function {
    std::list<double> endpoints;
    std::list<std::pair<double,double>> equations;
    piecewise_linear_function() = default;
    piecewise_linear_function(double r) : endpoints(1, r) {}
    /* modify our left end to include the approximation which is provided
     * on argument.
     * Warning: this is made constant time by using splice(), so that the
     * argument is destroyed ! */
    piecewise_linear_function& merge_left(piecewise_linear_function & o);
    /* guess what... */
    piecewise_linear_function& merge_right(piecewise_linear_function & o);
};

class piecewise_linear_approximator {
    cxx_double_poly const & f;
    cxx_double_poly f1;
    std::vector<double> f_roots;
    std::vector<double> f1_roots;
    double scale;
    std::vector<double> roots_off_course(cxx_double_poly const& uv, bool divide_root=false, double r = 0);
    piecewise_linear_function expand_at_root(double r);
    piecewise_linear_function C0_from_points(std::list<double> const & r);
    /* This assumes that the interval [i0,i1] is free of roots of the
     * polynomial f */
    piecewise_linear_function fill_gap(double i0, double i1);
    public:
    piecewise_linear_approximator(cxx_double_poly const & f, double scale);
    piecewise_linear_function logapprox(double i0, double i1);
};
