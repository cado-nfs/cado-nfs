#ifndef CADO_LOGAPPROX_HPP
#define CADO_LOGAPPROX_HPP

#include <list>     // for list
#include <utility>  // for pair
#include <vector>   // for vector

#include "polynomial.hpp"

struct piecewise_linear_function {
    std::vector<double> endpoints;
    std::vector<std::pair<double,double>> equations;
    bool has_precision_issues = false;
    piecewise_linear_function() = default;

    /* this is a precursor to piecewise_linear_function */
    struct precursor {
        std::list<double> endpoints;
        std::list<std::pair<double,double>> equations;
        bool has_precision_issues = false;
        precursor() = default;
        explicit precursor(double r) : endpoints(1, r) {}
        precursor(precursor const &) = default;
        precursor(precursor &&) = default;
        precursor & operator=(precursor const &) = default;
        precursor & operator=(precursor &&) = default;
        ~precursor() = default;
        /* modify our left end to include the approximation which is provided
         * on argument.
         * Warning: this is made constant time by using splice(), so that the
         * argument is destroyed ! */
        precursor& merge_left(precursor & o);
        /* guess what... */
        precursor& merge_right(precursor & o);
    };

    explicit piecewise_linear_function(precursor const & p)
        : endpoints(p.endpoints.begin(), p.endpoints.end())
        , equations(p.equations.begin(), p.equations.end())
        , has_precision_issues(p.has_precision_issues)
    {}
};

template<typename T>
class piecewise_linear_approximator {
    polynomial<T> f;
    polynomial<T> f1;
    polynomial<T> f2;
    std::vector<T> f_roots;
    std::vector<T> f1_roots;
    T scale;
    std::vector<T> roots_off_course(polynomial<T> const& uv) const;
    piecewise_linear_function::precursor expand_at_root(T r) const;
    piecewise_linear_function::precursor C0_from_points(std::list<T> const & r) const;
    /* This assumes that the interval [i0,i1] is free of roots of the
     * polynomial f */
    piecewise_linear_function::precursor fill_gap(T i0, T i1) const;
    public:
    piecewise_linear_approximator(polynomial<T> const & f, T scale);
    /* exceptions will only be thrown if the three conditions below hold:
     *  - (super)::report_precision_issues_on_double
     *  - T == double
     *  - allow_exceptions
     * otherwise the presence of precision issues is only reported in the
     * has_precision_issue flag of the result.
     */
    piecewise_linear_function logapprox(T i0, T i1) const;
};

extern template class piecewise_linear_approximator<double>;
extern template class piecewise_linear_approximator<long double>;

#endif	/* CADO_LOGAPPROX_HPP */
