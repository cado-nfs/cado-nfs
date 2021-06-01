#include "cado.h" // IWYU pragma: keep

#include <cmath>             // for exp
#include <cstdio>            // for fprintf, stderr
#include <algorithm>         // for sort, max
#include <initializer_list>  // for initializer_list
#include <iterator>          // for begin, end
#include <limits>            // for numeric_limits
#include "double_poly.h"  // cxx_double_poly
#include "logapprox.hpp"
#include "macros.h"          // for ASSERT_ALWAYS, UNLIKELY


#define xxxDEBUG_LOGAPPROX

/* This should be 1/Jmax */
#define SMALLEST_MEANINGFUL_LOGAPPROX_INTERVAL  (1.0/262144)

/* modify our left end to include the approximation which is provided
 * on argument.
 * Warning: this is made constant time by using splice(), so that the
 * argument is destroyed ! */
piecewise_linear_function& piecewise_linear_function::merge_left(piecewise_linear_function & o) {/*{{{*/
    ASSERT_ALWAYS(endpoints.front() == o.endpoints.back());
    endpoints.pop_front();
    endpoints.splice(endpoints.begin(), o.endpoints);
    equations.splice(equations.begin(), o.equations);
    return *this;
}/*}}}*/
piecewise_linear_function& piecewise_linear_function::merge_right(piecewise_linear_function & o) {/*{{{*/
    ASSERT_ALWAYS(endpoints.back() == o.endpoints.front());
    endpoints.pop_back();
    endpoints.splice(endpoints.end(), o.endpoints);
    equations.splice(equations.end(), o.equations);
    return *this;
}/*}}}*/

piecewise_linear_approximator::piecewise_linear_approximator(cxx_double_poly const & f, double scale) : f(f), scale(scale) {/*{{{*/
    f_roots.assign(f->deg, 0);
    unsigned int n = double_poly_compute_all_roots(f_roots.data(), f);
    f_roots.erase(f_roots.begin() + n, f_roots.end());
    sort(begin(f_roots), end(f_roots));

    if (f->deg == 1) return;

    double_poly_derivative(f1, f);
    f1_roots.assign(f1->deg, 0);
    n = double_poly_compute_all_roots(f1_roots.data(), f1);
    f1_roots.erase(f1_roots.begin() + n, f1_roots.end());
    sort(begin(f1_roots), end(f1_roots));
}/*}}}*/

std::vector<double> piecewise_linear_approximator::roots_off_course(cxx_double_poly const& uv, bool divide_root, double r)/*{{{*/
{
    std::vector<double> res;
    for(double m : { std::exp(scale), std::exp(-scale) }) {
        cxx_double_poly d;
        double_poly_set(d, f);
        double_poly_submul_double(d, uv, m);
        if (divide_root)
            double_poly_div_linear(d, d, r);
        std::vector<double> roots(d->deg, 0);
        unsigned int n = double_poly_compute_all_roots(roots.data(), d);
        res.insert(res.end(), roots.begin(), roots.begin() + n);
    }
    return res;
}/*}}}*/

piecewise_linear_function piecewise_linear_approximator::expand_at_root(double r)/*{{{*/
{
    double v = double_poly_eval(f1, r);
    double u = -r*v;
    cxx_double_poly uv(1);
    uv->coeff[0] = u;
    uv->coeff[1] = v;
    double_poly_cleandeg(uv, 1);
    double r0 = std::numeric_limits<double>::lowest();
    double r1 = std::numeric_limits<double>::max();

    if (!v)
        /* inflection points cause real trouble */
        fprintf(stderr, "# Warning: logapprox encountered a (quasi-) inflection point. We cannot fulfill our requirements\n");


    for(double x : roots_off_course(uv, true, r)) {
        if (x < r && x > r0) r0 = x;
        if (x > r && x < r1) r1 = x;
    }

    /* This should not happen, except maybe near inflection points */
    if (UNLIKELY(r0 == std::numeric_limits<double>::lowest())) {
        fprintf(stderr, "# Warning: using last resort bailout code in logapprox\n");
        r0 = r - SMALLEST_MEANINGFUL_LOGAPPROX_INTERVAL;
    }
    if (UNLIKELY(r1 == std::numeric_limits<double>::max())) {
        fprintf(stderr, "# Warning: using last resort bailout code in logapprox\n");
        r1 = r + SMALLEST_MEANINGFUL_LOGAPPROX_INTERVAL;
    }

    piecewise_linear_function res;
    res.endpoints.push_back(r0);
    res.endpoints.push_back(r1);
    res.equations.push_back(std::make_pair(u,v));
    return res;
}/*}}}*/

piecewise_linear_function piecewise_linear_approximator::C0_from_points(std::list<double> const & r)/*{{{*/
{
    piecewise_linear_function res;
    res.endpoints = r;
    auto it = r.begin();
    double a = *it++;
    double fa = double_poly_eval(f, a);
    for(; it != r.end() ; it++) {
        double b = *it;
        double fb = double_poly_eval(f, b);
        double u = (b*fa-a*fb)/(b-a);
        double v = (fb-fa)/(b-a);
        res.equations.push_back(std::make_pair(u,v));
        a = b;
        fa = fb;
    }
    return res;
}/*}}}*/

/* This function is where most of the stuff happens. We do our best to
 * find a piecewise linear function which stays within the required error
 * margin. It is not entirely obvious how one should proceed to obtain an
 * optimal result.
 */
/* This assumes that the interval [i0,i1] is free of roots of the
 * polynomial f */
piecewise_linear_function piecewise_linear_approximator::fill_gap(double i0, double i1) {/*{{{*/
    piecewise_linear_function todo, done;
    todo = C0_from_points(std::list<double>({i0,i1}));
    done = piecewise_linear_function(i0);
    int next_noderivativeroot=0;
    for( ; !todo.equations.empty() ; ) {
#ifdef DEBUG_LOGAPPROX
        printf("Done %zu pieces, %zu to go\n",
                done.equations.size(),
                todo.equations.size());
#endif
        double r0 = done.endpoints.back();
        /* This one really should be r0 anyway, and it's slightly
         * boring to have to deal with it.
         */
        todo.endpoints.pop_front();
        double r1 = todo.endpoints.front();
        /* Arrange so that we don't emit linear approximations on
         * subintervals that are absurdly close to the existing ones.
         */
        double r0s = r0 + SMALLEST_MEANINGFUL_LOGAPPROX_INTERVAL;
        double r1s = r1 - SMALLEST_MEANINGFUL_LOGAPPROX_INTERVAL;
#ifdef DEBUG_LOGAPPROX
        printf("Checking interval %f, %f\n", r0,r1);
#endif
        std::pair<double, double> uv=todo.equations.front();
        todo.equations.pop_front();

        cxx_double_poly guv(1);
        guv->coeff[0] = uv.first;
        guv->coeff[1] = uv.second;
        double_poly_cleandeg(guv,1);
        std::vector<double> roots;
        for(double r : roots_off_course(guv)) {
            if (r >= r0s && r <= r1s)
                roots.push_back(r);
        }
        sort(begin(roots), end(roots));
        /* intersect and sort */
#ifdef DEBUG_LOGAPPROX
        printf("Checking %f+x*%f as an estimator to f(x) in interval %f, %f. Roots are:\n", uv.first, uv.second, r0, r1);
        for(double x : roots) printf(" %f\n", x);
#endif
        if (roots.empty()) {
            done.endpoints.push_back(r1);
            done.equations.push_back(uv);
#ifdef DEBUG_LOGAPPROX
            printf("No conflict, that makes the %zu-th piece\n", done.equations.size());
#endif
            if (next_noderivativeroot)
                next_noderivativeroot--;
        } else {
            std::list<double> newsplits;
#ifdef DEBUG_LOGAPPROX
            printf("Conflict\n");
#endif
            if (next_noderivativeroot == 0) {
#ifdef DEBUG_LOGAPPROX
                printf("Haven't tried derivative roots on this one yet\n");
#endif
                /* Try to split at the roots of the derivative if we have any */
                newsplits.push_back(r0);
                for(auto r: f1_roots) {
                    if (r <= r0s) continue;
                    if (r >= r1s) break;
                    newsplits.push_back(r);
                }
                newsplits.push_back(r1);
                next_noderivativeroot = newsplits.size() - 1;
            } else {
                /* either we insert a midpoint which is the earliest
                 * of the off-course roots we identified, or we
                 * simply pick the middle of the segment. We opt for
                 * the latter, here.
                 *
                 * (for the former, we may use the boolean to mark
                 * the fact that we already know there's no root in
                 * the new interval).
                 */
                newsplits = std::list<double>({r0,(r0+r1)/2,r1});
                next_noderivativeroot--;
                next_noderivativeroot+=2;
            }
#ifdef DEBUG_LOGAPPROX
            printf("Set of new stop points is:\n");
            for(auto x : newsplits) printf(" %f\n", x);
#endif

            piecewise_linear_function G = C0_from_points(newsplits);
            G.endpoints.pop_front();
            todo.endpoints.pop_front();
            todo.endpoints.splice(todo.endpoints.begin(), G.endpoints);
            todo.equations.splice(todo.equations.begin(), G.equations);
            /* need to add r0 */
            todo.endpoints.push_front(done.endpoints.back());
        }
    }
    return done;
}/*}}}*/

piecewise_linear_function piecewise_linear_approximator::logapprox(double i0, double i1)/*{{{*/
{
    if (f->deg == 1) {
        /* the general code does not work well for linear functions.
         * Well, of course, a linear approximation to a linear function
         * is easy to find.
         */
        piecewise_linear_function res(i0);
        res.endpoints.push_back(i1);
        ASSERT_ALWAYS(f_roots.size() == 1);
        double r = f_roots.front();
        double v = double_poly_eval(f1, r);
        double u = -r*v;
        res.equations.push_back(std::make_pair(u,v));
        return res;
    }
    piecewise_linear_function done(i0);
    for(auto r: f_roots) {
        if (r < i0) continue;
        if (r > i1) break;
        piecewise_linear_function s = expand_at_root(r);
        double u1 = done.endpoints.back();
        double v0 = s.endpoints.front();
        if (v0 > u1 + SMALLEST_MEANINGFUL_LOGAPPROX_INTERVAL) {
            piecewise_linear_function G = fill_gap(u1, v0);
            done.merge_right(G);
        }
        /* important: otherwise we'll have overlapping intervals (and
         * a failing assert, too).
         */
        s.endpoints.front() = done.endpoints.back();
        done.merge_right(s);
    }
    double u1 = done.endpoints.back();
    if (u1 + SMALLEST_MEANINGFUL_LOGAPPROX_INTERVAL < i1) {
        piecewise_linear_function G = fill_gap(u1, i1);
        done.merge_right(G);
    }
    if (done.endpoints.back() >= i1)
        done.endpoints.back() = i1;

    return done;
}/*}}}*/
