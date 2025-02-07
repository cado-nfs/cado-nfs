#include "cado.h" // IWYU pragma: keep

#include <cfenv>
#include <cmath>             // for exp
#include <cstdio>            // for fprintf, stderr
#include <algorithm>         // for sort, max
#include <initializer_list>  // for initializer_list
#include <iterator>          // for begin, end
#include <limits>            // for numeric_limits
#include <vector>
#include <utility>
#include <list>

#include "fmt/core.h"

#include "logapprox.hpp"
#include "macros.h"          // for ASSERT_ALWAYS, UNLIKELY
#include "runtime_numeric_cast.hpp"
#include "polynomial.hpp"


#define xxxDEBUG_LOGAPPROX

/* This should be 1/Jmax */
#define SMALLEST_MEANINGFUL_LOGAPPROX_INTERVAL  (1.0/262144)

static std::vector<double> roots_helper(polynomial<double> const & f)
{
    /* This uses polynomial<T>, possibly for T=long double, for the roots
     * computation. T=long double is more accurate, but unfortunately it's
     * also considerably slower, to the point that it's probably barely
     * usable. The impact is per-special q and not per-sieve region, so it's
     * not *that* bad, but still, it's troubling enough (about 8 extra
     * seconds on the test case of #30107).
     */
    typedef double T;
    polynomial<T> lf;
    for(int i = 0 ; i <= f.degree() ; i++) lf[i] = f[i];
    auto v = lf.roots();
    std::vector<double> res(v.begin(), v.end());

    sort(begin(res), end(res));
    return res;
}

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

piecewise_linear_approximator::piecewise_linear_approximator(polynomial<double> const & f, double scale)
    : f(f)
    , f1(f.derivative())
    , f_roots(roots_helper(f))
    , f1_roots(roots_helper(f1))
    , scale(scale)
{ }

std::vector<double> piecewise_linear_approximator::roots_off_course(polynomial<double> const& uv, bool divide_root, double r) const/*{{{*/
{
    std::vector<double> res;
    for(double const m : { std::exp(scale), std::exp(-scale) }) {
        polynomial<double> d = f;
        d.submul(uv, m);
        if (divide_root)
            d = d.div_linear(r);
        auto roots = roots_helper(d);
        res.insert(res.end(), roots.begin(), roots.end());
    }
    return res;
}/*}}}*/

piecewise_linear_function piecewise_linear_approximator::expand_at_root(double r) const /*{{{*/
{
    double const v = f1(r);
    double const u = -r*v;
    polynomial<double> uv { u, v };
    double r0 = std::numeric_limits<double>::lowest();
    double r1 = std::numeric_limits<double>::max();

    if (v == 0)
        /* inflection points cause real trouble */
        fprintf(stderr, "# Warning: logapprox encountered a (quasi-) inflection point. We cannot fulfill our requirements\n");


    for(double const x : roots_off_course(uv, true, r)) {
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
    res.equations.emplace_back(u, v);
    return res;
}/*}}}*/

piecewise_linear_function piecewise_linear_approximator::C0_from_points(std::list<double> const & r) const /*{{{*/
{
    /* This computes the chords that connect the different points in the
     * list r. Each has equation u+v*x. r must have size at least two
     * (well, a size one list will yield no chord and that's it)
     */
    piecewise_linear_function res;
    res.endpoints = r;
    ASSERT_ALWAYS(!r.empty());
    auto it = r.begin();
    double a = *it++;
    double fa = f(a);
    for(; it != r.end() ; it++) {
        double const b = *it;
        double const fb = f(b);
        double const u = (b*fa-a*fb)/(b-a);
        double const v = (fb-fa)/(b-a);
        res.equations.emplace_back(u, v);
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
piecewise_linear_function piecewise_linear_approximator::fill_gap(double i0, double i1) const {/*{{{*/
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
        double const r0 = done.endpoints.back();
        /* This one really should be r0 anyway, and it's slightly
         * boring to have to deal with it.
         */
        todo.endpoints.pop_front();
        double const r1 = todo.endpoints.front();
        /* Arrange so that we don't emit linear approximations on
         * subintervals that are absurdly close to the existing ones.
         */
        double const r0s = r0 + SMALLEST_MEANINGFUL_LOGAPPROX_INTERVAL;
        double const r1s = r1 - SMALLEST_MEANINGFUL_LOGAPPROX_INTERVAL;
#ifdef DEBUG_LOGAPPROX
        printf("Checking interval %f, %f\n", r0,r1);
#endif
        std::pair<double, double> const uv=todo.equations.front();
        todo.equations.pop_front();

        polynomial<double> guv { uv.first, uv.second };
        std::vector<double> roots;
        for(double const r : roots_off_course(guv)) {
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
                next_noderivativeroot = runtime_numeric_cast<int>(newsplits.size()) - 1;
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

piecewise_linear_function piecewise_linear_approximator::logapprox(double i0, double i1) const /*{{{*/
{
    if (f.degree() == 1) {
        /* the general code does not work well for linear functions.
         * Well, of course, a linear approximation to a linear function
         * is easy to find.
         */
        piecewise_linear_function res(i0);
        res.endpoints.push_back(i1);
        ASSERT_ALWAYS(f_roots.size() == 1);
        double const r = f_roots.front();
        double const v = f1(r);
        double const u = -r*v;
        res.equations.emplace_back(u, v);
        return res;
    }
#ifdef DEBUG_LOGAPPROX
    fmt::print("Computing lognorm approximation for {} over interval [{},{}]\n",
            f, i0, i1);
#endif
    piecewise_linear_function done(i0);
    for(auto r: f_roots) {
        if (r < i0) continue;
        if (r > i1) break;
        piecewise_linear_function s = expand_at_root(r);
        double const u1 = done.endpoints.back();
        double const v0 = s.endpoints.front();
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
    double const u1 = done.endpoints.back();
    if (u1 + SMALLEST_MEANINGFUL_LOGAPPROX_INTERVAL < i1) {
        piecewise_linear_function G = fill_gap(u1, i1);
        done.merge_right(G);
    }
    done.endpoints.back() = std::min(done.endpoints.back(), i1);

    return done;
}/*}}}*/
