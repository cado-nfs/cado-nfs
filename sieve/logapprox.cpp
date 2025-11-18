#include "cado.h" // IWYU pragma: keep

#include <cmath>             // for exp
#include <cstdio>            // for fprintf, stderr

#include <algorithm>         // for sort, max
#include <array>
#include <limits>            // for numeric_limits
#include <list>
#include <utility>
#include <vector>

#include "fmt/base.h"

#include "logapprox.hpp"
#include "macros.h"          // for ASSERT_ALWAYS, UNLIKELY
#include "runtime_numeric_cast.hpp"
#include "polynomial.hpp"
#include "cado_math_aux.hpp"


#define xxxDEBUG_LOGAPPROX

/* This should be 1/Jmax */
#define SMALLEST_MEANINGFUL_LOGAPPROX_INTERVAL  (1.0/262144)


/* modify our left end to include the approximation which is provided
 * on argument.
 * Warning: this is made constant time by using splice(), so that the
 * argument is destroyed ! */
piecewise_linear_function::precursor& piecewise_linear_function::precursor::merge_left(piecewise_linear_function::precursor & o) {/*{{{*/
    ASSERT_ALWAYS(endpoints.front() == o.endpoints.back());
    endpoints.pop_front();
    endpoints.splice(endpoints.begin(), o.endpoints);
    equations.splice(equations.begin(), o.equations);
    if (o.has_precision_issues) has_precision_issues = true;
    return *this;
}/*}}}*/
piecewise_linear_function::precursor& piecewise_linear_function::precursor::merge_right(piecewise_linear_function::precursor & o) {/*{{{*/
    ASSERT_ALWAYS(endpoints.back() == o.endpoints.front());
    endpoints.pop_back();
    endpoints.splice(endpoints.end(), o.endpoints);
    equations.splice(equations.end(), o.equations);
    if (o.has_precision_issues) has_precision_issues = true;
    return *this;
}/*}}}*/

template<typename T>
piecewise_linear_approximator<T>::piecewise_linear_approximator(polynomial<T> const & f, T scale)
    : f(f)
    , f1(f.derivative())
    , f2(f1.derivative())
    , f_roots(f.template roots<T>())
    , f1_roots(f1.template roots<T>())
    , scale(scale)
{ }

template<typename T>
std::vector<T> piecewise_linear_approximator<T>::roots_off_course(polynomial<T> const& uv) const/*{{{*/
{
    std::vector<T> res;
    for(T const m : { std::exp(scale), std::exp(-scale) }) {
        polynomial<T> d = f;
        d.submul(uv, m);
        auto roots = d.template roots<T>();
        res.insert(res.end(), roots.begin(), roots.end());
    }
    return res;
}/*}}}*/

template<typename T>
piecewise_linear_function::precursor piecewise_linear_approximator<T>::expand_at_root(T r) const /*{{{*/
{
    const bool verbose = false;
    T const v = f1(r);
    T const u = -r*v;
    const polynomial<T> uv { -r*v, v };
    using namespace cado_math_aux;

    piecewise_linear_function::precursor res;
    res.equations.emplace_back(u, v);

    const T ulp = cado_math_aux::ulp(r);
    const T before = f(r-ulp);
    const T after = f(r+ulp);

    /* signs of the first and second order derivatives */
    const int D1 = sgn(v);
    const int D2 = sgn(f2(r));

    T r0 = std::numeric_limits<T>::lowest();
    T r1 = std::numeric_limits<T>::max();

    if (sgn(before) == sgn(after) || D1 != sgn(after) || v == 0) {
        if (!res.has_precision_issues && verbose)
            fmt::print(stderr,
                    "# Warning: we have precision issues around root {};"
                    " making the validity interval very short\n", r);
        res.has_precision_issues = true;
        r0 = r - SMALLEST_MEANINGFUL_LOGAPPROX_INTERVAL;
        r1 = r + SMALLEST_MEANINGFUL_LOGAPPROX_INTERVAL;
        /* the rest of the code will follow smoothly and keep r0 and r1
         * unchanged */
    }

    /* The approximation that is most endangered by tiny changes _at this
     * root_ is fences[1^(D1>0)^(D2>0)] above r, and [(D1>0)^(D2>0)]
     * below r.  However, long term, it might be that the other one
     * breaks first!
     */
    std::array<polynomial<T>, 2> fences {
        uv * std::exp(-scale),
        uv * std::exp(scale),
    };

#if 0
    /* This is visibly too much to ask. There are cases where what
     * happens at the infinitesimal level is far from satisfactory, and
     * yet the cuts that we compute with the root finding give us pretty
     * reasonable results.
     */
    if (fences[(D1>0)^(D2>0)](r-ulp) * D2 < before * D2) {
        if (!res.has_precision_issues && verbose)
            fmt::print("# root {} leads to precision issues;"
                    " we might want to start over with larger precision."
                    " For safety, we're not letting this interval expand"
                    " for longer than one ulp\n", r);
        res.has_precision_issues = true;
        r0 = r - SMALLEST_MEANINGFUL_LOGAPPROX_INTERVAL;
    }

    if (fences[(D1>0)^(D2>0)^1](r+ulp) * D2 < after * D2) {
        if (!res.has_precision_issues && verbose)
            fmt::print("# root {} leads to precision issues;"
                    " we might want to start over with larger precision."
                    " For safety, we're not letting this interval expand"
                    " for longer than one ulp\n", r);
        res.has_precision_issues = true;
        r1 = r + SMALLEST_MEANINGFUL_LOGAPPROX_INTERVAL;
    }
#endif

    /* So the general reasoning that consists in letting the root
     * computation guide our choices of r0 and r1 is of course the most
     * important one.
     */
    for(int i = 0 ; i < 2 ; i++) {
        auto const & approx = fences[i];
        /* Here we do not want to divide by x-r, by fear of precision
         * issues. I'm not totally sure of how legitimate it is, though.
         * It means that there's going to be a parasite root at r, which
         * is a bit annoying to deal with.
         */

        /* Not _all_ crossings of the approximations with f make sense.
         * There are some that we expect, based on our understanding of
         * the first and second order derivatives at r. If we do
         * encounter these, fine. If we don't then everything else we get
         * is likely to be spurious and then, by all means, we want to
         * disregard these crossings and shrink the validity interval as
         * much as we can (it doesn't have to be one ulp:
         * SMALLEST_MEANINGFUL_LOGAPPROX_INTERVAL can be used for that).
         */
        for(T const x : (f - approx).template roots<T>()) {
            const bool best_below = x < r && x > r0;
            const bool best_above = x > r && x < r1;
            if (!(best_above || best_below))
                continue;
            const bool away = i == 1 && D1*f1(x) > D1*approx[1];
            const bool near = i == 0 && D1*f1(x) < D1*approx[1];
            if (!(away || near))
                continue;
            if (best_below) {
                if (verbose)
                    fmt::print("# root {}: approximation {} is invalid below"
                                " {} [slope {} convexity {}] (r - {} ulps)"
                                " escapes towards {}\n",
                                r, "-+"[i], x, D1, D2, (r-x)/ulp,
                                away ? "infinity" : "0");
                r0 = x;
            } else if (best_above) {
                if (verbose)
                    fmt::print("# root {}: approximation {} is invalid above"
                                " {} [slope {} convexity {}] (r + {} ulps)"
                                " escapes towards {}\n",
                                r, "-+"[i], x, D1, D2, (x-r)/ulp,
                                away ? "infinity" : "0");
                r1 = x;
            }
        }
    }

    if (r0 == std::numeric_limits<T>::lowest()) {
        if (!res.has_precision_issues && verbose)
            fmt::print("# root {} leads to precision issues;"
                    " we might want to start over with larger precision."
                    " For safety, we're not letting this interval expand"
                    " for longer than one ulp\n", r);
        res.has_precision_issues = true;
        r0 = r - SMALLEST_MEANINGFUL_LOGAPPROX_INTERVAL;
    }

    if (r1 == std::numeric_limits<T>::max()) {
        if (!res.has_precision_issues && verbose)
            fmt::print("# root {} leads to precision issues;"
                    " we might want to start over with larger precision."
                    " For safety, we're not letting this interval expand"
                    " for longer than one ulp\n", r);
        res.has_precision_issues = true;
        r1 = r + SMALLEST_MEANINGFUL_LOGAPPROX_INTERVAL;
    }

    /* Here's a final sanity check. We want to be sure that we're indeed
     * in the correct interval.
     */
    unsigned int ns0 = 0;
    for(T x ; r0 < r - SMALLEST_MEANINGFUL_LOGAPPROX_INTERVAL ; r0 = x, ns0++) {
        x = (r+r0) / 2;
        T fx = f(x);
        if (D1 * fences[1](x) <= D1 * fx && D1 * fx <= D1 * fences[0](x))
            break;
    }
    if (ns0)
        fmt::print("# shortened approximation interval below {} by 2^-{}\n",
                r, ns0);

    unsigned int ns1 = 0;
    for(T x ; r1 > r + SMALLEST_MEANINGFUL_LOGAPPROX_INTERVAL ; r1 = x, ns1++) {
        x = (r+r1) / 2;
        T fx = f(x);
        if (D1 * fences[0](x) <= D1 * fx && D1 * fx <= D1 * fences[1](x))
            break;
    }
    if (ns1)
        fmt::print("# shortened approximation interval above {} by 2^-{}\n",
                r, ns1);




    res.endpoints.push_back(r0);
    res.endpoints.push_back(r1);
    return res;
}/*}}}*/

template<typename T>
piecewise_linear_function::precursor piecewise_linear_approximator<T>::C0_from_points(std::list<T> const & r) const /*{{{*/
{
    /* This computes the chords that connect the different points in the
     * list r. Each has equation u+v*x. r must have size at least two
     * (well, a size one list will yield no chord and that's it)
     */
    ASSERT_ALWAYS(!r.empty());
    piecewise_linear_function::precursor res;
    for(auto x : r) res.endpoints.push_back(double(x));
    auto it = r.begin();
    T a = *it++;
    T fa = f(a);
    for(; it != r.end() ; it++) {
        T const b = *it;
        T const fb = f(b);
        T const u = (b*fa-a*fb)/(b-a);
        T const v = (fb-fa)/(b-a);
        res.equations.emplace_back(double(u), double(v));
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
template<typename T>
piecewise_linear_function::precursor piecewise_linear_approximator<T>::fill_gap(T i0, T i1) const {/*{{{*/
    piecewise_linear_function::precursor todo, done;
    todo = C0_from_points(std::list<T>({i0,i1}));
    done = piecewise_linear_function::precursor(i0);
    int next_noderivativeroot=0;
    for( ; !todo.equations.empty() ; ) {
#ifdef DEBUG_LOGAPPROX
        printf("Done %zu pieces, %zu to go\n",
                done.equations.size(),
                todo.equations.size());
#endif
        T const r0 = done.endpoints.back();
        /* This one really should be r0 anyway, and it's slightly
         * boring to have to deal with it.
         */
        todo.endpoints.pop_front();
        T const r1 = todo.endpoints.front();
        /* Arrange so that we don't emit linear approximations on
         * subintervals that are absurdly close to the existing ones.
         */
        T const r0s = r0 + SMALLEST_MEANINGFUL_LOGAPPROX_INTERVAL;
        T const r1s = r1 - SMALLEST_MEANINGFUL_LOGAPPROX_INTERVAL;
#ifdef DEBUG_LOGAPPROX
        printf("Checking interval %f, %f\n", r0,r1);
#endif
        std::pair<T, T> const uv=todo.equations.front();
        todo.equations.pop_front();

        const polynomial<T> guv { uv.first, uv.second };
        std::vector<T> roots;
        for(T const r : roots_off_course(guv)) {
            if (r >= r0s && r <= r1s)
                roots.push_back(r);
        }
        sort(begin(roots), end(roots));
        /* intersect and sort */
#ifdef DEBUG_LOGAPPROX
        printf("Checking %f+x*%f as an estimator to f(x) in interval %f, %f. Roots are:\n", uv.first, uv.second, r0, r1);
        for(T x : roots) printf(" %f\n", x);
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
            std::list<T> newsplits;
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
                newsplits = std::list<T>({r0,(r0+r1)/2,r1});
                next_noderivativeroot--;
                next_noderivativeroot+=2;
            }
#ifdef DEBUG_LOGAPPROX
            printf("Set of new stop points is:\n");
            for(auto x : newsplits) printf(" %f\n", x);
#endif

            piecewise_linear_function::precursor G = C0_from_points(newsplits);
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

template<typename T>
piecewise_linear_function piecewise_linear_approximator<T>::logapprox(T i0, T i1) const /*{{{*/
{
    if (f.degree() == 1) {
        /* the general code does not work well for linear functions.
         * Well, of course, a linear approximation to a linear function
         * is easy to find.
         */
        piecewise_linear_function::precursor res(i0);
        res.endpoints.push_back(i1);
        ASSERT_ALWAYS(f_roots.size() == 1);
        T const r = f_roots.front();
        T const v = f1(r);
        T const u = -r*v;
        res.equations.emplace_back(u, v);
        return piecewise_linear_function(res);
    }
#ifdef DEBUG_LOGAPPROX
    fmt::print("Computing lognorm approximation for {} over interval [{},{}]\n",
            f, i0, i1);
#endif
    piecewise_linear_function::precursor done(i0);
    for(auto r: f_roots) {
        if (r < i0) continue;
        if (r > i1) break;
        piecewise_linear_function::precursor s = expand_at_root(r);
        T const u1 = done.endpoints.back();
        T const v0 = s.endpoints.front();
        if (v0 > u1 + SMALLEST_MEANINGFUL_LOGAPPROX_INTERVAL) {
            piecewise_linear_function::precursor G = fill_gap(u1, v0);
            done.merge_right(G);
        }
        /* important: otherwise we'll have overlapping intervals (and
         * a failing assert, too).
         */
        s.endpoints.front() = done.endpoints.back();
        done.merge_right(s);
    }
    T const u1 = done.endpoints.back();
    if (u1 + SMALLEST_MEANINGFUL_LOGAPPROX_INTERVAL < i1) {
        piecewise_linear_function::precursor G = fill_gap(u1, i1);
        done.merge_right(G);
    }
    done.endpoints.back() = std::min(done.endpoints.back(), double(i1));

    return piecewise_linear_function(done);
}/*}}}*/

template class piecewise_linear_approximator<double>;
template class piecewise_linear_approximator<long double>;
