#ifndef CADO_UTILS_POLYNOMIAL_HPP
#define CADO_UTILS_POLYNOMIAL_HPP

#include "cado_config.h"        // IWYU pragma: keep

/* This defines polynomial over arbitrary types, provided that these have
 * standard operator overloads defined. A priori we want to instantiate
 * these with float, double, and long double. But in principle it should
 * be possible to use this type more generically
 */

#include <cctype>
#include <cmath>
#include <cstddef>

#include <algorithm>
#include <initializer_list>
#include <ios>
#include <istream>
#include <ostream>
#include <sstream>
#include <string>      // for string
#include <type_traits>
#include <utility>
#include <vector>
#include <functional>


#include <gmp.h>
#include "fmt/base.h"
#include "fmt/ostream.h"

#include "macros.h"
#include "runtime_numeric_cast.hpp"
#include "number_context.hpp"
#include "mpz_poly.h"
#include "cxx_mpz.hpp"
#include "cado_math_aux.hpp"
#include "cado_expression_parser.hpp"
#include "coeff_proxy.hpp"
#include "cado_type_traits.hpp"
#include "named_proxy.hpp"
#include "extra_complex_overloads.hpp"

#ifdef HAVE_MPFR
#include "cxx_mpfr.hpp"
#endif

#ifdef HAVE_MPC
#include "cxx_mpc.hpp"
#endif

namespace polynomial_details {
template<typename T>
struct polynomial;
}

/* coeff_proxy is sort of a doomed approach, because it can't be used
 * transparently when the coefficient type is a template. The root cause
 * is that template argument deduction ignores implicit conversions in
 * that case. All we can do is try to provide overloads for some easy
 * cases that force the coeff_proxy transparent layer, but that's by no
 * means a definitive solution.
 */
#if 0
#define CADO_COEFF_PROXY_OVERLOAD(OP)					\
    template<typename T>						\
    static inline auto OP(					\
        T const & a,							\
        cado_details::coeff_proxy<polynomial_details::polynomial<T>> const & b) \
        -> decltype(OP(a, T(b)))                                        \
    {									\
        return OP(a, T(b));						\
    }									\
    template<typename T>						\
    static inline auto OP(					        \
        cado_details::coeff_proxy<polynomial_details::polynomial<T>> const & a, \
        T const & b)							\
        -> decltype(OP(T(a), b))                                        \
    {									\
        return OP(T(a), b);						\
    }									\
    template<typename T>						\
    static inline auto OP(					        \
        cado_details::coeff_proxy<polynomial_details::polynomial<T>> const & a, \
        cado_details::coeff_proxy<polynomial_details::polynomial<T>> const & b) \
        -> decltype(OP(T(a), T(b)))                                        \
    {									\
        return OP(T(a), T(b));						\
    }

CADO_COEFF_PROXY_OVERLOAD(operator*)
CADO_COEFF_PROXY_OVERLOAD(operator+)
CADO_COEFF_PROXY_OVERLOAD(operator-)
CADO_COEFF_PROXY_OVERLOAD(operator/)
CADO_COEFF_PROXY_OVERLOAD(operator==)
CADO_COEFF_PROXY_OVERLOAD(operator!=)
CADO_COEFF_PROXY_OVERLOAD(operator<=)
CADO_COEFF_PROXY_OVERLOAD(operator>=)
CADO_COEFF_PROXY_OVERLOAD(operator<)
CADO_COEFF_PROXY_OVERLOAD(operator>)
#else
#define CADO_COEFF_PROXY_OVERLOAD_META(CLASS, OP)			\
    template<typename T, typename U>					\
    static inline auto OP(					\
        U const & a,							\
        cado_details::CLASS<polynomial_details::polynomial<T>> const & b) \
        -> decltype(OP(a, T(b)))                                        \
    {									\
        return OP(a, T(b));						\
    }									\
    template<typename T, typename U>					\
    static inline auto OP(					        \
        cado_details::CLASS<polynomial_details::polynomial<T>> const & a, \
        U const & b)							\
        -> decltype(OP(T(a), b))                                        \
    {									\
        return OP(T(a), b);						\
    }									\
    template<typename T, typename U>					\
    static inline auto OP(					        \
        cado_details::CLASS<polynomial_details::polynomial<T>> const & a, \
        cado_details::CLASS<polynomial_details::polynomial<U>> const & b) \
        -> decltype(OP(T(a), U(b)))                                       \
    {									\
        return OP(T(a), U(b));						\
    }

#define CADO_COEFF_PROXY_OVERLOAD(OP)				        \
        CADO_COEFF_PROXY_OVERLOAD_META(coeff_proxy, OP)                 \
        CADO_COEFF_PROXY_OVERLOAD_META(const_coeff_proxy, OP)

CADO_COEFF_PROXY_OVERLOAD(operator*)
CADO_COEFF_PROXY_OVERLOAD(operator+)
CADO_COEFF_PROXY_OVERLOAD(operator-)
CADO_COEFF_PROXY_OVERLOAD(operator/)
CADO_COEFF_PROXY_OVERLOAD(operator==)
CADO_COEFF_PROXY_OVERLOAD(operator!=)
CADO_COEFF_PROXY_OVERLOAD(operator<=)
CADO_COEFF_PROXY_OVERLOAD(operator>=)
CADO_COEFF_PROXY_OVERLOAD(operator<)
CADO_COEFF_PROXY_OVERLOAD(operator>)
#endif

namespace polynomial_details {
template<typename T>
struct polynomial;

/* forward-declare the template. We need it in order to be able to
 * declare it as a friend of the polynomial<>struct
 */
template<typename T> std::istream& operator>>(std::istream& in, polynomial_details::polynomial<T> & F);
template<typename T> std::ostream& operator<<(std::ostream& o, polynomial_details::polynomial<T> const & f);
template<typename T> std::istream& operator>>(std::istream& in, cado::named_proxy<polynomial_details::polynomial<T> &> const & F);
template<typename T> std::ostream& operator<<(std::ostream& o, cado::named_proxy<polynomial_details::polynomial<T> const &> const & f);

/* Things that we could add (all are in the mpz_poly interface):
 *
 * - {mul,div}_xi
 * - more advanced operator*
 * - evaluation on polynomials (either with extra code or with an OR in
 *   the enable_if, that might work too).
 * - is_monomial_multiple / is_monomial
 * - is_monic / makemonic
 * - translation, rotation
 * - infinity norm
 * - homography
 * - discriminant
 * - coproduct_tree / prod / polyfromroots / multievaluate / interpolate
 *
 * and for integral coefficients:
 * - modular interface, including mod_f_mod_mpz, etc.
 * - divexact
 * - to_monic
 * - pseudodiv
 * - gcd, xgcd, content
 * - is_square_free
 * - polynomial factorization mod p ; factor_and_lift
 * - I don't understand the mpz_poly_base thing...
 *
 */

template<typename CoefficientType, typename PointType>
struct eval_type {
    static constexpr bool up = cado_math_aux::is_coercible<CoefficientType, PointType>::value;
    static constexpr bool down = cado_math_aux::is_strictly_coercible<PointType, CoefficientType>::value;
    using type = std::conditional_t<up, PointType, CoefficientType>;
};


template<typename CoefficientType, typename PointType>
using eval_type_t = typename eval_type<CoefficientType, PointType>::type;

/* each instance of a polynomial object with coefficient type T inherits
 * from number_context<T>. In most cases this is a trivial
 * layer with the identity function as operator(). With cxx_mpfr
 * however, which has variable precision, we want to attach a precision
 * to the polynomial, even before it has any coefficient. This will make
 * it possible to create coefficients such as 0 or 1 in the polynomial
 * with the correct precision.
 */

using cado::number_context;

template<typename T>
struct polynomial : public number_context<T>
{
    number_context<T> const & ctx() const { return *this; }
    number_context<T> & ctx() { return *this; }

    static constexpr int number_of_variables = 1;

    private:
    std::vector<T> coeffs;

    public:

    using coefficient_type = T;

    int degree() const {
        return runtime_numeric_cast<int>(coeffs.size())-1;
    }

    int valuation() const {
        for(int i = 0 ; i <= degree() ; i++)
            if (coeffs[i] != 0)
                return i;
        return INT_MAX;
    }

    unsigned int size() const { return coeffs.size(); }

    T lc() const { ASSERT_ALWAYS(!coeffs.empty()); return coeffs.back(); }

    polynomial() = default;
    explicit polynomial(number_context<T> const & ctx)
        : number_context<T>(ctx)
    {}
    ~polynomial() = default;
    polynomial(polynomial const&) = default;
    polynomial(polynomial &&) = default;
    polynomial& operator=(polynomial const&) = default;
    polynomial& operator=(polynomial &&) = default;

    /* all instantations love each other */
    template<typename U> friend struct polynomial;
    template<typename U>
        explicit polynomial(polynomial<U> const & a, number_context<T> const & tr)
        requires (!std::is_same_v<U, T>)
        : number_context<T>(tr)
        , coeffs { a.coeffs.begin(), a.coeffs.end() }
    {}
    template<typename U>
        explicit polynomial(polynomial<U> const & a)
        requires (!std::is_same_v<U, T>)
        : polynomial(a, number_context<T>(a.ctx()))
    {}

    void set_zero() { coeffs.clear(); }

    void set_xi(unsigned int i) {
        coeffs.assign((i+1), ctx()(0));
        coeffs[i] = ctx()(1);
    }

    polynomial& operator=(T v) {
        coeffs.clear();
        ctx() = number_context<T>(v);
        (*this)[0] = v;
        return *this;
    }
    explicit polynomial(T v)
        : number_context<T>(v)
        , coeffs(1, v)
    {
        cleandeg();
    }
    polynomial(std::initializer_list<T> l) 
        : coeffs(l.begin(), l.end())
    {
        if (l.begin() != l.end())
            ctx() = number_context<T>(*l.begin());
    }

    explicit operator cxx_mpz_poly() const {
        cxx_mpz_poly res;
        for(unsigned int i = 0 ; i < size() ; i++) {
            /* coeffs[i] may be float, double, or long double */
            mpz_poly_setcoeff(res, i, cado_math_aux::mpz_from<T>(coeffs[i]));
        }
        return res;
    }

    // NOLINTNEXTLINE(hicpp-explicit-conversions)
    explicit polynomial(std::string const & e, number_context<T> const & tr = {})
        : polynomial(tr)
    {
        std::istringstream is(e);
        if (!(operator>>(is, *this)))
            throw cado::parse_error();
    }
    // NOLINTNEXTLINE(hicpp-explicit-conversions)
    explicit polynomial(const char * e, number_context<T> const & tr = {})
        : polynomial(std::string(e), tr)
        {}

    private:
    void cleandeg(int deg) {
        ASSERT_ALWAYS(deg >= -1);
        ASSERT_ALWAYS((unsigned int) (deg + 1) <= size());
        for( ; deg >= 0 && coeffs[deg] == 0 ; deg--);
        coeffs.erase(coeffs.begin() + (deg + 1), coeffs.end());
    }
    void cleandeg() { cleandeg(degree()); }

    public:

    friend struct cado_details::coeff_proxy<polynomial>;
    friend struct cado_details::const_coeff_proxy<polynomial>;

    cado_details::coeff_proxy<polynomial> operator[](unsigned int i)
    { return { *this, i }; }
    cado_details::const_coeff_proxy<polynomial> operator[](unsigned int i) const
    { return { *this, i }; }

    /************** {{{ evaluation **************/
    /* {{{ evaluation at a point */

    private:

    template<typename U>
    number_context<U>
        number_context_for_evaluation(U const & x) const
        requires cado_math_aux::is_coercible_v<T, U>
    {
        return number_context<U>(x);
    }
    template<typename U>
    number_context<T> const &
        number_context_for_evaluation(U const &) const
        requires cado_math_aux::is_strictly_coercible_v<U, T>
    {
        /* We need to use the precision of the coefficients of the
         * polynomial, and not the default precision!
         * See also comment above. We call eval_with_reference<U>, but
         * eval_type_t<T, U> is T */
        return ctx();
    }


    /* Note that here, we do not ask that U is
     * eval_type_t<T, U>. This is because if T is larger (say T is
     * cxx_mpfr and U is int), we have obvious potential for better code
     * if we don't cast x to type T right away.
     */
    template<typename U>
    eval_type_t<T, U> eval_with_reference(
            number_context<eval_type_t<T, U>> const & tr,
            U const & x,
            std::function<void(eval_type_t<T, U> const &)> const & h
                = [](auto){}) const
    {
        /* h is used as a way to handle to coefficients that are produced
         * along the computation. They happen to be the coefficients of
         * the quotient polynomial, so div_linear actually uses them.
         */
#ifdef  __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"
#endif
        using cado_math_aux::fma;

        if (degree() < 0) return tr(0);

        T const * f = coeffs.data();

        auto s = tr(f[degree()]);
        auto xx = tr(x);

        int k = degree();
        using cado_math_aux::fma;
        switch(k-1) {
            case 8: h(s); s = fma(s, xx, tr(f[8])); no_break();
            case 7: h(s); s = fma(s, xx, tr(f[7])); no_break();
            case 6: h(s); s = fma(s, xx, tr(f[6])); no_break();
            case 5: h(s); s = fma(s, xx, tr(f[5])); no_break();
            case 4: h(s); s = fma(s, xx, tr(f[4])); no_break();
            case 3: h(s); s = fma(s, xx, tr(f[3])); no_break();
            case 2: h(s); s = fma(s, xx, tr(f[2])); no_break();
            case 1: h(s); s = fma(s, xx, tr(f[1])); no_break();
            case 0: h(s); s = fma(s, xx, tr(f[0])); no_break();
            case -1: break;
            default:
                    for ( ; k-- ; ) {
                        h(s);
                        s = fma(s, xx, tr(f[k]));
                    }
        }
        return s;
#ifdef  __GNUC__
#pragma GCC diagnostic pop
#endif
    }

    public:

    template<typename U> eval_type_t<T, U> eval(U const & x) const
    {
        /* Note that when T<U, we need to use the precision of the
         * coefficients of the polynomial, and not the default precision!
         * See also comment above. In that case, we call
         * eval_with_reference<U>, but eval_type_t<T, U> is T */
        return eval_with_reference(number_context_for_evaluation(x), x);
    }
    template<typename U> eval_type_t<T, U> operator()(U const & x) const
    {
        return eval(x);
    }

    /* }}} */

    /* {{{ Evaluate the homogenous polynomial induced by f at the pair
     * (x,y). That is, compute the sum f[i]*x^i*y^(n-i), where n is
     * degree(f)
     */
    template<typename U>
    eval_type_t<T, U> eval_with_reference(
            number_context<eval_type_t<T, U>> const & tr,
            U const & x,
            U const & y) const
    {
        using cado_math_aux::addmul;

        if (degree() < 0) return tr(0);

        T const * f = coeffs.data();

        auto s = tr(f[degree()]);
        auto py = tr(y);

        int k = degree();
        switch(k-1) {
            case 2: s *= x; addmul(s, py, tr(f[2])); py *= y;
                    no_break();
            case 1: s *= x; addmul(s, py, tr(f[1])); py *= y;
                    no_break();
            case 0: s *= x; addmul(s, py, tr(f[0])); py *= y;
                    no_break();
            case -1: break;
            default:
                    for ( ; k-- ; ) {
                        s *= x;
                        addmul(s, py, tr(f[k]));
                        /* the last multiplication is useless */
                        if (!k) break;
                        py *= y;
                    }
        }
        return s;
    }

    template<typename U>
    eval_type_t<T, U> eval(U const & x, U const & y) const
    {
        /* We have to assume that x and y have the same precision.
         * Otherwise we have no reason to choose instead of the other in
         * order to determine the working precision
         */
        return eval_with_reference(number_context_for_evaluation(x), x, y);
    }

    template<typename U>
    eval_type_t<T, U> operator()(U const & x, U const & y) const
    { return eval(x, y); }
    /* }}} */

    /* uses arbitrary precision */
    template<typename U>
    eval_type_t<T, U> eval_safe(U const & x) const
        requires cado_math_aux::is_real_v<eval_type_t<T, U>>
    {
        using cado_math_aux::exact_form;
        using cado_math_aux::ldexp;

        auto tr = number_context_for_evaluation(x);

        T const * f = coeffs.data();
        const int deg = degree();
        cxx_mpz xm; int xe; exact_form(xm, xe, x);
        /* We want to evaluate at xz * 2^xe */
        cxx_mpz vm; int ve; exact_form(vm, ve, f[deg]);
        for(int k = deg ; k-- ; ) {
            /* multiply by x = xm*2^xe */
            vm *= xm;
            ve += xe;
            /* add f[k] = fm*2^fe */
            cxx_mpz fm; int fe; exact_form(fm, fe, f[k]);
            if (fe < ve) {
                mpz_mul_2exp (vm, vm, ve - fe);
                ve = fe;
            } else {
                mpz_mul_2exp (fm, fm, fe - ve);
            }
            mpz_add (vm, vm, fm);
        }
        return ldexp(tr(vm), ve);
    }

    /* return the quotient of the division by x-r */
    private:
    template<typename U>
    std::pair<polynomial<eval_type_t<T, U>>, eval_type_t<T, U>>
    div_qr_xminusr_with_reference(
            number_context<eval_type_t<T, U>> const & tr,
            U const & r) const
    {
        /* we really want to set q to zero, however we want to import the
         * number context from e */
        using E = eval_type_t<T, U>;
        polynomial<E> q(tr);
        if (coeffs.empty())
            return { q, tr(0) };
        q.coeffs.assign(coeffs.size() - 1, tr(0));
        size_t k = coeffs.size()-1;
        std::function<void(E const &)> reverse_inserter = [&](E const & u) {
            ASSERT_ALWAYS(k);
            q[--k] = u;
        };
        auto ev = eval_with_reference(tr, r, reverse_inserter);
        ASSERT_ALWAYS(k == 0);
        return { q, ev };
    }
    public:
    template<typename U>
    std::pair<polynomial<eval_type_t<T, U>>, eval_type_t<T, U>>
    div_qr_xminusr(U const & r) const
    {
        return div_qr_xminusr_with_reference(
                number_context_for_evaluation(r),
                r);
    }
    template<typename U>
    polynomial<eval_type_t<T, U>> div_q_xminusr(U const & r) const
    {
        auto const [ q, e ] = div_qr_xminusr(r);
        return q;
    }
    /* }}} */


    /************** {{{ (real) root finding **************/

    /* care is taken to provide the functionality to compute roots in a
     * domain which is not restricted to the coefficient domain. To this
     * end, we can pass a number_context object.
     */
    private:

    /* {{{ Return a bound on the positive roots of p.
     * Assume the leading coefficient of p is positive, then for a
     * positive root r we have
     * p[d]*r^d + ... + p[1]*r + p[0] = 0 thus
     * p[d]*r^d <= -p[d-1]*r^(d-1) - ... - p[1]*r - p[0]
     * <= max(-p[d-1],0)*r^(d-1) + ... + max(-p[1],0)*r + max(-p[0],0)
     * thus q(r) <= 0 where q is the degree-d polynomial formed from p as follows:
     * q[d] = p[d]
     * q[i] = p[i] if p[i] < 0, and 0 otherwise for 0 <= i < d.
     * Since q has a unique positive root, say r0, and q(r) < 0 iff r < r0,
     * then positive roots of p are bounded by r0.
     *
     * More generally, if s in {-1,+1} is such that s*p[d] > 0, and we're
     * looking for a bound on roots of sign t in {-1,+1} (thus a bound on
     * the positive roots of p(tx), we have:
     *
     * s*p[d]*(tr)^d <= -s*p[d-1]*t*(tr)^(d-1) ... - s*p[1]*t^(d-1)*(tr) - s*p[0]*t^d
     * So if we let q[d] = s*p[d]
     * and q[i] = -min(s*p[i]*t^(d-i), 0)
     * We then have q(tr) > 0 for the bound r we're after.
     *
     * The question of what to store in q[i] then opens the question of
     * deciding whether
     *  -s * p[i] * t^(d-i) < 0
     *  (-1) * s * t^(d-i) * sgn(p[i]) < 0
     *
     * We can ignore the case p[i] == 0 since no matter what we do, we'll
     * end up setting q[i] = 0. So this simplifies as
     *
     *  1 + (lc() < 0) + (d-i) & (b < 0) + (p[i] < 0) is odd
     *  (lc() < 0) ^ (d-i) & (b < 0) ^ (p[i] < 0) == 0
     *  (lc() < 0) ^ (p[i] < 0) == (d-i) & (b < 0)
     *
     */

    T bound_positive_roots() const
    {
        const int d = degree();
        const int s = lc() < 0;
        polynomial q(ctx());
        q.coeffs.assign(coeffs.size(), ctx()(0));
        for(int i = 0 ; i < d ; i++) {
            T v = (s /* ^ (negative & (d-i)) */) ? -coeffs[i] : coeffs[i];
            /* simplifies to v==(-1)^s*coeffs[i] < 0 if negative == 0 */
            if (v < 0)
                q[i] = v;
        }
        q[degree()] = cado_math_aux::abs(lc());
        T b = 1;
        for( ; q.eval(b) < 0 ; b = b + b) ;
        return /* negative ? -b : */ b;
    }
    /* }}} */

    /* {{{ findroot_dichotomy
     * fall-back code to find a real root by dichotomy. Since this is an
     * internal function, both ends a and b are assumed to be within the
     * same number context (i.e. if they're both cxx_mpfr's, they have
     * the same precision).
     */
    template<typename U>
    U findroot_dichotomy(U a, U b, int sa) const
    {
        /* it would be an error to instantiate this with U distinct from
         * the real evaluation type */
        static_assert(std::is_same_v<U, eval_type_t<T, U>>);
        static_assert(cado_math_aux::is_real_v<U>);
        U s;
        for(;;) {
            /* The number context carries over from a and b to s */
            s = (a + b) / 2;
            cado_math_aux::do_not_outsmart_me(s);
            if (s == a || s == b) return s;
            using cado_math_aux::sgn;
            if (sgn(eval(s)) * sa > 0)
                a = s;
            else
                b = s;
        }
        return s;
    }
    /* }}} */

    /* {{{ findroot_falseposition
     *
     * assuming g(a)*g(b) < 0, and g has a single root in [a, b],
     * refines that root by the weighted false position method
     * Assumes sa is of same sign as g(a).
     *
     * The code is written with the case a<b in mind, but it may also be
     * called with b<a, in which case we need to adapt a few little
     * things.
     *
     * The appropriate number context is carried by a0 and b0
     */
    template<typename U>
    U findroot_falseposition(U const & a0, U const & b0, U const & pa0) const
    {
        /* it would be an error to instantiate this with U distinct from
         * the real evaluation type */
        static_assert(std::is_same_v<U, eval_type_t<T, U>>);
        static_assert(cado_math_aux::is_real_v<U>);
        int side=0;
        U a=a0, b=b0;

        if (a == b)
            return a;

        int sigma = cado_math_aux::sgn(b-a);

        ASSERT_ALWAYS(sigma*a < sigma*b);

        U pa=pa0, pb = eval(b);

        for(;;) {
            U s = (a*pb-b*pa)/(pb-pa);
            U middle = (a + b) / 2;

            cado_math_aux::do_not_outsmart_me(s);
            cado_math_aux::do_not_outsmart_me(middle);

            /* It may happen that because of overflow, (a*pb-b*pa)/(pb-pa)
             * reaches s==a or s==b too early. If it so happens that we're
             * doing this, while the middle cut doesn't behave this way, use
             * the middle cut instead.
             *
             * Note that almost by design, this countermeasure also cancels
             * some of the benefit of the false position method.
             */
            const bool escapes_range = sigma*s < sigma*a || sigma*s > sigma*b;
            const bool hits_bounds = s == a || s == b;
            const bool middle_cut_is_nice = !(middle == a || middle == b);
            if (escapes_range || (hits_bounds && middle_cut_is_nice))
                s = middle;
            if (s == a || s == b) return s;
            U ps = eval(s);
            using cado_math_aux::sgn;
            if (sgn(ps) * sgn(pa) > 0) {
                a = s; pa = ps;
                if (side==1) pb /= 2;
                side=1;
            } else {
                b = s; pb = ps;
                if (side==-1) pa /= 2;
                side=-1;
            }
            if (cado_math_aux::isnan(b)) {
                return findroot_dichotomy(a0, b0, sgn(pa0));
            }
        }
    }
    /* }}} */

    /* {{{ Descartes' rule of sign, and the derivative sign changes
     * knowing the positive sign changes of the derivative of *this given
     * in v , as well as a bound on the positive roots of *this, store in
     * v the positive roots of *this.  v is clobbered.
     */
    template<typename U>
    void positive_roots_from_derivative_sign_changes(std::vector<U> & v, U bound) const
    {
        /* U carries the number_context, and it would be an error to
         * instantiate this with U distinct from the real evaluation type
         *
         * all items in v are with respect to the same context.
         */
        static_assert(std::is_same_v<U, eval_type_t<T, U>>);
        static_assert(cado_math_aux::is_real_v<U>);
        number_context<U> tr(bound);

        using cado_math_aux::sgn;

        ASSERT_ALWAYS(degree() >= 1);
        if (degree() == 1) {
            /* A linear polynomial has at most one root.
             *
             * We want strictly positive roots here, so we must not
             * consider the case coeffs[0] == 0. On the other hand, the
             * bound counts.
             */
            const int s = sgn(coeffs[0]);
            if (s && s * eval(bound) <= 0)
                v.assign(1, - tr(coeffs[0]) / tr(coeffs[1]));
        } else {
            U a = tr(0);
            U va = tr(coeffs[0]);
            v.push_back(bound);
            size_t m = 0;
            /* If f(a)*f'(a+epsilon) > 0, we won't find a
             * root in the interval [a,b) (that is, until the sign of f'
             * changes).
             *
             * If f(a)*f'(a+epsilon) < 0, we may.
             *  - If we do, then f(b)*f'(b-epsilon) > 0, and
             *    f(b)*f'(b+epsilon) < 0.
             *  - If we don't, then f(b)*f'(b-epsilon) is still < 0,
             *    and then f(b)*f'(b+epsilon) > 0, so we can skip the
             *    next interval.
             */
            /* if coeffs[1] == 0, we may replace by coeffs[2] */
            /* if coeffs[0] == 0, we have a root at zero which doesn't
             * count as positive
             */
            using cado_math_aux::sgn;
            int c01 = sgn(coeffs[0]) * sgn(coeffs[1] ? coeffs[1] : (sgn(bound) * coeffs[2]));
            bool no_chance = c01 * sgn(bound) > 0 || coeffs[0] == 0;
            for(size_t i = 0 ; i < v.size() ; i++) {
                U b = v[i];
                U vb = eval(b);
                if (no_chance) {
                    no_chance = false;
                } else if (sgn(va) * sgn(vb) < 0) {
                    v[m++] = findroot_falseposition(a, b, va);
                } else {
                    /* we're in the case where this interval _looked_
                     * promising, and yet had no root. The next one
                     * certainly won't work.
                     */
                    no_chance = true;
                }
                a = b;
                va = vb;
            }
            v.erase(v.begin() + m, v.end());
        }
    }
    /* }}} */

    /* {{{ Descartes' rule of sign: wrapping it up */
    template<typename U>
    std::vector<U> positive_roots_inner(U bound) const
    {
        /* U carries the number_context, and it would be an error to
         * instantiate this with U distinct from the real evaluation type
         */
        static_assert(std::is_same_v<U, eval_type_t<T, U>>);
        static_assert(cado_math_aux::is_real_v<U>);

        const int d = degree();

        /* The roots of the zero polynomial are ill-defined. Bomb out */
        ASSERT_ALWAYS(d>=0);

        /* Handle constant polynomials separately */
        if (d == 0)
            return {}; /* Constant non-zero poly -> no roots */

        std::vector<polynomial> dg;     /* derivatives of *this */
        dg.reserve(d);
        dg.push_back(*this);
        for(int k = 1 ; k < d ; k++)
            dg.push_back(dg.back().derivative());

        /* work from the most derived polynomial, down to f.
         *
         * dg[d-1] is the (d-1)-th derivative, which is a linear
         * polynomial.
         */
        std::vector<U> res;
        for (int k = d; k-- ; )
            dg[k].positive_roots_from_derivative_sign_changes(res, bound);

        return res;
    }
    /* }}} */

    public:

    /* public interfaces for real root finding.
     *
     * We posit that there is no point in asking for integral roots of a
     * floating point polynomials, so we only allow *roots<U> when T is
     * coercible to U (and U is a real type). An optional number context
     * can be passed.
     *
     * To call these functions on a polynomial f with coefficient type T
     * in order to compute roots in a (larger) domain U, you have two choices.
     *  - either call f.template<U> roots(). For example if f is
     *  polynomial<int>, calling f.template<double> roots() can make
     *  sense.
     *  - or call f.roots(tr), where tr is an object of type
     *  number_context<U>. This is the preferred option when eval_type<U>
     *  is cxx_mpfr.
     *
     */
    /* {{{ positive_roots */
    template<typename U = T>
    std::vector<U> positive_roots(number_context<U> const & tr = {}) const
        requires (
        cado_math_aux::is_real_v<U> &&
        cado_math_aux::is_coercible_v<T, U>)
    {
        return positive_roots_inner(tr(bound_positive_roots()));
    }
    /* }}} */
    /* {{{ roots (positive and negative) */
    template<typename U>
    std::vector<U> roots(number_context<U> const & tr = {}) const
        requires (
        cado_math_aux::is_real_v<U> &&
        cado_math_aux::is_coercible_v<T, U>)
    {
        if (degree() == -1) return {};

        auto w = mirror().positive_roots(tr);
        for(auto & x : w) x = -x;

        std::reverse(w.begin(), w.end());

        if (coeffs[0] == 0)
            w.push_back(tr(0));

        const auto positive = positive_roots(tr);

        w.insert(w.end(), positive.begin(), positive.end());

        return w;
    }
    /* }}} */
    /* {{{ positive roots in ]0,s]. We can't use positive_roots_inner
     * because we want to relax the type requirement on the bound
     * argument. */
    template<typename U, typename B>
    std::vector<U> positive_roots(B const & bound, number_context<U> const & tr = {}) const
        requires (
        cado_math_aux::is_real_v<U> &&
        cado_math_aux::is_coercible_v<T, U>)
    {
        return positive_roots_inner(tr(bound));
    }
    /* }}} */
    /* {{{ roots in [-s,s]. */
    template<typename U, typename B>
    std::vector<U> roots(B const & bound, number_context<U> const & tr = {}) const
        requires (
        cado_math_aux::is_real_v<U> &&
        cado_math_aux::is_coercible_v<T, U>)
    {
        if (degree() == -1) return {};

        auto w = mirror().positive_roots_inner(tr(bound));
        for(auto & x : w) x = -x;

        std::reverse(w.begin(), w.end());

        if (coeffs[0] == 0)
            w.push_back(tr(0));

        const auto positive = positive_roots_inner(tr(bound));

        w.insert(w.end(), positive.begin(), positive.end());

        return w;
    }
    /* }}} */

    /* }}} */

    /************** {{{ (complex) root finding **************/
    private:
    double easy_root_lower_bound() const
        requires cado_math_aux::is_complex_v<T>
    {
        /* computes a lower bound on the moduli of the non-zero complex
         * roots of a polynomial p(x). N(x) is a polynomial whose i_th
         * coefficient is the modulus of the i_th coeffcient of p(x) but
         * whose constant term is negative. The lower bound is the
         * (unique) positive root x of N(x)
         *
         * Proof: assume x is a root. Then f(x) = fn x^n + ... + f0 = 0,
         * so |f0| <= \sum_{i=1}^n |f_i| |x|^i
         *
         * so N(|x|) >= 0, with N(t) = |f_n| t^n + ... + |f_1| |t| - |f_0|.
         *
         * There is only one positive root t0 of N, so we must have |x|
         * >= t0.
         */

        using cado_math_aux::abs;
        ASSERT_ALWAYS(!coeffs.empty());
        polynomial<double> N;
        size_t v = valuation();
        N.coeffs.assign(size() - v, 0);
        for(size_t i = 1 ; v + i < size() ; i++)
            N.coeffs[i] = double(abs(coeffs[v + i]));
        N.coeffs[0] = -double(abs(coeffs[v]));
        auto p = N.positive_roots();
        ASSERT_ALWAYS(p.size() == 1);
        return p[0];
    }

    double easy_root_upper_bound() const
        requires cado_math_aux::is_complex_v<T>
    {
        return 1 / reciprocal().easy_root_lower_bound();
    }

    double cauchy_upper_bound() const
        requires cado_math_aux::is_complex_v<T>
    {
        using cado_math_aux::abs;
        ASSERT_ALWAYS(!coeffs.empty());
        double vn = double(abs(coeffs.back()));
        double z = 0;
        for(size_t i = 0 ; i + 1 < size() ; i++)
            z = std::max(z, double(abs(coeffs[i])) / vn);
        z += 1;
        return z;
    }

    double cauchy_lower_bound() const
        requires cado_math_aux::is_complex_v<T>
    {
        using cado_math_aux::abs;
        ASSERT_ALWAYS(!coeffs.empty());
        return 1 / reciprocal().cauchy_upper_bound();
    }

    double lagrange_upper_bound() const
        requires cado_math_aux::is_complex_v<T>
    {
        using cado_math_aux::abs;
        ASSERT_ALWAYS(!coeffs.empty());
        double vn = double(abs(coeffs.back()));
        double z = 0;
        for(size_t i = 0 ; i + 1 < size() ; i++)
            z += double(abs(coeffs[i])) / vn;
        z = std::max(1.0, z);
        return z;
    }

    double lagrange_lower_bound() const
        requires cado_math_aux::is_complex_v<T>
    {
        using cado_math_aux::abs;
        ASSERT_ALWAYS(!coeffs.empty());
        return 1 / reciprocal().lagrange_upper_bound();
    }

    double weird_lower_bound_probably_wrong() const
        requires cado_math_aux::is_complex_v<T>
    {
        using cado_math_aux::abs;
        using cado_math_aux::log;
        using cado_math_aux::exp;

        const auto n = degree();

        polynomial<double> N;
        N.coeffs.reserve(size());
        for(auto const & c : coeffs)
            N.coeffs.emplace_back(double(cado_math_aux::abs(c)));
        N.coeffs[n] = -N.coeffs[n];

        /* compute upper estimate of bound: assume all the
         * middle terms of N are zero.
         *
         * FIXME: this looks just... wrong. And anyway, it's very
         * misguided.
         */

        double xmax = exp((log(-N[n]) - log(N[0])) / n);

        /* if ignoring the nonlinear terms of N produces
           a smaller root, use that instead */

        if (N[n - 1] != 0)
            xmax = std::min(xmax, -N[n] / N[n - 1]);

        /* chop the interval (0, x) until until x is about
           to make norms(x) change sign */

        double x;

        N = N.reciprocal();

        do {
            x = xmax;
            xmax = x / 10;
        } while (N(xmax) > 0.0);

        /* do Newton iteration until x converges to two decimal places */
        double dx = x;
        auto dN = N.derivative();

        while (abs(dx / x) > 0.005) {
            dx = N(x) / dN(x);
            x -= dx;
        }

        return x;
    }
    public:
    /* The "easy" bound is best, although its computation is in fact more
     * expensive than the others. The Cauchy bound is usually second
     * best, and is trivial to compute.
     */
    double lower_bound_complex_roots() const
        requires cado_math_aux::is_complex_v<T>
    {
        return easy_root_lower_bound();
    }
    double upper_bound_complex_roots() const
        requires cado_math_aux::is_complex_v<T>
    {
        return easy_root_upper_bound();
    }
    /* }}} */

    template<typename U>
    std::vector<U> roots(number_context<U> const & tr = {}) const
        requires (
        cado_math_aux::is_complex_v<U> &&
        cado_math_aux::is_coercible_v<T, U>)
    {
        /* Jenkins-Traub solver */

    }



    public:

    polynomial derivative() const
    {
        if (degree() <= 0) return {};
        polynomial df(ctx());
        df.coeffs.reserve(degree() - 1);
        for(int i = 1 ; i <= degree() ; i++)
            df.coeffs.push_back(coeffs[i] * i);
        return df;
    }

    polynomial operator*(T const & a) const
    {
        polynomial h(ctx());
        h.coeffs.reserve(coeffs.size());
        for(auto const & x : coeffs)
            h.coeffs.push_back(x * a);
        return h;
    }

    polynomial pow(unsigned long n) const
    {
        if (n == 0)
            return { ctx()(1) };
        if (degree() < 0)
            return {};

        unsigned long k = ((~0UL) >> 1) + 1;
        for (; k > n; k >>= 1)
            ;
        auto B = *this;
        for (; k >>= 1;) {
            B *= B;
            if (n & k)
                B *= *this;
        }
        return B;
    }

    polynomial operator/(T const & a) const
    {
        polynomial h(ctx());
        h.coeffs.reserve(coeffs.size());
        for(auto const & x : coeffs)
            h.coeffs.push_back(x / a);
        return h;
    }

    polynomial operator-() const
    {
        polynomial h(ctx());
        h.coeffs.reserve(coeffs.size());
        for(auto const & x : coeffs)
            h.coeffs.push_back(-x);
        return h;
    }

    polynomial operator*(polynomial const & g) const
    {
        polynomial const & f = *this;
        if (f == 0 || g == 0)
            return {};
        polynomial h(ctx());
        h.coeffs.assign(f.size() + g.size() - 1, ctx()(0));
        for(unsigned int i = 0 ; i < f.size() ; i++)
            for(unsigned int j = 0 ; j < g.size(); j++)
                h.coeffs[i+j] += f.coeffs[i] * g.coeffs[j];
        return h;
    }

    polynomial operator+(polynomial const & g) const
    {
        polynomial const & f = *this;
        polynomial h(ctx());
        h.coeffs.assign(std::max(f.size(), g.size()), ctx()(0));
        unsigned int i = 0;
        for( ; i < f.size() && i < g.size() ; i++)
            h.coeffs[i] = f.coeffs[i] + g.coeffs[i];
        for( ; i < f.size() ; i++)
            h.coeffs[i] = f.coeffs[i];
        for( ; i < g.size() ; i++)
            h.coeffs[i] = g.coeffs[i];
        h.cleandeg();
        return h;
    }

    polynomial operator-(polynomial const & g) const
    {
        polynomial const & f = *this;
        polynomial h(ctx());
        h.coeffs.assign(std::max(f.degree(), g.degree()) + 1, ctx()(0));
        unsigned int i = 0;
        for( ; i < f.size() && i < g.size() ; i++)
            h.coeffs[i] = f.coeffs[i] - g.coeffs[i];
        for( ; i < f.size() ; i++)
            h.coeffs[i] = f.coeffs[i];
        for( ; i < g.size() ; i++)
            h.coeffs[i] = -g.coeffs[i];
        h.cleandeg();
        return h;
    }

    /* all compound operators are done lazily, at least for now */
    polynomial& operator*=(T const & a) {
        return (*this) = (*this) * a;
    }

    polynomial& operator/=(T const & a) {
        return (*this) = (*this) / a;
    }

    polynomial& operator+=(polynomial const & g) {
        return (*this) = (*this) + g;
    }

    polynomial& operator-=(polynomial const & g) {
        return (*this) = (*this) - g;
    }

    polynomial& operator*=(polynomial const & g) {
        return (*this) = (*this) * g;
    }

    polynomial& addmul(polynomial const & a, polynomial const & b)
    {
        return (*this) += a*b;
    }

    polynomial& submul(polynomial const & a, polynomial const & b)
    {
        return (*this) -= a*b;
    }

    polynomial& addmul(polynomial const & a, T const & v)
    {
        return (*this) += a*v;
    }

    polynomial& submul(polynomial const & a, T const & v)
    {
        return (*this) -= a*v;
    }

    polynomial reciprocal() const
    {
        polynomial h(ctx());
        h.coeffs.reserve(coeffs.size());
        for(unsigned int i = 0 ; i < size() ; i++)
            h[size()-1-i] = coeffs[i];
        return h;
    }

    polynomial mirror() const
    {
        /* return f(-x) */
        polynomial m;
        m.coeffs.reserve(coeffs.size());
        int s = 1;
        for(auto const & c : coeffs) {
            m.coeffs.push_back(c * s);
            s = -s;
        }
        return m;
    }

    /* This computes u = scale^deg(f) * f(x/scale) (i.e., scale by
     * 1/scale). It is not the same as double_poly_scale, deleted in commit
     * c480fe82174a9de96e1cd35b2317fdf0de3678ab
     */
    polynomial inverse_scale(T const & scale) const
    {
        unsigned int sz = size();
        if (!sz)
            return {};
        polynomial h(ctx());
        h.coeffs.reserve(sz);
        h[--sz] = lc();
        for(T s = scale ; sz-- ; s *= scale) h[sz] = coeffs[sz] * s;
        return h;
    }

    std::string print(std::string const& var = "x") const
    {
        std::ostringstream os;
        if (degree() < 0) os << "0";
        for(int i = 0 ; i <= degree() ; i++) {
            T const & fi = coeffs[i];
            int const r = (fi > 0) - (fi < 0);
            if (r == 0) continue;
            /* We do not want to write "+-" ; in most cases, checking the
             * sign (r) is enough to know whether the printed
             * representation will start by - or not. Alas, for complex
             * types, it's more subtle.
             */
            auto s = fmt::format("{}", fi);
            ASSERT_ALWAYS(!s.empty());
            if (!os.str().empty() && s.front() != '-')
                os << "+";
            if (i == 0) {
                os << s;
            } else {
                if (fi == -1) {
                    os << "-";
                } else if (fi != 1) {
                    os << s << "*";
                }
                os << var;
                if (i > 1) os << "^" << i;
            }

        }
        return os.str();
    }

    explicit polynomial(cxx_mpz_poly const & f, number_context<T> const & tr = {})
        : number_context<T>(tr)
    {
        coeffs.assign(f.degree() + 1, ctx()(0));
        for(int i = 0 ; i <= f.degree() ; i++) {
            coeffs[i] = ctx()(mpz_poly_coeff_const(f, i));
        }
    }

    private:
    friend std::istream& operator>><T>(std::istream& in, cado::named_proxy<polynomial &> const & F);

    struct parser_traits {
        std::string x;
        explicit parser_traits(std::string x) : x(std::move(x)) {}
        struct unexpected_literal : public cado::parse_error {
            std::string msg;
            explicit unexpected_literal(std::string const & v)
                : msg(std::string("unexpected literal " + v))
            {}
            const char *what() const noexcept override {
                return msg.c_str();
            }
        };
        static constexpr const int accept_literals = 1;
        using type = polynomial;
        using number_type = T;
        void add(polynomial & c, polynomial const & a, polynomial const & b) {
            c = a + b;
        }
        void sub(polynomial & c, polynomial const & a, polynomial const & b) {
            c = a - b;
        }
        void neg(polynomial & c, polynomial const & a) {
            c = -a;
        }
        void mul(polynomial & c, polynomial const & a, polynomial const & b) {
            c = a * b;
        }
        void pow_ui(polynomial & c, polynomial const & a, unsigned int e) {
            c = a.pow(e);
        }
        void swap(polynomial & a, polynomial & b) {
            std::swap(a, b);
        }
        void set(polynomial & c, T const & z) {
            c = z;
        }
        void set_literal_power(polynomial & a, std::string const & v, unsigned long e) {
            if (v == x)
                a.set_xi(e);
            else
                throw unexpected_literal(v);
        }
    };
    public:

    cado::named_proxy<polynomial &> named(std::string const & x) {
        return { *this, x };
    }
    cado::named_proxy<polynomial const &> named(std::string const & x) const {
        return { *this, x };
    }

    polynomial(std::string const &, std::string const & var);

    private:
    int spaceship(polynomial const & f) const
    {
        int r = (degree() > f.degree()) - (f.degree() > degree());
        if (r) return r;
        for(int i = 0 ; i <= degree() ; i++) {
            T v = (*this)[i];
            r = (v > f[i]) - (f[i] > v);
            if (r) return r;
        }
        return 0;
    }
    public:
    int operator<=>(polynomial const & f) const { return spaceship(f); }
    bool operator==(polynomial const & f) const { return spaceship(f) == 0; }
    // bool operator!=(polynomial const & f) const { return !operator==(f); }
    bool operator<(polynomial const & f) const { return spaceship(f) < 0; }
    bool operator<=(polynomial const & f) const { return spaceship(f) <= 0; }
    bool operator>(polynomial const & f) const { return spaceship(f) > 0; }
    bool operator>=(polynomial const & f) const { return spaceship(f) >= 0; }
    template<typename U>
    bool operator!=(U v) const { return !operator==(v); }
    template<typename U>
        bool
    operator==(U v) const
    requires std::is_convertible_v<U, T>
    {
        return (degree() < 0 && v == 0) || (degree() == 0 && coeffs[0] == v);
    }

    bool has_nan() const {
        auto p = [](T const & x) { return cado_math_aux::isnan(x); };
        return std::any_of(coeffs.begin(), coeffs.end(), p);
    }

    bool has_inf() const {
        auto p = [](T const & x) { return cado_math_aux::isinf(x); };
        return std::any_of(coeffs.begin(), coeffs.end(), p);
    }

    private:

    /*
     * Compute the pseudo division of a and b such that
     *  lc(b)^(deg(a) - deg(b) + 1) * a = b * q + r with deg(r) < deg(b).
     *  See Henri Cohen, "A Course in Computational Algebraic Number Theory",
     *  for more information.
     *
     * Assume that deg(a) >= deg(b) and b is not the zero polynomial.
     *
     * Return true on success, false otherwise.
     *
     * No overlap allowed.
     */
    bool pseudo_division(polynomial * q, polynomial & r,
            polynomial const & b) const
    {
        polynomial const & a = *this;
        ASSERT(a.degree() >= b.degree());
        ASSERT(b.degree() != -1);

        int const m = a.degree();
        int const n = b.degree();
        T d = b.lc();
        int e = m - n + 1;
        polynomial s(ctx());

        if (q) *q = 0;

        r = a;

        while (r.degree() >= n) {
            s = ctx()(0);
            s[r.degree() - n] = r.lc();

            if (q) {
                *q *= d;
                *q += s;
            }
            r *= d;
            s *= b;
            int const nrdeg = r.degree() - 1;
            r -= s;
            /* We enforce this because the subtraction may miss the
             * cancellation of the leading term due to rounding.
             */
            if (r.degree() > nrdeg)
                r.cleandeg(nrdeg);
            e--;

            /* Cancellations can happen, here. We do have a test case
             * that triggers r===0, in which case we'll probably need to
             * resort to exact computations instead
             */
            if (r.has_nan() || r.has_inf() || r.degree() < 0)
                return false;
        }

        ASSERT(e >= 0);

        d = cado_math_aux::pow(d, static_cast<T>(e));
        if (q) *q *= d;
        r *= d;

        return true;
    }

    bool pseudo_remainder(polynomial & r, polynomial const & b) const
    {
        return pseudo_division(nullptr, r, b);
    }

    public:

    /* XXX caveat: This is not a proper resultant implementation. It
     * implicitly relies on the assumption that the inputs are integer
     * polynomials (but with floating-point representation).
     */
    T resultant(polynomial<T> const & q) const
    requires cado_math_aux::is_real_v<T>
    {
        polynomial const & p = *this;
        if (p.degree() < 0 || q.degree() < 0)
            return 0;

        polynomial a = p;
        polynomial b = q;
        polynomial r(ctx());

        int s = 1;
        int d;

        int pseudo_div = 1;

        T g = ctx()(1);
        T h = ctx()(1);

        if (a.degree() < b.degree()) {
            std::swap(a, b);

            if ((a.degree() % 2) == 1 && (b.degree() % 2) == 1)
                s = -1;
        }

        while (b.degree() > 0) {
            // TODO: verify if it is necessary.
            d = a.degree() - b.degree();

            if ((a.degree() % 2) == 1 && (b.degree() % 2) == 1)
                s = -s;

            pseudo_div = a.pseudo_remainder(r, b);
            if (!pseudo_div)
                break;

            a = b;

            ASSERT(d >= 0);

            b = r / (g * cado_math_aux::pow(h, d));

            g = a.lc();

            ASSERT(d != 0 || h == 1);

            h = cado_math_aux::pow(h, static_cast<T>(d - 1));
            h = cado_math_aux::pow(g, static_cast<T>(d)) / h;
        }

        if (pseudo_div) {
            ASSERT(a.degree() > 0);

            //For now, if b.degree() == -1, pseudo_div == 0.
            if (b.degree() == -1) {
                ASSERT(0);
            } else {

                ASSERT(a.degree() >= 0);

                h = cado_math_aux::pow(
                        static_cast<T>(b[0]), a.degree()) / cado_math_aux::pow(h, (a.degree() - 1));
                h *= s;
            }
        } else {
            // we encountered cancellations, so we need to resort to
            // exact arithmetic.
            // XXX This is only meant to be used for very special cases
            // of polynomials with integer coefficiens to start with, and
            // taking integer values. There doesn't seem to be a good
            // rationale for doing this kind of indiscriminate fallback.
            // TODO: use last version of a and b in pseudo_division.
            cxx_mpz val_z;
            cxx_mpz_poly pz(p);
            cxx_mpz_poly qz(q);
            mpz_poly_resultant(val_z, pz, qz);
            return cado_math_aux::mpz_get<T>(val_z);
        }
        return h;
    }

    T resultant(polynomial<T> const & q) const
    requires std::is_same_v<T, cxx_mpz>
    {
        cxx_mpz val_z;
        cxx_mpz_poly pz(*this);
        cxx_mpz_poly qz(q);
        mpz_poly_resultant(val_z, pz, qz);
        return val_z;
    }

};


static_assert(std::is_same_v<decltype(polynomial<int>{}(double())), double>);
static_assert(std::is_same_v<decltype(polynomial<int>{}(int())), int>);
static_assert(std::is_same_v<decltype(polynomial<cxx_mpz>{}(int())), cxx_mpz>);
static_assert(std::is_same_v<decltype(polynomial<cxx_mpz>{}(double())), double>);

// same idea. not a reason to pull cxx_mpfr.hpp or cxx_mpc.hpp though.
#ifdef HAVE_MPFR
static_assert(std::is_same_v<decltype(polynomial<cxx_mpz>{}(cxx_mpfr())), cxx_mpfr>);
#endif
// static_assert(std::is_same<decltype(polynomial<double>{}(cxx_mpc())), cxx_mpc>::value);


#ifdef HAVE_MPFR
/* This has to be a specific instantiation! */
template<>
inline polynomial<cxx_mpfr>::polynomial(cxx_mpz_poly const & f, number_context<cxx_mpfr> const & tr)
    : number_context<cxx_mpfr>(tr)
{
    coeffs.assign(f.degree() + 1, ctx()(0));
    for(int i = 0 ; i <= f.degree() ; i++)
        mpfr_set_z(coeffs[i], mpz_poly_coeff_const(f, i), MPFR_RNDN);
}
#endif

template<typename T>
std::istream& operator>>(std::istream& in, cado::named_proxy<polynomial<T> &> const & F)
{
    std::string line;
    for(;;in.get()) {
        int const c = in.peek();
        if (in.eof() || !isspace(c)) break;
    }
    if (!getline(in, line)) return in;
    std::istringstream is(line);

    using traits_type = typename polynomial<T>::parser_traits;
    cado_expression_parser<traits_type> P(F.x());
    P.tokenize(is);

    try {
        F.c = P.parse(number_context<T>(F.c));
    } catch (cado::parse_error const & p) {
        in.setstate(std::ios_base::failbit);
        return in;
    }

    return in;
}

template<typename T>
std::istream& operator>>(std::istream& in, polynomial<T> & F)
{
    return operator>>(in, F.named("x"));
}

/* printing needs a way to specify the variables... */
template<typename T>
inline std::ostream& operator<<(std::ostream& o, cado::named_proxy<polynomial<T> const &> const & f)
{
    return o << f.c.print(f.x());
}

/* we do have a default behaviour, though */
template<typename T>
inline std::ostream& operator<<(std::ostream& o, polynomial<T> const & f)
{
    return operator<<(o, f.named("x"));
}

} /* namespace polynomial_details */




template<typename T>
using polynomial = polynomial_details::polynomial<T>;


namespace fmt {
    template<typename T>
    struct formatter<polynomial<T>>: ostream_formatter {};
}


#endif	/* CADO_UTILS_POLYNOMIAL_HPP */
