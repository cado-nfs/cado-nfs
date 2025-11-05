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
#include <climits>

#include <algorithm>
#include <compare>
#include <initializer_list>
#include <ios>
#include <limits>
#include <istream>
#include <ostream>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <functional>
#include <tuple>
#include <array>
#include <stdexcept>
#include <memory>


#include <gmp.h>
#include "fmt/base.h"
#include "fmt/ostream.h"

#include "gmp_aux.h"
#include "macros.h"
#include "runtime_numeric_cast.hpp"
#include "number_context.hpp"
#include "mpz_poly.h"
#include "cxx_mpz.hpp"
#include "cado_math_aux.hpp"
#include "cado_addsubmul.hpp"
#include "cado_expression_parser.hpp"
#include "coeff_proxy.hpp"
#include "cado_type_traits.hpp"
#include "named_proxy.hpp"
#include "extra_complex_overloads.hpp"
#include "cado_constants.hpp"

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
#define CADO_COEFF_PROXY_OVERLOAD_META_X(CLASS1, CLASS2, OP)	\
    template<typename T, typename U>					\
    static inline auto OP(					        \
        cado_details::CLASS1<polynomial_details::polynomial<T>> const & a, \
        cado_details::CLASS2<polynomial_details::polynomial<U>> const & b) \
        -> decltype(OP(T(a), U(b)))                                       \
    {									\
        return OP(T(a), U(b));						\
    }
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
    CADO_COEFF_PROXY_OVERLOAD_META_X(CLASS, CLASS, OP)

#define CADO_COEFF_PROXY_OVERLOAD(OP)				        \
        CADO_COEFF_PROXY_OVERLOAD_META(coeff_proxy, OP)                 \
        CADO_COEFF_PROXY_OVERLOAD_META(const_coeff_proxy, OP)           \
        CADO_COEFF_PROXY_OVERLOAD_META_X(coeff_proxy, const_coeff_proxy, OP) \
        CADO_COEFF_PROXY_OVERLOAD_META_X(const_coeff_proxy, coeff_proxy, OP)

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
        polynomial(polynomial<U> const & a, number_context<T> const & tr)
        : number_context<T>(tr)
    {
        coeffs.reserve(a.coeffs.size());
        for(auto const & c : a.coeffs)
            coeffs.emplace_back(tr(c));
    }
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
        for(auto & c : coeffs)
            c = ctx()(c);
        cleandeg();
    }

    polynomial(std::vector<T> const & l, number_context<T> const & tr)
        : number_context<T>(tr)
        , coeffs(l.begin(), l.end())
    {
        for(auto & c : coeffs)
            c = ctx()(c);
        cleandeg();
    }

    explicit polynomial(std::vector<T> const & l)
        : coeffs(l.begin(), l.end())
    {
        if (l.begin() != l.end())
            ctx() = number_context<T>(*l.begin());
        for(auto & c : coeffs)
            c = ctx()(c);
        cleandeg();
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
    void chop() { coeffs.pop_back(); }

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
        number_context_for_evaluation(number_context<U> const & x) const
        requires cado_math_aux::is_coercible_v<T, U>
    {
        return x;
    }
    template<typename U>
    number_context<T> const &
        number_context_for_evaluation(number_context<U> const &) const
        requires cado_math_aux::is_strictly_coercible_v<U, T>
    {
        /* We need to use the precision of the coefficients of the
         * polynomial, and not the default precision!
         * See also comment above. We call eval_with_reference<U>, but
         * eval_type_t<T, U> is T */
        return ctx();
    }
    template<typename U>
    number_context<eval_type_t<T, U>>
        number_context_for_evaluation(U const & x) const
        {
            return number_context_for_evaluation(number_context<U>(x));
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

    /************** {{{ shift **************/
    /* compute the shift f(x+a) of a polynomial in soft-O(M(n))
     * we need two important auxiliary functions */

    /* Algorithm, broadly.
     *
     * Let f be an input polynomial of degree 2h-1. Compute the division
     * f = q * (x-a)^h + r, so that f(x+a) = q(x+a) * x^h + r(x+a).
     * By recursion on q and r, the result follows easily. Most of the
     * fine details are in the handling of the input degree and the
     * degree of the divisor.
     *
     * Let y=x-a. Let us reason on the recursion tree, and what divisors
     * are needed for the recursion.
     *
     * On an n-coefficient input, at depth i>=0, the tree nodes have either
     * h_i or 1+h_i coefficients, with h_i=floor(n/2^i).
     * We precompute the divisors y^{h_{i+1}} and y^{h_{i+1}+1}
     * 
     * note that two distinct situations may occur:
     *  - if h_i is even, h_i = 2*h_{i+1} and 1+h_i = h_{i+1} + (1+h_{i+1})
     *  - if h_i is odd,  h_i = h_{i+1} + (1+h_{i+1}) and 1+h_i = 2*(1+h_{i+1})
     * 
     * when we have h_i coefficients and h_i is even, we divide by
     * y^{h_{i+1}}, leaving us with two results, each having
     * h_{i+1} coefficients.
     * 
     * when we have 1+h_i coeficients and h_i is odd, we divide by
     * y^{1+h_{i+1}}, leaving us with two results, each
     * having 1+h_{i+1} coefficients.
     * 
     * otherwise we may divide by either divisor, and we'll have one
     * result with h_{i+1} coefficients, and the other with 1+h_{i+1}
     * coefficients.
     *
     */
    private:
    /* our first auxiliary function computes the ladder of divisors from
     * bottom to top, [ (y^{h_i}, y^{1+h_i}) ]. We compute them for i>0,
     * since we don't need them at depth 0. In practice, the input n of
     * this function is one half of the original input length, which
     * gives the same result.
     */
    template<typename Poly>
        static std::vector<std::array<Poly, 2>>
        power_ladder(Poly const & y, size_t n)
        {
            ASSERT_ALWAYS(n);
            int const k = nbits(n);
            std::array<Poly, 2> D { Poly(y.ctx()(1)), y };
            std::vector<decltype(D)> res;
            res.reserve(k+1);
            res.push_back(D);
            for(size_t m = 1 << (k-1) ; m ; m>>=1) {
                if (n & m) {
                    D[0] = D[0] * D[1];
                    D[1] = D[1] * D[1];
                } else {
                    D[1] = D[0] * D[1];
                    D[0] = D[0] * D[0];
                }
                res.push_back(D);
            }
            return res;
        }

    /* our second auxiliary function works down the ladder recursively
     * (from top to bottom), dealing with inputs with h_i or 1+h_i
     * coefficients at depth i.
     */
    template<typename U>
        static polynomial<U> shift_recursion(
                polynomial<U> const & u,
                size_t n,
                std::vector<std::array<polynomial<U>, 2>> const & ladder,
                number_context<U> const & tr,
                int depth = 0)
        {
            /* a polynomial that is constant or zero is easy to shift! */
            if (n <= 1)
                return u;
            ASSERT_ALWAYS(0 <= depth);
            /* within the recursion, degree drops are possible (in
             * remainders only)
             */
            ASSERT_ALWAYS(u.size() <= n);
            auto const & [ d0, d1 ] = ladder[ladder.size()-1-depth];
            size_t m0 = d0.size() - 1;
            size_t m1 = d1.size() - 1;
            ASSERT_ALWAYS(n == 2*m0 || n == (m0 + m1) || n == 2*m1);
            /* if n is 2*m0, we must divide by d0
             * if n is 2*m1, we must divide by d1
             * otherwise we can take either.  */
            auto const & divisor = (n == 2 * m0) ? d0 : d1;
            auto [ q, r ] = u.div_qr(divisor);
            size_t m = divisor.size() - 1;
            q = shift_recursion(q, n - m, ladder, tr, depth+1);
            r = shift_recursion(r, m, ladder, tr, depth+1);
            polynomial<U> res(tr);
            res.coeffs.assign(u.size(), tr(0));
            ASSERT_ALWAYS(q.size() + m == u.size());
            std::move(r.coeffs.begin(), r.coeffs.end(), res.coeffs.begin());
            std::move(q.coeffs.begin(), q.coeffs.end(), res.coeffs.begin() + m);
            /* composition never changes the degree, so cleandeg() is not
             * needed */
            return res;
        }

    public:

    /* compute f(x+a) */
    template<typename U>
        polynomial<eval_type_t<T, U>>
        shift(U const & a) const
        {
            using E = eval_type_t<T, U>;
            auto tr = number_context_for_evaluation(a);
            if (degree() <= 0)
                return { *this, tr };
            polynomial<E> x_minus_a { -tr(a), tr(1) };
            auto ladder = power_ladder(x_minus_a, (degree()+1)/2);
            return shift_recursion(polynomial<E>(*this), degree()+1, ladder, tr);
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

    /************** {{{ multievaluation **************/

    public:
    // {{{ product tree on a set of points
    struct tree {
        polynomial f;
        std::unique_ptr<tree> l {};
        std::unique_ptr<tree> r {};
    };
    template<typename U>
    static std::unique_ptr<tree> product_tree(std::vector<U> const & points, number_context<T> const & tr)
        requires cado_math_aux::is_coercible_v<U, T>
    {
        std::vector<std::unique_ptr<tree>> L;
        L.reserve(points.size());
        for(auto const & z : points) {
            polynomial xz { -tr(z), tr(1) };
            L.emplace_back(std::make_unique<tree>(xz));
        }
        for( ; L.size() > 1 ; ) {
            auto it = L.begin();
            size_t j = 0;
            for( ; j + 1 < L.size() ; ) {
                auto & l = L[j++];
                auto & r = L[j++];
                *it++ = std::make_unique<tree>(l->f * r->f,
                                               std::move(l), std::move(r));
            }
            if (j < L.size())
                *it++ = std::move(L[j]);
            L.erase(it, L.end());
        }
        return std::move(L[0]);
    }
    // }}}

    private:
    void multieval(std::vector<T> & it, tree const & A) const
    {
        ASSERT_ALWAYS(A.f.degree() > 0);
        auto fmod = div_r(A.f);
        ASSERT_ALWAYS(bool(A.l) == bool(A.r));
        if (A.l) {
            fmod.multieval(it, *A.l);
            fmod.multieval(it, *A.r);
        } else {
            ASSERT_ALWAYS(fmod.degree() <= 0);
            it.push_back(fmod[0]);
        }
    }

    public:
    template<typename U>
    std::vector<U>
    multieval(typename polynomial<U>::tree const & A, number_context<U> const & tr) const
    {
        std::vector<U> res;
        res.reserve(A.f.degree());
        polynomial<U>(*this, tr).multieval(res, A);
        return res;
    }

    // {{{ multi-evaluation at a set of points (easy case: just f)

    /* Evaluate f at the given points, and return results interpreted by
     * the number context tr.
     * The set of points must have size exactly degree(f). 
     * If f is constant (or zero), an empty vector is returned.
     */
    template<typename U>
        std::vector<eval_type_t<T, U>>
        multieval(std::vector<U> const & points, number_context<eval_type_t<T, U>> const & tr) const
        {
            using E = eval_type_t<T, U>;
            if (degree() <= 0) {
                ASSERT_ALWAYS(points.empty());
                return {};
            }
            auto A = polynomial<E>::product_tree(points, tr);
            return multieval<U>(A, tr);
        }
    // }}}

    /* }}} */

    /************** {{{ (complex) root finding: bounds **************/
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

    public:
    /* The "easy" bound is best, although its computation is in fact more
     * expensive than the others. The Cauchy bound is usually second
     * best, and is trivial to compute.
     *
     * Note that this is mostly relevant when the mean of the roots is
     * zero.
     */
    double lower_bound_complex_roots() const
        requires cado_math_aux::is_complex_v<T>
    {
        if (degree() < 0)
            throw std::domain_error("roots() is undefined on zero polynomial");
        else if (degree() == 0)
            return 0;
        return easy_root_lower_bound();
    }
    double upper_bound_complex_roots() const
        requires cado_math_aux::is_complex_v<T>
    {
        if (degree() < 0)
            throw std::domain_error("roots() is undefined on zero polynomial");
        else if (degree() == 0)
            return 0;
        return easy_root_upper_bound();
    }
    std::tuple<T, double, double>
    annulus_complex_roots() const
        requires cado_math_aux::is_complex_v<T>
    {
        if (degree() < 0)
            throw std::domain_error("roots() is undefined on zero polynomial");
        else if (degree() == 0)
            return { 0, 0, 0 };
        const int n = degree();
        const T mean = -coeffs[n-1]/coeffs[n]/n;
        auto g = shift(mean);
        const double lo = g.lower_bound_complex_roots();
        const double hi = g.upper_bound_complex_roots();
        return { mean, lo, hi };
    }
    /* }}} */

    /************** {{{ (complex) root finding: Aberth method **************/
    private:
    template<typename U>
    std::vector<eval_type_t<T, U>>
        complex_roots_inner_notreal(number_context<U> const & tr = {}) const
        requires cado_math_aux::is_complex_v<eval_type_t<T, U>>
    {
        using E = eval_type_t<T, U>;
        using Er = decltype(E().real());

        auto Rtr = tr.real();

        std::vector<E> z;
        int const n = degree();
        z.reserve(n);

        auto [ mean, lo, hi ] = polynomial<E>(*this, tr).annulus_complex_roots();
        /* pick n starting points in the annulus */
        for(int i = 0 ; i < n ; i++) {
            using cado_math_aux::pi_v;
            using cado_math_aux::exp;
            E theta { Rtr(0), 2 * pi_v(Rtr) * i / n + Rtr(0.125) };
            z.push_back(exp(theta) * Rtr(lo + i * (hi - lo) / n));
        }

        /* now do the iteration (Aberth-Ehrlich method) */
        for(int iter = 0 ; ; iter++) {
            auto const A = polynomial<E>::product_tree(z, tr);
            auto const f_zk = this->multieval(*A, tr);
            auto const df_zk = derivative().multieval(*A, tr);
            auto const & P = A->f;
            auto const dP = P.derivative();
            auto const ddP = dP.derivative();
            auto const dP_zk = dP.multieval(*A, tr);
            auto const ddP_zk = ddP.multieval(*A, tr);
            int best_bits = INT_MAX;
            for(int i = 0 ; i < n ; i++) {
                auto const a = ddP_zk[i] / dP_zk[i] / 2;
                auto const w = f_zk[i] / (df_zk[i] - f_zk[i] * a);
                /* Weierstrass-Durand-Kerner just does */
                // auto const w = f_zk[i] / dP_zk[i];
                // however we do encounter more frequent convergence
                // failures with WDK than with A-E
                auto const aw = double(cado_math_aux::abs(w));
                auto const az = double(cado_math_aux::abs(z[i]));
                if (aw != 0)
                    best_bits = std::min(best_bits, -std::ilogb(std::abs(aw/az)));
                z[i] -= w;
            }
            if constexpr (std::is_floating_point_v<Er>) {
                if (best_bits >= std::numeric_limits<Er>::digits)
                    break;
#ifdef HAVE_MPC
            } else if constexpr (std::is_same_v<E, cxx_mpc>) {
                if (best_bits >= tr.prec)
                    break;
#endif
            } else {
                /* it should be a static_assert, really. But we still
                 * have a few gcc-10, and those do not correcly elide the
                 * else branch of winning "if constexpr"
                 */
                ASSERT_ALWAYS(false);
            }
            if (iter >= 200) {
                /* I'm basically happy with this code, and it seems to
                 * behave as it should. However I would not feel
                 * comfortable asserting that it does so for every input.
                 * Maybe we want to add some leeway in the low bits. How
                 * much, I don't know.
                 */
                if (std::is_same_v<Er, long double> && cado_math_aux::valgrind_long_double_hopeless()) {
                    fmt::print("# infinite loop in complex root finding ; assuming this is valgrind's fault\n");
                    break;
                } else {
                    throw std::runtime_error("infinite loop in complex root finding");
                }
            }
        }

        return z;
    }

    public:
    template<typename U>
    std::vector<eval_type_t<T, U>>
        roots(number_context<U> const & tr = {}) const
        requires cado_math_aux::is_complex_v<eval_type_t<T, U>>
    {
        using E = eval_type_t<T, U>;
        using Er = decltype(E().real());
        auto Rtr = tr.real();

        std::vector<E> z;
        int const n = degree();
        z.reserve(n);

        if constexpr (cado_math_aux::is_real_v<T> || cado_math_aux::is_integral_v<T>) {
            /* If our input polynomial is of real type, we must make sure
             * that we identify the real roots. It's more painful than it
             * should. In fact, it seems easier to fire up our real root
             * finding code first.
             */
            auto fr = polynomial<Er>(*this, Rtr);
            auto real_roots = fr.roots(Rtr);
            for(auto & r : real_roots)
                z.emplace_back(tr(r));
            if (real_roots.size() < (size_t) degree()) {
                for(auto const & r : real_roots)
                    fr = fr.div_q_xminusr(r);
                for(auto const & r : fr.complex_roots_inner_notreal(tr))
                    z.emplace_back(r);
            }
            return z;
        } else {
            return complex_roots_inner_notreal(tr);
        }
    }
    /* }}} */

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

    /* TODO: templatize */
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

    /* all compound operators are done lazily, at least for now.
     * There would be place to use fma/fms, though. */
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
    /* {{{ parsing */
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
        void pow(polynomial & c, polynomial const & a, unsigned int e) {
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

    /* }}} */
    public:

    cado::named_proxy<polynomial &> named(std::string const & x) {
        return { *this, x };
    }
    cado::named_proxy<polynomial const &> named(std::string const & x) const {
        return { *this, x };
    }

    polynomial(std::string const &, std::string const & var);

    public:
    /* {{{ comparisons and predicates */
    decltype(T() <=> T()) operator<=>(polynomial const & b) const
    {
        using cmp_t = decltype(T() <=> T());
        polynomial<T> const & a(*this);
        if constexpr (std::is_same_v<cmp_t, std::partial_ordering>) {
            if (a.has_nan() || b.has_nan())
                return std::partial_ordering::unordered;
        }
        for (int i = std::max(a.degree(), b.degree()); i >= 0; i--) {
            if (i > b.degree())
                return a.coeffs[i] <=> 0;
            if (i > a.degree())
                return 0 <=> b.coeffs[i];
            if (auto r = a.coeffs[i] <=> b.coeffs[i]; r != 0)
                return r;
        }
        /* We return std::strong_ordering::equial, but the return type is
         * cmp_t, which might be std::partial_ordering. Since there's
         * implicit cast from the former to the latter, we get
         * std::partial_ordering::equivalent whenever we write
         * std::strong_ordering::equal
         */
        return std::strong_ordering::equal;
    }
    bool operator==(polynomial const & f) const {
        /* the degree() check is only a very mild short-circuit for the
         * different degree case: it just removes the nan check, and one
         * coefficient lookup.
         */
        if (degree() != f.degree())
            return false;
        return (*this <=> f) == 0;
    }

    bool operator==(T v) const {
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

    bool is_monic() const {
        return !coeffs.empty() && lc() == 1;
    }

    /* }}} */

    private:

    /* {{{ Compute the pseudo division of a and b such that
     *  lc(b)^(deg(a) - deg(b) + 1) * a = b * q + r with deg(r) < deg(b).
     * See Algorithm 3.1.2 in Cohen, "A Course in Computational Algebraic
     * Number Theory". This is only used in the integer resultant
     * algorithm.
     *
     * If b is monic, this is just Euclidean division.
     *
     * Assume that deg(a) >= deg(b) and b is not the zero polynomial.
     *
     * Return true on success, false otherwise.
     *
     * No overlap allowed.
     */
    private:
    template<typename U>
    polynomial<eval_type_t<T,U>> pseudo_division(
            polynomial<eval_type_t<T,U>> * q,
            polynomial<U> const & b) const
    requires cado_math_aux::is_integral_v<eval_type_t<T, U>>
    {
        ASSERT_ALWAYS(q != this);
        ASSERT(degree() >= b.degree());
        ASSERT(b.degree() != -1);

        using E = eval_type_t<T, U>;
        static_assert(std::is_empty_v<number_context<E>>);

        polynomial<E> r(*this);

        int const m = r.degree();
        int const n = b.degree();
        auto const lb = b.lc();
        int e = m - n + 1;

        if (q) *q = 0;

        while (r.degree() >= n) {
            /* d^(m-n+1-e)*a = b*q + r */
            int const j = r.degree()-n;
            /* we want to subtract lc(r)*x^j*b
             */
            // s = tr(0);
            // s[r.degree() - n] = tr(r.lc());

            /* note that lc(b)*r - lc(r)*b*x^(deg(r)-deg(b)) = lb*r-b*s
             * has degree at most deg(b)-1, so:
             *
             * lb^(m-n+1-(e-1))*a = lb*b*q + lb*r
             *                    = lb*b*q + b*s + (lb*r - b*s)
             *                    = b*(lb*q + s) + (lb*r - b*s)
             */

            if (q)
                /* replace q by lb*q + s */
                (*q)[j] += lb * r.lc();
            
            /* replace r by lb*r - b*s == lc(b)*r - lc(r)*b*x^(deg(r)-deg(b)) */
            auto const lr = r.lc();
            r.chop();
            r *= lb;
            for(int i = 0 ; i < n ; i++)
                r[i+j] -= lr * b[i];
            r.cleandeg();
            e--;
        }
        ASSERT(e >= 0);

        auto const d_e = cado_math_aux::pow(lb, e);
        if (q) *q *= d_e;
        r *= d_e;

        return r;
    }
    /* }}} */

    public:

    /* {{{ divisions for integral polynomials. We only want divisions
     * that are exact.  Those are not provided readily by the
     * pseudo-division algorithm above. Essentially, it means that we're
     * interested in monic divisors, although it's not a necessity.
     *
     * This returns the quotient and remainder IF the division is exact.
     * Otherwise, {0,0} is returned and the exact flag is set to false.
     */
    template<typename U>
    std::pair<
        polynomial<eval_type_t<T, U>>,
        polynomial<eval_type_t<T, U>>
        >
        div_qr(polynomial<U> const & b, bool & exact) const
        requires cado_math_aux::is_integral_v<eval_type_t<T, U>>
        {
            using E = eval_type_t<T, U>;
            static_assert(std::is_empty_v<number_context<E>>);

            polynomial<E> r = *this;
            polynomial<E> q;

            auto const & lb = b.lc();

            for ( ; r.degree() >= b.degree() ; ) {
                int const j = r.degree() - b.degree();
                if (r.lc() % lb) {
                    exact = false;
                    q = 0;
                    r = 0;
                    return { q, r };
                }
                q[j] = cado_math_aux::divexact(r.lc(), lb);
                /* do not bother updating r.lc() */
                r.chop();
                for (int i = 0 ; i < b.degree() ; i++)
                    r[i+j] -= q[j] * b[i];
                r.cleandeg();
            }
            exact = true;
            return { q, r };
        }

    template<typename U>
    polynomial<eval_type_t<T, U>>
        div_q(polynomial<U> const & b, bool & exact) const
        requires cado_math_aux::is_integral_v<eval_type_t<T, U>>
        {
            return div_qr(b, exact).first;
        }

    template<typename U>
    polynomial<eval_type_t<T, U>>
        div_r(polynomial<U> const & b, bool & exact) const
        requires cado_math_aux::is_integral_v<eval_type_t<T, U>>
        {
            return div_qr(b, exact).second;
        }

    /* it is also possible to just ignore the exact flag. Note that
     * inexact results will always return (0,0), which is in many cases a
     * good way to detect. (if a.div_qr = (0,0) and a!=0, it means that
     * the division was inexact).
     */
    template<typename U>
    std::pair<
        polynomial<eval_type_t<T, U>>,
        polynomial<eval_type_t<T, U>>
        >
        div_qr(polynomial<U> const & b) const
        requires cado_math_aux::is_integral_v<eval_type_t<T, U>>
        {
            bool e;
            return div_qr(b, e);
        }
    template<typename U>
        polynomial<eval_type_t<T, U>>
        div_q(polynomial<U> const & b) const
        requires cado_math_aux::is_integral_v<eval_type_t<T, U>>
        {
            bool e;
            return div_q(b, e);
        }
    template<typename U>
        polynomial<eval_type_t<T, U>>
        div_r(polynomial<U> const & b) const
        requires cado_math_aux::is_integral_v<eval_type_t<T, U>>
        {
            bool e;
            return div_r(b, e);
        }
    /* }}} */

    /* {{{ division of non-integral polynomials. */
    template<typename U>
    std::pair<
        polynomial<eval_type_t<T, U>>,
        polynomial<eval_type_t<T, U>>
        >
        div_qr(polynomial<U> const & b) const
        requires (!cado_math_aux::is_integral_v<eval_type_t<T, U>>)
        {
            using E = eval_type_t<T, U>;
            auto const tr = number_context_for_evaluation(b.ctx());

            polynomial<E> r { *this, tr };
            polynomial<E> q;

            if (degree() < b.degree())
                return { q, r };

            auto const ilb = tr(1) / b.lc();
            using cado_math_aux::submul;

            for (int k = degree() - b.degree(); k >= 0; k--) {
                q[k] = r[k + b.degree()] * ilb;
                /* do not bother updating r.lc() */
                for (int j = b.degree() + k - 1; j >= k; j--)
                    submul(r.coeffs[j], q.coeffs[k], b.coeffs[j-k]);
            }
            r.cleandeg(b.degree()-1);
            return { q, r };
        }

    template<typename U>
    polynomial<eval_type_t<T, U>>
        div_q(polynomial<U> const & b) const
        requires (!cado_math_aux::is_integral_v<eval_type_t<T, U>>)
        {
            return div_qr(b).first;
        }

    template<typename U>
    polynomial<eval_type_t<T, U>>
        div_r(polynomial<U> const & b) const
        requires (!cado_math_aux::is_integral_v<eval_type_t<T, U>>)
        {
            return div_qr(b).second;
        }
    /* }}} */

    public:

    /* {{{ resultant of integral polynomials
     *
     * See Algorithm 3.3.7 in Cohen, "A Course in Computational Algebraic
     * Number Theory".
     *
     * Using this algorithm on floating point types would be wrong (see
     * §3.3.3 in Cohen). If we really want to compute resultants with
     * floating point types, we should compute the determinant of the
     * Sylvester matrix.
     *
     * Note that while it is definitely possible to compute resultants of
     * two polynomial<int> values with this code, we're almost certain to
     * hit overflows (that is, the computation will be modulo whatever
     * power of two is appropriate)
     */
    template<typename U>
    eval_type_t<T, U> resultant(polynomial<U> const & q) const
    requires cado_math_aux::is_integral_v<eval_type_t<T, U>>
    {
        using E = eval_type_t<T, U>;
        static_assert(std::is_empty_v<number_context<E>>);
        polynomial const & p = *this;
        if (p.degree() < 0 || q.degree() < 0)
            return 0;

        polynomial<E> a = p;
        polynomial<E> b = q;
        polynomial<E> r;

        int s = 1;

        using cado_math_aux::pow;

        T g = 1;
        T h = 1;

        /* we might want to divide out by the content if there is one */

        if (a.degree() < b.degree()) {
            std::swap(a, b);

            if (a.degree() & b.degree() & 1)
                s = -1;
        }

        while (b.degree() > 0) {
            int d = a.degree() - b.degree();

            if (a.degree() & b.degree() & 1)
                s = -s;

            r = a.pseudo_division(nullptr, b);

            a = b;

            b = r / (g * pow(h, d));

            g = a.lc();

            ASSERT(d != 0 || h == 1);

            h = pow(g, d) / pow(h, d - 1);
        }

        ASSERT(a.degree() > 0);
        if (b.degree() < 0)
            return 0;

        return s * pow(b[0], a.degree()) / pow(h, (a.degree() - 1));
    }
    /* }}} */
};


static_assert(std::is_same_v<decltype(polynomial<int>{}(double())), double>);
static_assert(std::is_same_v<decltype(polynomial<int>{}(int())), int>);
static_assert(std::is_same_v<decltype(polynomial<cxx_mpz>{}(int())), cxx_mpz>);
static_assert(std::is_same_v<decltype(polynomial<cxx_mpz>{}(double())), double>);
static_assert(std::is_same_v<eval_type_t<long double, std::complex<long double> >, std::complex<long double>>);

#ifdef HAVE_MPFR
static_assert(std::is_same_v<decltype(polynomial<cxx_mpz>{}(cxx_mpfr())), cxx_mpfr>);
#endif
#ifdef HAVE_MPC
static_assert(std::is_same_v<decltype(polynomial<double>{}(cxx_mpc())), cxx_mpc>);
#endif


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

/* of course it's also okay to print non-const references */
template<typename T>
inline std::ostream& operator<<(std::ostream& o, cado::named_proxy<polynomial<T> &> const & f)
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
    template<typename T>
    struct formatter<cado::named_proxy<polynomial<T>&>>: ostream_formatter {};
    template<typename T>
    struct formatter<cado::named_proxy<polynomial<T> const &>>: ostream_formatter {};
}


#endif	/* CADO_UTILS_POLYNOMIAL_HPP */
