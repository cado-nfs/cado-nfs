#ifndef CADO_UTILS_POLYNOMIAL_HPP
#define CADO_UTILS_POLYNOMIAL_HPP

/* This defines polynomial over arbitrary types, provided that these have
 * standard operator overloads defined. A priori we want to instantiate
 * these with float, double, and long double. But in principle it shoud
 * be possible to use this type more generically
 */

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <ios>
#include <istream>
#include <ostream>
#include <sstream>
#include <string>      // for string
#include <type_traits>
#include <utility>
#include <vector>

#include <gmp.h>
#include "fmt/base.h"
#include "fmt/ostream.h"

#include "macros.h"
#include "runtime_numeric_cast.hpp"
#include "mpz_poly.h"
#include "cxx_mpz.hpp"
#include "cado_math_aux.hpp"
#include "cado_expression_parser.hpp"
#include "coeff_proxy.hpp"

namespace polynomial_details {
    template<typename U>
    class named_proxy {
        static_assert(std::is_reference_v<U>, "U must be a reference");
        using V = std::remove_reference_t<U>;
        using Vnc = std::remove_const_t<V>;
        using nc = named_proxy<Vnc &>;
        static constexpr const bool is_c = std::is_const_v<V>;
        public:
        U c;
        std::string x;
        named_proxy(U c, std::string x)
            : c(c), x(std::move(x))
        {}

        template<typename W = U>
        explicit named_proxy(W const & c)
                requires(std::is_same_v<W, nc>)
        : c(c.c), x(c.x) {}
    };
}

template<typename T>
struct polynomial;

/* forward-declare the template. We need it in order to be able to
 * declare it as a friend of the polynomial<>struct
 */
template<typename T> std::istream& operator>>(std::istream& in, polynomial_details::named_proxy<polynomial<T> &> const & F);

template<typename T>
struct polynomial {
    typedef T coefficient_type;
    private:
    std::vector<T> coeffs;

    public:
    int degree() const {
        return runtime_numeric_cast<int>(coeffs.size())-1;
    }
    unsigned int size() const { return coeffs.size(); }

    T lc() const { ASSERT_ALWAYS(!coeffs.empty()); return coeffs.back(); }

    polynomial() = default;
    ~polynomial() = default;
    polynomial(polynomial const&) = default;
    polynomial(polynomial &&) = default;
    polynomial& operator=(polynomial const&) = default;
    polynomial& operator=(polynomial &&) = default;

    /* all instantations love each other */
    template<typename U> friend struct polynomial;
    template<typename U>
        polynomial(polynomial<U> const & a)
        requires (!std::is_same_v<U, T>)
        : coeffs { a.coeffs.begin(), a.coeffs.end() }
    {}

    void set_zero() { coeffs.clear(); }
    void set_xi(unsigned int i) {
        coeffs.assign((i+1), T(0));
        coeffs[i] = 1;
    }
    /*
     * use P[i] += v instead
    polynomial& add_xi(unsigned int i, T v = 1) {
        if (coeffs.size() <= i)
            coeffs.insert(coeffs.end(), i + 1 - coeffs.size(), T(0));
        coeffs[i] += v;
        cleandeg();
        return *this;
    }
    */

    polynomial& operator=(T v) {
        coeffs.clear();
        (*this)[0] = v;
        return *this;
    }
    explicit polynomial(T v) : coeffs(1, v) { cleandeg(); }
    polynomial(std::initializer_list<T> l) : coeffs(l.begin(), l.end()) {}

    explicit operator cxx_mpz_poly() const {
        cxx_mpz_poly res;
        for(unsigned int i = 0 ; i < size() ; i++) {
            /* coeffs[i] may be float, double, or long double */
            mpz_poly_setcoeff(res, i, cado_math_aux::mpz_from<T>(coeffs[i]));
        }
        return res;
    }

    // NOLINTNEXTLINE(hicpp-explicit-conversions)
    polynomial(std::string const & e) : polynomial() {
        if (!(std::istringstream(e) >> *this))
            throw cado_expression_parser_details::parse_error();
    }
    // NOLINTNEXTLINE(hicpp-explicit-conversions)
    polynomial(const char * e) : polynomial(std::string(e)) {}

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


    T eval(T x) const {
        T const * f = coeffs.data();
        const int deg = degree();
        switch(deg) {
            case -1: return 0;
            case 0: return f[0];
            case 1: return f[0]+x*f[1];
            case 2: return f[0]+x*(f[1]+x*f[2]);
            case 3: return f[0]+x*(f[1]+x*(f[2]+x*f[3]));
            case 4: return f[0]+x*(f[1]+x*(f[2]+x*(f[3]+x*f[4])));
            case 5: return f[0]+x*(f[1]+x*(f[2]+x*(f[3]+x*(f[4]+x*f[5]))));
            case 6: return f[0]+x*(f[1]+x*(f[2]+x*(f[3]+x*(f[4]+x*(f[5]+x*f[6])))));
            case 7: return f[0]+x*(f[1]+x*(f[2]+x*(f[3]+x*(f[4]+x*(f[5]+x*(f[6]+x*f[7]))))));
            case 8: return f[0]+x*(f[1]+x*(f[2]+x*(f[3]+x*(f[4]+x*(f[5]+x*(f[6]+x*(f[7]+x*f[8])))))));
            case 9: return f[0]+x*(f[1]+x*(f[2]+x*(f[3]+x*(f[4]+x*(f[5]+x*(f[6]+x*(f[7]+x*(f[8]+x*f[9]))))))));
            default:
                    T r = f[deg];
                    for (int k = deg; k-- ; r = r * x + f[k]);
                    return r;
        }
    }
    T operator()(T x) const { return eval(x); }

    /* Evaluate the homogenous polynomial induced by f at the pair
     * (x,y). That is, compute the sum f[i]*x^i*y^(n-i), where n is
     * degree(f)
     */
    T eval(T x, T y) const
    {
        T const * f = coeffs.data();

        switch (size()) {
            case 0: return 0;
            case 1: return f[0];
            case 2: return y*f[0]+x*f[1];
            case 3: return y*y*f[0]+x*(y*f[1]+x*f[2]);
            default:
                    T s = 0;
                    T px = 1;
                    for(unsigned int k = 0 ; k < size() ; k++) {
                        s = s * y + f[k] * px;
                        px = px * x;
                    }
                    return s;
        }
    }
    T operator()(T x, T y) const { return eval(x, y); }

    /* uses arbitrary precision */
    T eval_safe(T x) const
    {
        T const * f = coeffs.data();
        const int deg = degree();
        cxx_mpz xm; int xe; exact_form(xm, xe, x);
        /* We want to evaluate at xz * 2^xe */
        cxx_mpz vm; int ve; exact_form(vm, ve, f[deg]);
        for(int k = deg ; k-- ; ) {
            /* multiply by x = xm*2^xe */
            mpz_mul (vm, vm, xm);
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
        T r = cado_math_aux::mpz_get<T> (vm);
        return ldexp(r, ve);
    }

    T findroot_dichotomy(T a, T b) const {
        return findroot_dichotomy(a, b, eval(a));
    }

    T findroot_dichotomy(T a, T b, T sa) const
    {
        T s;
        for(;;) {
            s = (a + b) * 0.5;
            cado_math_aux::do_not_outsmart_me(s);
            if (s == a || s == b) return s;
            using namespace cado_math_aux;
            if (sgn(eval(s)) * sgn(sa) > 0)
                a = s;
            else
                b = s;
        }
        return s;
    }

    polynomial derivative() const
    {
        if (degree() <= 0) return {};
        polynomial df;
        df.coeffs.reserve(degree() - 1);
        for(int i = 1 ; i <= degree() ; i++)
            df.coeffs.push_back(coeffs[i] * i);
        return df;
    }

    polynomial operator*(T a) const
    {
        polynomial h;
        h.coeffs.reserve(coeffs.size());
        for(auto const & x : coeffs)
            h.coeffs.push_back(x * a);
        return h;
    }

    polynomial operator/(T a) const
    {
        polynomial h;
        h.coeffs.reserve(coeffs.size());
        for(auto const & x : coeffs)
            h.coeffs.push_back(x / a);
        return h;
    }

    polynomial operator-() const
    {
        polynomial h;
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
        polynomial h;
        h.coeffs.assign(f.size() + g.size() - 1, 0);
        for(unsigned int i = 0 ; i < f.size() ; i++)
            for(unsigned int j = 0 ; j < g.size(); j++)
                h.coeffs[i+j] += f[i] * g[j];
        return h;
    }

    polynomial operator+(polynomial const & g) const
    {
        polynomial const & f = *this;
        polynomial h;
        h.coeffs.assign(std::max(f.size(), g.size()), 0);
        unsigned int i = 0;
        for( ; i < f.size() && i < g.size() ; i++)
            h.coeffs[i] = f[i] + g[i];
        for( ; i < f.size() ; i++)
            h.coeffs[i] = f[i];
        for( ; i < g.size() ; i++)
            h.coeffs[i] = g[i];
        h.cleandeg();
        return h;
    }

    polynomial operator-(polynomial const & g) const
    {
        polynomial const & f = *this;
        polynomial h;
        h.coeffs.assign(std::max(f.degree(), g.degree()) + 1, 0);
        unsigned int i = 0;
        for( ; i < f.size() && i < g.size() ; i++)
            h.coeffs[i] = f[i] - g[i];
        for( ; i < f.size() ; i++)
            h.coeffs[i] = f[i];
        for( ; i < g.size() ; i++)
            h.coeffs[i] = -g[i];
        h.cleandeg();
        return h;
    }

    /* all compound operators are done lazily, at least for now */
    polynomial& operator*=(T a) {
        return (*this) = (*this) * a;
    }

    polynomial& operator/=(T a) {
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

    polynomial& addmul(polynomial const & a, polynomial & b)
    {
        return (*this) += a*b;
    }

    polynomial& submul(polynomial const & a, polynomial & b)
    {
        return (*this) -= a*b;
    }

    polynomial& addmul(polynomial const & a, T v)
    {
        return (*this) += a*v;
    }

    polynomial& submul(polynomial const & a, T v)
    {
        return (*this) -= a*v;
    }

    polynomial reciprocal() const
    {
        polynomial h;
        h.coeffs.reserve(coeffs.size());
        for(unsigned int i = 0 ; i < size() ; i++)
            h[size()-1-i] = coeffs[i];
        return h;
    }

    /* This computes u = scale^deg(f) * f(x/scale)
     * not the same as double_poly_scale, deleted in commit
     * c480fe82174a9de96e1cd35b2317fdf0de3678ab
     */
    polynomial reverse_scale(double scale) const
    {
        unsigned int sz = size();
        if (!sz)
            return {};
        polynomial h;
        h.coeffs.reserve(sz);
        h[--sz] = lc();
        for(T s = scale ; sz-- ; s *= scale) h[sz] = coeffs[sz] * s;
        return h;
    }

    /* Return a bound on the positive roots of p.
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

    T bound_positive_roots(bool negative = false) const
    {
        const int d = degree();
        const int s = lc() < 0;
        polynomial q;
        q.coeffs.assign(coeffs.size(), 0);
        for(int i = 0 ; i < d ; i++) {
            T v = (s ^ (negative & (d-i))) ? -coeffs[i] : coeffs[i];
            /* simplifies to v==(-1)^s*coeffs[i] < 0 if negative == 0 */
            if (v < 0)
                q[i] = v;
        }
        q[degree()] = abs(lc());
        T b = 1;
        for( ; q.eval(b) < 0 ; b = b + b) ;
        return negative ? -b : b;
    }

    std::vector<T> positive_roots() const {
        return positive_roots(bound_positive_roots());
    }
    std::vector<T> negative_roots() const {
        return positive_roots(bound_positive_roots(true));
    }

    private:

    /* assuming g(a)*g(b) < 0, and g has a single root in [a, b],
     * refines that root by the weighted false position method
     * Assumes sa is of same sign as g(a).
     *
     * The code is written with the case a<b in mind, but it may also be
     * called with b<a, in which case we need to adapt a few little
     * things.
     */
    T findroot_falseposition(T a, T b, T pa) const
    {
        T pb;
        int side=0;
        T a0=a, b0=b, pa0=pa;

        if (a == b)
            return a;

        using namespace cado_math_aux;
        int sigma = cado_math_aux::sgn(b-a);

        ASSERT_ALWAYS(sigma*a < sigma*b);

        pb = eval(b);

        for(;;) {
            T s, middle;

            s = (a*pb-b*pa)/(pb-pa);
            middle = (a + b) * 0.5;

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
            T ps = eval(s);
            using namespace cado_math_aux;
            if (sgn(ps) * sgn(pa) > 0) {
                a = s; pa = ps;
                if (side==1) pb /= 2;
                side=1;
            } else {
                b = s; pb = ps;
                if (side==-1) pa /= 2;
                side=-1;
            }
            if (std::isnan(b)) {
                return findroot_dichotomy(a0, b0, pa0);
            }
        }
    }

    /* knowing the positive sign changes of the derivative of *this given
     * in v , as well as a bound on the positive roots of *this, store in
     * v the positive roots of *this.  v is clobbered.
     */
    void positive_roots_from_derivative_sign_changes(std::vector<T> & v, T bound)
    {
        using namespace cado_math_aux;
        if (degree() <= 0) {
            /* A constant polynomial has no sign changes */
            v.clear();
        } else if (degree() == 1) {
            /* A linear polynomial has at most one root.
             *
             * We want strictly positive roots here, so we must not
             * consider the case coeffs[0] == 0. On the other hand, the
             * bound counts.
             */
            const int s = sgn(coeffs[0]);
            if (s && s * eval(bound) <= 0)
                v.assign(1, - coeffs[0] / coeffs[1]);
        } else {
            T a = 0;
            T va = coeffs[0];
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
            T c01 = sgn(coeffs[0]) * sgn(coeffs[1] ? coeffs[1] : (sgn(bound) * coeffs[2]));
            bool no_chance = c01 * sgn(bound) > 0 || coeffs[0] == 0;
            for(size_t i = 0 ; i < v.size() ; i++) {
                T b = v[i];
                T vb = eval(b);
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
    public:

    std::vector<T> positive_roots(T bound) const
    {
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

        /* work from the most derived polynomial, down to f */
        std::vector<T> res;
        for (int k = d; k-- ; )
            dg[k].positive_roots_from_derivative_sign_changes(res, bound);

        return res;
    }

    std::vector<T> roots() const
    {
        if (degree() == -1) return {};

        auto w = negative_roots();

        if (coeffs[0] == 0)
            w.push_back(0);

        const auto positive = positive_roots();

        w.insert(w.end(), positive.begin(), positive.end());
        std::sort(w.begin(), w.end());

        return w;
    }


    /* divide by x-r */
    polynomial div_linear(T r) const
    {
        T const * f = coeffs.data();
        const int deg = degree();
        if (deg < 0)
            return {};
        T u = f[deg];
        polynomial q;
        q.coeffs.reserve(coeffs.size() - 1);
        for (int k = deg ; k-- ; ) {
            T const c = f[k];
            q[k] = u;
            u = u * r + c;
        }
        /* u is f(r), we might as well return it in some way. */
        return q;
    }

    std::string print(std::string const& var = "x") const
    {
        std::ostringstream os;
        if (degree() < 0) os << "0";
        for(int i = 0 ; i <= degree() ; i++) {
            T const & fi = coeffs[i];
            int const r = (fi > 0) - (fi < 0);
            if (r == 0) continue;
            if (r > 0 && !os.str().empty())
                os << "+";
            if (i == 0) {
                os << fi;
            } else {
                if (fi == -1) {
                    os << "-";
                } else if (fi != 1) {
                    os << fi << "*";
                }
                os << var;
                if (i > 1) os << "^" << i;
            }

        }
        return os.str();
    }

    explicit polynomial(cxx_mpz_poly const & f)
    {
        coeffs.assign(f.degree() + 1, 0);
        for(int i = 0 ; i <= f.degree() ; i++)
            coeffs[i] = cado_math_aux::mpz_get<T>(mpz_poly_coeff_const(f, i));
    }

    private:
    friend std::istream& operator>><T>(std::istream& in, polynomial_details::named_proxy<polynomial &> const & F);

    struct parser_traits {
        std::string x;
        explicit parser_traits(std::string x) : x(std::move(x)) {}
        struct unexpected_literal : public cado_expression_parser_details::parse_error {
            std::string msg;
            explicit unexpected_literal(std::string const & v)
                : msg(std::string("unexpected literal " + v))
            {}
            const char *what() const noexcept override {
                return msg.c_str();
            }
        };
        static constexpr const int accept_literals = 1;
        typedef polynomial type;
        typedef T number_type;
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
        void pow_ui(polynomial & c, polynomial const & a, T e) {
            c = a * e;
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

    polynomial_details::named_proxy<polynomial &> named(std::string const & x) {
        return { *this, x };
    }
    polynomial_details::named_proxy<polynomial const &> named(std::string const & x) const {
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
    bool operator!=(polynomial const & f) const { return !operator==(f); }
    bool operator<(polynomial const & f) const { return spaceship(f) < 0; }
    bool operator<=(polynomial const & f) const { return spaceship(f) <= 0; }
    bool operator>(polynomial const & f) const { return spaceship(f) > 0; }
    bool operator>=(polynomial const & f) const { return spaceship(f) >= 0; }
    bool operator!=(T v) const { return !operator==(v); }
    bool operator==(T v) const {
        return (degree() < 0 && v == 0) || (degree() == 0 && coeffs[0] == v);
    }

    bool has_nan() const {
        for(auto c: coeffs)
            if (std::isnan(c))
                return true;
        return false;
    }

    bool has_inf() const {
        for(auto c: coeffs)
            if (std::isinf(c))
                return true;
        return false;
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
            polynomial const & b)
    {
        polynomial const & a = *this;
        ASSERT(a.degree() >= b.degree());
        ASSERT(b.degree() != -1);

        int const m = a.degree();
        int const n = b.degree();
        T d = b.lc();
        int e = m - n + 1;
        polynomial s;

        if (q) *q = 0;

        r = a;

        while (r.degree() >= n) {
            s = 0;
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

        d = std::pow(d, (T) e);
        if (q) *q *= d;
        r *= d;

        return true;
    }

    bool pseudo_remainder(polynomial & r, polynomial const & b)
    {
        return pseudo_division(nullptr, r, b);
    }

    public:

    T resultant(polynomial const & q) const
    {
        polynomial const & p = *this;
        if (p.degree() < 0 || q.degree() < 0)
            return 0;

        polynomial a = p;
        polynomial b = q;
        polynomial r;

        int s = 1;
        int d;

        int pseudo_div = 1;

        T g = 1;
        T h = 1;

        if (a.degree() < b.degree()) {
            std::swap(a, b);

            if ((a.degree() % 2) == 1 && (b.degree() % 2) == 1)
                s = -1;
        }

        while (b.degree() > 0) {
            //TODO: verify if it is necessary.
            d = a.degree() - b.degree();

            if ((a.degree() % 2) == 1 && (b.degree() % 2) == 1)
                s = -s;

            pseudo_div = a.pseudo_remainder(r, b);
            if (!pseudo_div)
                break;

            a = b;

            ASSERT(d >= 0);

            b = r / (g * std::pow(h, d));

            g = a.lc();

            ASSERT(d != 0 || h == 1);

            h = std::pow(h, (T) (d - 1));
            h = std::pow(g, (T) d) / h;
        }

        if (pseudo_div) {
            ASSERT(a.degree() > 0);

            //For now, if b.degree() == -1, pseudo_div == 0.
            if (b.degree() == -1) {
                ASSERT(0);
            } else {

                ASSERT(a.degree() >= 0);

                h = std::pow(b[0], (T) a.degree()) / std::pow(h, (T) (a.degree() - 1));
                h *= s;
            }
        } else {
            // we encountered cancellations, so we need to resort to
            // exact arithmetic.
            // TODO: use last version of a and b in pseudo_division.
            cxx_mpz val_z;
            cxx_mpz_poly pz(p);
            cxx_mpz_poly qz(q);
            mpz_poly_resultant(val_z, pz, qz);
            return cado_math_aux::mpz_get<T>(val_z);
        }
        return h;
    }
};


template<typename T>
std::istream& operator>>(std::istream& in, polynomial_details::named_proxy<polynomial<T> &> const & F)
{
    std::string line;
    for(;;in.get()) {
        int const c = in.peek();
        if (in.eof() || !isspace(c)) break;
    }
    if (!getline(in, line)) return in;
    std::istringstream is(line);

    typedef typename polynomial<T>::parser_traits traits_type;
    typedef cado_expression_parser<traits_type> poly_parser;
    poly_parser P(F.x);
    P.tokenize(is);

    try {
        F.c = P.parse();
    } catch (cado_expression_parser_details::parse_error const & p) {
        in.setstate(std::ios_base::failbit);
        return in;
    }

    return in;
}

template<typename T>
std::istream& operator>>(std::istream& in, polynomial<T> & F)
{
    return in >> F.named("x");
}

/* printing needs a way to specify the variables... */
template<typename T>
inline std::ostream& operator<<(std::ostream& o, polynomial_details::named_proxy<polynomial<T> const &> const & f)
{
    return o << f.c.print(f.x);
}

/* we do have a default behaviour, though */
template<typename T>
inline std::ostream& operator<<(std::ostream& o, polynomial<T> const & f)
{
    return o << f.named("x");
}

namespace fmt {
    template<typename T>
    struct formatter<polynomial<T>>: ostream_formatter {};
}


#endif	/* CADO_UTILS_POLYNOMIAL_HPP */
