#include "cado.h" // IWYU pragma: keep

#include <cstdint>
#include <cstdlib>

#include <stdexcept>

#include "cado_constants.hpp"
#include "cado_math_aux.hpp"
#ifdef HAVE_MPFR
#include "cxx_mpfr.hpp"
#endif

#include "macros.h"
#include "mpz_poly.h"
#include "number_context.hpp"
#include "polynomial.hpp"
#include "polyselect_norms.hpp"

/************************* norm and skewness *********************************/

/* The L2-(log)norm here as defined by Kleinjung is
 * log(1/2 sqrt(int(int((F(sx,y)/s^(d/2))^2, x=-1..1), y=-1..1))).
 * Since we only want to compare norms, we don't consider the log(1/2) term,
 * and compute only 1/2 log(int(int(...))) [here the 1/2 factor is important,
 * since it is added to the alpha root property term].
 *
 * In cado, we compute the circular L2-norm (the lognorm is 1/2*log(the
 * L2 norm)), where we integrate over the unit circle: the formula is
 *
 * I = int(int(F(r*cos(t)*s,r*sin(t))^2*r/s^d, r=0..1), t=0..2*Pi).
 *
 * Cf Remark 3.2 in Kleinjung paper, Math. of Comp., 2006.
 *
 * As it happens, the circular L2-norm is a quadratic form
 * in the coefficients of the input polynomial. We can use it to compute
 * the L2-lognorm at a very small cost.
 *
 * L2_lognorm() takes an mpz_poly_srcptr, while L2_lognorm() takes a
 * double_poly_srcptr
 *
 * Note that by homogeneity we have
 *
 * I = int(int(F(r*cos(t)*s,r*sin(t))^2*r/s^d, r=0..1), t=0..2*Pi)
 *   = int(int(F(cos(t)*s,sin(t))^2*r^(2*d+1)/s^d, r=0..1), t=0..2*Pi)
 *   = int(r^(2*d+1),r=0..1)*int(F(cos(t)*s,sin(t))^2/s^d, t=0..2*Pi)
 *   = 1/(2*d+2)*int(F(cos(t)*s^(1/2),sin(t)/s^(1/2))^2, t=0..2*Pi)
 *
 * Sagemath code computing the circular L2-norm (up to scaling)
    S.<a> = InfinitePolynomialRing(SR)
    def sum_d(d, a, force_poly=True):
        s = SR.var('s')
        t = SR.var('t')
        x = s^(1/2)*cos(t)
        y = s^(-1/2)*sin(t)
        F = sum(a[i]*x^i*y^(d-i) for i in range(d+1))
        if force_poly:
            I = sum([a[i]*a[j]*integrate(x^(i+j)*y^(2*d-i-j),(t,0,2*pi))/(2*d+2)
                     for i in range(d+1) for j in range(d+1)])
        else:
            I = integrate(F^2,(t,0,2*pi))/(2*d+2)
        return I

    for d in range(7):
        print(sum_d(d, a, False).expand().collect(pi))
    
 * e.g. for d==6 and s==1 we have:
 *
    n = 231 * (a6 * a6 + a0 * a0) + 42.0 * (a6 * a4 + a2 * a0)
       + 21 * (a5 * a5 + a1 * a1) + 7.0 * (a4 * a4 + a2 * a2)
       + 14 * (a6 * a2 + a5 * a3 + a4 * a0 + a3 * a1)
       + 10 * (a6 * a0 + a5 * a1 + a4 * a2) + 5.0 * a3 * a3;
     n = n * pi/7168

 * Note that we can factor out stuff that depends only on the degree,
 * since this adds up to a constant offset in the lognorm. So the
 * expression below is just as useful.
 *
   for d in range(7):
       print(((sum_d(d, a, False)*2^(2*d)*(d+1)).subs(s=1).expand().collect(pi)/pi))
 *
 * This expression I* is 2^(2*d)/pi*(d+1) times the similar expression
 * I above. (Actually I* has nontrivial content for all d>=1, but by
 * keeping this stray 2 we have something that casts to Z even in the
 * degenerate case d=0). We can write I*, for example for d=6 and s==1, as:
 *
 *   924*(a_0^2 + a_6^2)
 * + 84*(a_1^2 + 2*a_0*a_2 + 2*a_4*a_6 + a_5^2)
 * + 28*(a_2^2 + 2*a_1*a_3 + 2*a_2*a_6 + a_4^2 + 2*a_0*a_4 + 2*a_3*a_5)
 * + 20*(a_3^2 + 2*a_2*a_4 + 2*a_1*a_5 + 2*a_0*a_6)
 *
 * As a matter of fact, the coefficients have a simple expression. Let
 * k,d be such that 0<=k<=d, and:
 *
 *      W = lambda d, k: integrate(cos(t)^k*sin(t)^(d-k),(t,0,2*pi))
 *      w(d,k) = 2^(2*d-1)/pi * W(2d,2k)
 * 
 * It is easy to see that:
 *  - W is non-zero only for k,d even 
 *  - W is symmetric with respect to k <--> d-k
 *  - w(d,0) = 2*(2*d-1)/d*w(d-1,0) = \binom{2d}{d}
 *  - w(d,k) = 2*(2k-1)/d*w(d-1,k-1) = w(d,0)*\binom{d}{k}/\binom{2d}{2k}
 *
 * The following code prints I and I*
 *
    W = lambda d, k: integrate(cos(t)^k*sin(t)^(d-k),(t,0,2*pi))
    w = lambda d, k: W(2*d,2*k)/pi*2^(d-1)
    w = lambda d, k: binomial(2*d,d) * binomial(d,k) / binomial(2*d,2*k)
    for d in range(1,7):
        for k in range(d+1):
            assert w(d,k) == W(2*d,2*k)/pi*2^(2*d-1)
    H = lambda PR,d: sum([ZZ(w(d,i))*PR.gen()^(2*i) for i in range(d+1)])

    def Istar(d, a):
        S = a[0].parent()
        SP = S['x']
        hh = H(SP,d)
        f = SP([a[i] for i in range(d+1)])
        return S((hh*f^2)[2*d])

    # sanity check
    d = 6
    I = sum_d(d, a)
    # see https://github.com/sagemath/sage/issues/41656 for the ugly
    # quirk. It's not easy to mix-and-match SR and InfinitePolynomialRing
    assert I.map_coefficients(lambda c: c.subs(s=1)) == pi/(2^(2*d)*(d+1)) * Istar(d, a)
    b = [a[i]*s^(i-d/2) for i in range(d+1)]
    assert I == pi/(2^(2*d)*(d+1)) * Istar(d, b)


    for d in range(1,7):
        hh = H(S['x'],d)
        print(f"d={d}:")
        print(f"\tH_{d}(x)={hh}")
        Is = Istar(d, a)
        print(f"\tI*=({Is})")
        g = gcd(Is.coefficients())
        print(f"\tI=pi/{2^(2*d)*(d+1)/g}*({Is/g})")

 * 
 * and the derivative of the L2 norm with respect to s can also be
 * computed with similar code.
 *

    def Istar_diff(d, a, force_poly=True):
        ckrange = lambda k: range(max(2*k-d,0),min(d,2*k)+1)
        ck = lambda k: sum([a[i]*a[2*k-i] for i in ckrange(k)])
        tt = sum([w(d, k)*ck(k)*(2*k-d)*s^(2*k-d-1) for k in range(d+1)])
        if force_poly:
            return tt
        else:
            return SR(tt)

    # sanity check
    d = 6
    b = [a[i]*s^(i-d/2) for i in range(d+1)]
    assert Istar(d, b).map_coefficients(lambda c:c.diff(s)) == Istar_diff(d, a)
        
 * The quadratic form that is defined on the input polynomial, and its
 * derivative, are obtained with the following code:
 *
    def Q(f):
        # f must be a polynomial
        return FractionField(QQ['x'])((sum_d(f.degree(),f)/pi).subs(s=x))
    def dQ(f):
        # f must be a polynomial
        d = f.degree()
        return FractionField(QQ['x'])(SR(sum_d_diff(d,f)).subs(s=x^2)*x^(-d-1)/(d+1)/2^(2*d))

 */

#define L2_DMAX 16
const uint32_t coeffs_integral[L2_DMAX][L2_DMAX/2] = 
{
    /* [[ binomial(2*d,d) * binomial(d,k) / binomial(2*d,2*k) /((k|(d&1))?1:2)
     *    for k in range((d+1)//2, d+1)]
     *   for d in range(16)]
     */
    { 1 },
    { 2 },
    { 1, 6 },
    { 4, 20 },
    { 3, 10, 70 },
    { 12, 28, 252 },
    { 10, 28, 84, 924 },
    { 40, 72, 264, 3432 },
    { 35, 90, 198, 858, 12870 },
    { 140, 220, 572, 2860, 48620 },
    { 126, 308, 572, 1716, 9724, 184756 },
    { 504, 728, 1560, 5304, 33592, 705432 },
    { 462, 1092, 1820, 4420, 16796, 117572, 2704156 },
    { 1848, 2520, 4760, 12920, 54264, 416024, 10400600 },
    { 1716, 3960, 6120, 12920, 38760, 178296, 1485800, 40116600 },
    { 6864, 8976, 15504, 36176, 118864, 594320, 5348880, 155117520 }
};

template<typename T>
static T L2_lognorm(polynomial<T> const & p, T const & s)
{
    const int d = p.degree();
    ASSERT_ALWAYS(1 <= p.degree());
    ASSERT_ALWAYS(d <= L2_DMAX);
    /* a rewritten version of the sum above:
       def sum_d(d, a):
           tt = 0
           invs = 1/s
           s2 = s * s
           is2 = invs * invs
           u0 = 1 if (d % 2 == 0) else invs
           u1 = 1 if (d % 2 == 0) else s
           p0 = d//2
           p1 = d - p0
           for k in range(d//2+1):
               t0 = 0
               t1 = 0
               for i in range(1, p0-k+1):
                   t0 += a[p0-k-i] * a[p0-k+i]
                   t1 += a[p1+k-i] * a[p1+k+i]
               t0 = (2 * t0 + a[p0-k] * a[p0-k]) * u0
               t1 = (2 * t1 + a[p1+k] * a[p1+k]) * u1
               u0 *= is2
               u1 *= s2
               tt += coeff2(d, k)*(t0 + t1)
           return tt*pi/2^(2*d)/(d+1)
       def check2(sumfunc):
           var('r,s,t,y')
           R.<x> = PolynomialRing(ZZ)
           S.<a> = InfinitePolynomialRing(R)
           for d in range(2,10):
               f = SR(sum(a[i]*x^i for i in range(d+1)))
               F = expand(f(x=x/y)*y^d)
               F = F.subs(x=s^(1/2)*r*cos(t),y=r/s^(1/2)*sin(t))
               v = integrate(integrate(F^2*r,(r,0,1)),(t,0,2*pi))
               print(v - sumfunc(d, a))
       check2(sum_d)
     */

    const T * a = p.get_coeff_data();
    auto const & tr = cado::number_context<T>(a[0]);

    const T invs = tr(1) / tr(s);
    const T s2 = tr(s) * tr(s);
    const T is2 = invs * invs;
    T tt = tr(0);
    T u0 = (d & 1) ? invs : tr(1);
    T u1 = (d & 1) ? tr(s) : tr(1);
    const int p0 = d / 2;
    const int p1 = d - p0;
    for(int k = 0 ; k <= d/2 ; k++) {
        T t0 = tr(0);
        T t1 = tr(0);
        for(int i = 1 ; i <= d/2 - k ; i++) {
            t0 += a[p0-k-i] * a[p0-k+i];
            t1 += a[p1+k-i] * a[p1+k+i];
        }
        t0 = (2 * t0 + a[p0-k] * a[p0-k]) * u0;
        t1 = (2 * t1 + a[p1+k] * a[p1+k]) * u1;
        u0 *= is2;
        u1 *= s2;
        tt += tr(coeffs_integral[d][k])*(t0 + t1);
    }
    tt = cado_math_aux::ldexp(tt * cado_math_aux::pi_v(tr) / (d+1), -2*d);
    if (cado_math_aux::isnan (tt) || cado_math_aux::isinf (tt))
        throw std::overflow_error("");
    return cado_math_aux::log(tt) / 2;
}

/* The L2 lognorm is (the image by a monotonous function of) a polynomial
 * in s. We first form this polynomial, and then compute its roots, which
 * will allow us to compute the best skewness that minimizes the L2
 * lognorm.
 *
 * We thus want to minimize the polynomial given by sum_d above, but
 * without the pi/2^(2*d-1)/(d+1) factor, which we obviously don't care
 * about because it's constant. Thus we're after the zeros of
 * (the numerator of) its derivative.
 *
 * sum_d actually a polynomial that is either odd or even depending on
 * whether d is itself odd or even. Its expression is reasonably simple,
 * and given by the python code sum_d_diff near the top of this file.
 *
 * So if dP is the polynomial of degree d that is computed by sum_d_diff
 * (and by the C code below), the minima s of the lognorm are such that
 * s^2 is a zero of dP. We thus want to find the positive zeros of dP.
 */

template<typename T>
static polynomial<T> L2_skewness_derivative_numerator(polynomial<T> const & p)
{
    const int d = p.degree();
    ASSERT_ALWAYS(1 <= p.degree());
    ASSERT_ALWAYS(d <= L2_DMAX);

    const T * a = p.get_coeff_data();
    auto const & tr = cado::number_context<T>(a[0]);

    polynomial<T> dP(tr);

    // double_poly_realloc(dP, d+1);
    for(int k = 0 ; k <= d ; k++) {
        T t0 = tr(0);
        for(int i = 1 ; i <= k && i <= d - k ; i++)
            t0 += a[k-i] * a[k+i];
        t0 = (2 * t0 + a[k] * a[k]);
        const int j = (ABS(d-2*k)-(d&1))/2;
        t0 = t0 * coeffs_integral[d][j]*(2*k-d);
        dP[k] = t0; // ldexp(t0/(d+1), -2*d);
    }
    if (dP.has_inf() || dP.has_nan())
        throw std::overflow_error("");

    return dP;
}

/* Do the same just for the numerator of the L2 skewness. Again, this is
 * a polynomial in s^2 !!!
    def sum_xd(d, a):
        sq = lambda k: a[k]^2
        cross = lambda k: 2*sum([a[k-i]*a[k+i] for i in range(1, min(d-k,k)+1)])
        ck = lambda k: (sq(k)+cross(k))
        tt = 0
        for l in range(d+1):
            tt += coeff2(d, d//2-l)*ck(l)*s^(2*l)/(1 if 2*l == d else 2)
        return tt*pi/2^(2*d-1)/(d+1)
    check2(lambda d,a:sum_xd(d,a)/s^d)
 */
template<typename T>
static polynomial<T> L2_skewness_numerator(polynomial<T> const & p)
{
    const int d = p.degree();
    ASSERT_ALWAYS(1 <= p.degree());
    ASSERT_ALWAYS(d <= L2_DMAX);

    const T * a = p.get_coeff_data();
    auto const & tr = cado::number_context<T>(a[0]);

    polynomial<T> dP(tr);

    // double_poly_realloc(dP, d+1);
    for(int k = 0 ; k <= d ; k++) {
        T t0 = tr(0);
        for(int i = 1 ; i <= k && i <= d - k ; i++)
            t0 += a[k-i] * a[k+i];
        t0 = (2 * t0 + a[k] * a[k]);
        const int j = (ABS(d-2*k)-(d&1))/2;
        t0 = t0 * coeffs_integral[d][j];
        if (2*k != d)
            t0 /= 2;
        dP[k] = t0; // ldexp(t0/(d+1), -(2*d-1));
    }
    // cleandeg() should be automatic with polynomial<T>
    // double_poly_cleandeg(dP, d);
    if (dP.has_inf() || dP.has_nan())
        throw std::overflow_error("");

    return dP;
}

/* return the skewness giving the best lognorm sum for two polynomials.
 * We do so by computing the formal expression of the quadratic form (up
 * to a factor that does not vary, see the commented out ldexp calls
 * above), and then we compute the roots.
 *
 * We're interested in log(Q(f))+log(Q(g)) = log(Q(f)Q(g)), so we want to
 * minimize dQ(f)*Q(g)+Q(f)*dQ(g). Pay attention to the subtleties with
 * the square roots.
 */

template<typename T>
static T L2_combined_skewness2(polynomial<T> const & F, polynomial<T> const & G)
{
    auto QF = L2_skewness_numerator(F);
    auto dQF = L2_skewness_derivative_numerator(F);
    auto QG = L2_skewness_numerator(G);
    auto dQG = L2_skewness_derivative_numerator(G);

    auto S = dQF * QG + QF * dQG;

    if (S.has_inf() || S.has_nan())
        ASSERT_ALWAYS(0);

    auto roots = S.positive_roots();

    /* some pathological examples for gfpn have integer roots, and this
     * wreaks havoc.
     */

    ASSERT_ALWAYS(!roots.empty());

    auto best_s = cado_math_aux::sqrt(roots[0]);
    auto l = L2_lognorm(F, best_s) + L2_lognorm(G, best_s);

    for(size_t i = 1 ; i < roots.size() ; i++) {
        const auto s = cado_math_aux::sqrt(roots[i]);
        const auto li = L2_lognorm(F, s) + L2_lognorm(G, s);
        if (li < l) {
            best_s = s;
            l = li;
        }
    }

    return best_s;
}

template<typename T>
static T
L2_skewness(polynomial<T> const & P)
{
    auto dP = L2_skewness_derivative_numerator(P);

    /* some pathological examples for gfpn have integer roots, and this
     * wreaks havoc.
     */
    auto roots = dP.positive_roots();

    /* We often have a single zero, but not always. Here's an example
     * with three. It typically happens when we have pairs of roots that
     * are almost symmetrical along the axes. I'm not sure I see a
     * generalization of this phenomenon to more than two extrema, actually.
     * 3520*x^4 - 14848*x^3 - 26038239232*x^2 + 3120*x + 110
     */
    ASSERT_ALWAYS(!roots.empty());

    T best_s = cado_math_aux::sqrt(roots[0]);
    T l = L2_lognorm(P, best_s);
    for(size_t i = 1 ; i < roots.size() ; i++) {
        const auto s = cado_math_aux::sqrt(roots[i]);
        const auto li = L2_lognorm(P, s);
        if (li < l) {
            best_s = s;
            l = li;
        }
    }

    return best_s;
}

double L2_lognorm (mpz_poly_srcptr f, double s)
{
    try {
        return L2_lognorm(polynomial<double>(f), s);
    } catch (std::overflow_error const &) {
        try {
            const long double sl = s;
            const auto l = L2_lognorm(polynomial<long double>(f), sl);
            return static_cast<double>(l);
#ifdef HAVE_MPFR
        } catch (std::overflow_error const &) {
            const cxx_mpfr sl = s;
            const auto l = L2_lognorm(polynomial<cxx_mpfr>(f), sl);
            return static_cast<double>(l);
#else
        } catch (...) {
            throw;
#endif
        }
    }
}

double L2_skewness (mpz_poly_srcptr f)
{
    try {
        return L2_skewness(polynomial<double>(f));
    } catch (std::overflow_error const &) {
        try {
            const auto s = L2_skewness(polynomial<long double>(f));
            return static_cast<double>(s);
#ifdef HAVE_MPFR
        } catch (std::overflow_error const &) {
            const auto s = L2_skewness(polynomial<cxx_mpfr>(f));
            return static_cast<double>(s);
#else
        } catch (...) {
            throw;
#endif
        }
    }
}

double L2_skew_lognorm (mpz_poly_srcptr f)
{
    try {
        const auto F = polynomial<double>(f);
        return L2_lognorm (F, L2_skewness (F));
    } catch (std::overflow_error const &) {
        try {
            const auto F = polynomial<long double>(f);
            const auto l = L2_lognorm (F, L2_skewness (F));
            return static_cast<double>(l);
#ifdef HAVE_MPFR
        } catch (std::overflow_error const &) {
            const auto F = polynomial<cxx_mpfr>(f);
            const auto l = L2_lognorm (F, L2_skewness (F));
            return static_cast<double>(l);
#else
        } catch (...) {
            throw;
#endif
        }
    }
}

double L2_combined_skewness2 (mpz_poly_srcptr f, mpz_poly_srcptr g)
{
    try {
        const polynomial<double> F(f);
        const polynomial<double> G(g);
        return L2_combined_skewness2(F, G);
    } catch (std::overflow_error const &) {
        try {
            const polynomial<long double> F(f);
            const polynomial<long double> G(g);
            return static_cast<double>(L2_combined_skewness2(F, G));
#ifdef HAVE_MPFR
        } catch (std::overflow_error const &) {
            const auto F = polynomial<cxx_mpfr>(f);
            const auto G = polynomial<cxx_mpfr>(g);
            return static_cast<double>(L2_combined_skewness2(F, G));
#else
        } catch (...) {
            throw;
#endif
        }
    }
}
