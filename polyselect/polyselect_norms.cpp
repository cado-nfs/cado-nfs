#include "cado.h" // IWYU pragma: keep

#include <cfloat>
#include <cstdint>
#include <cstdlib>
#include <cmath>

#include "double_poly.h"
#include "macros.h"
#include "polyselect_norms.h"

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

static double
L2_lognorm_d (double_poly_srcptr p, double s)
{
    const double * a = p->coeff;
    const int d = p->deg;
    ASSERT_ALWAYS(1 <= p->deg);
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

    const double invs = 1.0 / s;
    const double s2 = s * s;
    const double is2 = invs * invs;
    double tt = 0;
    double u0 = (d & 1) ? invs : 1;
    double u1 = (d & 1) ? s : 1;
    const int p0 = d / 2;
    const int p1 = d - p0;
    for(int k = 0 ; k <= d/2 ; k++) {
        double t0 = 0;
        double t1 = 0;
        for(int i = 1 ; i <= d/2 - k ; i++) {
            t0 += a[p0-k-i] * a[p0-k+i];
            t1 += a[p1+k-i] * a[p1+k+i];
        }
        t0 = (2 * t0 + a[p0-k] * a[p0-k]) * u0;
        t1 = (2 * t1 + a[p1+k] * a[p1+k]) * u1;
        u0 *= is2;
        u1 *= s2;
        tt += (coeffs_integral[d][k])*(t0 + t1);
    }
    tt = ldexp(tt * M_PI / (double) (d+1), -2*d);
    if (std::isnan (tt) || std::isinf (tt))
        tt = DBL_MAX;
    return log(tt) / 2;
}

double
L2_lognorm (mpz_poly_srcptr f, double s)
{
    /* we used to have code doing this computation over the integers as
     * well. It's easy to add, given the 20-line code above.
     */
    double res;
    double_poly a;
    double_poly_init(a, f->deg);
    double_poly_set_mpz_poly (a, f);

    res = L2_lognorm_d (a, s);

    double_poly_clear(a);
    return res;
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

static void L2_skewness_derivative_numerator(double_poly_ptr dP, double_poly_srcptr p)
{
    const double * a = p->coeff;
    const int d = p->deg;
    ASSERT_ALWAYS(1 <= p->deg);
    ASSERT_ALWAYS(d <= L2_DMAX);
    double_poly_realloc(dP, d+1);
    for(int k = 0 ; k <= d ; k++) {
        double t0 = 0;
        for(int i = 1 ; i <= k && i <= d - k ; i++)
            t0 += a[k-i] * a[k+i];
        t0 = (2 * t0 + a[k] * a[k]);
        const int j = (ABS(d-2*k)-(d&1))/2;
        t0 = t0 * coeffs_integral[d][j]*(2*k-d);
        dP->coeff[k] = t0; // ldexp(t0/(d+1), -2*d);
    }
    double_poly_cleandeg(dP, d);
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
static void L2_skewness_numerator(double_poly_ptr dP, double_poly_srcptr p)
{
    const double * a = p->coeff;
    const int d = p->deg;
    ASSERT_ALWAYS(1 <= p->deg);
    ASSERT_ALWAYS(d <= L2_DMAX);
    double_poly_realloc(dP, d+1);
    for(int k = 0 ; k <= d ; k++) {
        double t0 = 0;
        for(int i = 1 ; i <= k && i <= d - k ; i++)
            t0 += a[k-i] * a[k+i];
        t0 = (2 * t0 + a[k] * a[k]);
        const int j = (ABS(d-2*k)-(d&1))/2;
        t0 = t0 * coeffs_integral[d][j];
        if (2*k != d)
            t0 /= 2;
        dP->coeff[k] = t0; // ldexp(t0/(d+1), -(2*d-1));
    }
    double_poly_cleandeg(dP, d);
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
double
L2_combined_skewness2 (mpz_poly_srcptr f, mpz_poly_srcptr g)
{
    double_poly dQF, QF, F;
    double_poly dQG, QG, G;
    double_poly_init(F, -1);
    double_poly_init(QF, -1);
    double_poly_init(dQF, -1);
    double_poly_init(G, -1);
    double_poly_init(QG, -1);
    double_poly_init(dQG, -1);
    double_poly_set_mpz_poly (F, f);
    double_poly_set_mpz_poly (G, g);
    L2_skewness_numerator(QF, F);
    L2_skewness_derivative_numerator(dQF, F);
    L2_skewness_numerator(QG, G);
    L2_skewness_derivative_numerator(dQG, G);
    double_poly t, S;
    double_poly_init(t, -1);
    double_poly_init(S, -1);
    double_poly_mul(S, dQF, QG);
    double_poly_mul(t, QF, dQG);
    double_poly_add(S, S, t);

    double * roots = (double *) malloc(S->deg * sizeof(double));

    const double B = double_poly_bound_roots(S);
    /* some pathological examples for gfpn have integer roots, and this
     * wreaks havoc.
     */
    const unsigned int nroots = double_poly_compute_roots(roots, S, B + 1);

    double best_s = sqrt(roots[0]);
    if (nroots > 1) {
        double l = L2_lognorm_d(F, best_s) + L2_lognorm_d(G, best_s);
        for(unsigned int i = 1 ; i < nroots ; i++) {
            const double s = sqrt(roots[i]);
            const double li = L2_lognorm_d(F, s) + L2_lognorm_d(G, s);
            if (li < l) {
                best_s = s;
                l = li;
            }
        }
    }

    free(roots);

    double_poly_clear(S);
    double_poly_clear(t);
    double_poly_clear(F);
    double_poly_clear(QF);
    double_poly_clear(dQF);
    double_poly_clear(G);
    double_poly_clear(QG);
    double_poly_clear(dQG);

    return best_s;
}

double L2_skewness (mpz_poly_srcptr f)
{
    const int d = f->deg;
    double_poly dP, P;
    double_poly_init(P, -1);
    double_poly_init(dP, -1);
    double_poly_set_mpz_poly (P, f);
    L2_skewness_derivative_numerator(dP, P);

    double * roots = (double *) malloc(d * sizeof(double));

    const double B = double_poly_bound_roots(dP);
    /* some pathological examples for gfpn have integer roots, and this
     * wreaks havoc.
     */
    const unsigned int nroots = double_poly_compute_roots(roots, dP, B + 1);
    /* We often have a single zero, but not always. Here's an example
     * with three. It typically happens when we have pairs of roots that
     * are almost symmetrical along the axes. I'm not sure I see a
     * generalization of this phenomenon to more than two extrema, actually.
     * 3520*x^4 - 14848*x^3 - 26038239232*x^2 + 3120*x + 110
     */
    ASSERT_ALWAYS(nroots >= 1);

    double best_s = sqrt(roots[0]);
    if (nroots > 1) {
        double l = L2_lognorm_d(P, best_s);
        for(unsigned int i = 1 ; i < nroots ; i++) {
            const double s = sqrt(roots[i]);
            const double li = L2_lognorm_d(P, s);
            if (li < l) {
                best_s = s;
                l = li;
            }
        }
    }


    free(roots);
    double_poly_clear(P);
    double_poly_clear(dP);

    return best_s;
}

double L2_skew_lognorm (mpz_poly_srcptr f)
{
  return L2_lognorm (f, L2_skewness (f));
}
