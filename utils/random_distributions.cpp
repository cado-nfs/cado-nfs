#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <climits>
#include <cstring>

#include <utility>

#include <gmp.h>

#include "gmp_aux.h"
#include "random_distributions.hpp"
#include "macros.h"

#ifndef M_PI
#error "Uh, on a posix platform, we should have M_PI defined"
#endif

/* random picking */
double random_uniform(gmp_randstate_t rstate)
{
    /* the constant is 2^-53 */
    // return gmp_urandomm_ui(rstate, 1UL<<53) * 1.11022302462515654042363166809E-16;
    return gmp_urandomm_ui(rstate, ULONG_MAX) / (double) ULONG_MAX;
    /*
    mpf_t x;
    mpf_init2(x, 53);
    mpf_urandomb(x, rstate, 53);
    double y = mpf_get_d(x);
    mpf_clear(x);
    return y;
    */
}

double random_normal_standard(cxx_gmp_randstate & rstate)
{
    /* Box-Muller transform */
#ifdef  ALLOW_NON_REENTRANT_random_normal_standard
    static int last = 0;
    static double vlast = 0;
    if (last) { --last; return vlast; }
#endif
    const double u = random_uniform(rstate);
    const double v = random_uniform(rstate);
    const double rho = sqrt(-2*log(u));
    const double theta = 2*M_PI*v;
    const double x = rho*cos(theta);
#ifdef ALLOW_NON_REENTRANT_random_normal_standard
    const double y = rho*sin(theta);
    vlast=y;
    last++;
#endif
    return x;
}

double random_normal(cxx_gmp_randstate & rstate, double mean, double sdev)
{
    return mean + sdev * random_normal_standard(rstate);
}

/* return the expected maximum among n picks for a normal low with given
 * mean and sdev */
double extreme_normal(double n, double mean, double sdev)
{
    /* See Cramer, Mathematical methods of statistics, p. 575 (7th printing). */
    const double x = sqrt(2*log(n))+(log(log(n))+log(4*M_PI)-2*0.5772)/(2*sqrt(2*log(n)));
    return mean + x * sdev;
}

double random_normal_constrained(cxx_gmp_randstate & rstate, double mean, double sdev, double a, double b)
{
    /* By cutting the tail, we're doing a real heresy. This increases the
     * average and standard deviation significantly. See the function
     * below. In fact, the normal approximation will never be useful in
     * our case.
     */
    for(;;) {
        const double x = round(random_normal(rstate, mean, sdev));
        if (x >= a && x < b) return x;
    }
}

/* given a probability mass function which gives the gaussian with mean
 * and sdev given my mx[0] and mx[1], but truncated to the interval
 * [a,b[, return the mean and sdev of the resulting distribution.
 *
 * This is just an illustration, which can be used to witness how the
 * normal approximation can end up being catastrophic if we're more
 * Poisson-like.
 */
void accuracy_of_normal_approximation_to_binomial(double * my, const double *mx, unsigned long a, unsigned long b)
{
    /* Let e(x) = 1/sqrt(2*pi)*exp(-x^2/2).
     * We have:
     *  e'(x) = -x e(x).
     * \int_{-\infty}^{+\infty} e(x) = 1
     * \int_{-\infty}^{+\infty} xe(x) = 0
     * \int_{-\infty}^{+\infty} x^2e(x) = 1
     *
     * and more generally:
     * \int_a^b e(x) = 1/2*(erf(b/sqrt(2))-erf(a/sqrt(2))) = S0(a,b)
     * \int_a^b xe(x) = e(a) - e(b)                        = S1(a,b)
     * \int_a^b x^2e(x) = a*e(a)-b*e(b)+S0(a,b)            = S2(a,b)
     *
     * Let now e*(x) = 1/s*e((x-m)/s), pdf of a gaussian with mean and sdev
     * equal to m and s. Let x*=(x-m)/s, a*=(a-m)/s, and b*=(b-m)/s.
     * So that dx = s d{x*} ; note that e*(x)dx = e(x*)d{x*}.
     *
     * We have:
     * M0(a,b) = \int_a^b e*(x) dx
     *         = \int_{a*}^{b*}e(x*)d{x*}
     *         = S0(a*,b*)
     * M1(a,b) = \int_a^b x e*(x) dx
     *         = \int_{a*}^{b*} (m+s*x*) e(x*)d{x*}
     *         = m*S0(a*,b*) + s*S1(a*,b*)
     * M2(a,b) = \int_a^b x^2 e*(x) dx
     *         = \int_{a*}^{b*} (m+s*x*)^2 e(x*) d{x*}
     *         = m^2*S0(a*,b*) + 2*m*s*S1(a*,b*) + s^2*S2(a*,b*)
     * 
     * when scaled, we get:
     *
     * M0 = 1
     * M1 = m + s * (S1/S0)(a*,b*)
     * M2 = m^2 + 2*m*s * (S1/S0)(a*,b*) + s^2 * (S2/S0)(a*,b*)
     * sdev = s * sqrt(((S2-S1^2)/S0)(a*,b*))
     */
    const double m = mx[0];
    const double s = mx[1];
    const double as = ((double) a-m)/s;
    const double bs = ((double) b-m)/s;
    const double eas = exp(-as*as/2)/sqrt(2*M_PI);
    const double ebs = exp(-bs*bs/2)/sqrt(2*M_PI);
    const double S0 = (erf(bs/sqrt(2))-erf(as/sqrt(2)))/2;
    const double S1 = eas - ebs;
    const double S2 = as*eas - bs*ebs + S0;
    /*
       double M0 = s * S0;
       double M1 = s * (m*S0 + s*S1);
       double M2 = s * (m^2*S0 + 2*m*s*S1 + s^2*S2);
       */
    // double M0 = 1;
    // double M1 = (m + s*S1/S0);
    // double M2 = (m^2 + 2*m*s*S1/S0 + s^2*S2/S0);
    my[0] = m + s * S1/S0;
    my[1] = s * sqrt((S2-S1*S1)/S0);
}

double random_poisson(cxx_gmp_randstate & rstate, double lambda)
{
    /* "method PA" from "The Computer Generation of Poisson Random
     * Variables" by A. C. Atkinson, Journal of the Royal Statistical
     * Society Series C (Applied Statistics) Vol. 28, No. 1. (1979),
     * pages 29-35.
     */
    if (lambda < 10) {
        return random_uniform(rstate)*2*lambda;
    }
    const double c = 0.767 - 3.36/lambda;
    const double beta = M_PI/sqrt(3.0*lambda);
    const double alpha = beta*lambda;
    const double k = log(c) - lambda - log(beta);

    for(;;) {
        const double u = random_uniform(rstate);
        const double x = (alpha - log((1.0 - u)/u))/beta;
        const int n = (int) floor(x + 0.5);
        if (n < 0)
            continue;
        const double v = random_uniform(rstate);
        const double y = alpha - beta*x;
        const double temp = 1.0 + exp(y);
        const double lhs = y + log(v/(temp*temp));
        const double rhs = k + n*log(lambda) - lgamma(n-1);
        if (lhs <= rhs)
            return n;
    }
}

/* This is the random variable associated to the *size* of the sample */
double random_binomial(cxx_gmp_randstate & rstate, unsigned long n, double p)
{
    /*
     * This first way of doing things is appropriate when mean \pm 3
     * times sdev is good.
     */
    const double mean = (double) n * p;
    const double sdev = sqrt((double) n * p * (1-p));
    if (0 <= mean - 3 * sdev && mean + 3 * sdev <= (double) n) {
        return random_normal_constrained(rstate, mean, sdev, 0, (double) n);
    }
    /* otherwise we'll return the Poisson approximation, which does not
     * care much about the standard deviation, but matches relatively
     * well as far as our application is concerned. */

    double r;
    for( ; (r = random_poisson(rstate, mean)) >= (double) n ; ) ;
    return r;
}

void punched_interval::recycle(node_t && c, pool_t & pool)
{
    if (!c) return;
    // fprintf(stderr, "STOW %p\n", c);
    /* enqueue both children to the free pool */
    recycle(std::move(c->left), pool);
    recycle(std::move(c->right), pool);
    /* also store the count */
    c->has_left = 1 + (pool ? pool->has_left : 0);
    c->left = std::move(pool);
    pool = std::move(c);
}

punched_interval::node_t punched_interval::alloc(pool_t & pool, double b0, double b1)
{
    node_t x;
    if (pool) {
        x = std::move(pool);
        pool = std::move(x->left);
    } else {
        x.reset(new punched_interval);
    }

    x->b0 = b0;
    x->b1 = b1;
    x->has_left = 0;
    x->holes = 0;
    ASSERT_ALWAYS(!x->left);
    ASSERT_ALWAYS(!x->right);
    return x;
}

void punched_interval::punch_inner(pool_t & pool, double x0, double x1)
{
    /* This function is misleading. It is only ever called on full
     * intervals, so that we always have the following...
     */
    ASSERT_ALWAYS(holes == 0);
    ASSERT_ALWAYS(!has_left);
    ASSERT_ALWAYS(!left);
    ASSERT_ALWAYS(!right);

    holes += x1 - x0;
    left = alloc(pool, b0, x0);
    has_left=1;
    right = alloc(pool, x1, b1);

}

unsigned long punched_interval::pick_inner(pool_t & pool,
        matrix_column_distribution const & D,
        double x)
{
    /* x should be within [c->b0, c->b1 - c->holes] */
    ASSERT_ALWAYS(x >= b0);
    ASSERT_ALWAYS(x + holes < b1);
    if (!has_left) {
        /* no holes ! */
        double const r = D.qrev(x);
        unsigned long i;
        if (r < 0) {
            i = 0;
        } else {
            i = floor(r);
        }
        const double x0 = D.q(double(i));
        const double x1 = D.q(double(i) + 1);
        punch_inner(pool, x0, x1);
        return i;
    }
    /* try to correct x with all left holes */
    double xc = x + left->holes;
    if (xc < left->b1) {
        double const h = left->holes;
        const unsigned long i = left->pick_inner(pool, D, x);
        holes += left->holes - h;
        return i;
    } else {
        /* modify x. It's more than just xc ! */
        xc += right->b0 - left->b1;
        double const h = right->holes;
        const unsigned long i = right->pick_inner(pool, D, xc);
        holes += right->holes - h;
        return i;
    }
}

unsigned long punched_interval::pick(pool_t & pool,
        matrix_column_distribution const & D,
        cxx_gmp_randstate & rstate)
{
    double const x = random_uniform(rstate) * (b1 - holes);
    return pick_inner(pool, D, x);
}


void punched_interval::print_rec(FILE * f) const
{
    if (!has_left) return;
    left->print_rec(f);
    fprintf(f, "(\e[31m%.4f...%.4f\e[0m)...", left->b1, right->b0);
    right->print_rec(f);
}

void punched_interval::print(FILE * f) const
{
    fprintf(f, "%.4f...", b0);
    print_rec(f);
    fprintf(f, "%.4f\n", b1);
}
