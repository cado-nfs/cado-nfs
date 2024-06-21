#include "cado.h"
#include <math.h>
#include "random_distributions.hpp"

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

double random_normal_standard(gmp_randstate_t rstate)
{
    /* Box-Muller transform */
#ifdef  ALLOW_NON_REENTRANT_random_normal_standard
    static int last = 0;
    static double vlast = 0;
    if (last) { --last; return vlast; }
#endif
    double u = random_uniform(rstate);
    double v = random_uniform(rstate);
    double rho = sqrt(-2*log(u));
    double theta = 2*M_PI*v;
    double x = rho*cos(theta);
#ifdef ALLOW_NON_REENTRANT_random_normal_standard
    double y = rho*sin(theta);
    vlast=y;
    last++;
#endif
    return x;
}

double random_normal(gmp_randstate_t rstate, double mean, double sdev)
{
    return mean + sdev * random_normal_standard(rstate);
}

/* return the expected maximum among n picks for a normal low with given
 * mean and sdev */
double extreme_normal(double n, double mean, double sdev)
{
    /* See Cramer, Mathematical methods of statistics, p. 575 (7th printing). */
    double x = sqrt(2*log(n))+(log(log(n))+log(4*M_PI)-2*0.5772)/(2*sqrt(2*log(n)));
    return mean + x * sdev;
}

double random_normal_constrained(gmp_randstate_t rstate, double mean, double sdev, double a, double b)
{
    /* By cutting the tail, we're doing a real heresy. This increases the
     * average and standard deviation significantly. See the function
     * below. In fact, the normal approximation will never be useful in
     * our case.
     */
    for(;;) {
        double x = round(random_normal(rstate, mean, sdev));
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
void accuracy_of_normal_approximation_to_binomial(double * my, double *mx, unsigned long a, unsigned long b)
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
    double m = mx[0];
    double s = mx[1];
    double as = (a-m)/s;
    double bs = (b-m)/s;
    double eas = exp(-as*as/2)/sqrt(2*M_PI);
    double ebs = exp(-bs*bs/2)/sqrt(2*M_PI);
    double S0 = (erf(bs/sqrt(2))-erf(as/sqrt(2)))/2;
    double S1 = eas - ebs;
    double S2 = as*eas - bs*ebs + S0;
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

double random_poisson(gmp_randstate_t rstate, double lambda)
{
    /* "method PA" from "The Computer Generation of Poisson Random
     * Variables" by A. C. Atkinson, Journal of the Royal Statistical
     * Society Series C (Applied Statistics) Vol. 28, No. 1. (1979),
     * pages 29-35.
     */
    if (lambda < 10) {
        return random_uniform(rstate)*2*lambda;
    }
    double c = 0.767 - 3.36/lambda;
    double beta = M_PI/sqrt(3.0*lambda);
    double alpha = beta*lambda;
    double k = log(c) - lambda - log(beta);

    for(;;) {
        double u = random_uniform(rstate);
        double x = (alpha - log((1.0 - u)/u))/beta;
        int n = (int) floor(x + 0.5);
        if (n < 0)
            continue;
        double v = random_uniform(rstate);
        double y = alpha - beta*x;
        double temp = 1.0 + exp(y);
        double lhs = y + log(v/(temp*temp));
        double rhs = k + n*log(lambda) - lgamma(n-1);
        if (lhs <= rhs)
            return n;
    }
}

/* This is the random variable associated to the *size* of the sample */
double random_binomial(gmp_randstate_t rstate, unsigned long n, double p)
{
    /*
     * This first way of doing things is appropriate when mean \pm 3
     * times sdev is good.
     */
    double mean = n * p;
    double sdev = sqrt(n * p * (1-p));
    if (0 <= mean - 3 * sdev && mean + 3 * sdev <= n) {
        return random_normal_constrained(rstate, mean, sdev, 0, n);
    }
    /* otherwise we'll return the Poisson approximation, which does not
     * care much about the standard deviation, but matches relatively
     * well as far as our application is concerned. */

    double r;
    for( ; (r = random_poisson(rstate, mean)) >= n ; ) ;
    return r;
}


