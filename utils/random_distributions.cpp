#include "cado.h"
#include <stdlib.h>
#include <math.h>
#include "random_distributions.hpp"
#include "macros.h"

/* random picking */
double random_uniform(cxx_gmp_randstate & rstate)
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

double random_normal(cxx_gmp_randstate & rstate, double mean, double sdev)
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

double random_normal_constrained(cxx_gmp_randstate & rstate, double mean, double sdev, double a, double b)
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
double random_binomial(cxx_gmp_randstate & rstate, unsigned long n, double p)
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

void punched_interval_free(punched_interval_ptr c, punched_interval_ptr * pool)
{
    if (!c) return;
    // fprintf(stderr, "STOW %p\n", c);
    /* enqueue both children to the free pool */
    punched_interval_free(c->left, pool);
    punched_interval_free(c->right, pool);
    c->left = *pool;
    /* also store the count */
    c->has_left = 1 + ((*pool) ? (*pool)->has_left : 0);
    *pool = c;
}

punched_interval_ptr punched_interval_alloc(punched_interval_ptr * pool, double b0, double b1)
{
    punched_interval_ptr x;
    if (*pool) {
        x = *pool;
        // fprintf(stderr, "REUSE %p\n", x);
        *pool = x->left;
    } else {
        x = (struct punched_interval_s *) malloc(sizeof(struct punched_interval_s));
        // fprintf(stderr, "ALLOC %p\n", x);
    }
    memset(x, 0, sizeof(struct punched_interval_s));

    x->b0 = b0;
    x->b1 = b1;
    x->has_left = 0;
    x->holes = 0;
    x->left = nullptr;
    x->right = nullptr;
    return x;
}

void punched_interval_free_pool(punched_interval_ptr * pool)
{
    for(punched_interval_ptr q = *pool, v ; q ; q = v) {
        v = q->left;
        // fprintf(stderr, "FREE %p\n", q);
        free(q);
    }
    *pool = nullptr;
}

void punched_interval_pre_free_pool(punched_interval_ptr * pool, int max, int print)
{
    if (!*pool) return;
    if ((*pool)->has_left < 2 * max) return;
    if (print) {
        fprintf(stderr, "Reducing punched_interval pool from size %d to %d\n",
                (*pool)->has_left, max);
    }
    punched_interval_ptr q = * pool;
    int size = (*pool)->has_left;
    for(int i = 0 ; q->has_left >= max ; i++) {
        ASSERT_ALWAYS(q->left);
        ASSERT_ALWAYS(q->has_left == size - i);
        punched_interval_ptr nq = q->left;
        free(q);
        q = nq;
    }
    *pool = q;
}


static void punched_interval_punch_inner(punched_interval_ptr * pool, punched_interval_ptr c, double x0, double x1)
{
    /* This function is misleading. It is only ever called on full
     * intervals, so that we always have the following...
     */
    ASSERT_ALWAYS(c->holes == 0);
    ASSERT_ALWAYS(!c->has_left);
    ASSERT_ALWAYS(!c->left);
    ASSERT_ALWAYS(!c->right);

    c->holes += x1 - x0;
    c->left = punched_interval_alloc(pool, c->b0, x0);
    c->has_left=1;
    c->right = punched_interval_alloc(pool, x1, c->b1);

    /* A "longer" version exists here, but I doubt it's correct.
     * The code below has branches that are not callable, and whether
     * they make any sense is really not clear at all. I should probably
     * get rid of them.
     */
#if 0
    c->holes += x1 - x0;
    if (!c->left) {
        c->left = punched_interval_alloc(pool, c->b0, x0);
    } else {
        punched_interval_set_full(c->left, c->b0, x0);
    }
    c->has_left=1;
    if (!c->right) {
        c->right = punched_interval_alloc(pool, x1, c->b1);
    } else {
        punched_interval_set_full(c->right, x1, c->b1);
    }
#endif
}

static unsigned long punched_interval_pick_inner(punched_interval_ptr * pool, punched_interval_ptr c,
        matrix_column_distribution const & D,
        double x)
{
    /* x should be within [c->b0, c->b1 - c->holes] */
    ASSERT_ALWAYS(x >= c->b0);
    ASSERT_ALWAYS(x + c->holes < c->b1);
    if (!c->has_left) {
        /* no holes ! */
        double const r = D.qrev(x);
        unsigned long i;
        if (r < 0) {
            i = 0;
        } else {
            i = floor(r);
        }
        double x0 = D.q(i);
        double x1 = D.q(i + 1);
        punched_interval_punch_inner(pool, c, x0, x1);
        return i;
    }
    /* try to correct x with all left holes */
    double xc = x + c->left->holes;
    if (xc < c->left->b1) {
        double const h = c->left->holes;
        unsigned long i = punched_interval_pick_inner(pool, c->left, D, x);
        c->holes += c->left->holes - h;
        return i;
    } else {
        /* modify x. It's more than just xc ! */
        xc += c->right->b0 - c->left->b1;
        double const h = c->right->holes;
        unsigned long i = punched_interval_pick_inner(pool, c->right, D, xc);
        c->holes += c->right->holes - h;
        return i;
    }
}

unsigned long punched_interval_pick(punched_interval_ptr * pool, punched_interval_ptr c,
        matrix_column_distribution const & D,
        cxx_gmp_randstate & rstate)
{
    double const x = random_uniform(rstate) * (c->b1 - c->holes);
    return punched_interval_pick_inner(pool, c, D, x);
}


void punched_interval_print_rec(FILE * f, punched_interval_ptr c)
{
    if (!c->has_left) return;
    punched_interval_print_rec(f, c->left);
    fprintf(f, "(\e[31m%.4f...%.4f\e[0m)...", c->left->b1, c->right->b0);
    punched_interval_print_rec(f, c->right);
}

void punched_interval_print(FILE * f, punched_interval_ptr c)
{
    fprintf(f, "%.4f...", c->b0);
    punched_interval_print_rec(f, c);
    fprintf(f, "%.4f\n", c->b1);
}
