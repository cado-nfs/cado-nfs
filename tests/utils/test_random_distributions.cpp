#include "cado.h"
#include "tests_common.h"
#include "random_distributions.hpp"
#include "macros.h"
#include <math.h>

/*  Some tests */
void test_random_normal_standard(gmp_randstate_t rstate, unsigned long N, int logscale_report)
{
    double s = 0, ss = 0;
    for(unsigned long i = 0, l = logscale_report ; l <= N ; l *= logscale_report) {
        for(  ; i < l ; i++) {
            double x = random_normal_standard(rstate);
            s += x;
            ss += x * x;
        }
        double m = s / l;
        double sd = sqrt(ss / l - m * m);
        fprintf(stderr, "%s: after %lu picks, mean=%.3f sdev=%.3f\n",
                __func__,
                l, m, sd);

        if (l >= 1024) {
            ASSERT_ALWAYS(abs(m) <= 0.1);
            ASSERT_ALWAYS(abs(sd-1) <= 0.1);
            fprintf(stderr, "%s: checks passed\n", __func__);
        }
    }
}

void test_random_normal(gmp_randstate_t rstate, double xm, double xs, unsigned long N, int logscale_report)
{
    double s = 0, ss = 0;
    for(unsigned long i = 0, l = logscale_report ; l <= N ; l *= logscale_report) {
        for(  ; i < l ; i++) {
            double x = random_normal(rstate, xm, xs);
            s += x;
            ss += x * x;
        }
        double m = s / l;
        double sd = sqrt(ss / l - m * m);
        fprintf(stderr, "%s: after %lu picks, mean=%.3f sdev=%.3f\n",
                __func__,
                l, m, sd);
        if (l >= 1024) {
            ASSERT_ALWAYS(abs(m-xm) <= 0.1 * xs);
            ASSERT_ALWAYS(abs(sd-xs) <= 0.1 * xs);
            fprintf(stderr, "%s: checks passed\n", __func__);
        }
    }
}

void test_random_normal_constrained(gmp_randstate_t rstate, double xm, double xs, unsigned long a, unsigned long b, unsigned long N, int logscale_report)
{
    double s = 0, ss = 0;
    double mmx[2]={xm, xs}, mmy[2];
    accuracy_of_normal_approximation_to_binomial(mmy, mmx, a, b);
    fprintf(stderr, "%s: want (%.3f,%.3f), expect instead (%.3f,%.3f) when truncating to [%lu,%lu]\n",
            __func__,
            mmx[0], mmx[1],
            mmy[0], mmy[1],
            a, b);
    for(unsigned long i = 0, l = logscale_report ; l <= N ; l *= logscale_report) {
        for(  ; i < l ; i++) {
            double x = random_normal_constrained(rstate, xm, xs, a, b);
            s += x;
            ss += x * x;
        }
        double m = s / l;
        double sd = sqrt(ss / l - m * m);
        fprintf(stderr, "%s: after %lu picks, mean=%.3f sdev=%.3f\n",
                __func__,
                l, m, sd);
        if (l >= 1024) {
            ASSERT_ALWAYS(abs(m-mmy[0]) <= 0.1 * mmy[1]);
            ASSERT_ALWAYS(abs(sd-mmy[1]) <= 0.1 * mmy[1]);
            fprintf(stderr, "%s: checks passed\n", __func__);
        }
    }
}

/* Do random picks of a poissn distributed random variable with mean xm,
 * and pdf Pr(X=k)=xm^k/k!*exp(-xm), but restrict the output to values
 * that are at most n
 */
void test_random_poisson(gmp_randstate_t rstate, double xm, unsigned long n, unsigned long N, int logscale_report)
{
    double s = 0, ss = 0;
    for(unsigned long i = 0, l = logscale_report ; l <= N ; l *= logscale_report) {
        for(  ; i < l ; i++) {
            unsigned long x;
            for( ; (x = random_poisson(rstate, xm)) >= n ; ) ;
            s += x;
            ss += x * x;
        }
        double m = s / l;
        double sd = sqrt(ss / l - m * m);
        fprintf(stderr, "%s: after %lu picks, mean=%.3f sdev=%.3f\n",
                __func__,
                l, m, sd);
        if (l >= 1024) {
            ASSERT_ALWAYS(abs(m/xm-1) <= 0.1);
            ASSERT_ALWAYS(abs(sd/sqrt(xm)-1) <= 0.1);
            fprintf(stderr, "%s: checks passed\n", __func__);
        }
    }
}

int main(int argc, const char * argv[])
{
  unsigned long iter = 20000;

  tests_common_cmdline (&argc, &argv, PARSE_SEED | PARSE_ITER);
  tests_common_get_iter (&iter);


  gmp_randstate_t rstate;
  gmp_randinit_default(rstate);

  test_random_normal_standard(rstate, iter, 2);
  test_random_normal(rstate, 4242.0, 1717.0, iter, 2);
  test_random_normal_constrained(rstate, 4242.0, 1717.0, 0, ULONG_MAX, iter, 2);//F->ncols);
  test_random_poisson(rstate, 4242.0, ULONG_MAX, iter, 2);

  gmp_randclear(rstate);
  tests_common_clear ();
  exit (EXIT_SUCCESS);
}
