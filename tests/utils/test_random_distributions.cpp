#include "cado.h"
#include "macros.h"
#include "random_distributions.hpp"
#include "tests_common.h"
#include <algorithm>
#include <cmath>
#include <vector>

/* Some tests */
void test_random_normal_standard(gmp_randstate_t rstate, unsigned long N,
                                 int logscale_report)
{
    double s = 0, ss = 0;
    for (unsigned long i = 0, l = logscale_report; l <= N;
         l *= logscale_report) {
        for (; i < l; i++) {
            double const x = random_normal_standard(rstate);
            s += x;
            ss += x * x;
        }
        double const m = s / double(l);
        double const sd = sqrt(ss / double(l) - m * m);
        fprintf(stderr, "%s: after %lu picks, mean=%.3f sdev=%.3f\n", __func__,
                l, m, sd);

        if (l >= 1024) {
            ASSERT_ALWAYS(abs(m) <= 0.1);
            ASSERT_ALWAYS(abs(sd - 1) <= 0.1);
            fprintf(stderr, "%s: checks passed\n", __func__);
        }
    }
}

void test_random_normal(gmp_randstate_t rstate, double xm, double xs,
                        unsigned long N, int logscale_report)
{
    double s = 0, ss = 0;
    for (unsigned long i = 0, l = logscale_report; l <= N;
         l *= logscale_report) {
        for (; i < l; i++) {
            double const x = random_normal(rstate, xm, xs);
            s += x;
            ss += x * x;
        }
        double const m = s / double(l);
        double const sd = sqrt(ss / double(l) - m * m);
        fprintf(stderr, "%s: after %lu picks, mean=%.3f sdev=%.3f\n", __func__,
                l, m, sd);
        if (l >= 1024) {
            ASSERT_ALWAYS(abs(m - xm) <= 0.1 * xs);
            ASSERT_ALWAYS(abs(sd - xs) <= 0.1 * xs);
            fprintf(stderr, "%s: checks passed\n", __func__);
        }
    }
}

void test_random_normal_constrained(gmp_randstate_t rstate, double xm,
                                    double xs, unsigned long a, unsigned long b,
                                    unsigned long N, int logscale_report)
{
    double s = 0, ss = 0;
    double mmx[2] = {xm, xs}, mmy[2];
    accuracy_of_normal_approximation_to_binomial(mmy, mmx, a, b);
    fprintf(stderr,
            "%s: want (%.3f,%.3f), expect instead (%.3f,%.3f) when truncating "
            "to [%lu,%lu]\n",
            __func__, mmx[0], mmx[1], mmy[0], mmy[1], a, b);
    for (unsigned long i = 0, l = logscale_report; l <= N;
         l *= logscale_report) {
        for (; i < l; i++) {
            double const x =
                random_normal_constrained(rstate, xm, xs, double(a), double(b));
            s += x;
            ss += x * x;
        }
        double const m = s / double(l);
        double const sd = sqrt(ss / double(l) - m * m);
        fprintf(stderr, "%s: after %lu picks, mean=%.3f sdev=%.3f\n", __func__,
                l, m, sd);
        if (l >= 1024) {
            ASSERT_ALWAYS(abs(m - mmy[0]) <= 0.1 * mmy[1]);
            ASSERT_ALWAYS(abs(sd - mmy[1]) <= 0.1 * mmy[1]);
            fprintf(stderr, "%s: checks passed\n", __func__);
        }
    }
}

/* Do random picks of a poissn distributed random variable with mean xm,
 * and pdf Pr(X=k)=xm^k/k!*exp(-xm), but restrict the output to values
 * that are at most n
 */
void test_random_poisson(gmp_randstate_t rstate, double xm, unsigned long n,
                         unsigned long N, int logscale_report)
{
    double s = 0, ss = 0;
    for (unsigned long i = 0, l = logscale_report; l <= N;
         l *= logscale_report) {
        for (; i < l; i++) {
            unsigned long x;
            for (; (x = (unsigned long)random_poisson(rstate, xm)) >= n;)
                ;
            s += double(x);
            ss += double(x) * double(x);
        }
        double const m = s / double(l);
        double const sd = sqrt(ss / double(l) - m * m);
        fprintf(stderr, "%s: after %lu picks, mean=%.3f sdev=%.3f\n", __func__,
                l, m, sd);
        if (l >= 1024) {
            ASSERT_ALWAYS(abs(m / xm - 1) <= 0.1);
            ASSERT_ALWAYS(abs(sd / sqrt(xm) - 1) <= 0.1);
            fprintf(stderr, "%s: checks passed\n", __func__);
        }
    }
}

struct cdf_pow_description {
    size_t N;
    double e;
};

double cdf_pow(void const * f, double x)
{
    auto const * C = static_cast<const struct cdf_pow_description *>(f);
    return pow(x / double(C->N), C->e);
}
double cdf_invpow(void const * f, double x)
{
    auto const * C = static_cast<const struct cdf_pow_description *>(f);
    return pow(x, 1.0 / C->e) * double(C->N);
}

void test_arbitrary(gmp_randstate_ptr rstate, unsigned long iter)
{
    for (int ee = 1; ee < 6; ee++) {
        struct cdf_pow_description C {
            8192, double(ee)
        };

        size_t K = ceil(sqrt(double(C.N)));
        for (; K * K < C.N; K++)
            ;

        std::vector<int> bins(K, 0);

        if (iter < K) {
            fprintf(stderr, "too few iterations, check skipped");
            continue;
        }

        for (unsigned long j = 0; j < iter / K; j++) {
            punched_interval_ptr pool = nullptr;
            punched_interval_ptr range = punched_interval_alloc(&pool, 0, 1);
            for (unsigned long i = 0; i < K; i++) {
                unsigned long const a = punched_interval_pick(
                    &pool, range, (double (*)(void const *, double))cdf_pow,
                    (double (*)(void const *, double))cdf_invpow,
                    (void const *)&C, rstate);
                ASSERT_ALWAYS(a < C.N);
                bins[a / K]++;
            }
            /* just for fun, print the contents after the first round */
            if (j == 0)
                punched_interval_print(stdout, range);
            punched_interval_free(range, &pool);
            punched_interval_free_pool(&pool);
        }

        unsigned long const npicks = (iter / K) * K;
        double d = 0;

        double s0 = 0, s1 = 0, s2 = 0;
        for (size_t i = 0; i < K; i++) {
            double const nd = cdf_pow(&C, double(std::min((i + 1) * K, C.N)));
            int const got = bins[i];
            double const expect = double(npicks) * (nd - d);
            // printf("%zu %d %.1f\n", i, got, expect);
            s0++;
            s1 += (got - expect);
            s2 += (got - expect) * (got - expect);
            d = nd;
        }
        double const mean = s1 / s0,
                     sdev = sqrt(s2 / s0 - (s1 / s0) * (s1 / s0));
        printf("[[e=%.1f]] After %lu rounds of filling %zu bins with %zu "
               "random picks (which means %.1f in each bin on average), "
               "comparison with expectation: mean=%.1f, sdev=%.1f\n",
               C.e, iter / K, K, K, double(npicks) / double(K), mean, sdev);

        ASSERT_ALWAYS(mean < 10);
        ASSERT_ALWAYS(sdev < (sqrt(iter) * 10));
        printf("checks passed\n");
    }
}

int main(int argc, char const * argv[])
{
    unsigned long iter = 20000;

    tests_common_cmdline(&argc, &argv, PARSE_SEED | PARSE_ITER);
    tests_common_get_iter(&iter);

    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);

    test_random_normal_standard(rstate, iter, 2);
    test_random_normal(rstate, 4242.0, 1717.0, iter, 2);
    test_random_normal_constrained(rstate, 4242.0, 1717.0, 0, ULONG_MAX, iter,
                                   2); // F->ncols);
    test_random_poisson(rstate, 4242.0, ULONG_MAX, iter, 2);

    test_arbitrary(rstate, iter);

    gmp_randclear(rstate);
    tests_common_clear();
    return EXIT_SUCCESS;
}
