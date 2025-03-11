#include "cado.h" // IWYU pragma: keep

#include <cstddef>
#include <cmath>
#include <climits>
#include <cstdlib>

#include <algorithm>
#include <utility>
#include <vector>

#include "fmt/base.h"

#include "gmp_aux.h"
#include "macros.h"
#include "random_distributions.hpp"
#include "tests_common.h"

/* Some tests */
static void test_random_normal_standard(cxx_gmp_randstate & rstate,
                                        unsigned long N, int logscale_report)
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
        fmt::print(stderr, "{}: after {} picks, mean={:.3f} sdev={:.3f}\n",
                   __func__, l, m, sd);

        if (l >= 1024) {
            ASSERT_ALWAYS(std::abs(m) <= 0.1);
            ASSERT_ALWAYS(std::abs(sd - 1) <= 0.1);
            fmt::print(stderr, "{}: checks passed\n", __func__);
        }
    }
}

static void test_random_normal(cxx_gmp_randstate & rstate, double xm, double xs,
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
        fmt::print(stderr, "{}: after {} picks, mean={:.3f} sdev={:.3f}\n", __func__,
                l, m, sd);
        if (l >= 1024) {
            ASSERT_ALWAYS(abs(m - xm) <= 0.1 * xs);
            ASSERT_ALWAYS(abs(sd - xs) <= 0.1 * xs);
            fmt::print(stderr, "{}: checks passed\n", __func__);
        }
    }
}

static void test_random_normal_constrained(cxx_gmp_randstate & rstate,
                                           double xm, double xs,
                                           unsigned long a, unsigned long b,
                                           unsigned long N, int logscale_report)
{
    double s = 0, ss = 0;
    double mmx[2] = {xm, xs}, mmy[2];
    accuracy_of_normal_approximation_to_binomial(mmy, mmx, a, b);
    fmt::print(stderr,
            "{}: want ({:.3f},{:.3f}), expect instead ({:.3f},{:.3f}) when truncating "
            "to [{},{}]\n",
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
        fmt::print(stderr, "{}: after {} picks, mean={:.3f} sdev={:.3f}\n", __func__,
                l, m, sd);
        if (l >= 1024) {
            ASSERT_ALWAYS(abs(m - mmy[0]) <= 0.1 * mmy[1]);
            ASSERT_ALWAYS(abs(sd - mmy[1]) <= 0.1 * mmy[1]);
            fmt::print(stderr, "{}: checks passed\n", __func__);
        }
    }
}

/* Do random picks of a poissn distributed random variable with mean xm,
 * and pdf Pr(X=k)=xm^k/k!*exp(-xm), but restrict the output to values
 * that are at most n
 */
static void test_random_poisson(cxx_gmp_randstate & rstate, double xm,
                                unsigned long n, unsigned long N,
                                int logscale_report)
{
    double s = 0, ss = 0;
    ASSERT_ALWAYS(logscale_report > 1);
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
        fmt::print(stderr, "{}: after {} picks, mean={:.3f} sdev={:.3f}\n", __func__,
                l, m, sd);
        if (l >= 1024) {
            ASSERT_ALWAYS(abs(m / xm - 1) <= 0.1);
            ASSERT_ALWAYS(abs(sd / sqrt(xm) - 1) <= 0.1);
            fmt::print(stderr, "{}: checks passed\n", __func__);
        }
    }
}

struct test_column_distribution : public matrix_column_distribution {
    size_t N;
    double e;
    test_column_distribution(size_t N, double e)
        : N(N)
        , e(e)
    {
    }
    double q(double x) const override { return pow(x / double(N), e); }
    double qrev(double x) const override
    {
        return pow(x, 1.0 / e) * double(N);
    }
};

static void test_arbitrary(cxx_gmp_randstate & rstate, unsigned long iter)
{
    for (int ee = 1; ee < 6; ee++) {
        test_column_distribution const C(8192, ee);

        size_t K = ceil(sqrt(double(C.N)));
        for (; K * K < C.N; K++)
            ;

        std::vector<int> bins(K, 0);

        if (iter < K) {
            fmt::print(stderr, "too few iterations, check skipped");
            continue;
        }

        for (unsigned long j = 0; j < iter / K; j++) {
            punched_interval::pool_t pool;
            auto range = punched_interval::alloc(pool, 0, 1);
            for (unsigned long i = 0; i < K; i++) {
                unsigned long const a = range->pick(pool, C, rstate);
                ASSERT_ALWAYS(a < C.N);
                bins[a / K]++;
            }
            /* just for fun, print the contents after the first round */
            if (j == 0)
                range->print(stdout);
            punched_interval::recycle(std::move(range), pool);
        }

        unsigned long const npicks = (iter / K) * K;
        double d = 0;

        double s0 = 0, s1 = 0, s2 = 0;
        for (size_t i = 0; i < K; i++) {
            double const nd = C.q(double(std::min((i + 1) * K, C.N)));
            int const got = bins[i];
            double const expect = double(npicks) * (nd - d);
            // fmt::print("{} {} {:.1f}\n", i, got, expect);
            s0++;
            s1 += (got - expect);
            s2 += (got - expect) * (got - expect);
            d = nd;
        }
        double const mean = s1 / s0,
                     sdev = sqrt(s2 / s0 - (s1 / s0) * (s1 / s0));
        fmt::print(
                "[[e={:.1f}]] After {} rounds of filling {} bins with {} "
                "random picks (which means {:.1f} in each bin on average), "
                "comparison with expectation: mean={:.1f}, sdev={:.1f}\n",
                C.e, iter / K, K, K, double(npicks) / double(K), mean, sdev);

        ASSERT_ALWAYS(mean < 10);
        ASSERT_ALWAYS(sdev < (sqrt(iter) * 10));
        fmt::print("checks passed\n");
    }
}

int main(int argc, char const * argv[])
{
    unsigned long iter = 20000;

    tests_common_cmdline(&argc, &argv, PARSE_SEED | PARSE_ITER);
    tests_common_get_iter(&iter);

    cxx_gmp_randstate rstate;

    test_random_normal_standard(rstate, iter, 2);
    test_random_normal(rstate, 4242.0, 1717.0, iter, 2);
    test_random_normal_constrained(rstate, 4242.0, 1717.0, 0, ULONG_MAX, iter,
                                   2); // F->ncols);
    test_random_poisson(rstate, 4242.0, ULONG_MAX, iter, 2);

    test_arbitrary(rstate, iter);

    tests_common_clear();
    return EXIT_SUCCESS;
}
