#include "cado.h" // IWYU pragma: keep

#include <cmath>
#include <cstdlib>

#include <algorithm>
#include <vector>
#include <utility>

#include "fmt/base.h"

#include "misc.h"

#include "decomp.hpp"
#include "gen_decomp.hpp"
#include "tab_decomp.hpp"
#include "runtime_numeric_cast.hpp"
#include "gmp_aux.h"
#include "integer_partitions.hpp"
#include "macros.h"
#include "polynomial.hpp"
#include "utils_cxx.hpp"

/* We want to estimate number of products of primes:
 * l[0]-bit * l[1]-bit * ... * l[k-1]-bit
 * with product of exactly mfb bits, and each prime >= lim
 *
 * Two options exist.
 *
 *  - The first one is to generate primes (or numbers that follow the
 *    distribution of primes) in each of the bit ranges,
 *    determine how often this leads to a product with the exact
 *    requested bit size. This gives us a ratio. Apply the ratio to the
 *    complete population of choices of primes. This approach is
 *    probabilistic and requires a certain number of random picks in
 *    order to give something useful.
 *
 *  - The second approach is to work with the generating series of the
 *    product. This works best if replace bit sizes by ranges of number
 *    of digits in base gamma, for some gamma \in (1,2]. Smaller values
 *    of gamma (or equivalently larger values of log(beta)/log(gamma))
 *    lead to more accurate but more expensive calculations.
 *
 * In our experiments, it is reasonably easy to make these estimations
 * agree within 5%, with PREC = 100000 in the firt test, and gamma=1.005
 * in the second test.
 */

/* We need two utility functions. The first one computes the prime counts
 * of geometrically increasing intervals, up to some length. The returned
 * vector B has size n, and is such that B[i] is the number of primes in
 * the interval [gamma^i, gamma^(i+1)).
 */
static std::vector<double> prime_counts_in_geometric_intervals(double gamma, size_t n, double lim=0)
{
    std::vector<double> B;
    double p0 = 1;
    /* A simple way to gain a factor of two: nprimes_interval uses the
     * log-integral computation, which is expensive. And it's a
     * primitive, so we compute Li(p1)-Li(p0) for all intervals.
     *
     * Our code arbitrarily defines Li(x)=0 for all x<=1. So it suffices
     * to compute Li(p)-Li(1) for all the boundary values, and subtract
     * what we get.
     */
    double Z0 = nprimes_interval(1, std::max(lim, p0));
    for (unsigned int i = 0; i < n; i++) {
        const double p1 = p0 * gamma;
        const double Z1 = nprimes_interval(1, std::max(lim, p1));
        B.push_back(Z1 - Z0);
        Z0 = Z1;
        p0 = p1;
    }
    return B;
}

/* This takes into account that we might have repeated sizes in the
 * decomposition, for which the order obviously does not matter. So we
 * compute the product of the k!'s for each size that is repeated k times */
static double multiset_correction(std::vector<unsigned int> const & q)
{
    double res = 1;
    for(unsigned int i = 0 ; i < q.size() ; i++) {
        unsigned int j;
        for(j = 1 ; j <= i && q[i] == q[i-j] ; j++);
        res *= j;
    }
    return res;
}


/* This is the first method. */
static double psi_probabilistic(std::vector<double> const & by_interval,
        double mfb0,
        double mfb1,
        std::vector<unsigned int> const & l,
        unsigned long lim, size_t number_of_trials)
{
    cxx_gmp_randstate rstate;

    double S = 1;
    for (auto c : l) {
        ASSERT_ALWAYS(c <= by_interval.size());
        S *= by_interval[c-1];
    }
    S /= multiset_correction(l);
    if (S == 0)
        return 0;

    const double rmin = ldexp(1.0, int(mfb0));
    const double rmax = ldexp(1.0, int(mfb1));

    size_t ok = 0;
    size_t n_above_lim = 0;

    for(size_t tot = 0 ; tot < number_of_trials ; tot++) {
        /* see how often a product with these sizes falls in the correct
         * bit range
         */
        double r = 1;
        for (auto const c : l) {
            const double p = random_along_prime_distribution(c, rstate);
            if (p < (double)lim) {
                r = 0;
                break;
            }
            r *= p;
        }
        if (r == 0)
            continue;
        n_above_lim++;
        ok += (rmin <= r && r <= rmax);
    }
    if (!n_above_lim)
        return 0;
    return (double) ok * S / (double) n_above_lim;
}

/* This is the second method. We give a rough estimate of the number of
 * numbers between beta^n0 and beta^n1 that are products of q.size()
 * prime numbers between
 * beta^(q[i]-1) and beta^(q[i]), for 0<=i<q.size()
 *
 * computations are done based on basis gamma. Bringing gamma closer to 1
 * makes the computation more precise, but also more expensive.
 *
 * beta and gamma are not given explicitly, but via
 * scale=log(beta)/log(gamma).
 *
 * Note that the by_interval array must have been computed with the same
 * gamma as the one we want to work with here.
 */
static double psi_series(std::vector<double> const & by_interval,
        double n0, double n1,
        std::vector<unsigned int> const & q,
        double,
        double scale)
{
    using R = runtime_numeric_cast<unsigned int>;
    /* the scale below is the main driving factor of the computation cost.  */

    polynomial<double> S("1");
    unsigned int offset = 0;
    for(double c : q) {
        /* we're interested in the prime factors such that
         * c-1 <= log(p)/log(beta) < c
         *
         * (and also loglim <= log(p)/log(beta))
         *
         * find d0 and d1 such that this range is a subset of the range
         * of primes such that
         * d0-1 <= log(p) / log(gamma) < d1
         *
         * so we need d0-1 <= (c-1) * scale < d0
         * and d1-1 < c*scale <= d1
         */
        const unsigned int d0 = R(lround((c-1) * scale) + 1);
        const unsigned int d1 = R(lround(ceil(c * scale)));
        polynomial<double> C;
        offset += d0;
        /* The first count we need is the number of primes such that
         * d0-1 <= log(p) / log(gamma) < d0
         * gamma^(d0-1) <= p < gamma^d0
         */
        for(unsigned int d = std::max(d0, 1U) ; d < d1 ; d++)
            C[d - d0] += by_interval[d-1];
        S *= C;
        // fmt::print("{} -> {}-{}\n", c, d0, d1);
        // fmt::print("{} -> {}-{} {}\n", c, d0, d1, C);
    }
    const unsigned int m0 = R(lround(n0 * scale)+1);
    const unsigned int m1 = R(lround(ceil(n1 * scale)));
    double res = 0;
    for(unsigned int m = m0 ; m < m1 ; m++) {
        if (m >= offset)
            res += S[m-offset];
    }
    return res / multiset_correction(q);
}

struct psi_backend_base {
    using R = runtime_numeric_cast<unsigned int>;
    unsigned int mfb;
    unsigned long lim;
    double loglim;
    unsigned int max_factors;
    psi_backend_base(unsigned int mfb, unsigned long lim)
        : mfb(mfb)
        , lim(lim)
        , loglim(log2(double(lim)))
        , max_factors(R(lround(mfb / loglim)))
    {}
};

struct psi_backend_probabilistic {
    psi_backend_base const & psi0;
    static constexpr const size_t ntrials = 100000;
    using R = runtime_numeric_cast<unsigned int>;
    std::vector<double> B;
    explicit psi_backend_probabilistic(psi_backend_base const & psi0)
        : psi0(psi0)
        , B(prime_counts_in_geometric_intervals(2,
                    psi0.mfb + psi0.max_factors + 1,
                    double(psi0.lim)))
    {
    }
    double operator()(std::vector<unsigned int> const & q) const
    {
        return psi_probabilistic(B, psi0.mfb-1, psi0.mfb, q, psi0.lim, ntrials);
    }
};


struct psi_backend_series {
    psi_backend_base const & psi0;
    static constexpr const double gamma = 1.005;
    using R = runtime_numeric_cast<unsigned int>;
    double scale;
    unsigned int max_gamma;
    std::vector<double> B;
    explicit psi_backend_series(psi_backend_base const & psi0)
        : psi0(psi0)
        , scale(log(2) / log(gamma))
        , max_gamma(R(lround(ceil(scale * (psi0.mfb + psi0.max_factors)))) + 1)
        , B(prime_counts_in_geometric_intervals(gamma,
                    max_gamma,
                    double(psi0.lim)))
    {
    }
    double operator()(std::vector<unsigned int> const & q) const
    {
        return psi_series(B, psi0.mfb-1, psi0.mfb, q, psi0.loglim, scale);
    }
};


template<typename T>
static tabular_decomp generate_all_decomp(unsigned int mfb, unsigned long lim)
{
    tabular_decomp res;

    const psi_backend_base psi0(mfb, lim);
    const T psi(psi0);

    for(unsigned int nfactors = 2 ; nfactors <= psi0.max_factors ; nfactors++) {
        /* generate partitions of mfb + k with exactly nfactors summands,
         * for offset values k such that 0 <= k < nfactors.
         * This is because if we increase by one bit the range of all
         * factors but one, then there's still room for the product to
         * have mfb bits.
         */
        const auto imin = (unsigned int) (psi0.loglim+1);
        for(unsigned int k = 0 ; k < nfactors ; k++) {
            auto P = integer_partitions_in_k_parts(mfb + k, nfactors, imin);
            for(auto const & q : P) {
                decomp D { psi(q), q };
                if (D.nb_elem < 1)
                    continue;
                res.emplace_back(std::move(D));
            }
        }
    }
    return res;
}

void generate_all_decomp_compare(unsigned int mfb, unsigned int lim)
{
    const psi_backend_base psi0(mfb, lim);
    const psi_backend_probabilistic psi1(psi0);
    const psi_backend_series psi2(psi0);

    for(unsigned int nfactors = 2 ; nfactors <= psi0.max_factors ; nfactors++) {
        /* generate partitions of mfb + k with exactly nfactors summands,
         * for offset values k such that 0 <= k < nfactors.
         * This is because if we increase by one bit the range of all
         * factors but one, then there's still room for the product to
         * have mfb bits.
         */
        const auto imin = (unsigned int) (psi0.loglim+1);
        for(unsigned int k = 0 ; k < nfactors ; k++) {
            auto P = integer_partitions_in_k_parts(mfb + k, nfactors, imin);
            for(auto const & q : P) {
                auto E1 = psi1(q);
                auto E2 = psi2(q);
                if (E1 < 1 && E2 < 1)
                    continue;
                fmt::print("{} {:6.3e} {:6.3e} {:.3f}\n", join(q, " "), E1, E2, E2 / E1);

            }
        }
    }
}

tabular_decomp generate_all_decomp(unsigned int mfb, unsigned long lim)
{
    return generate_all_decomp<psi_backend_series>(mfb, lim);
}
