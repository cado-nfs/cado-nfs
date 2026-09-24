#ifndef CADO_TESTS_BENCH_ARITH_VS_ARITHXX_HPP
#define CADO_TESTS_BENCH_ARITH_VS_ARITHXX_HPP

/* Kernel benchmark of the old arith layer against arithxx.
 *
 * Each kernel is a dependent chain of operations, run the same way on
 * both sides from the same starting values, so that the final values
 * (the "check" values) must agree. See bench_arith_vs_arithxx.cpp.
 */

#include <cstddef>
#include <cstdint>

#include <chrono>
#include <string>
#include <vector>

#include "cxx_mpz.hpp"

struct kbench_iters {
    size_t mul, pow, inv, prime;
};

struct kbench_result {
    std::string op;
    double ns;
    uint64_t chk;
};

typedef std::vector<kbench_result> kbench_results;

/* best of 5 runs of f(), which does n iterations; in ns per iteration */
template<typename F>
double kbench_time(size_t n, F && f)
{
    double best = 1e300;
    for (int rep = 0; rep < 5; rep++) {
        auto const t0 = std::chrono::steady_clock::now();
        f();
        auto const t1 = std::chrono::steady_clock::now();
        double const ns =
            std::chrono::duration<double, std::nano>(t1 - t0).count() /
            double(n);
        if (ns < best)
            best = ns;
    }
    return best;
}

/* the exponents used by the pow kernels, as unsigned long so that both
 * sides use the same ones on 32-bit machines too */
static constexpr unsigned long kbench_pow_exponent = (unsigned long) 0xd1b54a32d192ed03ULL;
static constexpr unsigned long kbench_pow2_exponent = (unsigned long) 0x9e3779b97f4a7c15ULL;

/* the sequence of exponents of the 2pow kernel. It does not depend on
 * the results, so that the loop times 2pow alone (getting a residue out
 * of mod_mpz_new allocates, for instance) */
static inline unsigned long kbench_next_exponent(unsigned long e)
{
    return e * (unsigned long) 6364136223846793005ULL
        + (unsigned long) 1442695040888963407ULL;
}

/* low 64 bits of a nonnegative integer */
uint64_t kbench_low64(mpz_srcptr z);

kbench_results bench_old_ul(cxx_mpz const & M, kbench_iters const & it);
kbench_results bench_old_15ul(cxx_mpz const & M, kbench_iters const & it);
kbench_results bench_old_2ul2(cxx_mpz const & M, kbench_iters const & it);
kbench_results bench_old_mpz(cxx_mpz const & M, kbench_iters const & it);

#endif /* CADO_TESTS_BENCH_ARITH_VS_ARITHXX_HPP */
