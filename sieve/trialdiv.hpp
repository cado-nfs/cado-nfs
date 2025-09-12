#ifndef CADO_TRIALDIV_HPP
#define CADO_TRIALDIV_HPP

#include <cstddef>   // for size_t
#include <cstdint>   // for uint64_t, SIZE_MAX
#include <climits>
#include <vector>    // for vector
#include "cado_compile_time_hacks.hpp"

struct cxx_mpz;

/* The maximum number of words in the numbers to be trial-divided.
   The $l$ value from the thesis text is equal to TRIALDIV_MAXLEN - 1 */
#ifndef TRIALDIV_MAXLEN
#define TRIALDIV_MAXLEN 8
#endif

typedef struct {
  unsigned long p;
  unsigned long w[TRIALDIV_MAXLEN - 1]; /* w[i] = w^{i+1} mod p */
  unsigned long pinv; /* pinv == 1/p (mod w) */
  unsigned long plim; /* plim = (w-1)/p */
} trialdiv_divisor_t;

#ifdef __cplusplus
extern "C" {
#endif

struct trialdiv_data : public std::vector<trialdiv_divisor_t>
{
    trialdiv_data() = default;
    explicit trialdiv_data(std::vector<unsigned long> const & primes, size_t skip = 0);

    /* A shortcoming that we face is that std::sqrt only becomes
     * constexpr with c++26. We have a dichotomy implementation in
     * cado_math_aux which should be good enough.
     *
     * (I don't quite understand what max_p is, to be honest. It's seldom
     * used. It was introduced in b04031eb1 and 0af05e731).
     */
    static constexpr unsigned long max_p = 
            (TRIALDIV_MAXLEN == 1) ?
                ULONG_MAX :
                (cado_math_aux::constant_sqrt(ULONG_MAX / (TRIALDIV_MAXLEN - 1)) - 1);

    /* TODO ; input primes are unsigned long, output primes are uint64_t,
     * it's a bit ridiculous.
     *
     * (rationale for output primes: they go in the factor_list, which
     * contains all sorts of primes, some of which may even exceed 32 bits).
     *
     * max_factors: give up after finding that many factors.
     *
     * The trial_divide(L, N, max) methods *appends* to L. (and stops at
     * max *new* factors if needed).
     */
    size_t trial_divide(std::vector<uint64_t> &, cxx_mpz & N, size_t max_factors = SIZE_MAX) const;
    std::vector<uint64_t> trial_divide(cxx_mpz & N, size_t max_factors = SIZE_MAX) const {
        std::vector<uint64_t> res;
        trial_divide(res, N, max_factors);
        return res;     /* copy elision */
    }
};

#ifdef __cplusplus
}
#endif


#endif	/* CADO_TRIALDIV_HPP */
