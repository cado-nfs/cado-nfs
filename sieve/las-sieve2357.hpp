#ifndef LAS_SIEVE2357_HPP
#define LAS_SIEVE2357_HPP

#include "cado_config.h"  // for HAVE_AVX2, HAVE_SSSE3

#include <cstddef>       // for size_t
#include <cstdint>       // for uint8_t
#ifdef HAVE_SSE2
/* The base type that we'll use is good ol sse2 __m128i, but we need
 * SSSE3 instructions to deal with them in the .cpp file */
#include <emmintrin.h>  // IWYU pragma: keep
#endif
#ifdef HAVE_AVX2
#include <immintrin.h>
#endif
#ifdef HAVE_ARM_NEON
#include <arm_neon.h>
#endif

#include "fb-types.h"     // for fbprime_t
struct where_am_I;


class sieve2357base {
public:
#if defined(HAVE_AVX2)
    typedef __m256i preferred_simd_type;
#elif defined(HAVE_SSSE3)
    typedef __m128i preferred_simd_type;
#elif defined(HAVE_ARM_NEON)
    typedef uint8x16_t preferred_simd_type;
#elif LONG_BIT == 64 /* FIXME: this is false on, e.g., MinGW. What's a better condition here? */
    typedef uint64_t preferred_simd_type;
#else
    typedef uint32_t preferred_simd_type;
#endif

    /* Define ordering on q which sieve2357 expect: first q=1, then powers of 2,
       then odd primes in increasing order */
    static inline bool order_lt(fbprime_t q1, fbprime_t q2) {
        if (q1 == 1 && q2 != 1)
          return true;
        if (q1 != 1 && q2 == 1)
          return false;
        if (q1 % 2 == 0 && q2 % 2 == 1)
          return true;
        if (q1 % 2 == 1 && q2 % 2 == 0)
          return false;
        return q1 < q2;
    }

    /* We hit sievearray[x] for x = idx + i*q with 0 <= x < arraylen */
    struct prime_t {
        fbprime_t q, idx;
        unsigned char logp;
    };

    /* The sieve() function can either write to the sieve array, overwriting any
       data that was previously there, or add to the sieve array.
       If sieve2357::sieve() is to write to a sieve array that contains only
       zeroes, then the former is faster (don't have to zero out the array and
       sieve() saves one read and add per word).
       If the sieve array contains non-zero data, then obviously update_add is
       required. */
    enum {
        update_set,
        update_add
    };
};

template <typename SIMDTYPE, typename ELEMTYPE>
class sieve2357: public sieve2357base {
private:
    static inline void
    sieve_odd_prime(SIMDTYPE * const result, const ELEMTYPE logp,
        const fbprime_t stride, const fbprime_t idx, const SIMDTYPE even_mask);
    static void
    big_loop_set(SIMDTYPE * __restrict__ sievearray, const SIMDTYPE * const __restrict__ pattern235,
        const SIMDTYPE * const __restrict__ pattern7, size_t l7);
    static void
    big_loop_add(SIMDTYPE * __restrict__ sievearray, const SIMDTYPE * const __restrict__ pattern235,
        const SIMDTYPE * const __restrict__ pattern7, size_t l7);
    static SIMDTYPE
    sieve2(const fbprime_t q, const fbprime_t idx, const uint8_t logp);
    static SIMDTYPE
    get_even_mask(const int skip_mod_2);

public:
    /* A predicate that tells whether a prime power q = p^k can be sieved by
       sieve2357 with a given SIMD and element data type */
    static inline bool can_sieve(const fbprime_t q)
    {
        const size_t N = sizeof(SIMDTYPE) / sizeof(ELEMTYPE);
        return q == 1 || (q % 2 == 0 && N % q  == 0) || q == 3 || q == 5
            || q == 7;
    }

    static void sieve(SIMDTYPE *sievearray, size_t arraylen, const prime_t *primes,
        int skip_mod_2, int update_operation, where_am_I & w);

};
#endif
