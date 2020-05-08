#include "cado.h" // IWYU pragma: keep
#include <math.h>

#include "las-arith.hpp"

/* utility : is_prime_power */
/* Assume q is a prime power p^k with k>=1, return p if k > 1, 0 otherwise. */
/* This is cheap enough */
fbprime_t is_prime_power(fbprime_t q)
{
    unsigned int maxk, k;
    uint32_t p;

    for (maxk = 0, p = q; p > 1; p /= 2, maxk++);
    for (k = maxk; k >= 2; k--) {
        p = (fbprime_t) (pow((double) q, 1.0 / (double) k) + 0.5);
        if (q % p == 0)
            return p;
    }
    return 0;
}


