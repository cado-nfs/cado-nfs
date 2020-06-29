#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#ifdef HAVE_MPIR
#include <mpir.h>
#else
#include <gmp.h>
#endif

/* 
 * If this program fails, it means that gmp's default random number
 * generator has changed, and this is likely to ruin our tests that are
 * dependent on its behaviour.
 */
gmp_randstate_t rstate;

static inline uint64_t long_random() {
#if ULONG_BITS == 64
    return gmp_urandomb_ui(rstate, 64);
#else
    mpz_t z;
    mpz_init(z);
    mpz_urandomb(z, rstate, 64);
    uint64_t res = mpz_get_ui (z) + (((uint64_t) mpz_getlimbn(z,1)) << 32);
    mpz_clear(z);
    return res;
#endif
}

int main()
{
    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, 0);
    if (long_random() != UINT64_C(0xd41c91186caf806b)) exit(EXIT_FAILURE);
    if (long_random() != UINT64_C(0x45558c7335696741)) exit(EXIT_FAILURE);
    if (long_random() != UINT64_C(0x71096848fde90ec7)) exit(EXIT_FAILURE);
    if (long_random() != UINT64_C(0x7b34411325e1217a)) exit(EXIT_FAILURE);
    gmp_randclear(rstate);
    exit(EXIT_SUCCESS);
}
