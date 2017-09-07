#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#ifdef HAVE_MPIR
#include <mpir.h>
#else
#include <gmp.h>
#endif

/* We expect that this program prints exactly:
 *
 * d41c91186caf806b_45558c7335696741_71096848fde90ec7_7b34411325e1217a
 *
 * If it doesn't, then it means that gmp's default random number
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
    printf("%" PRIx64 "_", long_random());
    printf("%" PRIx64 "_", long_random());
    printf("%" PRIx64 "_", long_random());
    printf("%" PRIx64 "\n", long_random());
    gmp_randclear(rstate);
}
