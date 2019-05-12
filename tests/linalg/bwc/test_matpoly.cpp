#include "cado.h"
#include "lingen-matpoly.hpp"
#include "lingen-matpoly-ft.hpp"
#include "utils.h"

void one_test(cxx_mpz const & p, unsigned int m, unsigned int len1, unsigned int len2, gmp_randstate_t rstate)
{
    abfield ab;
    abfield_init(ab);
    abfield_specify(ab, MPFQ_PRIME_MPZ, p);

    matpoly P(ab, 1, m, len1);
    matpoly Q(ab, m, 1, len2);
    matpoly R0;
    matpoly R1;
    matpoly M0;
    matpoly M1;

    P.fill_random(len1, rstate);
    Q.fill_random(len2, rstate);

    R0.mul(P, Q);
    matpoly_mul_caching(R1, P, Q, NULL);

    /* segfault ? */
    M0.mp(P, Q);
    matpoly_mp_caching(M1, P, Q, NULL);

    ASSERT_ALWAYS(R0.cmp(R1) == 0);
    ASSERT_ALWAYS(M0.cmp(M1) == 0);

    abfield_clear(ab);
}

int main(int argc, char * argv[])
{
    cxx_mpz p;
    unsigned long seed = 1;
    gmp_randstate_t rstate;

    for(argv++,argc--;argc;argv++,argc--) {
        if (argc >= 2 && strcmp(*argv, "--prime") == 0) {
            argv++,argc--;
            mpz_set_str(p, *argv, 0);
            continue;
        }
        if (argc >= 2 && strcmp(*argv, "--seed") == 0) {
            argv++,argc--;
            seed = atoi(*argv);
            continue;
        }
        fprintf(stderr, "Unexpected argument: %s\n", *argv);
        exit(EXIT_FAILURE);
    }

    if (mpz_cmp_ui(p, 0) == 0) {
        fprintf(stderr, "--prime is mandatory\n");
        exit(EXIT_FAILURE);
    }

    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, seed);


    one_test(p, 4, 1000, 600, rstate);


    gmp_randclear(rstate);
}
