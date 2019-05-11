#include "cado.h"
#include "lingen-matpoly.h"
#include "lingen-matpoly-ft.h"
#include "utils.h"

void one_test(cxx_mpz const & p, unsigned int m, unsigned int len1, unsigned int len2, gmp_randstate_t rstate)
{
    abfield ab;
    abfield_init(ab);
    abfield_specify(ab, MPFQ_PRIME_MPZ, p);
    matpoly P, Q, R0, R1, M0, M1;

    matpoly_init(ab, P, 1, m, len1);
    matpoly_init(ab, Q, m, 1, len2);
    matpoly_init(ab, R0, 0, 0, 0);
    matpoly_init(ab, R1, 0, 0, 0);
    matpoly_init(ab, M0, 0, 0, 0);
    matpoly_init(ab, M1, 0, 0, 0);

    matpoly_fill_random(ab, P, len1, rstate);
    matpoly_fill_random(ab, Q, len2, rstate);

    matpoly_mul(ab, R0, P, Q);
    matpoly_mul_caching(ab, R1, P, Q, NULL);

    /* segfault ? */
    matpoly_mp(ab, M0, P, Q);
    matpoly_mp_caching(ab, M1, P, Q, NULL);

    ASSERT_ALWAYS(matpoly_cmp(ab, R0, R1) == 0);
    ASSERT_ALWAYS(matpoly_cmp(ab, M0, M1) == 0);

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
