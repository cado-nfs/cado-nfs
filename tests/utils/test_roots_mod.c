#include "cado.h" // IWYU pragma: keep

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gmp.h>

#include "getprime.h"
#include "gmp_aux.h"
#include "mpz_poly.h"
#include "rootfinder.h"

int roots_mod_uint64(uint64_t * r, uint64_t a, int d, uint64_t p,
                     gmp_randstate_ptr rstate);

/* sort the roots r[0], ..., r[n-1] in increasing order */
static void sort_roots(uint64_t * r, unsigned int n)
{
    uint64_t t;

    for (unsigned int i = 1; i < n; i++) {
        t = r[i];
        unsigned int j;
        for (j = i; j > 0 && r[j - 1] > t; j--)
            r[j] = r[j - 1];
        r[j] = t;
    }
}

unsigned int roots_mod_mpz(uint64_t * r, uint64_t a, int d, uint64_t p,
                           gmp_randstate_ptr rstate)
{
    mpz_poly F;
    mpz_poly_init(F, d);
    mpz_poly_set_xi(F, d);
    mpz_set_uint64(mpz_poly_coeff(F, 0), p - a);
    unsigned int n = mpz_poly_roots_uint64(r, F, p, rstate);
    mpz_poly_clear(F);

    sort_roots(r, n);
    return n;
}

int main(int argc, char const * argv[])
{
    uint64_t *r1, *r2;
    unsigned long a, p, d;
    unsigned long minp = 3, maxp = 10000, mina = 0, maxa = 100, mind = 1,
                  maxd = 10;
    int check = 1;

    if (argc > 1 && strcmp(argv[1], "-nc") == 0) {
        printf("Checking results disabled.\n");
        check = 0;
        argc--;
        argv++;
    }
    if (argc > 1)
        minp = strtoul(argv[1], NULL, 10);
    if (argc > 2)
        maxp = strtoul(argv[2], NULL, 10);
    if (argc > 3)
        mina = strtoul(argv[3], NULL, 10);
    if (argc > 4)
        maxa = strtoul(argv[4], NULL, 10);
    if (argc > 5)
        mind = strtoul(argv[5], NULL, 10);
    if (argc > 6)
        maxd = strtoul(argv[6], NULL, 10);

    printf("minp = %lu, maxp = %lu, mina = %lu, maxa = %lu, mind = %lu, maxd = "
           "%lu\n",
           minp, maxp, mina, maxa, mind, maxd);

    r1 = malloc(sizeof(uint64_t) * maxd);
    r2 = malloc(sizeof(uint64_t) * maxd);

    prime_info pi;
    prime_info_init(pi);
    for (p = 2; p < minp; p = getprime_mt(pi))
        ;

    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);

#ifdef __COVERITY__
    __coverity_mark_pointee_as_sanitized__(&maxp, LOOP_BOUND);
    __coverity_mark_pointee_as_sanitized__(&maxa, LOOP_BOUND);
    __coverity_mark_pointee_as_sanitized__(&maxd, LOOP_BOUND);
#endif

    while (p <= maxp) {
        for (a = mina; a <= maxa && a < p; a++) {
            for (d = mind; d <= maxd; d++) {
                unsigned int n1 = roots_mod_uint64(r1, a, (int)d, p, rstate);
                if (check) {
                    unsigned int n2 = roots_mod_mpz(r2, a, (int)d, p, rstate);
                    if (n1 != n2) {
                        fprintf(
                            stderr,
                            "Error: for a=%lu, d=%lu, p=%lu, roots_mod_uint64()"
                            " reports %d roots, roots_mod_mpz() reports %d\n",
                            a, d, p, n1, n2);
                        exit(EXIT_FAILURE);
                    }
                    unsigned int i;
                    for (i = 0; i < n1 && r1[i] == r2[i]; i++)
                        ;
                    if (i != n1) {
                        fprintf(
                            stderr,
                            "Error: for a=%lu, d=%lu, p=%lu, roots_mod_uint64()"
                            " reports roots: ",
                            a, d, p);
                        for (i = 0; i < n1; i++)
                            fprintf(stderr, "%" PRIu64 " ", r1[i]);
                        fprintf(stderr, ", roots_mod_mpz() reports ");
                        for (i = 0; i < n2; i++)
                            fprintf(stderr, "%" PRIu64 " ", r2[i]);
                        fprintf(stderr, "\n");
                        exit(EXIT_FAILURE);
                    }
                }
            }
        }
        p = getprime_mt(pi);
    }
    gmp_randclear(rstate);
    prime_info_clear(pi);
    free(r1);
    free(r2);
    exit(EXIT_SUCCESS);
}
