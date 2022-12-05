#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include "polyselect_norms.h"
#include "polyselect_alpha.h"
#include "auxiliary.h" /* for common routines with polyselect.c */
#include "area.h"
#include "cado_poly.h"
#include "mpz_poly.h"

static void usage_and_die(char *argv0) {
    fprintf(stderr, "usage: %s [-area a] [-I nnn] [-Bf b] [-Bg c] [-skew s] [-B b] poly j k\n", argv0);
    fprintf(stderr, "  apply rotation f += (j*x+k)*g to poly.\n");
    fprintf(stderr, "  poly: filename of polynomial\n");
    fprintf(stderr, "  j,k : integers\n");
    exit(1);
}

int main(int argc, char **argv) {
    cado_poly cpoly;
    long j, k;
    int I = 0;
    mpz_t b, m;
    double skew = 0.0;

    mpz_init(b);
    mpz_init(m);
    while (argc >= 2 && argv[1][0] == '-')
      {
        if (strcmp (argv[1], "-area") == 0)
          {
            area = atof (argv [2]);
            argv += 2;
            argc -= 2;
          }
        else if (strcmp (argv[1], "-Bf") == 0)
          {
            bound_f = atof (argv [2]);
            argv += 2;
            argc -= 2;
          }
        else if (strcmp (argv[1], "-I") == 0)
          {
            I = atoi (argv [2]);
            argv += 2;
            argc -= 2;
          }
        else if (strcmp (argv[1], "-Bg") == 0)
          {
            bound_g = atof (argv [2]);
            argv += 2;
            argc -= 2;
          }
        else if (strcmp (argv[1], "-skew") == 0)
          {
            skew = atof (argv [2]);
            argv += 2;
            argc -= 2;
          }
        else if (strcmp (argv[1], "-B") == 0)
          {
            set_alpha_bound (strtoul (argv [2], NULL, 10));
            argv += 2;
            argc -= 2;
          }
        else
          break;
      }
    if (argc != 4)
        usage_and_die(argv[0]);

    if (I != 0)
      area = bound_f * pow (2.0, (double) (2 * I - 1));

    cado_poly_init (cpoly);
    if (!cado_poly_read(cpoly, argv[1])) {
        fprintf(stderr, "Problem when reading file %s\n", argv[1]);
        usage_and_die(argv[0]);
    }
    j = strtol(argv[2], NULL, 10);
    k = strtol(argv[3], NULL, 10);

    if (skew != 0.0)
      cpoly->skew = skew; /* command-line overrides skewness given in file (if any) */

    /* If skewness is not given in file nor in command line, compute it. */
    if (cpoly->skew == 0.0)
      cpoly->skew = L2_combined_skewness2 (cpoly->pols[0], cpoly->pols[1],
                                          SKEWNESS_DEFAULT_PREC);

    mpz_set (b, cpoly->pols[RAT_SIDE]->coeff[1]);
    mpz_set (m, cpoly->pols[RAT_SIDE]->coeff[0]);
    rotate_aux (cpoly->pols[ALG_SIDE]->coeff, b, m, 0, k, 0);
    rotate_aux (cpoly->pols[ALG_SIDE]->coeff, b, m, 0, j, 1);
    mpz_clear(b);
    mpz_clear(m);

    print_cadopoly_extra (stdout, cpoly, argc, argv, 0);
    return 0;
} 
