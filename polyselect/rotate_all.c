#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "auxiliary.h" /* for common routines with polyselect.c */

extern int MAX_k;

static void
usage_and_die (char *argv0)
{
  fprintf (stderr, "usage: %s [-v] poly kmax\n", argv0);
  fprintf (stderr, "  apply rotation f += (j*x+k)*g to poly.\n");
  fprintf (stderr, "  poly: filename of polynomial\n");
  fprintf (stderr, "  j,k : integers\n");
  exit (1);
}

int
main (int argc, char **argv)
{
    cado_poly poly;
    long kmax, jmin, kmin;
    unsigned long alim = 2000;
    mpz_t b, m;
    int argc0 = argc, verbose = 0;
    char **argv0 = argv;

    while (argc >= 2 && strcmp (argv[1], "-v") == 0)
      {
        argv ++;
        argc --;
        verbose ++;
      }

    mpz_init(b);
    mpz_init(m);
    if (argc != 3)
        usage_and_die (argv0[0]);
    cado_poly_init (poly);
    if (!cado_poly_read(poly, argv[1])) {
        fprintf(stderr, "Problem when reading file %s\n", argv[1]);
        usage_and_die (argv0[0]);
    }
    kmax = strtol(argv[2], NULL, 10);
    MAX_k = kmax;

    poly->skew = L2_skewness (poly->f, poly->degree, SKEWNESS_DEFAULT_PREC,
                           DEFAULT_L2_METHOD);
    if (verbose)
      print_poly (stdout, poly, argc0, argv0, 0, 1);
    else
      printf ("Initial skewness=%1.2f, alpha=%1.2f\n", poly->skew,
              get_alpha (poly->f, poly->degree, ALPHA_BOUND));
    optimize (poly->f, poly->degree, poly->g, verbose, 1);
    poly->skew = L2_skewness (poly->f, poly->degree, SKEWNESS_DEFAULT_PREC,
                           DEFAULT_L2_METHOD);
    if (verbose)
      print_poly (stdout, poly, argc0, argv0, 0, 1);
    else
      printf ("After norm optimization, skewness=%1.2f, alpha=%1.2f\n",
              poly->skew, get_alpha (poly->f, poly->degree, ALPHA_BOUND));

    mpz_set (b, poly->g[1]);
    mpz_neg (m, poly->g[0]);
    rotate (poly->f, poly->degree, alim, m, b, &jmin, &kmin, 0, verbose,
            DEFAULT_L2_METHOD);
    mpz_set (poly->g[1], b);
    mpz_neg (poly->g[0], m);
    /* optimize again, but only translation */
    fprint_polynomial (stdout, poly->g, poly->degreeg);
    optimize (poly->f, poly->degree, poly->g, verbose, 0);
    poly->skew = L2_skewness (poly->f, poly->degree, SKEWNESS_DEFAULT_PREC,
                           DEFAULT_L2_METHOD);
    mpz_clear(b);
    mpz_clear(m);

    print_poly (stdout, poly, argc0, argv0, 0, 1);
    return 0;
} 
