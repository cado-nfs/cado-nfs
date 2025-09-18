#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <gmp.h>
#include "polyselect_norms.h"
#include "polyselect_alpha.h"
#include "auxiliary.h" /* for common routines with polyselect.c */
#include "area.h"
#include "cado_poly.h"
#include "mpz_poly.h"
#include "params.h"

/* This simplistic binary just computes the rotation f+k*g, starting from a
 * polynomial pair (f,g), and a rotation polynomial k (given least
 * coefficient first).
 */

static void usage(const char *argv, const char *missing, param_list_ptr pl)
{
    fprintf(stderr, "usage: %s [parameters] <polynomial file> [<rotation coefficients> ...]\n", argv);
    if (missing)
    {
        fprintf(stderr, "\nError: missing or invalid parameter \"-%s\"\n",
                missing);
    }
    param_list_print_usage(pl, argv, stderr);
    param_list_clear(pl);
    exit(EXIT_FAILURE);
}

void declare_usage(param_list_ptr pl)
{
    param_list_decl_usage(pl, "area", "sieving area (bound on a,b?) used for the computation of MurphyE");
    param_list_decl_usage(pl, "I", "width of the I,J rectangle used for the computation of MurphyE");
    param_list_decl_usage(pl, "Bf", "smoothness bound on side 1, used for the computation of MurphyE");
    param_list_decl_usage(pl, "Bg", "smoothness bound on side 0, used for the computation of MurphyE");
    param_list_decl_usage(pl, "skew", "skewness (default to automatic)\n");
    param_list_decl_usage(pl, "B", "alpha bound used for the computation of alpha");
}


int main(int argc, char const * argv[])
{
    char const ** argv0 = argv;
    int argc0 = argc;
    cado_poly cpoly;
    int I = 0;
    double skew = 0.0;
    param_list pl;
    mpz_poly rot;
    mpz_t tmp;

    mpz_init(tmp);
    mpz_poly_init(rot, -1);

    const char * polyfilename = NULL;

    param_list_init(pl);

    declare_usage(pl);

    argv++, argc--;
    for (int wild = 0; argc; )
    {
        if (param_list_update_cmdline(pl, &argc, &argv))
            continue;
        if (wild == 0) {
            polyfilename = *argv;
            argc--,argv++,wild++;
        } else if (isdigit(argv[0][0]) || (argv[0][0] == '-' && isdigit(argv[0][1]))) {
            int ok = mpz_set_str(tmp, *argv, 0) == 0;
            if (!ok) {
                fprintf(stderr, "Error, %s is not an integer\n", argv[0]);
                usage(argv0[0], NULL, pl);
            }
            mpz_poly_setcoeff(rot, wild-1, tmp);
            argc--,argv++,wild++;
        } else {
            fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
            usage(argv0[0], NULL, pl);
        }
    }

    param_list_parse_double(pl, "area", &area);
    param_list_parse_double(pl, "Bf", &bound_f);
    param_list_parse_double(pl, "Bg", &bound_g);
    param_list_parse_double(pl, "skew", &skew);
    param_list_parse_int(pl, "I", &I);
    {
        unsigned long B;
        if (param_list_parse_ulong(pl, "B", &B))
            set_alpha_bound (B);
    }

    if (param_list_warn_unused(pl))
        usage (argv0[0], NULL, pl);

    if (I != 0)
      area = bound_f * pow (2.0, (double) (2 * I - 1));

    cado_poly_init (cpoly);
    if (!cado_poly_read(cpoly, polyfilename)) {
        fprintf(stderr, "Problem when reading file %s\n", argv[1]);
        usage(argv[0], NULL, pl);
    }

    if (skew != 0.0)
      cpoly->skew = skew; /* command-line overrides skewness given in file (if any) */

    /* If skewness is not given in file nor in command line, compute it. */
    if (cpoly->skew == 0.0)
      cpoly->skew = L2_combined_skewness2 (cpoly->pols[0], cpoly->pols[1]);

    /* Well, it's really very simple. */

    mpz_poly_mul(rot, rot, cpoly->pols[RAT_SIDE]);
    mpz_poly_add(cpoly->pols[ALG_SIDE], cpoly->pols[ALG_SIDE], rot);

    print_cadopoly_extra (stdout, cpoly, argc0, argv0, 0);

    cado_poly_clear(cpoly);

    param_list_clear(pl);
    mpz_clear(tmp);
    mpz_poly_clear(rot);

    return 0;
} 
