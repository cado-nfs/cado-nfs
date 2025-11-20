/* Schirokauer maps

   Given a list of a,b pairs, compute the corresponding SMs.

   Each line of input can also be a polynomial, with the following format:
     p deg c0 c1 ... c_deg
   where 'p' is just the letter p and the rest are integers.

   */

#include "cado.h" // IWYU pragma: keep
#include <cstdio>
#include <cstdlib>

#include <gmp.h>

#include "cado_poly.h"
#include "macros.h"
#include "mpz_poly.h"
#include "params.h"
#include "sm_utils.hpp"
#include "timing.h"
#include "verbose.h"

static void my_sm(char const * outfile, char const * infile,
                  std::vector<sm_side_info> const & sm_info, int nb_polys)
{
    FILE * in;
    in = fopen(infile, "r");
    if (in == NULL) {
        fprintf(stderr, "Error: could not open %s for reading\n", infile);
        exit(EXIT_FAILURE);
    }
    FILE * out = stdout;
    if (outfile != NULL) {
        out = fopen(outfile, "w");
        if (out == NULL) {
            fprintf(stderr, "Error: could not open %s for writing\n", outfile);
            exit(EXIT_FAILURE);
        }
    }

    char buf[1024];
    cxx_mpz_poly pol, smpol;
    int maxdeg = 0;
    for (int side = 0; side < nb_polys; side++)
        maxdeg = MAX(maxdeg, sm_info[side].f->deg);
    while (fgets(buf, 1024, in)) {
        if (buf[0] == '#')
            continue;
        if (buf[0] == 'p') {
            // we read a polynomial
            // buf = "p deg c0 c1 ... c_deg"
            int deg;
            char * tmp = buf + 2;
            sscanf(tmp, "%d", &deg);
#ifdef __COVERITY__
            __coverity_mark_pointee_as_sanitized__(&deg, LOOP_BOUND);
#endif
            mpz_t zbuf;
            mpz_init(zbuf);
            for (int i = 0; i <= deg; i++) {
                for (++tmp; *tmp != ' '; tmp++)
                    ;
                gmp_sscanf(tmp, "%Zd", zbuf);
                mpz_poly_setcoeff(pol, i, zbuf);
            }
            mpz_clear(zbuf);
            if (0) {
                fprintf(stderr, "Poly read: ");
                mpz_poly_fprintf(stderr, pol);
                fprintf(stderr, "\n");
            }
        } else {
            // we read a relation
            cxx_mpz a, b;
            int const ret = gmp_sscanf(buf, "%Zd,%Zd:", (mpz_ptr)a, (mpz_ptr)b);
            ASSERT_ALWAYS(ret == 2);
            mpz_poly_set_mpz_ab(pol, a, b);
        }
        for (int side = 0; side < nb_polys; ++side) {
            sm_info[side].compute_piecewise(smpol, pol);
            print_sm(out, sm_info[side], smpol);
            if (side == 0 && sm_info[0].nsm > 0 && sm_info[1].nsm > 0)
                fprintf(out, " ");
        }
        fprintf(out, "\n");
    }

    if (out != NULL)
        fclose(out);
    fclose(in);
}

static void declare_usage(param_list pl)
{
    param_list_decl_usage(pl, "poly", "(required) poly file");
    param_list_decl_usage(pl, "inp",
                          "(required) input file containing relations");
    param_list_decl_usage(pl, "out", "output file");
    param_list_decl_usage(pl, "ell", "(required) group order");
    param_list_decl_usage(pl, "sm-mode", "SM mode (see sm-portability.h)");
    verbose_decl_usage(pl);
}

static void usage(char const * argv, char const * missing, param_list pl)
{
    if (missing) {
        fprintf(stderr, "\nError: missing or invalid parameter \"-%s\"\n",
                missing);
    }
    param_list_print_usage(pl, argv, stderr);
    exit(EXIT_FAILURE);
}

/* -------------------------------------------------------------------------- */

// coverity[root_function]
int main(int argc, char const * argv[])
{
    char const * argv0 = argv[0];

    char const * polyfile = NULL;
    char const * infile = NULL;
    char const * outfile = NULL;

    param_list pl;
    cado_poly cpoly;

    mpz_t ell, ell2;
    double t0;

    /* read params */
    param_list_init(pl);
    declare_usage(pl);

    if (argc == 1)
        usage(argv[0], NULL, pl);

    argc--, argv++;
    for (; argc;) {
        if (param_list_update_cmdline(pl, &argc, &argv)) {
            continue;
        }
        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        usage(argv0, NULL, pl);
    }

    /* Read poly filename from command line */
    if ((polyfile = param_list_lookup_string(pl, "poly")) == NULL) {
        fprintf(stderr, "Error: parameter -poly is mandatory\n");
        param_list_print_usage(pl, argv0, stderr);
        exit(EXIT_FAILURE);
    }

    /* Read purged filename from command line */
    if ((infile = param_list_lookup_string(pl, "inp")) == NULL) {
        fprintf(stderr, "Error: parameter -inp is mandatory\n");
        param_list_print_usage(pl, argv0, stderr);
        exit(EXIT_FAILURE);
    }

    /* Read outfile filename from command line ; defaults to stdout. */
    outfile = param_list_lookup_string(pl, "out");

    /* Read ell from command line (assuming radix 10) */
    mpz_init(ell);
    if (!param_list_parse_mpz(pl, "ell", ell)) {
        fprintf(stderr, "Error: parameter -ell is mandatory\n");
        param_list_print_usage(pl, argv0, stderr);
        exit(EXIT_FAILURE);
    }

    /* Init polynomial */
    cado_poly_init(cpoly);
    cado_poly_read(cpoly, polyfile);

    std::vector<mpz_poly_srcptr> F(cpoly->nb_polys, NULL);

    for (int side = 0; side < cpoly->nb_polys; side++)
        F[side] = cpoly->pols[side];

    char const * sm_mode_string = param_list_lookup_string(pl, "sm-mode");

    if (param_list_warn_unused(pl))
        usage(argv0, NULL, pl);
    verbose_interpret_parameters(pl);
    param_list_print_command_line(stdout, pl);

    mpz_init(ell2);
    mpz_mul(ell2, ell, ell);

    std::vector<sm_side_info> sm_info;

    for (int side = 0; side < cpoly->nb_polys; side++) {
        sm_info.emplace_back(F[side], ell, 0);
        sm_info[side].set_mode(sm_mode_string);
    }

    for (int side = 0; side < cpoly->nb_polys; side++) {
        fprintf(stdout, "\n# Polynomial on side %d:\nF[%d] = ", side, side);
        mpz_poly_fprintf(stdout, F[side]);

        printf("# SM info on side %d:\n", side);
        sm_info[side].print(stdout);

        fflush(stdout);
    }

    t0 = seconds();

    my_sm(outfile, infile, sm_info, cpoly->nb_polys);

    fprintf(stdout, "\n# sm completed in %2.2lf seconds\n", seconds() - t0);
    fflush(stdout);

    mpz_clear(ell);
    mpz_clear(ell2);
    cado_poly_clear(cpoly);
    param_list_clear(pl);

    return 0;
}
