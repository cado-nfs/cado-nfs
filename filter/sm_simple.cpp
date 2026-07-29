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

#include "cado_poly.hpp"
#include "macros.h"
#include "mpz_poly.h"
#include "params.hpp"
#include "sm_utils.hpp"
#include "timing.h"
#include "verbose.hpp"

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

static void declare_usage(cxx_param_list& pl)
{
    pl.declare_usage("poly", "(required) poly file");
    pl.declare_usage("inp",
                          "(required) input file containing relations");
    pl.declare_usage("out", "output file");
    pl.declare_usage("ell", "(required) group order");
    pl.declare_usage("sm-mode", "SM mode (see sm-portability.h)");
    verbose_decl_usage(pl);
}

/* -------------------------------------------------------------------------- */

// coverity[root_function]
int main(int argc, char const * argv[])
{
    char const * polyfile = NULL;
    char const * infile = NULL;
    char const * outfile = NULL;

    cxx_param_list pl;
    cxx_cado_poly cpoly;

    mpz_t ell, ell2;
    double t0;

    /* read params */
    declare_usage(pl);

    pl.process_command_line(argc, argv, false);

    /* Read poly filename from command line */
    if ((polyfile = pl.lookup_old("poly")) == NULL)
        pl.fail("Error: parameter -poly is mandatory\n");

    /* Read purged filename from command line */
    if ((infile = pl.lookup_old("inp")) == NULL)
        pl.fail("Error: parameter -inp is mandatory\n");

    /* Read outfile filename from command line ; defaults to stdout. */
    outfile = pl.lookup_old("out");

    /* Read ell from command line (assuming radix 10) */
    mpz_init(ell);
    pl.parse_mandatory("ell", ell);

    /* Init polynomial */
    cpoly.read(polyfile);

    std::vector<mpz_poly_srcptr> F(cpoly.nsides(), NULL);

    for (int side = 0; side < cpoly.nsides(); side++)
        F[side] = cpoly[side];

    char const * sm_mode_string = pl.lookup_old("sm-mode");

    if (pl.warn_unused())
        pl.fail("Unused parameters are given");

    verbose_interpret_parameters(pl);
    pl.print_command_line(stdout);

    mpz_init(ell2);
    mpz_mul(ell2, ell, ell);

    std::vector<sm_side_info> sm_info;

    for (int side = 0; side < cpoly.nsides(); side++) {
        sm_info.emplace_back(F[side], ell, 0);
        sm_info[side].set_mode(sm_mode_string);
    }

    for (int side = 0; side < cpoly.nsides(); side++) {
        fprintf(stdout, "\n# Polynomial on side %d:\nF[%d] = ", side, side);
        mpz_poly_fprintf(stdout, F[side]);

        printf("# SM info on side %d:\n", side);
        sm_info[side].print(stdout);

        fflush(stdout);
    }

    t0 = seconds();

    my_sm(outfile, infile, sm_info, cpoly.nsides());

    fprintf(stdout, "\n# sm completed in %2.2lf seconds\n", seconds() - t0);
    fflush(stdout);

    mpz_clear(ell);
    mpz_clear(ell2);

    return 0;
}
