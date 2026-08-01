#include "cado.h" // IWYU pragma: keep
#include <cstdio>
#include <cstdlib>
#include <gmp.h>
#include "gfpkdlpolyselect.hpp"
#include "params.hpp"
#include "cxx_mpz.hpp"
#include "verbose.hpp"             // verbose_output_print

static void
declare_usage(cxx_param_list & pl)
{
    pl.declare_usage("p", "defines the prime p in GF(p^n)");
    pl.declare_usage("n", "defines the degree n in GF(p^n), atm n=2 is mandatory");
    pl.declare_usage("ell", "defines the prime order subgroup of GF(p^n), ell divides Phi_n(p)");
    pl.declare_usage("mnfs", "defines the number of number fields on second side ('g' side or 'linear' side), should be >= 1, and, for the moment, <= 2");
    pl.declare_usage("out",
            "filename where to write selected polynomials");
    verbose_decl_usage(pl);
}

int main (int argc, char const * argv[])
{
    cxx_mpz p, ell;
    unsigned int n = 0;
    unsigned int mnfs = 0;

    /* read params */
    cxx_param_list pl;
    declare_usage(pl);

    pl.process_command_line(argc, argv, false);

    if (!pl.parse("p", p)) // fill mpz_t p with p if a right value is given
        pl.fail("missing parameter: p");
    if (!pl.parse("ell", ell)){
        // fill mpz_t ell with ell if a right value is given
        ell = 42; // not used anyway...
    }

    if (!pl.parse("n", n) || (n != 2))
        pl.fail("missing parameter: n");

    if (!pl.parse("mnfs", mnfs) || (mnfs < 1) || (mnfs > 2))
        mnfs = 1; // let's default to 1.

    const char * out = pl.lookup_old("out");
    verbose_interpret_parameters(pl);

    /* check unused and print command line */
    if (pl.warn_unused())
        pl.fail("Unused parameters are given");
    pl.print_command_line(stdout);

    // gfpk_print_params(n,p,ell);

    int return_code = gfpkdlpolyselect( n, p, ell, mnfs, out);

    if (!return_code){
      fprintf(stderr, "no polynomial found.\n");
      return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
