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
    param_list_decl_usage(pl, "p", "defines the prime p in GF(p^n)");
    param_list_decl_usage(pl, "n", "defines the degree n in GF(p^n), atm n=2 is mandatory");
    param_list_decl_usage(pl, "ell", "defines the prime order subgroup of GF(p^n), ell divides Phi_n(p)");
    param_list_decl_usage(pl, "mnfs", "defines the number of number fields on second side (\'g\' side or \'linear\' side), should be >= 1, and, for the moment, <= 2");
    param_list_decl_usage(pl, "out",
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

    param_list_process_command_line(pl, &argc, &argv, false);

    if (!param_list_parse_mpz(pl, "p", p)) // fill mpz_t p with p if a right value is given
        pl.fail("missing parameter: p");
    if (!param_list_parse_mpz(pl, "ell", ell)){
        // fill mpz_t ell with ell if a right value is given
        ell = 42; // not used anyway...
    }

    if (!param_list_parse_uint(pl, "n", &n) || (n != 2))
        pl.fail("missing parameter: n");

    if (!param_list_parse_uint(pl, "mnfs", &mnfs) || (mnfs < 1) || (mnfs > 2))
        mnfs = 1; // let's default to 1.

    const char *out = param_list_lookup_string (pl, "out");
    verbose_interpret_parameters(pl);

    /* check unused and print command line */
    if (param_list_warn_unused(pl))
        pl.fail("Unused parameters are given");
    param_list_print_command_line (stdout, pl);

    // gfpk_print_params(n,p,ell);

    int return_code = gfpkdlpolyselect( n, p, ell, mnfs, out);

    if (!return_code){
      fprintf(stderr, "no polynomial found.\n");
      return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
