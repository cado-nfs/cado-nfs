#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <memory>

#include <gmp.h>

#include "gmp_aux.h"
#include "badideals.hpp"
#include "cxx_mpz.hpp"
#include "mpz_poly.h"
#include "params.h"
#include "cado_poly.h"
#include "sm_utils.hpp"

using namespace std;

/* This program is intended to replicate exactly the behaviour of the
 * scripts/badideals.mag program of old.
 *
 * Now, we acknowledge several bizarre things in the form of the
 * .badideals and .badidealinfo file.
 *
 */
/*
  A line in badidealinfo has the form:
    p k r side v0 v1 ... vs
  If means that if (a/b) mod p^k is r (with the usual convention for
  projective roots, see below), then the s corresponding columns for the
  given side should be filled with the values v0, v1, ... vs.

  More precisely, if vi is positive, then this is indeed the value, but if vi
  is negative, then the value (e - |vi|) should be put in the column, where e
  is the valuation of p in the norm.

  Remarks:
  - there should be a line of the form
     p,(r mod p):side: s
    in the .badideals file, in order to "declare" the appropriate number of
    columns for this (p,r).
  - the badidealinfo is supposed to cover all the cases, but not necessarily
    in a simple way (the set of congruences might involve some "mod p" and
    some "mod p^2" rules, for instance).

  Projective roots:
  If we are in the case where (b/a) == 0 mod p, we write "p^k + (1/r)"
  instead of r.
  E.g.
    2 3 10 1 -2 2
  means that we are dealing with the case (a:b) = (1:2) mod 2^3, and that in
  that case, we should write (e-2) and 2 in the two corresponding columns on
  the side 1.

  p and r are read in basis 10 in the badidealsinfo file.
 */

/*
 * For reference, here is a test case:
 *
 * > clear; DEBUG:=1; load "scripts/badideals.mag";
 * Loading "scripts/badideals.mag"
 * 2,2:0: 2
 * d,5:0: 2
 * 2,1:1: 2
 * 2,2:1: 2
 * 3,1:1: 2
 * # bad ideals for poly[0]=240*x^4 + 4327846*x^3 - 11463289949*x^2 -
 * 48524924823332*x + 99623823815957215
 * 2 2 4 0 1 -1
 * 2 2 6 0 -1 1
 * 13 2 18 0 1 1
 * 13 2 31 0 1 1
 * 13 2 44 0 1 1
 * 13 2 57 0 1 1
 * 13 2 70 0 1 1
 * 13 2 83 0 1 1
 * 13 2 109 0 1 1
 * 13 2 122 0 1 1
 * 13 2 135 0 1 1
 * 13 2 148 0 1 1
 * 13 2 161 0 1 1
 * 13 2 5 0 1 -1
 * 13 2 96 0 -1 1
 * # bad ideals for poly[1]=840*x^5 - 16210*x^4 + 2610983*x^3 -
 * 2560484269*x^2 - 34656366009*x + 253976285986533
 * 2 1 1 1 1 -1
 * 2 2 4 1 -1 1
 * 2 2 6 1 1 -1
 * 3 2 7 1 1 1
 * 3 2 1 1 1 -1
 * 3 2 4 1 -1 1
 *
 *
 * The bizarre points are the following:
 *
 *  - That we have both a .badideals and a .badidealinfo file is very
 *    odd.
 *    Actually maybe this info would be just as well saved in the .roots file,
 *    perhaps.
 *  - The "side" info in the .badidealinfo file is weird.
 *  - The 13 lines above are inelegant. A more concise way would be:
 *        13 2 5 0 1 -1
 *        13 2 96 0 -1 1
 *        13 1 5 0 1 1
 *    But that would be non commutative, so this would be a clear change
 *    in semantics.
 *  - The file format, for now, is not ready to handle larger-degree
 *    prime ideals.
 */ 

static void badideals_declare_usage(cxx_param_list & pl)/*{{{*/
{
    param_list_decl_usage(pl, "badideals", "badideals file");
    param_list_decl_usage(pl, "badidealinfo", "badidealinfo file");
    param_list_decl_usage(pl, "polystr", "polynomial (string)");
    param_list_decl_usage(pl, "poly", "polynomial file");
    param_list_decl_usage(pl, "seed", "random seed");
    param_list_decl_usage(pl, "ell", "ell (for computing default number of maps ; not used for bad ideals)");
}/*}}}*/

static void usage(param_list_ptr pl, char const ** argv, const char * msg = nullptr)/*{{{*/
{
    param_list_print_usage(pl, argv[0], stderr);
    if (msg) {
        fprintf(stderr, "%s\n", msg);
    }
    exit(EXIT_FAILURE);
}/*}}}*/

// coverity[root_function]
int main(int argc, char const * argv[])
{
    char const ** original_argv;

    setvbuf(stderr, nullptr, _IONBF, 0);
    setvbuf(stdout, nullptr, _IONBF, 0);

    cxx_param_list pl;

    badideals_declare_usage(pl);
    param_list_configure_alias(pl, "polystr", "f");

    original_argv = argv;

    argv++,argc--;
    /* switches, if any. See below */
    /* aliases, if any. See below */

    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        /* Do perhaps some other things on the argument that haven't
         * been eaten at all. Like check whether it is a valid file to
         * source in order to get more options. See
         * param_list_read_stream and param_list_read_file for that. */
        fprintf (stderr, "Unknown option: %s\n", argv[0]);
        usage(pl, original_argv);
    }

    cxx_gmp_randstate state;

    unsigned long seed = 1;
    if (param_list_parse_ulong(pl, "seed", &seed)) {
        gmp_randseed_ui(state, seed);
    }

    const char * tmp;
    if ((tmp = param_list_lookup_string(pl, "polystr")) != nullptr) {
        int const side = 0;
        cxx_mpz_poly f(tmp);

        vector<badideal> const badideals = badideals_for_polynomial(f, side);
        cout << "--- .badideals data ---\n";
        for(auto const & b : badideals)
            b.print_dot_badideals_file(cout, side);

        cout << "--- .badidealinfo data ---\n";
        cout << "# bad ideals for poly"<<side<<"=" << f.print_poly("x") << "\n";
        for(auto const & b : badideals)
            b.print_dot_badidealinfo_file(cout, side);
    } else if ((tmp = param_list_lookup_string(pl, "poly")) != nullptr) {
        cado_poly cpoly;
        cado_poly_init(cpoly);
        cado_poly_read(cpoly, tmp);

        /* We're no longer using this functionality, but it's still
         * present.
         */
        const char * fbname = param_list_lookup_string(pl, "badideals");
        const char * fbiname = param_list_lookup_string(pl, "badidealinfo");

        std::unique_ptr<std::ofstream> fb;
        std::unique_ptr<std::ofstream> fbi;

        if (fbname)
            fb = std::make_unique<std::ofstream>(fbname);
        if (fbname)
            fbi = std::make_unique<std::ofstream>(fbiname);

        for(int side = 0 ; side < cpoly->nb_polys ; side++) {
            cxx_mpz_poly f(cpoly->pols[side]);
            if (f->deg == 1) continue;
            vector<badideal> const badideals = badideals_for_polynomial(f, side);
            if (fb) {
                for(auto const & b : badideals) {
                    b.print_dot_badideals_file(*fb, side);
                }
            }
            if (fbi) {
                *fbi << "# bad ideals for poly"<<side<<"=" << f.print_poly("x") << "\n";
                for(auto const & b : badideals) {
                    b.print_dot_badidealinfo_file(*fbi, side);
                }
            }
        }
        cxx_mpz ell;
        if (param_list_parse_mpz(pl, "ell", ell)) {
            for(int side = 0 ; side < cpoly->nb_polys ; side++) {
                sm_side_info const sm(cpoly->pols[side], ell, false);
                cout << "# nmaps" << side << " " << sm.nsm << "\n";
            }
        }
        cado_poly_clear(cpoly);
    } else {
        usage(pl, original_argv, "-poly or -polystr are mandatory");
    }
}
