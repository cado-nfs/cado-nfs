#include "cado.h" // IWYU pragma: keep
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cstdint>      // for uint64_t
#include <algorithm>     // for copy, max
#include <set>
#include <string>        // for string
#include <vector>        // for vector
#include <sstream>
#include <iostream>
#include <map>
#include <gmp.h>         // for gmp_randclear, gmp_randinit_default, gmp_ran...
#include <fmt/format.h>
#include "cxx_mpz.hpp"
#include "params.h"     // param_list
#include "cado_poly.h"  // cado_poly
#include "verbose.h"    // verbose_decl_usage
#include "typedefs.h"   // index_t 
#include "renumber.hpp" // renumber_t
#include "timing.h"     // seconds wct_seconds
#include "misc.h"     // size_disp
#include "indexed_relation.hpp"

/* This program takes (from stdin) an indexed relation (typically as
 * found in the .purged.gz file, or in the .dup1/ subdirectories after
 * dup2) and spells out the correspondence between indices and ideals in
 * sagemath format. The renumber table is needed, but it can be computed
 * on the fly.
 *
 * Do not forget to add the --dl flag if you have DLP data!
 *
 * The program output should be _almost_ accurate. In some cases though,
 * there is an expected (and fixable) glitch, which is that there is a
 * "combined" ideal J0J1 which exists on both sides, and is not correctly
 * printed.
 *
 * Except for the caveat above, the output should be pasteable into a
 * sagemath console, and asserts should all pass.
 *
 * Note that presently, this _only_ works with very small number fields,
 * because we let sage compute the maximal order (out of laziness,
 * somehow).
 *
 * HOW TO USE (EXAMPLE)
 * --------------------
 *
 * ./cado-nfs.py --workdir $PWD/work/ 90377629292003121684002147101760858109247336549001090677693
 * make explain_indexed_relation
 * $bindir/misc/explain_indexed_relation -poly work/c60.poly -renumber work/c60.renumber.gz < $(zcat work/c60.purged.gz | sort -R | head -n 1000) > check.sage
 * docker run --rm -t -v $PWD:/host sagemath/sagemath sage /host/check3.sage
 */

static void declare_usage(cxx_param_list & pl)
{
  param_list_decl_usage(pl, "poly", "input polynomial file");
  param_list_decl_usage(pl, "renumber", "input file for renumbering table ; exclusive with --build");
  param_list_decl_usage(pl, "lpbs", "large primes bounds (comma-separated list, for --build only)");
  param_list_decl_usage(pl, "build", "build the renumbering table on the fly, instead of loading it (requires --lpbs)");
  param_list_decl_usage(pl, "dl", "interpret as DL-related data.");
  verbose_decl_usage(pl);
}

static void
usage (cxx_param_list & pl, char *argv0)
{
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
}


int
main (int argc, char *argv[])
{
    int build = 0;
    int for_dl = 0;
    char *argv0 = argv[0];
    cxx_cado_poly cpoly;

    cxx_param_list pl;
    declare_usage(pl);
    renumber_t::builder_declare_usage(pl);

    param_list_configure_switch(pl, "build", &build);
    param_list_configure_switch(pl, "dl", &for_dl);

    argv++, argc--;
    if (argc == 0)
      usage (pl, argv0);

    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv))
            continue;
        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        usage (pl, argv0);
    }
    /* print command-line arguments */
    verbose_interpret_parameters(pl);
    param_list_print_command_line (stdout, pl);
    fflush(stdout);

    const char *polyfilename = param_list_lookup_string(pl, "poly");
    const char *renumberfilename = param_list_lookup_string(pl, "renumber");

    renumber_t::builder_lookup_parameters(pl);

    if (polyfilename == NULL)
    {
      fprintf (stderr, "Error, missing -poly command line argument\n");
      usage (pl, argv0);
    }
    if (renumberfilename == NULL && !build) {
      fprintf (stderr, "Error, missing -renumber command line argument\n");
      usage (pl, argv0);
    }
    if (renumberfilename != NULL && build) {
      fprintf (stderr, "Error, --build and -renumber are exclusive\n");
      usage (pl, argv0);
    }
    if (!param_list_lookup_string(pl, "lpbs") && build) {
      fprintf (stderr, "Error, --build requires -lpbs\n");
      usage (pl, argv0);
    }
    if (param_list_lookup_string(pl, "lpbs") && !build) {
      fprintf (stderr, "Error, --lpbs is only valid with --build\n");
      usage (pl, argv0);
    }

    if (!cado_poly_read (cpoly, polyfilename))
    {
      fprintf (stderr, "Error reading polynomial file\n");
      exit (EXIT_FAILURE);
    }

    renumber_t tab(cpoly);

    if (build) {
        std::vector<unsigned int> lpb(tab.get_nb_polys(),0);
        param_list_parse_uint_list(pl, "lpbs", lpb.data(), tab.get_nb_polys(), ",");
        tab.set_lpb(lpb);
        tab.build(pl, for_dl);
    } else {
        tab.read_from_file(renumberfilename, for_dl);
        tab.recompute_debug_number_theoretic_stuff();
    }

    // sage preparser is woefully inefficient. A large source file will
    // easily takes many minutes to preparse
    fmt::print("ZP.<x> = PolynomialRing(Integers())\n");
    for(int side = 0 ; side < cpoly->nb_polys ; side++) {
        std::ostringstream os;
        os << cxx_mpz_poly(cpoly->pols[side]);
        fmt::print("K{0}.<alpha{0}>=NumberField({1})\n", side, os.str());
        fmt::print("OK{0}=K{0}.maximal_order()\n", side);
        fmt::print("J{0}=OK{0}.fractional_ideal(1,alpha{0})^-1\n", side);
        if (!for_dl)
            fmt::print("is_sq=lambda I:prod([v[1]%2==0 for v in I.factor()])\n");
        fmt::print("is_1=lambda x:x==1\n");

        fmt::print("def check(f, s, e):\n\tif not f(e):\n\t\traise AssertionError(\"Failed check on \"+s)\n");
    }

    std::set<index_t> printed;

    int line = 0;
    for(std::string s ; std::getline(std::cin, s) ;) {
        ++line;
        if (s.empty() || s[0] == '#')
            continue;

        std::istringstream is(s);
        indexed_relation rel;

        if (!(is >> rel)) {
            throw std::runtime_error(fmt::format("Parse error on line {}: {}\n", line, s));
        }

        fmt::print("print(\"{}\")\n", s);
        fmt::print("a={}; b={}\n", rel.az, rel.bz);

        std::vector<std::vector<std::string>> ideals_per_side(cpoly->nb_polys);

        for(int side = 0 ; side < cpoly->nb_polys ; side++) {
            /* In the DL case, the factorization has to include J, and
             * the thing on the left hand side is really the principal
             * ideal generated by a-b*alpha
             *
             * In the factorization case, the ideal J is forcibly
             * skipped, and does not appear in the factorization. It is
             * trivial anyway to retrieve its valuation from the number
             * of (a,b) pairs used during sqrt.
             *
             * Note that free relations, in any case, have degree zero in
             * alpha and therefore must not include J
             */
            if (for_dl || rel.b == 0) {
                fmt::print("ab{0}=(OK{0}.fractional_ideal({1}-{2}*alpha{0}))\n",
                        side, rel.az, rel.bz);
            } else {
                fmt::print("ab{0}=(OK{0}.fractional_ideal({1}-{2}*alpha{0})*J{0})\n",
                        side, rel.az, rel.bz);
            }
        }

        for(auto const & c : rel.data) {
            auto x = printed.find(c);
            auto it = tab.p_r_from_index(c);
            if (x == printed.end()) {
                printed.insert(c);
                auto s = tab.debug_data_sagemath(c);

                fmt::print("I{:x}={};", c, s);
                if (!tab.is_additional_column(c))
                    fmt::print(" check(is_prime,\"I{0:x}\",I{0:x})", c);
                fmt::print("\n");
            }
            ideals_per_side[it.side].push_back(fmt::format("I{:x}", c));
        }

        for(int side = 0 ; side < cpoly->nb_polys ; side++) {
            bool empty = true;
            std::ostringstream os;
            for(auto const & s : ideals_per_side[side]) {
                if (!empty) os << "*";
                empty = false;
                os << s;
            }
            if (empty) os << "OK" << side;
            if (for_dl) {
                std::cout << fmt::format("check(is_1, \"ab{0}<<{2},{3}>>\", ab{0}/{1})", side, os.str(), rel.az, rel.bz) << "\n";
            } else {
                std::cout << fmt::format("check(is_sq, \"ab{0}<<{2},{3}>>\", ab{0}/({1}))", side, os.str(), rel.az, rel.bz) << "\n";
            }
        }
    }


    return EXIT_SUCCESS;
}

