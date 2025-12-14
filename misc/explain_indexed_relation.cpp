#include "cado.h" // IWYU pragma: keep
#include <cstdio>
#include <cstdlib>
#include <set>
#include <string>        // for string
#include <vector>        // for vector
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <memory>

#include "fmt/base.h"
#include "fmt/format.h"

#include "mpz_poly.h"
#include "params.h"     // param_list
#include "cado_poly.h"  // cado_poly
#include "verbose.h"    // verbose_decl_usage
#include "typedefs.h"   // index_t 
#include "renumber.hpp" // renumber_t
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
    param_list_decl_usage(pl, "renumber", "input file for renumbering table ; exclusive with --build");
    param_list_decl_usage(pl, "relations", "explain indexed relations from this file.");
    verbose_decl_usage(pl);
}

static void
usage (cxx_param_list & pl, const char *argv0)
{
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
}

static std::string rewrite_carets(std::string const & s)
{
    std::string t;
    for(auto c : s) {
        if (c == '^') {
            t += "**";
        } else {
            t += c;
        }
    }
    return t;
}

struct command_line {
    int python = 0;
    int raw = 0;
    int dl = 0;
    int all = 0;
    int build = 0;
    int skip_ideal_checks = 0;
    const char *polyfilename = nullptr;
    const char *renumberfilename = nullptr;
    const char *relationsfilename = nullptr;
    std::vector<unsigned int> lpb;
    cxx_cado_poly cpoly;

    static void declare_usage(cxx_param_list & pl) {
        param_list_decl_usage(pl, "dl", "interpret as DL-related data.");
        param_list_decl_usage(pl, "python", "output python directly (skip the sage preparser).");
        param_list_decl_usage(pl, "raw", "output only machine readable contents.");
        param_list_decl_usage(pl, "all", "output code with the definition of all ideals.");
        param_list_decl_usage(pl, "build", "build the renumbering table on the fly, instead of loading it (requires --lpbs)");
        param_list_decl_usage(pl, "skip-ideal-checks", "do not include the sagemath checks for ideals.");
        param_list_decl_usage(pl, "lpbs", "large primes bounds (comma-separated list, for --build only)");
        param_list_decl_usage(pl, "poly", "input polynomial file");
    }

    void configure_switches(cxx_param_list & pl) {
        param_list_configure_switch(pl, "dl", &dl);
        param_list_configure_switch(pl, "python", &python);
        param_list_configure_switch(pl, "raw", &raw);
        param_list_configure_switch(pl, "all", &all);
        param_list_configure_switch(pl, "build", &build);
        param_list_configure_switch(pl, "skip-ideal-checks", &skip_ideal_checks);
    }

    void lookup_parameters(cxx_param_list & pl) {
        polyfilename = param_list_lookup_string(pl, "poly");
        renumberfilename = param_list_lookup_string(pl, "renumber");
        relationsfilename = param_list_lookup_string(pl, "relations");
        param_list_parse(pl, "lpbs", lpb);
        if (!cado_poly_read (cpoly, polyfilename))
        {
            fmt::print (stderr, "Error reading polynomial file\n");
            exit (EXIT_FAILURE);
        }

    }

    void check_inconsistencies(const char * argv0, cxx_param_list & pl) const {
        if (python && raw)
        {
            fmt::print (stderr, "Error, -python and -raw are incompatible\n");
            usage (pl, argv0);
        }

        if (relationsfilename && raw)
        {
            fmt::print (stderr, "Error, -relations and -raw are incompatible\n");
            usage (pl, argv0);
        }

        if (!all && raw)
        {
            fmt::print (stderr, "Error, -raw requires -all\n");
            usage (pl, argv0);
        }

        if (!polyfilename)
        {
            fmt::print (stderr, "Error, missing -poly command line argument\n");
            usage (pl, argv0);
        }
        if (!renumberfilename && !build) {
            fmt::print (stderr, "Error, missing -renumber command line argument\n");
            usage (pl, argv0);
        }
        if (renumberfilename && build) {
            fmt::print (stderr, "Error, --build and -renumber are exclusive\n");
            usage (pl, argv0);
        }
        if (lpb.empty() && build) {
            fmt::print (stderr, "Error, --build requires -lpbs\n");
            usage (pl, argv0);
        }
        if (!lpb.empty() && !build) {
            fmt::print (stderr, "Error, --lpbs is only valid with --build\n");
            usage (pl, argv0);
        }
        if (build && lpb.size() != (size_t) cpoly->nb_polys) {
            fmt::print (stderr, "Error, --build requires one lpb per side\n");
            usage (pl, argv0);
        }

    }
};


static void output_prologue(command_line const & cmdline)
{
    cxx_cado_poly const & cpoly(cmdline.cpoly);

    if (cmdline.python) {
        std::vector<std::pair<std::string, std::string>> const imports {
            //{ "sage.categories.category", "" },
            //{ "sage.categories.commutative_rings", "" },
            //{ "sage.categories.commutative_additive_groups", "" },
            { "sage.rings.polynomial.polynomial_element", "Polynomial" },
            { "sage.rings.polynomial.polynomial_ring_constructor", "PolynomialRing" },
            { "sage.rings.integer_ring", "ZZ" },
            { "sage.rings.integer", "Integer" },
            { "sage.rings.number_field.number_field", "NumberField" },
            { "sage.misc.misc_c", "prod" },
            // { "sage.repl.preparse", "preparse" },
        };
        for(auto const & i : imports)
            if (! i.second.empty()) {
                fmt::print("from {} import {}\n", i.first, i.second);
            } else {
                fmt::print("import {}\n", i.first);
            }
        fmt::print("ZP = PolynomialRing(ZZ, names=('x',)); x = ZP.gen()\n");
    } else {
        fmt::print("ZP.<x> = ZZ[]\n");
    }
    for(int side = 0 ; side < cpoly->nb_polys ; side++) {
        std::ostringstream os;
        os << cxx_mpz_poly(cpoly->pols[side]);
        if (cmdline.python) {
            fmt::print("K{0}=NumberField({1}, names=('alpha{0}',)); alpha{0}=K{0}.gen()\n",
                    side, rewrite_carets(os.str()));
        } else {
            fmt::print("K{0}.<alpha{0}>=NumberField({1})\n", side, os.str());
        }
        fmt::print("OK{0}=K{0}.maximal_order()\n", side);
        fmt::print("J{0}=OK{0}.fractional_ideal(1,alpha{0})**-1\n", side);
        if (!cmdline.dl)
            fmt::print("is_sq=lambda I:prod([v[1]%2==0 for v in I.factor()])\n");
        fmt::print("is_1=lambda x:x==1\n");
        fmt::print("is_prime_ideal=lambda x:x.is_prime()\n");

        fmt::print("def check(f, s, e):\n\tif not f(e):\n\t\traise AssertionError(\"Failed check on \"+s)\n");
    }
}

static void print_one_ideal(index_t c, renumber_t const & tab, command_line const & cmdline, bool printing_all)
{
        if (!cmdline.raw || !printing_all) {
            auto s = tab.debug_data_sagemath(c);

            if (tab.is_additional_column(c) && tab.has_merged_additional_column()) {
                /* The J0J1 ideal does not exist, really. We'll just
                 * report it as a comment in the text
                 */
                fmt::print("# I{:x}={}; # virtually the product of J0^-1 and J1^_1\n", c, s);
                return;
            }
            if (cmdline.python) {
                fmt::print("I{:x}={};", c, rewrite_carets(s));
            } else {
                fmt::print("I{:x}={};", c, s);
            }
            if (!cmdline.skip_ideal_checks && !tab.is_additional_column(c))
                fmt::print(" check(is_prime_ideal,\"I{0:x}\",I{0:x});", c);
            fmt::print("\n");
            if (printing_all)
                fmt::print(" all_ideals.append(I{:x})\n", c);
        } else {
            fmt::print("{}\n", tab.debug_data_machine_description(c));
        }
}

static std::set<index_t> print_all_ideals(renumber_t const & tab, command_line const & cmdline)
{
    std::set<index_t> printed;
    if (!cmdline.raw)
        fmt::print("all_ideals=[]\n");

    for(index_t c = 0 ; c < tab.size() ; ++c) {
        printed.insert(c);
        print_one_ideal(c, tab, cmdline, true);
    }

    return printed;
}

static
void examine_one_relation(renumber_t const & tab, std::set<index_t> & printed, indexed_relation & rel, command_line const & cmdline)
{
    fmt::print("a={}; b={}\n", rel.az, rel.bz);

    std::vector<std::vector<std::string>> ideals_per_side(cmdline.cpoly->nb_polys);

    for(int side = 0 ; side < cmdline.cpoly->nb_polys ; side++) {// {{{
        /* In the DL case, the factorization has to include J, and
         * the thing that we're factoring is really the principal
         * ideal generated by a-b*alpha. The only catch is that we
         * have a single additional column that reflects the presence
         * of the ideal J on both sides, and furthermore the exponent
         * is -1. So it's easier to put it in the denominator.
         *
         * In the factorization case, the ideal J is forcibly
         * skipped, and does not appear in the factorization. It is
         * trivial anyway to retrieve its valuation from the number
         * of (a,b) pairs used during sqrt.
         *
         * Note that free relations, in any case, have degree zero in
         * alpha and therefore must not include J
         */
        if (rel.b == 0) {
            fmt::print("ab{0}=(OK{0}.fractional_ideal({1}-{2}*alpha{0}))\n",
                    side, rel.az, rel.bz);
        } else {
            // there's something very fishy in the handling of
            // positional arguments with the following format.
            // Every once in a while, I get 'ab780' instead of a,
            // but _not_ when under gdb.
            // fmt::print("ab{0}=(OK{0}.fractional_ideal({1}-{2}*alpha{0})*J{0})\n",
            //         side, rel.az, rel.bz);
            auto gen = fmt::format("{}-{}*alpha{}", rel.az, rel.bz, side);
            auto ab = fmt::format("ab{}", side);
            auto I = fmt::format("OK{}.fractional_ideal({})", side, gen);
            auto J = fmt::format("J{}", side);

            fmt::print("{}=({}*{})\n", ab, I, J);
        }
    }// }}}

    for(auto const & c : rel.data) {
        auto x = printed.find(c);
        auto it = tab.p_r_from_index(c);
        if (x == printed.end()) {
            printed.insert(c);
            print_one_ideal(c, tab, cmdline, false);
        }
        if (tab.is_additional_column(c)) {
            if (tab.has_merged_additional_column()) {
                ideals_per_side[0].emplace_back("1");
                ideals_per_side[1].emplace_back("1");
            } else {
                ideals_per_side[it.side].emplace_back("1");
            }
        } else {
            ideals_per_side[it.side].emplace_back(fmt::format("I{:x}", c));
        }
    }

    for(int side = 0 ; side < cmdline.cpoly->nb_polys ; side++) {
        bool empty = true;
        std::ostringstream os;
        for(auto const & s : ideals_per_side[side]) {
            if (!empty) os << "*";
            empty = false;
            os << s;
        }
        if (empty) os << "OK" << side;
        if (cmdline.dl) {
            std::cout << fmt::format("check(is_1, \"ab{0}<<{2},{3}>>\", ab{0}/({1}))", side, os.str(), rel.az, rel.bz) << "\n";
        } else {
            std::cout << fmt::format("check(is_sq, \"ab{0}<<{2},{3}>>\", ab{0}/({1}))", side, os.str(), rel.az, rel.bz) << "\n";
        }
    }
}

int main(int argc, char const * argv[])
{
    command_line cmdline;
    const char *argv0 = argv[0];

    cxx_param_list pl;
    declare_usage(pl);
    renumber_t::builder_declare_usage(pl);
    command_line::declare_usage(pl);

    cmdline.configure_switches(pl);

    argv++, argc--;
    if (argc == 0)
      usage (pl, argv0);

    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv))
            continue;
        fmt::print(stderr, "Unhandled parameter {}\n", argv[0]);
        usage (pl, argv0);
    }
    /* print command-line arguments */
    verbose_interpret_parameters(pl);
    param_list_print_command_line (stdout, pl);
    fflush(stdout);

    cmdline.lookup_parameters(pl);
    renumber_t::builder_lookup_parameters(pl);

    cmdline.check_inconsistencies(argv0, pl);


    renumber_t tab(cmdline.cpoly);

    if (cmdline.build) {
        tab.set_lpb(cmdline.lpb);
        tab.build(pl, cmdline.dl);
    } else {
        tab.read_from_file(cmdline.renumberfilename, cmdline.dl);
        tab.recompute_debug_number_theoretic_stuff();
    }

    if (!cmdline.raw) {
        output_prologue(cmdline);
    }

    std::set<index_t> printed;

    if (cmdline.all)
        printed = print_all_ideals(tab, cmdline);

    if (cmdline.relationsfilename) {
        std::istream * cin;
        std::unique_ptr<std::ifstream> real_file;
        auto const & s = cmdline.relationsfilename;
        if (s == std::string("-")) {
            cin = &std::cin;
        } else {
            real_file = std::make_unique<std::ifstream>(s);
            cin = real_file.get();
        }

        int line = 0;
        for(std::string s ; std::getline(*cin, s) ;) {
            ++line;
            if (s.empty() || s[0] == '#')
                continue;

            std::istringstream is(s);
            indexed_relation rel;

            if (!(is >> rel)) {
                throw std::runtime_error(fmt::format("Parse error on line {}: {}\n", line, s));
            }

            fmt::print("print(\"{}\")\n", s);

            examine_one_relation(tab, printed, rel, cmdline);
        }
    }

    return EXIT_SUCCESS;
}

