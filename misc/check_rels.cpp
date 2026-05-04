/* check_rels: check the factorization of the norm on both alg and rat side.
   Can also, in option, check the primality of ideal and correct uncomplete
   relation or relation with non prime ideal */

#include "cado.h" // IWYU pragma: keep
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <array>
#include <iostream>

#include <gmp.h>
#include "fmt/base.h"
#include "fmt/ostream.h"

#include "arith/mod_ul.h"
#include "cado_poly.hpp"
#include "filelist.hpp"
#include "filter_io.hpp"
#include "getprime.h"
#include "gzip.h"
#include "macros.h"
#include "mpz_poly.h"
#include "params.hpp"
#include "relation-tools.h"
#include "timing.h"
#include "typedefs.h"
#include "verbose.hpp"
#include "fstream_maybe_compressed.hpp"

#define FACTOR_DOES_NOT_DIVIDE 1UL
#define FACTOR_NOT_PRIME 2UL
#define FACTORIZATION_NOT_COMPLETE 4UL
#define FACTOR_ABOVE_LPB 8UL
#define REL_FULLY_FIXED 16UL

uint64_t nrels_read = 0, nrels_ok = 0, nrels_err = 0, nrels_fullyfixed = 0,
         nrels_doesnotdivide = 0, nrels_notprime = 0, nrels_notcomplete = 0,
         nrels_abovelpb = 0, nrels_fullycompleted = 0, nrels_fixednotprime = 0;
cxx_cado_poly cpoly;
std::vector<unsigned long> lpb; // full bounds, with the 2^...
unsigned long lpb_max = 0;
int verbose = 0;
int abhexa = 0;
int check_primality = 0; /* By default no primality check */
int fix_it = 0;          /* By default, we just check the rels */

/* used for counting time in different processes */
timingstats_dict_t stats;

static bool isprime(unsigned long p)
{
    /* TODO: what is fastest between this and mpz_probab_prime_p, which
     * is used by ulong_isprime?
     */
    modulusul_t m;
    modul_initmod_ul(m, p);
    bool const b = modul_isprime(m);
    modul_clearmod(m);
    return b;
}

template <typename relation_type>
static void rel_add_prime(relation_type & rel, p_r_values_t p, exponent_t e,
                          int side_index)
{
    for (auto & pse: rel.primes) {
        if (pse.p == p && pse.side == side_index) {
            pse.e += e;
            return;
        }
    }
    rel.primes.push_back({.p = p, .e = e, .side = side_index});
}

template <typename relation_type>
static inline void
print_error_line(relation_type const & rel,
                 typename relation_type::prime_type const & pes,
                 std::array<cxx_mpz, 2> const & norm, unsigned long err_type,
                 int will_be_fixed)
{
    auto const & [p, e, side_index] = pes;
    auto side = rel.active_sides[side_index];
    char const * str = (will_be_fixed) ? "Warning" : "Error";

    if (err_type == FACTOR_DOES_NOT_DIVIDE) {
        fmt::print(stderr,
                   "#   Error, given factor {} with exponent {} does "
                   "not divide the norm on side {}\n",
                   p, e, side);
    } else if (err_type == FACTOR_NOT_PRIME) {
        fmt::print(stderr,
                   "#   {}, given factor {} is not prime on side "
                   "{}\n",
                   str, p, side);
    } else if (err_type == FACTOR_ABOVE_LPB) {
        fmt::print(stderr,
                   "#   {}, given factor {} is greater than lpb "
                   "(={}) on side {}\n",
                   str, p, lpb[side], side);
    } else if (err_type == FACTORIZATION_NOT_COMPLETE) {
        fmt::print(stderr,
                   "#   {}, factorization of the norm on side {} is not "
                   "complete ({} is missing)\n",
                   str, side, norm[side_index]);
    }
    fflush(stderr);
}

/* assuming a given prime in a relation is actually _not_ prime, correct
 * it, and possibly add a few new primes to the relation.
 */
template <typename relation_type>
static inline void
factor_nonprime_ideal(relation_type & rel,
                      typename relation_type::prime_type & pes)
{
    using prime_type = relation_type::prime_type;
    auto [p, e, side_index] = pes;
    auto side = side_index;
    for (p_r_values_t pr: prime_range(2)) {
        exponent_t e_pr_in_p = 0;
        // funny linter bug, here.
        // NOLINTNEXTLINE(bugprone-infinite-loop)
        for( ; p % pr == 0 ; p = p / pr)
            e_pr_in_p++;

        if (e_pr_in_p) {
            const prime_type newprime = {.p = pr, .e = (e * e_pr_in_p), .side = side};
            if (p == 1) {
                /* we had a perfect power */
                pes = newprime;
                return;
            } else {
                rel.primes.push_back(newprime);
            }
        }

        if (isprime(p)) {
            const prime_type newprime = {.p = p, .e = e, .side = side};
            pes = newprime;
            return;
        }
    }
    /* remaining p is prime. The exponent is the same as we originally had
     * */
}

static int both_equal_to_1(std::array<cxx_mpz, 2> const & norm)
{
    return norm[0] == 1 && norm[1] == 1;
}

/* return 0 if everything is ok (factorization, primality, and complete)
 * else the ith bit of the return value is set if
 *  i = 1: a factor does not divide the norms (the check is then stopped)
 *  i = 2: a factor is not prime
 *  i = 3: the factorization of a norm is not complete
 *  i = 4: a factor is above a lpb
 *  i = 5: (only with fix_it != 0) non-prime factors were successfully factored
/bin/bash: q: command not found
 *  i = 6: less than 2 sides are used
 */
template <typename relation_type>
static unsigned long process_one_relation(relation_type & rel)
{
    using prime_type = relation_type::prime_type;
    std::array<cxx_mpz, 2> norm = {1, 1};
    const int nsides = cpoly.nsides();
    ASSERT_ALWAYS(nsides == 1 || nsides == 2);
    unsigned long err = 0;

    auto end = [](unsigned long err) {
        if (verbose) {
            fmt::print(stderr, "#   end with error status 0x{:x}\n", err);
            fflush(stderr);
        }
        return err;
    };

    if (verbose) {
        fmt::print(stderr, "# relation {} with (a,b) = ({},{}):\n", rel.num,
                   rel.a, rel.b);
        fmt::print(stderr, "# {}\n", rel);
        fflush(stderr);
    }

    /* compute the norms */
    for (int side_index = 0; side_index < nsides; side_index++) {
        int side = rel.active_sides[side_index];
        auto const & f = cpoly[side];
        mpz_poly_homogeneous_eval_siui(norm[side_index], f, rel.a, rel.b);
        mpz_abs(norm[side_index], norm[side_index]);
        if (verbose) {
            fmt::print(stderr, "#   norm on side {} = {}\n", side,
                       norm[side_index]);
            fflush(stderr);
        }
    }

    /* check for correctness of the factorization of the norms */
    for (auto const & pes: rel.primes) {
        auto const & [p, e, side_index] = pes;
        /* otherwise the relation is bad */
        ASSERT_ALWAYS(p != 0); /* could reveal a problem in parsing */
        ASSERT_ALWAYS(e > 0);  /* non positive exponent is not possible */
        for (int j = 0; j < e; j++) {
            if (!mpz_divisible_ui_p(norm[side_index], p)) {
                err |= FACTOR_DOES_NOT_DIVIDE;
                if (verbose != 0)
                    print_error_line(rel, pes, norm, FACTOR_DOES_NOT_DIVIDE,
                                     fix_it);
            } else {
                mpz_divexact_ui(norm[side_index], norm[side_index], p);
            }
        }
    }

    /* With an error in the factorization of the norm, no need to continue */
    if (err)
        return end(err);

    /* check primality of all ideals appearing in the relations */
    if (check_primality) {
        for (size_t i = 0; i < rel.primes.size(); i++) {
            auto & pes = rel.primes[i];
            if (!isprime(pes.p)) {
                err |= FACTOR_NOT_PRIME;
                if (verbose != 0)
                    print_error_line(rel, pes, norm, FACTOR_NOT_PRIME, fix_it);
                if (fix_it) /* if fix_it = 1, we factor it */
                    factor_nonprime_ideal(rel, pes);
            }
        }
    }

    /* Check that the product of the factors is equal to the norm. */
    for (int side_index = 0; side_index < nsides; side_index++) {
        if (norm[side_index] != 1) {
            err |= FACTORIZATION_NOT_COMPLETE;
            if (verbose) {
                prime_type const fake = {.p = 0, .e = 0, .side = side_index};
                print_error_line(rel, fake, norm, FACTORIZATION_NOT_COMPLETE,
                                 fix_it);
            }
        }
    }

    /* complete relation if it is asked and necessary */
    if (fix_it) {
        /* complete at least for primes up to 10000 (because of GGNFS and Msieve
         * that skip these primes) */
        for (auto p: prime_range(2, std::max(lpb_max, 10000UL))) {
            if (both_equal_to_1(norm))
                break;
            for (int side_index = 0; side_index < nsides;
                 side_index++) {
                int side = rel.active_sides[side_index];
                exponent_t e = 0;
                while (mpz_divisible_ui_p(norm[side_index], p)) {
                    e++;
                    mpz_divexact_ui(norm[side_index], norm[side_index], p);
                }
                if (e != 0) {
                    rel_add_prime(rel, p, e, side);
                    if (verbose != 0) {
                        fmt::print(
                            stderr,
                            "#   factorization of the norm on side {} is "
                            "not complete, need to add {}^{}\n",
                            side, p, e);
                    }
                }
            }
        }

        if (!both_equal_to_1(norm)) {
            fmt::print(stderr, "#   factorization of the norm is still not "
                               "complete on at least one side\n");
            fflush(stderr);
        }
    }

    /* check that ideals appearing in the relations are below the lpb. We
     * skip this check if all the primes dividing the norms are not known
     * (because we do not know if the missing primes are below or above the
     * lpbs). */
    if (both_equal_to_1(norm)) {
        for (auto const & pes: rel.primes) {
            auto const & [p, e, side_index] = pes;
            auto side = rel.active_sides[side_index];
            if (p > lpb[side]) {
                err |= FACTOR_ABOVE_LPB;
                if (verbose != 0)
                    print_error_line(rel, pes, norm, FACTOR_ABOVE_LPB, fix_it);
            }
        }
    }

    if (fix_it && (err & FACTOR_NOT_PRIME || err & FACTORIZATION_NOT_COMPLETE))
        if (both_equal_to_1(norm))
            err |= REL_FULLY_FIXED;

    return end(err);
}

/* Callback function called by filter_rels */

template <typename relation_type>
static void thread_callback(relation_type & rel, std::ostream * outfile)
{
    nrels_read++;
    int const ret = process_one_relation(rel);
    int is_printable;

    if (ret == 0) {
        nrels_ok++;
        is_printable = 1;
    } else if ((ret & REL_FULLY_FIXED) && !(ret & FACTOR_ABOVE_LPB)) {
        nrels_fullyfixed++;
        is_printable = 1;
    } else {
        nrels_err++;
        is_printable = 0;
    }

    if (ret & FACTOR_DOES_NOT_DIVIDE)
        nrels_doesnotdivide++;
    if (ret & FACTOR_NOT_PRIME) {
        if (ret & REL_FULLY_FIXED && is_printable)
            nrels_fixednotprime++;
        else if (!(ret & REL_FULLY_FIXED))
            nrels_notprime++;
    }
    if (ret & FACTORIZATION_NOT_COMPLETE) {
        if (ret & REL_FULLY_FIXED && is_printable)
            nrels_fullycompleted++;
        else if (!(ret & REL_FULLY_FIXED))
            nrels_notcomplete++;
    }
    if (ret & FACTOR_ABOVE_LPB)
        nrels_abovelpb++;

    if (outfile && is_printable)
        fmt::print(*outfile, "{}\n", rel);
}

static void declare_usage(cxx_param_list & pl)
{
    param_list_decl_usage(pl, "filelist",
                          "file containing a list of input files");
    param_list_decl_usage(pl, "basepath",
                          "path added to all files in filelist");
    param_list_decl_usage(pl, "poly", "polynomials file (mandatory)");
    param_list_decl_usage(pl, "abhexa",
                          "read and write a and b as hexa "
                          "(instead of decimal)");
    param_list_decl_usage(pl, "fixit", "Try to fix wrong relations");
    param_list_decl_usage(pl, "check_primality",
                          "check primality of "
                          "primes (default, no checking)");
    param_list_decl_usage(pl, "out",
                          "optional output file for correct or fixed "
                          "relations");
    param_list_decl_usage(pl, "lpb0", "large prime bound on side 0");
    param_list_decl_usage(pl, "lpb1", "large prime bound on side 1");
    param_list_decl_usage(pl, "lpbs",
                          "large primes bounds (comma-separated list) "
                          "(for MNFS)");
    param_list_decl_usage(pl, "v", "more verbose output");
    param_list_decl_usage(pl, "path_antebuffer", "path to antebuffer program");
}

template <typename relation_type>
static size_t filter(filelist const & input, std::ostream * outfile)
{
    return filter_rels<relation_type>(
        input.create_file_list(), nullptr, stats,
        [&](relation_type & rel) { thread_callback(rel, outfile); });
}

static size_t filter(filelist const & input, bool abhexa,
                     std::ostream * outfile)
{
    if (abhexa) {
        using relation_type = cado::relation_building_blocks::primes_block<
            prime_type_for_sieve_relations,
            cado::relation_building_blocks::ab_block<uint64_t, 16>>;
        return filter<relation_type>(input, outfile);
    } else {
        using relation_type = cado::relation_building_blocks::primes_block<
            prime_type_for_sieve_relations,
            cado::relation_building_blocks::ab_block<uint64_t, 10>>;
        return filter<relation_type>(input, outfile);
    }
}

int main(int argc, char const * argv[])
{
    cxx_param_list pl;

    declare_usage(pl);
    filelist::configure(pl);

    param_list_configure_switch(pl, "abhexa", &abhexa);
    param_list_configure_switch(pl, "fixit", &fix_it);
    param_list_configure_switch(pl, "v", &verbose);
    param_list_configure_switch(pl, "check_primality", &check_primality);
    cado::filter_io_details::configure(pl);

    pl.process_command_line(argc, argv, true);

    cado::filter_io_details::interpret_parameters(pl);
    verbose_interpret_parameters(pl);
    pl.print_command_line(stdout);

    fflush(stdout);

    filelist const input(pl, argc, argv);

    char const * polyfilename = param_list_lookup_string(pl, "poly");
    char const * outfilename = param_list_lookup_string(pl, "out");
    param_list_lookup_string(pl, "lpb0");
    param_list_lookup_string(pl, "lpb1");
    param_list_lookup_string(pl, "lpbs");

    if (pl.warn_unused())
        pl.fail("Unused parameters are given");

    if (!polyfilename)
        pl.fail("missing -poly command line argument\n");

    if (!cpoly.read(polyfilename))
        pl.fail("cannot read {}", polyfilename);

    pl.parse_per_side("lpb", lpb, cpoly.nsides(),
                      cado::params::copy_previous_side());
    for (auto & x: lpb)
        x = 1UL << x;
    lpb_max = *std::ranges::max_element(lpb);

    ofstream_maybe_compressed outfile;

    std::ostream * out = nullptr;

    if (outfilename && strcmp(outfilename, "-") == 0) {
        out = &std::cout;
        fmt::print("# Correct relations will be written to stdout\n");
        outfilename = nullptr;
    } else if (outfilename) {
        outfile.open(outfilename);
        ASSERT_ALWAYS(outfile.good());
        fmt::print("# Correct relations will be written to {}\n", outfilename);
        out = &outfile;
    }

    fmt::print("# Verbose output: {}\n", verbose ? "yes" : "no");
    fmt::print("# Correct wrong relations if possible: {}\n",
           fix_it ? "yes" : "no");
    fmt::print("# Check primality of ideal: {}\n", check_primality ? "yes" : "no");
    for (int i = 0; i < cpoly.nsides(); i++) {
        fmt::print("# On side {}, ", i);
        if (lpb[i] > 0)
            fmt::print("will check that norms of ideals are below {}\n", lpb[i]);
        else
            fmt::print("will not check the size of the norms of the ideals\n");
    }

    timingstats_dict_init(stats);

    filter(input, abhexa, out);

    fmt::print("Number of read relations: {}\n", nrels_read);
    fmt::print("Number of correct relations: {}\n", nrels_ok);
    if (fix_it) {
        fmt::print("Number of fixed relations: {}\n"
                   "   among which {} were fully completed\n"
                   "           and {} contained at least 1 non-prime factor\n",
                   nrels_fullyfixed, nrels_fullycompleted, nrels_fixednotprime);
        fmt::print("Number of wrong relations: {}\n", nrels_err);
        fmt::print("   among which {} had a factor not dividing the norm\n"
                   "           and {} could not be fully completed\n"
                   "           and {} contained an ideal larger than a lpb\n",
                   nrels_doesnotdivide, nrels_notcomplete, nrels_abovelpb);
    } else {
        fmt::print("Number of wrong relations: {}\n", nrels_err);
        fmt::print("   among which {} had a factor not dividing the norm\n"
                   "           and {} contained at least 1 non-prime factor\n"
                   "           and {} were not complete\n"
                   "           and {} contained an ideal larger than a lpb\n",
                   nrels_doesnotdivide, nrels_notprime, nrels_notcomplete,
                   nrels_abovelpb);
    }

    if (outfilename)
        outfile.close();

    timingstats_dict_add_mythread(stats, "main");
    timingstats_dict_disp(stats);
    timingstats_dict_clear(stats);

    return (nrels_err) ? EXIT_FAILURE : EXIT_SUCCESS;
}
