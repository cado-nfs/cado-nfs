#include "cado.h" // IWYU pragma: keep
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cstdint>      // for uint64_t
#include <algorithm>     // for copy, max
#include <string>        // for string
#include <vector>        // for vector
#include <gmp.h>         // for gmp_randclear, gmp_randinit_default, gmp_ran...
#include "params.h"     // param_list
#include "cado_poly.h"  // cado_poly
#include "verbose.h"    // verbose_decl_usage
#include "typedefs.h"   // index_t 
#include "renumber.hpp" // renumber_t
#include "timing.h"     // seconds wct_seconds
#include "misc.h"     // size_disp

static void declare_usage(cxx_param_list & pl)
{
  param_list_decl_usage(pl, "poly", "input polynomial file");
  param_list_decl_usage(pl, "renumber", "input file for renumbering table ; exclusive with --build");
  param_list_decl_usage(pl, "lpbs", "large primes bounds (comma-separated list, for --build only)");
  param_list_decl_usage(pl, "check", "check the renumbering table");
  param_list_decl_usage(pl, "build", "build the renumbering table on the fly, instead of loading it (requires --lpbs)");
  param_list_decl_usage(pl, "bench", "bench lookup performance in the renumbering table");
  param_list_decl_usage(pl, "quiet", "do not print the renumbering table contents");
  param_list_decl_usage(pl, "dl", "interpret as DL-related data.");
  verbose_decl_usage(pl);
}

static void
usage (cxx_param_list & pl, const char *argv0)
{
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
}

// coverity[root_function]
int main(int argc, char const * argv[])
{
    int check = 0;
    int build = 0;
    int quiet = 0;
    int bench = 0;
    int for_dl = 0;
    const char * argv0 = argv[0];
    cxx_cado_poly cpoly;

    cxx_param_list pl;
    declare_usage(pl);
    renumber_t::builder_declare_usage(pl);

    param_list_configure_switch(pl, "check", &check);
    param_list_configure_switch(pl, "build", &build);
    param_list_configure_switch(pl, "quiet", &quiet);
    param_list_configure_switch(pl, "bench", &bench);
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
        double const wtt = wct_seconds();
        double const tt = seconds();
        std::vector<unsigned int> lpb(tab.get_nb_polys(),0);
        param_list_parse_uint_list(pl, "lpbs", lpb.data(), tab.get_nb_polys(), ",");
        tab.set_lpb(lpb);
        tab.build(pl, for_dl);

        if (bench) {
            printf("# Build time: %.2f (%.2f on cpu)\n",
                    wct_seconds() -  wtt,
                    seconds() -  tt);
        }

    } else {
        double const wtt = wct_seconds();
        double const tt = seconds();

        tab.read_from_file(renumberfilename, for_dl);

        if (bench) {
            printf("# Read time: %.2f (%.2f on cpu)\n",
                    wct_seconds() -  wtt,
                    seconds() -  tt);
        }
    }

    if (bench) {
        char buf[16];
        printf("# memory size of the renumbering table: %s\n",
                size_disp(tab.get_memory_size(), buf));
    }


    if (bench) {
        double tt = seconds();
        std::vector<size_t> counts(tab.get_nb_polys(),0);
        for(auto const x : tab) {
            /*
            if (x.side == tab.get_rational_side())
                printf("%" PRpr " rat %d\n", x.p, x.side);
            else
                printf("%" PRpr " %" PRpr " %d\n", x.p, x.r, x.side);
                */
            counts[x.side]++;
        }
        tt = seconds() - tt;
        printf("# full table traversal (%" PRIu64 " entries):"
                " %.3g (time per 1000000 steps: %.3g)\n",
                tab.get_size(), tt, tt*1.0e6/tab.get_size());
        printf("# ideals per side:");
        for(auto x : counts)
            printf(" %zu", x);
        printf("\n");
    }

    if (!quiet) {
        for (index_t i = 0; i < tab.get_size() ; i++) {
            std::string const s = tab.debug_data(i);
            printf ("%s\n", s.c_str());
        }
    }

    /* Check for all indices if mapping i <--> (p,r,side) works
     * consistently both ways.
     */
    if (check) {
        uint64_t nerrors = 0;
        for (index_t i = 0; i < tab.get_size(); i++)
        {
            if (tab.is_additional_column(i) || tab.is_bad(i)) continue;
            renumber_t::p_r_side x = tab.p_r_from_index(i);
            index_t j = tab.index_from_p_r(x);
            if (i == j) {
                if (!quiet) printf("## %" PRid ": Ok\n", i);
            } else {
                printf("#### Error:");
                printf(" i=%" PRid " p=%" PRpr " r=%" PRpr " side=%d",
                        i, x.p, x.r, x.side);
                x = tab.p_r_from_index(j);
                printf(" --> i=%" PRid " p=%" PRpr " r=%" PRpr " side=%d",
                        j, x.p, x.r, x.side);
                j = tab.index_from_p_r(x);
                x = tab.p_r_from_index(j);
                printf(" --> i=%" PRid " p=%" PRpr " r=%" PRpr " side=%d",
                        j, x.p, x.r, x.side);
                if (i != j) printf(" --> ...");
                printf("\n");

                nerrors++;
            }
        }
        printf("Number of errors: %" PRIu64 "\n", nerrors);
        if (nerrors) return EXIT_FAILURE;
    }

    if (bench) {
        gmp_randstate_t rstate;
        gmp_randinit_default(rstate);
        double tt;
        volatile unsigned long sum_h MAYBE_UNUSED = 0;
        int const nlookups = 1000 * 1000;

        std::vector<renumber_t::p_r_side> sample;

        tt = seconds();
        for(int i = 0 ; i < nlookups ; i++) {
            index_t const h = gmp_urandomm_ui(rstate, tab.get_size());
            renumber_t::p_r_side const x = tab.p_r_from_index(h);
            if (!tab.is_additional_column(h))
                sample.push_back(x);
        }
        printf("# Time for %d random lookups"
                " (p_r_from_index, arbitrary primes): %.3g\n",
                nlookups,
                seconds() - tt);

        tt = seconds();
        std::vector<renumber_t::p_r_side> sample_cached;
        for(auto const & x : sample) {
            if (x.p >> RENUMBER_MAX_LOG_CACHED) continue;
            sample_cached.push_back(x);
        }

        int nlookups_cached = 0;
        for(auto const & x : sample_cached) {
            nlookups_cached += 1;
            index_t const h = tab.index_from_p_r(x);
            sum_h += h;
        }
        printf("# Time for %d random lookups"
                " (index_from_p_r, cached primes): %.3g\n",
                nlookups_cached,
                seconds() - tt);
        if (nlookups_cached) {
            printf("# (scaled time for %d random lookups"
                    " (index_from_p_r, cached primes): %.3g)\n",
                    nlookups,
                    (seconds() - tt) * nlookups / nlookups_cached);
        }

        tt = seconds();
        for(auto const & x : sample) {
            index_t const h = tab.index_from_p_r(x);
            sum_h += h;
        }
        printf("# Time for %d random lookups"
                " (index_from_p_r, arbitrary primes): %.3g\n",
                nlookups,
                seconds() - tt);


        gmp_randclear(rstate);
    }

    return EXIT_SUCCESS;
}

