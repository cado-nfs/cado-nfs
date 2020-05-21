#include "cado.h"
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "mod_ul.c"
#include "params.h"     // param_list
#include "cado_poly.h"  // cado_poly
#include "verbose.h"    // verbose_decl_usage
#include "typedefs.h"   // index_t 
#include "renumber.hpp" // renumber_t

static void declare_usage(cxx_param_list & pl)
{
  param_list_decl_usage(pl, "poly", "input polynomial file");
  param_list_decl_usage(pl, "renumber", "input file for renumbering table ; exclusive with --build");
  param_list_decl_usage(pl, "lpbs", "large primes bounds (comma-separated list, for --build only)");
  param_list_decl_usage(pl, "check", "check the renumbering table");
  param_list_decl_usage(pl, "build", "build the renumbering table on the fly, instead of loading it (requires --lpbs)");
  param_list_decl_usage(pl, "quiet", "do not print the renumbering table contents");
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
    int check = 0;
    int build = 0;
    int quiet = 0;
    char *argv0 = argv[0];
    cxx_cado_poly cpoly;

    cxx_param_list pl;
    declare_usage(pl);
    param_list_configure_switch(pl, "check", &check);
    param_list_configure_switch(pl, "build", &build);
    param_list_configure_switch(pl, "quiet", &quiet);

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
        param_list_parse_uint_list(pl, "lpbs", &lpb[0], tab.get_nb_polys(), ",");
        tab.set_lpb(lpb);
        tab.build();
    } else {
        tab.read_from_file(renumberfilename);
    }

    if (!quiet) {
        for (index_t i = 0; i < tab.get_size() ; i++) {
            std::string s = tab.debug_data(i);
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
    return EXIT_SUCCESS;
}

