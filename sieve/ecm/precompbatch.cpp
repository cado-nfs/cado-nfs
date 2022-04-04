#include "cado.h" // IWYU pragma: keep

#include <climits>     // for ULONG_MAX
#include <cstdio>      // for fprintf, stderr, NULL, fclose, fopen, stdout
#include <cstdlib>     // for exit, EXIT_FAILURE, EXIT_SUCCESS
#include "batch.hpp"    // for create_batch_file
#include "cado_poly.h"
#include "cxx_mpz.hpp"
#include "mpz_poly.h"
#include "verbose.h"    // verbose_decl_usage
#include "params.h"
#include "las-side-config.hpp"


static void declare_usage(cxx_param_list & pl)
{
    param_list_decl_usage(pl, "poly", "poly file");
    param_list_decl_usage(pl, "t", "number of threads");

    siever_side_config::declare_usage(pl);
    batch_side_config::declare_usage(pl);

    verbose_decl_usage(pl);
}

// coverity[root_function]
int
main (int argc, char *argv[])
{
  cxx_param_list pl;
  cxx_cado_poly cpoly;
  char *argv0 = argv[0];
  unsigned long nb_threads = 1;

  declare_usage(pl);

  argv++, argc--;
  for( ; argc ; ) {
      FILE *f;
      if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }

      /* Could also be a file */
      if ((f = fopen(argv[0], "r")) != NULL) {
          param_list_read_stream(pl, f, 0);
          fclose(f);
          argv++,argc--;
          continue;
      }

      fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
      param_list_print_usage(pl, argv0, stderr);
      exit (EXIT_FAILURE);
  }
  verbose_interpret_parameters(pl);
  param_list_print_command_line(stdout, pl);

  const char * filename;
  if ((filename = param_list_lookup_string(pl, "poly")) == NULL) {
      fprintf(stderr, "Error: parameter -poly is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }
  if (!cado_poly_read(cpoly, filename)) {
      fprintf (stderr, "Error reading polynomial file %s\n", filename);
      exit (EXIT_FAILURE);
  }

  param_list_parse_ulong(pl, "t"   , &nb_threads);

  int nsides = cpoly->nb_polys;

  std::vector<siever_side_config> sides;
  siever_side_config::parse(pl, sides, nsides, { "lim" });

  std::vector<batch_side_config> bsides;
  batch_side_config::parse(pl, bsides, nsides, { "batchlpb", "batchfile" });
  
  cxx_mpz batchP[2];

  double extra_time = 0;
  for (int side = 0; side < nsides; ++side) {
      siever_side_config & S(sides[side]);
      batch_side_config & bS(bsides[side]);
      create_batch_file (bS.batchfilename, batchP[side], S.lim,
              1UL << bS.batchlpb, cpoly->pols[side], stdout, nb_threads, extra_time);
  }

  return EXIT_SUCCESS;
}
