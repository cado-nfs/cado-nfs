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


static void declare_usage(param_list pl)
{
    param_list_decl_usage(pl, "poly", "poly file");
    param_list_decl_usage(pl, "t", "number of threads");
    param_list_decl_usage(pl, "lim0", "sieved bound on side 0");
    param_list_decl_usage(pl, "lim1", "sieved bound on side 1");
    param_list_decl_usage(pl, "batchlpb0", "batch bound on side 0");
    param_list_decl_usage(pl, "batchlpb1", "batch bound on side 1");
    param_list_decl_usage(pl, "batch0", "file of product of primes on side 0");
    param_list_decl_usage(pl, "batch1", "file of product of primes on side 1");

    verbose_decl_usage(pl);
}

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
  
  unsigned long lim[2] = {ULONG_MAX, ULONG_MAX};
  param_list_parse_ulong(pl, "lim0", &lim[0]);
  param_list_parse_ulong(pl, "lim1", &lim[1]);
  if (lim[0] == ULONG_MAX || lim[1] == ULONG_MAX) {
      fprintf(stderr,
              "Error: parameters lim[01] are mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  int batchlpb[2] = {0, 0};
  param_list_parse_int(pl, "batchlpb0", &batchlpb[0]);
  param_list_parse_int(pl, "batchlpb1", &batchlpb[1]);

  if (batchlpb[0] * batchlpb[1] == 0) {
      fprintf(stderr,
              "Error: parameters batchlpb[01] are mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  const char *batch_file[2];
  batch_file[0] = param_list_lookup_string (pl, "batch0");
  batch_file[1] = param_list_lookup_string (pl, "batch1");
  if (batch_file[0] == NULL || batch_file[1] == NULL) {
      fprintf(stderr, "Error: parameters batch[01] are mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  cxx_mpz batchP[2];

  double extra_time = 0;
  for (int side = 0; side < 2; ++side) {
      create_batch_file (batch_file[side], batchP[side], lim[side],
              1UL << batchlpb[side], cpoly->pols[side], stdout, nb_threads, extra_time);
  }

  return EXIT_SUCCESS;
}
