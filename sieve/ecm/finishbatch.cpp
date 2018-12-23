#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* for strcmp() */
#include <pthread.h>
#include <sstream>
#include "portability.h"
#include "utils.h"
#include "memusage.h"
#include "batch.hpp"


static void declare_usage(param_list pl)
{
    param_list_decl_usage(pl, "in", "input file");
    param_list_decl_usage(pl, "poly", "poly file");
    param_list_decl_usage(pl, "t", "number of threads");
    param_list_decl_usage(pl, "lim0", "sieved bound on side 0");
    param_list_decl_usage(pl, "lim1", "sieved bound on side 1");
    param_list_decl_usage(pl, "lpb0", "large prime bound on side 0");
    param_list_decl_usage(pl, "lpb1", "large prime bound on side 1");
    param_list_decl_usage(pl, "batchmfb0", "threshold after batch on side 0");
    param_list_decl_usage(pl, "batchmfb1", "threshold after batch on side 1");
    param_list_decl_usage(pl, "batchlpb0", "batch bound on side 0");
    param_list_decl_usage(pl, "batchlpb1", "batch bound on side 1");
    param_list_decl_usage(pl, "batch0", "file of product of primes on side 0");
    param_list_decl_usage(pl, "batch1", "file of product of primes on side 1");
    param_list_decl_usage(pl, "doecm", "finish with ECM [default = no]");

    verbose_decl_usage(pl);
}

int
main (int argc, char *argv[])
{
  cxx_param_list pl;
  cxx_cado_poly cpoly;
  char *argv0 = argv[0];
  unsigned long nb_threads = 1;
  int doecm = 0;

  declare_usage(pl);

  param_list_configure_switch(pl, "-doecm", &doecm);

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

  const char * infilename;
  if ((infilename = param_list_lookup_string(pl, "in")) == NULL) {
      fprintf(stderr, "Error: parameter -in is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }
  param_list_parse_ulong(pl, "t"   , &nb_threads);
  
  unsigned long lim[2] = {ULONG_MAX, ULONG_MAX};
  param_list_parse_ulong(pl, "lim0", &lim[0]);
  param_list_parse_ulong(pl, "lim1", &lim[1]);
  if (lim[0] == ULONG_MAX || lim[1] == ULONG_MAX) {
      fprintf(stderr, "Error: parameters lim[01] are mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  int lpb[2] = {0, 0};
  param_list_parse_int(pl, "lpb0", &lpb[0]);
  param_list_parse_int(pl, "lpb1", &lpb[1]);

  if (lpb[0] * lpb[1] == 0) {
      fprintf(stderr, "Error: parameters lpb[01] are mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }


  const char *batch_file[2];
  int batchlpb[2] = {0, 0};
  int batchmfb[2] = {0, 0};
  batch_file[0] = param_list_lookup_string (pl, "batch0");
  batch_file[1] = param_list_lookup_string (pl, "batch1");
  param_list_parse_int(pl, "batchlpb0", &(batchlpb[0]));
  param_list_parse_int(pl, "batchlpb1", &(batchlpb[1]));
  param_list_parse_int(pl, "batchmfb0", &(batchmfb[0]));
  param_list_parse_int(pl, "batchmfb1", &(batchmfb[1]));
  if (batchlpb[0] * batchlpb[1] * batchmfb[0] * batchmfb[1] == 0) {
      fprintf(stderr, "Error: parameters batchlpb[01] and batchmfb[01] are mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  std::array<cxx_mpz, 2> batchP;

  double extra_time = 0;

  for (int side = 0; side < 2; ++side) {
      create_batch_file (batch_file[side], batchP[side], lim[side],
              1UL << batchlpb[side], cpoly->pols[side], stdout, nb_threads,
              extra_time);
  }

  // Read list from the input file.
  cofac_list List;
  FILE * inp = fopen_maybe_compressed(infilename, "r");
  ASSERT_ALWAYS(inp);

#define MAX_SIZE 2048
  char str[MAX_SIZE];
  std::array<cxx_mpz, 2> norms;
  cxx_mpz q;
  mpz_set_ui(q, 1);
  // Create a fake special-q
  std::vector<uint64_t> empty;
  las_todo_entry doing(q, q, 0, empty);
  long a;
  unsigned long b;
  while (fgets(str, MAX_SIZE, inp)) {
      if (str[0] == '#') continue;
      gmp_sscanf(str, "%ld %lu %Zd %Zd\n", &a, &b, (mpz_ptr) norms[0], (mpz_ptr) norms[1]);
      List.emplace_back(a, b, norms, &doing);
  }
  fclose_maybe_compressed(inp, infilename);

  find_smooth(List, batchP, batchlpb, lpb, batchmfb, stdout, nb_threads, extra_time);
  
  if (doecm) {
      std::list<relation> smooth = factor(List, cpoly, batchlpb, lpb, stdout, nb_threads, extra_time);
        for(auto const & rel : smooth) {
            std::ostringstream os;
            os << rel << "\n";
            printf("%s\n", os.str().c_str());
        }
  } else {
      for (auto const & x : List) {
          gmp_printf("%" PRIi64 " %" PRIu64 " %Zd %Zd\n",
                  x.a, x.b,
                  (mpz_srcptr) x.cofactor[0],
                  (mpz_srcptr) x.cofactor[1]);
      }
  }

    const long peakmem = PeakMemusage();
    if (peakmem > 0)
        printf ("# PeakMemusage (MB) = %ld \n",
                peakmem >> 10);

  return EXIT_SUCCESS;
}
