#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* for strcmp() */
#include <pthread.h>
#include "portability.h"
#include "utils.h"
#include "batch.h"
#include "facul.hpp"

static void declare_usage(param_list pl)
{
    param_list_decl_usage(pl, "in", "input file");
    param_list_decl_usage(pl, "lpb0", "large prime bound on side 0");
    param_list_decl_usage(pl, "lpb1", "large prime bound on side 1");
    param_list_decl_usage(pl, "batchlpb0", "batch bound on side 0");
    param_list_decl_usage(pl, "batchlpb1", "batch bound on side 1");

    verbose_decl_usage(pl);
}

int
main (int argc, char *argv[])
{
  param_list pl;
  char *argv0 = argv[0];
  double st, wct;
  st = seconds();
  wct = wct_seconds();

  param_list_init(pl);
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

    const char * infilename;
  if ((infilename = param_list_lookup_string(pl, "in")) == NULL) {
      fprintf(stderr, "Error: parameter -in is mandatory\n");
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

  int batchlpb[2] = {0, 0};
  param_list_parse_int(pl, "batchlpb0", &(batchlpb[0]));
  param_list_parse_int(pl, "batchlpb1", &(batchlpb[1]));
  if (batchlpb[0] * batchlpb[1] == 0) {
      fprintf(stderr, "Error: parameters batchlpb[01] are mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  unsigned long B[2];
  B[0] = 1UL<<batchlpb[0];
  B[1] = 1UL<<batchlpb[1];

  facul_strategy_t* strategy0 = facul_make_strategy(
          B[0], lpb[0], 4*lpb[0],
          30, 0);
  strategy0->assume_prime_thresh = 0.0;
  facul_strategy_t* strategy1 = facul_make_strategy(
          B[1], lpb[1], 4*lpb[1],
          30, 0);
  strategy1->assume_prime_thresh = 0.0;

  // Read list from the input file.
  FILE * inp = fopen_maybe_compressed(infilename, "r");
  ASSERT_ALWAYS(inp);

#define MAX_SIZE 1024
  char str[MAX_SIZE];
  cxx_mpz A, R;
  long a;
  unsigned long b;
  unsigned long nrels = 0;
  while (fgets(str, MAX_SIZE, inp)) {
      if (str[0] == '#') continue;
      gmp_sscanf(str, "%ld %lu %Zd %Zd\n", &a, &b, mpz_ptr(R), mpz_ptr(A));
      std::vector<cxx_mpz> factorsR, factorsA;
  
      if (mpz_sizeinbase(R, 2) > (size_t)lpb[0]) {
          int found = facul(factorsR, R, strategy0);
          if (found <= 0)
              continue;
      }
      if (mpz_sizeinbase(A, 2) > (size_t)lpb[1]) {
          int found = facul(factorsA, A, strategy1);
          if (found <= 0)
              continue;
      }
      nrels++;
      gmp_printf("%ld, %lu\n", a, b);
  }
  printf("# Finish ECM found %lu rels in %.2f s (wct %2.f s)\n",
          nrels, seconds()-st, wct_seconds()-wct);
  
  fclose_maybe_compressed(inp, infilename);
  facul_clear_strategy(strategy0);
  facul_clear_strategy(strategy1);
  param_list_clear(pl);

  return EXIT_SUCCESS;
}
