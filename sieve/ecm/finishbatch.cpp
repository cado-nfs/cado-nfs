#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* for strcmp() */
#include <pthread.h>
#include "portability.h"
#include "utils.h"
#include "batch.h"


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
    param_list_decl_usage(pl, "doecm", "(switch) finish with ECM [default = no]");

    verbose_decl_usage(pl);
}

int
main (int argc, char *argv[])
{
  param_list pl;
  cado_poly cpoly;
  char *argv0 = argv[0];
  unsigned long nb_threads = 1;
  int doecm = 0;

  param_list_init(pl);
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
  cado_poly_init(cpoly);
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
  
  unsigned long lim[2] = {0, 0};
  param_list_parse_ulong(pl, "lim0", &lim[0]);
  param_list_parse_ulong(pl, "lim1", &lim[1]);

  int lpb[2] = {0, 0};
  param_list_parse_int(pl, "lpb0", &lpb[0]);
  param_list_parse_int(pl, "lpb1", &lpb[1]);

  if (lpb[0] * lpb[1] * lim[0] * lim[1] == 0) {
      fprintf(stderr, "Error: parameters lpb[01] and lim[01] are mandatory\n");
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

  mpz_t batchP[2];
  mpz_init (batchP[0]);
  mpz_init (batchP[1]);

  for (int side = 0; side < 2; ++side) {
      create_batch_file (batch_file[side], batchP[side], lim[side],
              1UL << batchlpb[side], cpoly->pols[side], stdout, nb_threads);
  }

  // Read list from the input file.
  cofac_list List;
  cofac_list_init(List);
  FILE * inp = fopen_maybe_compressed(infilename, "r");
  ASSERT_ALWAYS(inp);

#define MAX_SIZE 2048
  char str[MAX_SIZE];
  mpz_t A, R, q;
  mpz_init(A);
  mpz_init(R);
  mpz_init(q);
  mpz_set_ui(q, 1);
  long a;
  unsigned long b;
  while (fgets(str, MAX_SIZE, inp)) {
      if (str[0] == '#') continue;
      gmp_sscanf(str, "%ld %lu %Zd %Zd\n", &a, &b, R, A);
      cofac_list_add(List, a, b, R, A, 0, q);
  }
  fclose_maybe_compressed(inp, infilename);
  mpz_clear(A);
  mpz_clear(R);
  mpz_clear(q);

  mpz_t B[2], L[2], M[2];
  for(int side = 0 ; side < 2 ; side++) {
      mpz_init(B[side]);
      mpz_init(L[side]);
      mpz_init(M[side]);
      mpz_ui_pow_ui(B[side], 2, batchlpb[side]);
      mpz_ui_pow_ui(L[side], 2, lpb[side]);
      mpz_ui_pow_ui(M[side], 2, batchmfb[side]);
  }
  unsigned long nrels = find_smooth(List, batchP, B, L, M, stdout, nb_threads);
  
  if (doecm) {
      factor(List, nrels, cpoly, batchlpb, lpb, stdout, nb_threads);
  } else {
      for (unsigned long i = 0; i < nrels; i++) {
          uint32_t *perm = List->perm;
          gmp_printf("%" PRIi64 " %" PRIu64 " %Zd %Zd\n",
                  List->a[perm[i]], List->b[perm[i]],
                  List->R[perm[i]], List->A[perm[i]]);
      }
  }

  for(int side = 0 ; side < 2 ; side++) {
      mpz_clear(B[side]);
      mpz_clear(L[side]);
      mpz_clear(M[side]);
  }
  cofac_list_clear(List);

  cado_poly_clear(cpoly);
  mpz_clear (batchP[0]);
  mpz_clear (batchP[1]);
  param_list_clear(pl);

  return EXIT_SUCCESS;
}
