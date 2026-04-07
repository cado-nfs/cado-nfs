#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>
#include <vector>

#include <gmp.h>

#include "macros.h"
#include "params.h"
#include "timing.h"
#include "verbose.hpp"
#include "gzip.h"
#include "cxx_mpz.hpp"

#include "facul.hpp"
#include "facul_strategies.hpp"

static void declare_usage(cxx_param_list & pl)
{
    pl.declare_usage("in", "input file");
    pl.declare_usage("lpb0", "large prime bound on side 0");
    pl.declare_usage("lpb1", "large prime bound on side 1");
    pl.declare_usage("batchlpb0", "batch bound on side 0");
    pl.declare_usage("batchlpb1", "batch bound on side 1");

    verbose_decl_usage(pl);
}

// coverity[root_function]
int
main (int argc, char const *argv[])
{
  cxx_param_list pl;
  double st, wct;
  st = seconds();
  wct = wct_seconds();

  declare_usage(pl);

  param_list_process_command_line_and_extra_parameter_files(pl, &argc, &argv);

  verbose_interpret_parameters(pl);
  param_list_print_command_line(stdout, pl);

    const char * infilename;
  if ((infilename = param_list_lookup_string(pl, "in")) == nullptr)
      pl.fail("Error: parameter -in is mandatory\n");
  
  int lpb[2] = {0, 0};
  param_list_parse_int(pl, "lpb0", &lpb[0]);
  param_list_parse_int(pl, "lpb1", &lpb[1]);

  if (lpb[0] * lpb[1] == 0)
      pl.fail("Error: parameters lpb[01] are mandatory\n");

  int batchlpb[2] = {0, 0};
  param_list_parse_int(pl, "batchlpb0", &(batchlpb[0]));
  param_list_parse_int(pl, "batchlpb1", &(batchlpb[1]));
  if (batchlpb[0] * batchlpb[1] == 0)
      pl.fail("Error: parameters batchlpb[01] are mandatory\n");

  unsigned long B[2];
  B[0] = 1UL<<batchlpb[0];
  B[1] = 1UL<<batchlpb[1];

  /* keep it simple: use 30 curves in both cases */
  facul_strategy_oneside strategy0(B[0], lpb[0], 4*lpb[0], 30, 0);
  strategy0.BB = 0.0;

  facul_strategy_oneside strategy1(B[1], lpb[1], 4*lpb[1], 30, 0);
  strategy1.BB = 0.0;

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
  
      if (mpz_sizeinbase(R, 2) > (size_t)lpb[0]) {
          if (facul(R, strategy0).status != FACUL_SMOOTH)
              continue;
      }
      if (mpz_sizeinbase(A, 2) > (size_t)lpb[1]) {
          if (facul(A, strategy1).status != FACUL_SMOOTH)
              continue;
      }
      nrels++;
      gmp_printf("%ld, %lu\n", a, b);
  }
  printf("# Finish ECM found %lu rels in %.2f s (wct %.2f s)\n",
          nrels, seconds()-st, wct_seconds()-wct);
  
  fclose_maybe_compressed(inp, infilename);

  return EXIT_SUCCESS;
}
