#include "cado.h" // IWYU pragma: keep
#include <cinttypes>          // for PRIi64, PRIu64
#include <climits>            // for ULONG_MAX
#include <cstdio>             // for fprintf, stderr, printf, stdout, NULL
#include <cstdlib>            // for exit, EXIT_FAILURE, EXIT_SUCCESS
#include <array>               // for array, array<>::value_type
#include <cstdint>             // for uint64_t
#include <list>                // for list
#include <sstream>             // for operator<<, ostringstream, ostream
#include <string>              // for basic_string
#include <vector>              // for vector
#include <gmp.h>               // for gmp_sscanf, mpz_ptr, mpz_srcptr, gmp_p...
#include "batch.hpp"           // for cofac_candidate, cofac_list, create_ba...
#include "cado_poly.h"
#include "cxx_mpz.hpp"
#include "gmp_aux.h"
#include "gzip.h"       // fopen_maybe_compressed
#include "las-todo-entry.hpp"  // for las_todo_entry
#include "macros.h"            // for ASSERT_ALWAYS
#include "memusage.h"   // PeakMemusage
#include "mpz_poly.h"
#include "relation.hpp"        // for operator<<, relation
#include "verbose.h"    // verbose_decl_usage
#include "params.h"
#include "las-side-config.hpp"

static void declare_usage(cxx_param_list & pl)
{
    param_list_decl_usage(pl, "in", "input file");
    param_list_decl_usage(pl, "poly", "poly file");
    param_list_decl_usage(pl, "t", "number of threads");
    siever_side_config::declare_usage(pl);
    batch_side_config::declare_usage(pl);
    param_list_decl_usage(pl, "doecm", "finish with ECM [default = no]");
    param_list_decl_usage(pl, "dont_recomp_norm", "given integers are the norms themselves (w/ sq) [default = no]");
    param_list_decl_usage(pl, "ncurves", "number of curves to be used in ECM [default = 50]");

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
  int doecm = 0;
  int no_recomp_norm = 0;
  int ncurves = 50;

  declare_usage(pl);

  param_list_configure_switch(pl, "-doecm", &doecm);
  param_list_configure_switch(pl, "-dont_recomp_norm", &no_recomp_norm);

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
  
  int nsides = cpoly->nb_polys;

  std::vector<siever_side_config> sides;
  siever_side_config::parse(pl, sides, nsides, { "lpb", "mfb" });

  std::vector<batch_side_config> bsides;
  batch_side_config::parse(pl, bsides, nsides, { "batchlpb", "batchmfb" });

  std::vector<cxx_mpz> batchP;

  double extra_time = 0;

  auto lpb = siever_side_config::collect_lpb(sides);
  auto batchlpb = batch_side_config::collect_batchlpb(bsides);
  auto batchmfb = batch_side_config::collect_batchmfb(bsides);

  for (int side = 0; side < 2; ++side) {
      batchP.push_back(0);
      create_batch_file (bsides[side].batchfilename,
              batchP[side],
              sides[side].lim,
              1UL << bsides[side].batchlpb,
              cpoly->pols[side],
              stdout,
              nb_threads,
              extra_time);
  }

  // Read list from the input file.
  cofac_list List;
  FILE * inp = fopen_maybe_compressed(infilename, "r");
  ASSERT_ALWAYS(inp);

#define MAX_SIZE 2048
  char str[MAX_SIZE];
  cxx_mpz q;
  mpz_set_ui(q, 1);
  // Create a fake special-q
  std::vector<uint64_t> empty;
  las_todo_entry fake_q(q, q, 0, empty);

  // If the special-q info is present, we will use it. Otherwise, the
  // fake sq will be used everywhere. This list keeps in memory all the
  // special q encountered.
  std::list<las_todo_entry> list_q;
  list_q.push_back(fake_q);

  long a;
  unsigned long b;
  for(int lnum = 0; fgets(str, MAX_SIZE, inp) ; lnum++) {
      if (str[0] == '#') {
          cxx_mpz r;
          int side;
          int ret = gmp_sscanf(str, "# q = (%Zd, %Zd, %d)",
                  &q, &r, &side);
          if (ret == 3) {
              std::vector<uint64_t> primes;
              uint64_t qq = mpz_get_uint64(q);
              primes.push_back(qq);
              las_todo_entry this_q(q, r, side, primes);
              list_q.push_back(this_q);
          }
          continue;
      }
      std::vector<cxx_mpz> norms(nsides, 0);
      std::istringstream ss(str);
      ss >> a >> b;
      for(int side = 0 ; side < nsides ; side++)
          ss >> norms[side];
      if (!ss) {
          fprintf(stderr, "parse error at line %d in cofactor input file\n", lnum);
          exit(EXIT_FAILURE);
      }
      List.emplace_back(a, b, norms, &list_q.back());
  }
  fclose_maybe_compressed(inp, infilename);

  find_smooth(List, batchP, batchlpb, lpb, batchmfb, stdout, nb_threads, extra_time);
  
  if (doecm) {
      std::list<relation> smooth = factor(List, cpoly, batchlpb, lpb, ncurves, stdout, nb_threads, extra_time, !no_recomp_norm);
      for(auto const & rel : smooth) {
          std::ostringstream os;
          os << rel << "\n";
          printf("%s", os.str().c_str());
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
