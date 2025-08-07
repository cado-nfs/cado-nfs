#include "cado.h" // IWYU pragma: keep

#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstdint>

#include <list>
#include <sstream>
#include <utility>
#include <vector>

#include <gmp.h>
#include "fmt/base.h"

#include "batch.hpp"
#include "cado_poly.h"
#include "cxx_mpz.hpp"
#include "gzip.h"
#include "las-side-config.hpp"
#include "special-q.hpp"
#include "macros.h"
#include "memusage.h"
#include "mpz_poly.h"
#include "params.h"
#include "relation.hpp"
#include "verbose.h"

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
main (int argc, char const *argv[])
{
  cxx_param_list pl;
  cxx_cado_poly cpoly;
  char const * argv0 = argv[0];
  unsigned long nb_threads = 1;
  int doecm = 0;
  int no_recomp_norm = 0;
  int const ncurves = 50;

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
  
  int const nsides = cpoly->nb_polys;

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
  std::list<std::pair<special_q, std::list<cofac_candidate>>> List;
  FILE * inp = fopen_maybe_compressed(infilename, "r");
  ASSERT_ALWAYS(inp);

#define MAX_SIZE 2048
  char str[MAX_SIZE];

  // If the special-q info is present, we will use it. Otherwise, the
  // fake sq will be used everywhere. This list keeps in memory all the
  // special q encountered.
  // Create a fake special-q
  special_q sq(1, 1, 0);

  long a;
  unsigned long b;
  std::list<cofac_candidate> for_this_q;

  for(int lnum = 0; fgets(str, MAX_SIZE, inp) ; lnum++) {
      cxx_mpz q;
      if (str[0] == '#') {
          cxx_mpz r;
          int side;
          int const ret = gmp_sscanf(str, "# q = (%Zd, %Zd, %d)",
                  &q, &r, &side);
          if (ret == 3) {
              if (!for_this_q.empty())
                  List.emplace_back(sq, std::move(for_this_q));
              for_this_q.clear();
              /* make sure we don't try to factor q */
              const std::vector<uint64_t> primes { uint64_t(q) };
              sq = special_q(q, r, side, primes);
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
      for_this_q.emplace_back(a, b, norms);
  }
  if (!for_this_q.empty())
      List.emplace_back(sq, std::move(for_this_q));
  fclose_maybe_compressed(inp, infilename);

  find_smooth(List, batchP, batchlpb, lpb, batchmfb, stdout, nb_threads, extra_time);
  
  if (doecm) {
      for(auto const & lq : List) {
          fmt::print("# {}\n", lq.first);
          std::list<relation> const smooth = factor(lq.second, cpoly, lq.first, batchlpb, lpb, ncurves, stdout, nb_threads, extra_time, !no_recomp_norm);
          for(auto const & rel : smooth)
              fmt::print("{}\n", rel);
      }
  } else {
      for (auto const & lq : List) {
          for (auto const & x : lq.second) {
              fmt::print("{} {} {}\n",
                      x.a, x.b,
                      join(x.cofactor, " "));
          }
      }
  }

    const size_t peakmem = PeakMemusage();
    if (peakmem > 0)
        printf ("# PeakMemusage (MB) = %zu \n",
                peakmem >> 10);

  return EXIT_SUCCESS;
}
