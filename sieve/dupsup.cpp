#include "cado.h" // IWYU pragma: keep

#include <cstdlib>
#include <cstdio>
#include <climits>

#include <sstream>

#include "fmt/base.h"

#include "gzip.h"
#include "las-duplicate.hpp"
#include "las-info.hpp"
#include "las-siever-config.hpp"
#include "las-side-config.hpp"
#include "special-q.hpp"
#include "params.hpp"
#include "relation.hpp"
#include "sieve-methods.hpp"
#include "verbose.hpp"
#include "utils_cxx.hpp"

static constexpr bool reprint_special_q_at_each_rel = false;

static void *
dupsup (FILE *output, relation & rel, special_q const& doing, const int is_dupe)
{
  if (reprint_special_q_at_each_rel)
      fmt::print("# {}\n", doing);
  rel.print(output, is_dupe ? "# DUPE " : "");
  return nullptr;
}

/* If the line is a special-q comment, sets sq and rho and returns 1.
   Otherwise returns 0. */
static bool
read_sq_comment(special_q & doing, const char *line)
{
    std::istringstream is(line);
    if (!(is >> expect("# Sieving ")))
        return false;
    is >> doing;
    return is.good();
}

/* FIXME: we should really have static class functions like
 * las_info::declare_usage
 *
 * Well, except that indeed some parameters are used in general by
 * las_info, and not in this case.
 */
static void declare_usage(cxx_param_list & pl)
{
  pl.declare_usage("path_antebuffer", "path to antebuffer program");
  pl.declare_usage("poly", "polynomial file");
  pl.declare_usage("I",    "set sieving region to 2^I times J");
  pl.declare_usage("A",    "set sieving region to 2^A");
  pl.declare_usage("skew", "skewness");
  pl.declare_usage("dup-qmin", "lower limit of global q-range for 2-sided duplicate removal");
  pl.declare_usage("dup-qmax", "upper limit of global q-range for 2-sided duplicate removal");
  pl.declare_usage("sqside", "side of special-q (default=1)");
  /* those are typical from las invocations, we wish to keep them
   * accepted */
  pl.declare_usage("out",  "filename where relations are written, instead of stdout");
  siever_side_config::declare_usage(pl);
  pl.declare_usage("fb",   "(unused)");
  pl.declare_usage("fbc",  "(unused)");
  pl.declare_usage("q0",   "(unused)");
  pl.declare_usage("q1",   "(unused)");
  pl.declare_usage("nq",   "(unused)");
  pl.declare_usage("v",    "(unused)");
  verbose_decl_usage(pl);
  las_info::declare_usage(pl);
}

// coverity[root_function]
int
main (int argc, char const * argv[])
{
    cxx_param_list pl;
    declare_usage(pl);

    pl.configure_switch_old("-v", nullptr);

    pl.process_command_line(argc, argv, true);

    verbose_interpret_parameters(pl);
    pl.print_command_line(stdout);
    fflush(stdout);

    if (argc == 0)
      pl.fail("Error, provide freeform file names\n");

    pl.lookup("fb0");
    pl.lookup("fb1");
    pl.lookup("fbc");
    pl.lookup("q0");
    pl.lookup("q1");
    pl.lookup("nq");
    const char * outputname = pl.lookup_old("out");

    las_info las(pl, NFS{});

    las.prepare_sieve_shared_data(pl);

    if (pl.warn_unused())
      pl.fail("Unused parameters are given\n");

    FILE * output = stdout;
    if (outputname) {
	if (!(output = fopen_maybe_compressed(outputname, "w"))) {
	    fprintf(stderr, "Could not open %s for writing\n", outputname);
	    exit(EXIT_FAILURE);
	}
    }

    setvbuf(output, nullptr, _IOLBF, 0);

    special_q doing;

    for (int argi = 0; argi < argc; argi++) {
      FILE *f = fopen_maybe_compressed(argv[argi], "rb");
      if (f == nullptr) {
          perror(argv[argi]);
          abort();
      }
      for (; !feof(f) ;) {
        char line[1024];
        if (fgets(line, sizeof(line), f) == nullptr)
          break;

        if (read_sq_comment(doing, line)) {
            /* If qmin is not given, use lim on the (current) special-q
             * side by default.  This makes sense only if the relevant
             * fields have been filled from the command line.
             *
             * This is a kludge, really. If we have special-q's on two
             * sides, the only reliable way to go is to provide both
             * dup-qmin arguments. The default below kinda makes sense as
             * a shortcut for the special-q-on-only-one-side setting.
             *
             * Note that in particular, this mandates the use of the
             * -sync argument.
             */
            if (las.dupqmin[doing.side] == ULONG_MAX)
                las.dupqmin[doing.side] = las.config_pool.base.sides[doing.side].lim;

            continue;
        }

        relation rel;
        if (rel.parse(line)) {
            int const is_dupe = relation_is_duplicate(rel, doing, las);
            dupsup(output, rel, doing, is_dupe);
        }
      }
      fclose_maybe_compressed(f, argv[argi]);
    }
    
    if (outputname)
        fclose_maybe_compressed(output, outputname);

    return 0;
}
