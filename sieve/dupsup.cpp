#include "cado.h" // IWYU pragma: keep

#include <array>
#include <cstdlib>
#include <cstdio>
#include <climits> /* for CHAR_BIT */
#ifdef HAVE_MINGW
#include <fcntl.h>   /* for _O_BINARY */
#endif
#include <gmp.h>
#include "cxx_mpz.hpp"
#include "gzip.h"
#include "las-duplicate.hpp"
#include "las-info.hpp"
#include "las-siever-config.hpp"
#include "las-todo-entry.hpp"
#include "params.h"
#include "relation.hpp"
#include "verbose.h"


static void *
dupsup (FILE *output, relation & rel, las_todo_entry const& doing, const int is_dupe)
{
  if (0)
    gmp_fprintf (output, "# sq = %Zd, rho = %Zd, side = %d\n",
            (mpz_srcptr) doing.p,
            (mpz_srcptr) doing.r, doing.side);
  rel.print(output, is_dupe ? "# DUPE " : "");
  return NULL;
}

/* If the line is a special-q comment, sets sq and rho and returns 1.
   Otherwise returns 0. */
static int
read_sq_comment(las_todo_entry & doing, const char *line)
{
    int side;
    cxx_mpz p,r;
  if (gmp_sscanf(line, "# Sieving side-%d q=%Zd; rho=%Zd;",
              &side,
              (mpz_ptr) p,
              (mpz_ptr) r) == 3) {
      /* this is the way to go if we want proper initialization of all
       * fields */
      doing = las_todo_entry(p, r, side);
    return 1;
  }
  return 0;
}

/* FIXME: we should really have static class functions like
 * las_info::declare_usage
 *
 * Well, except that indeed some parameters are used in general by
 * las_info, and not in this case.
 */
static void declare_usage(cxx_param_list & pl)
{
  param_list_decl_usage(pl, "path_antebuffer", "path to antebuffer program");
  param_list_decl_usage(pl, "poly", "polynomial file");
  param_list_decl_usage(pl, "I",    "set sieving region to 2^I times J");
  param_list_decl_usage(pl, "A",    "set sieving region to 2^A");
  param_list_decl_usage(pl, "skew", "skewness");
  param_list_decl_usage(pl, "dup-qmin", "lower limit of global q-range for 2-sided duplicate removal");
  param_list_decl_usage(pl, "dup-qmax", "upper limit of global q-range for 2-sided duplicate removal");
  param_list_decl_usage(pl, "sqside", "side of special-q (default=1)");
  /* those are typical from las invocations, we wish to keep them
   * accepted */
  param_list_decl_usage(pl, "out",  "filename where relations are written, instead of stdout");
  siever_side_config::declare_usage(pl);
  param_list_decl_usage(pl, "fb",   "(unused)");
  param_list_decl_usage(pl, "fbc",  "(unused)");
  param_list_decl_usage(pl, "q0",   "(unused)");
  param_list_decl_usage(pl, "q1",   "(unused)");
  param_list_decl_usage(pl, "nq",   "(unused)");
  param_list_decl_usage(pl, "v",    "(unused)");
  verbose_decl_usage(pl);
  las_info::declare_usage(pl);
}

static void
usage (param_list_ptr pl, char const * argv0)
{
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
}

// coverity[root_function]
int
main (int argc, char const * argv[])
{
    const char * argv0 = argv[0];

    cxx_param_list pl;
    declare_usage(pl);
    argv++,argc--;

#ifdef HAVE_MINGW
    _fmode = _O_BINARY;     /* Binary open for all files */
#endif

    if (argc == 0)
      usage (pl, argv0);

    param_list_configure_switch(pl, "-v", NULL);

    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        /* Since we accept file names freeform, we decide to never abort
         * on unrecognized options */
        break;
        // fprintf (stderr, "Unknown option: %s\n", argv[0]);
        // abort();
    }
    /* print command-line arguments */
    verbose_interpret_parameters(pl);
    param_list_print_command_line (stdout, pl);
    fflush(stdout);

    if (argc == 0) {
      fprintf(stderr, "Error, provide freeform file names\n");
      usage(pl, argv0);
    }

    param_list_lookup_string(pl, "fb0");
    param_list_lookup_string(pl, "fb1");
    param_list_lookup_string(pl, "fbc");
    param_list_lookup_string(pl, "q0");
    param_list_lookup_string(pl, "q1");
    param_list_lookup_string(pl, "nq");
    const char * outputname = param_list_lookup_string(pl, "out");

    cxx_mpz sq, rho;

    las_info las(pl);

    las.prepare_sieve_shared_data(pl);

    if (param_list_warn_unused(pl))
    {
      fprintf(stderr, "Error, unused parameters are given\n");
      usage(pl, argv0);
    }

    FILE * output = stdout;
    if (outputname) {
	if (!(output = fopen_maybe_compressed(outputname, "w"))) {
	    fprintf(stderr, "Could not open %s for writing\n", outputname);
	    exit(EXIT_FAILURE);
	}
    }

    setvbuf(output, NULL, _IOLBF, 0);      /* mingw has no setlinebuf */

    las_todo_entry doing;

    for (int argi = 0; argi < argc; argi++) {
      FILE *f = fopen_maybe_compressed(argv[argi], "rb");
      if (f == NULL) {
          perror(argv[argi]);
          abort();
      }
      for (; !feof(f) ;) {
        char line[1024];
        if (fgets(line, sizeof(line), f) == NULL)
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
