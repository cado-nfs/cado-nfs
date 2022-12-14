#include "cado.h" // IWYU pragma: keep
#include <cstdio>
#include <cstdlib>
#include <fcntl.h>   /* for _O_BINARY */

#include "renumber.hpp" // renumber_t
#include "params.h" // param_list
#include "cado_poly.h"  // cado_poly
#include "gzip.h"       // fopen_maybe_compressed
#include "verbose.h"   // verbose_interpret_parameters
#include "macros.h"

char *argv0; /* = argv[0] */

static void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "in", "input file; if not given, read from stdin");
  param_list_decl_usage(pl, "renumber", "renumber table (mandatory)");
  param_list_decl_usage(pl, "inverse", "convert index to ideal");
  param_list_decl_usage(pl, "poly", "polynomial file (required for -inverse)");
  verbose_decl_usage(pl);
}

static void
usage (param_list pl, char *argv0)
{
  fprintf (stderr, "Convert ideal to index in the renumbering table\n");
  fprintf (stderr, "Ideals are read from file (or stdin), each line should have"
                   " the following format:\n    p r side\nwhere p and r are in "
                   "hexa and side is 0 or 1. For a rational side, r can be "
                   "anything.\n");
  fprintf (stderr, "If -inverse parameter is given, each line should have the "
                   "following format:\n    id\nwhere id is in hexa.\n\n");
  param_list_print_usage(pl, argv0, stderr);
  exit(EXIT_FAILURE);
}

int
main (int argc, char *argv[])
{
  argv0 = argv[0];
  cxx_cado_poly cpoly;
  int inverse = 0;

  cxx_param_list pl;
  declare_usage(pl);
  param_list_configure_switch(pl, "inverse", &inverse);
  argv++,argc--;

#ifdef HAVE_MINGW
  _fmode = _O_BINARY;     /* Binary open for all files */
#endif

  if (argc == 0)
    usage (pl, argv0);

  for( ; argc ; )
  {
    if (param_list_update_cmdline(pl, &argc, &argv))
      continue;
    fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
    usage (pl, argv0);
  }
  /* print command-line arguments */
  verbose_interpret_parameters(pl);
  param_list_print_command_line (stdout, pl);
  fflush(stdout);

  const char *polyfilename = param_list_lookup_string(pl, "poly");
  const char *renumberfilename = param_list_lookup_string(pl, "renumber");
  const char * infilename = param_list_lookup_string(pl, "in");

  if (param_list_warn_unused(pl))
  {
    fprintf(stderr, "Error, unused parameters are given\n");
    usage(pl, argv0);
  }

  if (renumberfilename == NULL)
  {
    fprintf (stderr, "Error, missing -renumber command line argument\n");
    usage(pl, argv0);
  }
  if (inverse && polyfilename == NULL)
  {
    fprintf (stderr, "Error, missing -poly command line argument\n");
    usage(pl, argv0);
  }

  if (inverse)
  {
    if (!cado_poly_read (poly, polyfilename))
    {
      fprintf (stderr, "Error reading polynomial file\n");
      exit (EXIT_FAILURE);
    }
  }

  FILE *infile = NULL;
  if (infilename == NULL)
    infile = stdin;
  else
  {
    infile = fopen_maybe_compressed (infilename, "r");
  }
  ASSERT_ALWAYS (infile != NULL);

  renumber_t renumber_table(renumberfilename);

  if (inverse) {
      index_t i;
      while (fscanf (infile, "%" SCNid "\n", &i) == 1)
      {
          renumber_t::p_r_side x = renumber_table.p_r_from_index(i);
          printf("CONVERT: index = %" PRid "   ===>   p = %" PRpr " r = %" PRpr
                  " side = %d\n", i, x.p, x.r, x.side);
      }
  } else {
      p_r_values_t p, r;
      int side;
      while (fscanf (infile, "%" SCNpr " %" SCNpr " %d\n", &p, &r, &side) == 3)
      {
          index_t i = renumber_table.index_from_p_r (p, r, side);
          printf("CONVERT: p = %" PRpr " r = %" PRpr " side = %d   ===>  "
                  "index = %" PRid "\n", p, r, side, i);
      }
  }
  ASSERT_ALWAYS (feof(infile));

  if (infilename != NULL)
    fclose_maybe_compressed (infile, infilename); 

  return 0;
}
