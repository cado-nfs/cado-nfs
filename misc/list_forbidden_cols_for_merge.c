#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>   /* for _O_BINARY */

#include "portability.h"
#include "utils_with_io.h"

char *argv0; /* = argv[0] */

renumber_t renumber_tab;
cado_poly poly;

static void *
compute_side_info (void * context_data, earlyparsed_relation_ptr rel)
{
  uint64_t *cols_side_info = (uint64_t *) context_data;

  uint64_t nonvoidside = 0; /* bit vector of which sides appear in the rel */

  for (unsigned int i = 0; i < rel->nb; i++)
  {
    int side = renumber_get_side_from_index (renumber_tab, rel->primes[i].h, poly);
    nonvoidside |= ((uint64_t) 1) << side;
  }

  for (unsigned int i = 0; i < rel->nb; i++)
    cols_side_info[rel->primes[i].h] |= nonvoidside;

  return NULL;
}

static void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "purgedfile", "file with relations from purge");
  param_list_decl_usage(pl, "renumber", "input file for renumbering table");
  param_list_decl_usage(pl, "poly", "input polynomial file");
  param_list_decl_usage(pl, "outfile", "output file");
  param_list_decl_usage(pl, "force-posix-threads", "(switch)");
  param_list_decl_usage(pl, "path_antebuffer", "path to antebuffer program");
  verbose_decl_usage(pl);
}

static void
usage (param_list pl, char *argv0)
{
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
}

int
main (int argc, char *argv[])
{
  argv0 = argv[0];

  param_list pl;
  param_list_init(pl);
  declare_usage(pl);
  argv++,argc--;

  param_list_configure_switch(pl, "force-posix-threads", &filter_rels_force_posix_threads);

#ifdef HAVE_MINGW
  _fmode = _O_BINARY;     /* Binary open for all files */
#endif

  if (argc == 0)
    usage (pl, argv0);

  for( ; argc ; )
  {
    if (param_list_update_cmdline(pl, &argc, &argv))
      continue;
    fprintf (stderr, "Unknown option: %s\n", argv[0]);
    abort();
  }
  /* print command-line arguments */
  verbose_interpret_parameters(pl);
  param_list_print_command_line (stdout, pl);
  fflush(stdout);

  const char * purgedfilename = param_list_lookup_string(pl, "purgedfile");
  const char * outfilename = param_list_lookup_string(pl, "outfile");
  const char * renumberfilename = param_list_lookup_string(pl, "renumber");
  const char * polyfilename = param_list_lookup_string(pl, "poly");
  const char * path_antebuffer = param_list_lookup_string(pl, "path_antebuffer");

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
  if (purgedfilename == NULL)
  {
    fprintf (stderr, "Error, missing -purgedfile command line argument\n");
    usage(pl, argv0);
  }
  if (outfilename == NULL)
  {
    fprintf (stderr, "Error, missing -outfile command line argument\n");
    usage(pl, argv0);
  }
  if (polyfilename == NULL)
  {
    fprintf (stderr, "Error, missing -poly command line argument\n");
    usage (pl, argv0);
  }

  set_antebuffer_path (argv0, path_antebuffer);

  cado_poly_init(poly);
  if (!cado_poly_read (poly, polyfilename))
  {
    fprintf (stderr, "Error reading polynomial file\n");
    exit (EXIT_FAILURE);
  }

  renumber_init_for_reading (renumber_tab);
  renumber_read_table (renumber_tab, renumberfilename);
  

  // malloc and memset to 0 info on all cols
  uint64_t *cols_side_info = NULL;
  cols_side_info = (uint64_t *) malloc (renumber_tab->size * sizeof (uint64_t));
  ASSERT_ALWAYS (cols_side_info != NULL);
  memset (cols_side_info, 0, renumber_tab->size * sizeof (uint64_t));

  char * filelist[2] = { (char *) purgedfilename, NULL};

  filter_rels(filelist, (filter_rels_callback_t) &compute_side_info,
          cols_side_info,
          EARLYPARSE_NEED_AB_HEXA | EARLYPARSE_NEED_INDEX, NULL, NULL);

  FILE *outfile = NULL;
  outfile = fopen_maybe_compressed (outfilename, "w");
  ASSERT_ALWAYS (outfile != NULL);

  for (index_t i = 0; i < renumber_tab->size; i++)
  {
    if (__builtin_popcountl(cols_side_info[i]) > 2)
      fprintf (outfile, "%" PRid "\n", i);
  }

  fclose_maybe_compressed (outfile, outfilename);

  free (cols_side_info);

  param_list_clear(pl);
  cado_poly_clear (poly);
  renumber_clear (renumber_tab);
  return 0;
}
