#include "cado.h"
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>  /* for _O_BINARY */
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <pthread.h>
#include <errno.h>
#include <pthread.h>
#ifdef HAVE_LIBGEN_H
#include <libgen.h>
#endif

#include "portability.h"
#include "utils_with_io.h"

cado_poly poly;

FILE *outfiles[NB_POLYS_MAX][NB_POLYS_MAX];
uint64_t nrels[NB_POLYS_MAX][NB_POLYS_MAX];

static void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "filelist", "file containing a list of input files");
  param_list_decl_usage(pl, "basepath", "path added to all file in filelist");
  param_list_decl_usage(pl, "outdir", "output directory");
  param_list_decl_usage(pl, "outprefix", "prefix for the output files");
  param_list_decl_usage(pl, "poly", "input polynomial file");
  param_list_decl_usage(pl, "renumber", "input file for renumbering table");
  verbose_decl_usage(pl);
}

static void
usage (param_list pl, char *argv0)
{
  param_list_print_usage(pl, argv0, stderr);
  exit(EXIT_FAILURE);
}

static void *
thread_split_rel (renumber_t tab, earlyparsed_relation_ptr rel)
{
  uint64_t nonvoidside = 0; /* bit vector of which sides appear in the rel */
  for (weight_t i = 0; i < rel->nb; i++)
  {
    index_t h = rel->primes[i].h;
    rel->primes[i].side = renumber_get_side_from_index (tab, h, poly);
    nonvoidside |= ((uint64_t) 1) << rel->primes[i].side;
  }

  for (int s1 = 0; s1 < poly->nb_polys; s1++)
  {
    for (int s2 = s1+1; s2 < poly->nb_polys; s2++)
    {
      if (nonvoidside & (((uint64_t) 1) << s1) &&
          nonvoidside & (((uint64_t) 1) << s2))
      {
        char buf[1 << 12], *p, *op;
        size_t t;
        unsigned int j; 

        p = d64toa16(buf, rel->a);
        *p++ = ',';
        p = u64toa16(p, rel->b);
        *p++ = ':';

        for (weight_t i = 0; i < rel->nb; i++)
        {
          int side = rel->primes[i].side;
          if (side == s1 || side == s2)
          {
            op = p;
            p = u64toa16(p, (uint64_t) rel->primes[i].h);
            *p++ = ',';
            t = p - op;
            for (j = (unsigned int) ((rel->primes[i].e) - 1); j--; p += t)
              memcpy(p, op, t);
          }
        }

        *(--p) = '\n';
        p[1] = 0;
        if (fputs(buf, outfiles[s1][s2]) == EOF)
        {
          perror("Error writing relation");
          abort();
        }
        nrels[s1][s2]++;
      }
    }
  }

  return NULL;
}

/*************************** main ********************************************/

int main(int argc, char **argv)
{
  char * argv0 = argv[0];
  param_list pl;
  renumber_t tab;

#ifdef HAVE_MINGW
    _fmode = _O_BINARY;		/* Binary open for all files */
#endif

  double wct0 = wct_seconds();

  param_list_init(pl);
  declare_usage(pl);
  argv++,argc--;

  if (argc == 0)
    usage (pl, argv0);

  for( ; argc ; ) {
    if (param_list_update_cmdline(pl, &argc, &argv))
      continue;
    /* Since we accept file names freeform, we decide to never abort
     * on unrecognized options */
    break;
  }
  /* print command-line arguments */
  verbose_interpret_parameters(pl);
  param_list_print_command_line (stdout, pl);
  fflush(stdout);

  const char *polyfilename = param_list_lookup_string(pl, "poly");
  const char *renumberfilename = param_list_lookup_string(pl, "renumber");
  const char * filelist = param_list_lookup_string(pl, "filelist");
  const char * basepath = param_list_lookup_string(pl, "basepath");
  const char * outprefix = param_list_lookup_string(pl, "outprefix");
  const char * outdir = param_list_lookup_string(pl, "outdir");

  if (param_list_warn_unused(pl))
  {
    fprintf(stderr, "Error, unused parameters are given\n");
    usage(pl, argv0);
  }

  if (polyfilename == NULL)
  {
    fprintf (stderr, "Error, missing -poly command line argument\n");
    usage (pl, argv0);
  }
  if (renumberfilename == NULL)
  {
    fprintf (stderr, "Error, missing -renumber command line argument\n");
    usage (pl, argv0);
  }
  if (outprefix == NULL)
  {
    fprintf (stderr, "Error, missing -outprefix command line argument\n");
    usage (pl, argv0);
  }
  if (outdir == NULL)
  {
    fprintf (stderr, "Error, missing -outdir command line argument\n");
    usage (pl, argv0);
  }
  if (basepath && !filelist)
  {
    fprintf(stderr, "Error, -basepath only valid with -filelist\n");
    usage(pl, argv0);
  }
  if ((filelist != NULL) + (argc != 0) != 1)
  {
    fprintf(stderr, "Error, provide either -filelist or freeform file names\n");
    usage(pl, argv0);
  }

  cado_poly_init(poly);
  if (!cado_poly_read (poly, polyfilename))
  {
    fprintf (stderr, "Error reading polynomial file\n");
    exit (EXIT_FAILURE);
  }

  renumber_init_for_reading (tab);
  renumber_read_table (tab, renumberfilename);

  for (unsigned int s1 = 0; s1 < tab->nb_polys; s1++)
    for (unsigned int s2 = s1+1; s2 < tab->nb_polys; s2++)
    {
      char *filename;
      int rc;
      rc = asprintf(&filename, "%s/%s.%d%d.rels",  outdir, outprefix, s1, s2);
      ASSERT_ALWAYS (rc >= 0);
      outfiles[s1][s2] = fopen_maybe_compressed(filename, "w");
      ASSERT_ALWAYS (outfiles[s1][s2]);
      free (filename);
    }

  char ** infiles = filelist ? filelist_from_file (basepath, filelist, 0) : argv;

  filter_rels(infiles, (filter_rels_callback_t) &thread_split_rel,
              (void *) tab, EARLYPARSE_NEED_INDEX | EARLYPARSE_NEED_AB_HEXA,
              NULL, NULL);

  for (unsigned int s1 = 0; s1 < tab->nb_polys; s1++)
    for (unsigned int s2 = s1+1; s2 < tab->nb_polys; s2++)
    {
      char *filename;
      int rc;
      rc = asprintf(&filename, "%s/%s.%d%d.rels",  outdir, outprefix, s1, s2);
      ASSERT_ALWAYS (rc >= 0);
      fclose_maybe_compressed (outfiles[s1][s2], filename);

      printf ("# %" PRIu64 " relations of type %d%d in %s\n", nrels[s1][s2], s1,
              s2, filename);
      free (filename);
    }

  /* print usage of time and memory */
  print_timing_and_memory(wct0);

  if (filelist)
    filelist_clear(infiles);
  renumber_clear (tab);
  cado_poly_clear (poly);
  param_list_clear(pl);
  return EXIT_SUCCESS;
}
