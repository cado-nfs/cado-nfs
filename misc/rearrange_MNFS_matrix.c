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

#define SORT_NOSORTING 0
#define SORT_SIDE_INCR 1
#define SORT_PRIME_INCR 2
#define SORT_SIDE_INCR_THEN_PRIME_INCR 3

typedef struct {
  index_t index;
  int side;
  uint64_t new_id;
} col_data_t ;

static inline int
cmp_prime_incr (const void *px, const void *py)
{
  col_data_t *x = (col_data_t *) px;
  col_data_t *y = (col_data_t *) py;

  return (x->index > y->index);
}

static inline int
cmp_side_incr_prime_incr (const void *px, const void *py)
{
  col_data_t *x = (col_data_t *) px;
  col_data_t *y = (col_data_t *) py;

  if (x->side == y->side)
    return (x->index > y->index);
  else
    return (x->side > y->side);
}

static void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "poly", "input polynomial file");
  param_list_decl_usage(pl, "renumber", "input file for renumbering table");
  param_list_decl_usage(pl, "ideals", "ideals file produced by replay");
  param_list_decl_usage(pl, "matrix", "matrix file produced by replay");
  param_list_decl_usage(pl, "sm", "SM file produced by sm");
  param_list_decl_usage(pl, "nsm", "number of SM on side 0,1,...");
  param_list_decl_usage(pl, "t", "number of threads (default 1)");
  param_list_decl_usage(pl, "new-ideals", "new ideals file");
  param_list_decl_usage(pl, "new-matrix", "new matrix file");
  param_list_decl_usage(pl, "new-sm", "new SM file");
  param_list_decl_usage(pl, "side-info", "info on sides of rows");
  param_list_decl_usage(pl, "sorting", "How to sort the colums (in [0..3]).\n"
            "                 0 :  no sorting, keep the same order as in input\n"
            "                 1 :  sort by increasing side, inside each side "
                                 "keep the same order as in input (default)\n"
            "                 2 :  sort by increasing prime value (in fact by "
                                      "index value in the renumbering table)\n"
            "                 3 :  sort by increasing side then by increassing "
                                                                 "prime value\n");
  verbose_decl_usage(pl);
}

static void
usage (param_list pl, char *argv0)
{
  param_list_print_usage(pl, argv0, stderr);
  exit(EXIT_FAILURE);
}

/*************************** main ********************************************/
int main(int argc, char **argv)
{
  char * argv0 = argv[0];
  param_list pl;
  renumber_t tab;
  uint64_t ncols = 0;
  cado_poly poly;

  unsigned int sorting = SORT_SIDE_INCR;

#ifdef HAVE_MINGW
    _fmode = _O_BINARY;  /* Binary open for all files */
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
    fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
    usage (pl, argv0);
  }
  /* print command-line arguments */
  verbose_interpret_parameters(pl);
  param_list_print_command_line (stdout, pl);
  fflush(stdout);

  const char *polyfilename = param_list_lookup_string(pl, "poly");
  const char *renumberfilename = param_list_lookup_string(pl, "renumber");
  const char *idealsfilename = param_list_lookup_string(pl, "ideals");
  const char *smfilename = param_list_lookup_string(pl, "sm");
  const char *matrixfilename = param_list_lookup_string(pl, "matrix");
  const char *newidealsfilename = param_list_lookup_string(pl, "new-ideals");
  const char *newmatrixfilename = param_list_lookup_string(pl, "new-matrix");
  const char *newsmfilename = param_list_lookup_string(pl, "new-sm");
  const char *sideinfofilename = param_list_lookup_string(pl, "side-info");

  param_list_parse_uint64 (pl, "ncols", &ncols);
  param_list_parse_uint (pl, "sorting", &sorting);

  int nsm[NB_POLYS_MAX];
  int nread = param_list_parse_int_list (pl, "nsm", nsm, NB_POLYS_MAX, ",");

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
  if (idealsfilename == NULL)
  {
    fprintf (stderr, "Error, missing -ideals command line argument\n");
    usage (pl, argv0);
  }
  if (matrixfilename == NULL)
  {
    fprintf (stderr, "Error, missing -matrix command line argument\n");
    usage (pl, argv0);
  }
  if (smfilename == NULL)
  {
    fprintf (stderr, "Error, missing -sm command line argument\n");
    usage (pl, argv0);
  }
  if (newidealsfilename == NULL)
  {
    fprintf (stderr, "Error, missing -new-ideals command line argument\n");
    usage (pl, argv0);
  }
  if (newmatrixfilename == NULL)
  {
    fprintf (stderr, "Error, missing -new-matrix command line argument\n");
    usage (pl, argv0);
  }
  if (newsmfilename == NULL)
  {
    fprintf (stderr, "Error, missing -new-sm command line argument\n");
    usage (pl, argv0);
  }
  if (sideinfofilename == NULL)
  {
    fprintf (stderr, "Error, missing -side-info command line argument\n");
    usage (pl, argv0);
  }
  if (sorting != SORT_NOSORTING && sorting != SORT_SIDE_INCR &&
      sorting != SORT_PRIME_INCR && sorting != SORT_SIDE_INCR_THEN_PRIME_INCR)
  {
    fprintf (stderr, "Error, wrong value for -sorting command line argument\n");
    usage (pl, argv0);
  }


  cado_poly_init(poly);
  if (!cado_poly_read (poly, polyfilename))
  {
    fprintf (stderr, "Error reading polynomial file\n");
    exit (EXIT_FAILURE);
  }

  int nsm_tot = 0;
  for(int side = 0; side < poly->nb_polys; side++)
  {
    if (side >= nread)
      nsm[side] = poly->pols[side]->deg;
    nsm_tot += nsm[side];
  }

  renumber_init_for_reading (tab);
  renumber_read_table (tab, renumberfilename);

  col_data_t *cols_data = NULL;

  /* Read the input file containing links between column number and index in the
   * renumbering table*/
  {
    FILE *inid = NULL;
    uint64_t i = 0, col;
    index_t h;

    inid = fopen_maybe_compressed (idealsfilename, "r");
    FATAL_ERROR_CHECK (inid == NULL, "Cannot open ideals file");

    if (fscanf (inid, "# %" SCNu64 "\n", &ncols) != 1)
    {
      fprintf(stderr, "Error while reading first line of %s\n", idealsfilename);
      abort();
    }

    /* malloc column data structure */
    cols_data = (col_data_t *) malloc (ncols * sizeof (col_data_t));
    ASSERT_ALWAYS (cols_data != NULL);

    stats_data_t stats; /* struct for printing progress */
    stats_init (stats, stdout, &i, nbits(ncols)-5, "Read", "lines", "", "lines");
    while (fscanf (inid, "%" SCNu64 " %" PRid "\n", &col, &h) == 2)
    {
      FATAL_ERROR_CHECK (col >= ncols, "Too big value of column number");
      FATAL_ERROR_CHECK (h >= tab->size, "Too big value of index");
      ASSERT_ALWAYS (col == i);

      cols_data[i].index = h;
      cols_data[i].side = renumber_get_side_from_index (tab, h, poly);

      i++;
      if (stats_test_progress (stats))
        stats_print_progress (stats, i, 0, 0, 0);
    }
    stats_print_progress (stats, i, 0, 0, 1);
    ASSERT_ALWAYS (feof(inid));
    ASSERT_ALWAYS (i == ncols);
    fclose_maybe_compressed (inid, idealsfilename);
  }

  /* Sort the column and print new ideal file */
  {
    FILE *outid = NULL;
    outid = fopen_maybe_compressed (newidealsfilename, "w");
    FATAL_ERROR_CHECK (outid == NULL, "Cannot open new-ideals file");
    fprintf (outid, "# %" PRIu64 "\n", ncols);

    /* do not sort, keep the same order as in input */
    if (sorting == SORT_NOSORTING)
    {
      for (uint64_t i = 0; i < ncols; i++)
      {
        cols_data[i].new_id = i;
        fprintf (outid, "%" PRIu64 " %" PRid "\n", cols_data[i].new_id,
                                                   cols_data[i].index);
      }
    }
    /* sort by increasing side, inside each side keep the same order as in input */
    else if (sorting == SORT_SIDE_INCR)
    {
      uint64_t next_id = 0;
      for (int side = 0; side < poly->nb_polys; side++)
      {
        printf ("side %d begins at columns %" PRIu64 "\n", side, next_id);
        for (uint64_t i = 0; i < ncols; i++)
          if (cols_data[i].side == side)
          {
            cols_data[i].new_id = next_id++;
            fprintf (outid, "%" PRIu64 " %" PRid "\n", cols_data[i].new_id,
                                                       cols_data[i].index);
          }
      }
      ASSERT_ALWAYS (next_id == ncols);
    }
    /* sort by increasing index in renumbering table */
    else if (sorting == SORT_PRIME_INCR)
    {
      col_data_t *sort_data = NULL;
      sort_data = (col_data_t *) malloc (ncols * sizeof (col_data_t));
      ASSERT_ALWAYS (sort_data != NULL);
      memcpy (sort_data, cols_data, ncols * sizeof (col_data_t));
      for (uint64_t i = 0; i < ncols; i++)
        sort_data[i].new_id = i;

      qsort (sort_data, ncols, sizeof (col_data_t), cmp_prime_incr);

      for (uint64_t i = 0; i < ncols; i++)
      {
        uint64_t old_id = sort_data[i].new_id;
        cols_data[old_id].new_id = i;
        fprintf (outid, "%" PRIu64 " %" PRid "\n", i, cols_data[old_id].index);
      }
      free (sort_data);
    }
    /* sort by increasing side then by increasing index in renumbering table */
    else if (sorting == SORT_SIDE_INCR_THEN_PRIME_INCR)
    {
      col_data_t *sort_data = NULL;
      sort_data = (col_data_t *) malloc (ncols * sizeof (col_data_t));
      ASSERT_ALWAYS (sort_data != NULL);
      memcpy (sort_data, cols_data, ncols * sizeof (col_data_t));
      for (uint64_t i = 0; i < ncols; i++)
        sort_data[i].new_id = i;

      qsort (sort_data, ncols, sizeof (col_data_t), cmp_side_incr_prime_incr);

      for (uint64_t i = 0; i < ncols; i++)
      {
        uint64_t old_id = sort_data[i].new_id;
        cols_data[old_id].new_id = i;
        fprintf (outid, "%" PRIu64 " %" PRid "\n", i, cols_data[old_id].index);
      }

      free (sort_data);
    }
    else /* should never happen */
      ASSERT_ALWAYS (0);

    fclose_maybe_compressed (outid, newidealsfilename);
  }

  FILE *inmat = NULL, *outmat = NULL;
  FILE *insm = NULL, *outsm = NULL;
  FILE *outsideinfo = NULL;

  outsideinfo = fopen_maybe_compressed (sideinfofilename, "w");
  FATAL_ERROR_CHECK (outsideinfo == NULL, "Cannot open new-ideals file");

  inmat = fopen_maybe_compressed (matrixfilename, "r");
  FATAL_ERROR_CHECK (inmat == NULL, "Cannot open matrix file");
  insm = fopen_maybe_compressed (smfilename, "r");
  FATAL_ERROR_CHECK (insm == NULL, "Cannot open sm file");

  outmat = fopen_maybe_compressed (newmatrixfilename, "w");
  FATAL_ERROR_CHECK (outmat == NULL, "Cannot open new-matrix file");
  outsm = fopen_maybe_compressed (newsmfilename, "w");
  FATAL_ERROR_CHECK (outsm == NULL, "Cannot open new-sm file");

  uint64_t nr, nr2, nc;
  int nsm_tot_read;
  char tmp[1024];
  if (fscanf (inmat, "%" SCNu64 " %" SCNu64 "\n", &nr, &nc) != 2)
  {
    fprintf(stderr, "Error while reading first line of %s\n", matrixfilename);
    abort();
  }
  ASSERT_ALWAYS (nc == ncols);
  fprintf (outmat, "%" PRIu64 " %" PRIu64 "\n", nr, nc);

  if (fscanf (insm, "%" SCNu64 " %d %s\n", &nr2, &nsm_tot_read, tmp) != 3)
  {
    fprintf(stderr, "Error while reading first line of %s\n", smfilename);
    abort();
  }
  ASSERT_ALWAYS (nr == nr2);
  ASSERT_ALWAYS (nsm_tot_read == nsm_tot);
  fprintf (outsm, "%" PRIu64 " %d %s\n", nr2, nsm_tot_read, tmp);

  uint64_t i = 0;
  stats_data_t stats; /* struct for printing progress */
  stats_init (stats, stdout, &i, nbits(nr)-5, "Read", "rows", "", "rows");
  unsigned int rw;
  while (fscanf (inmat, "%u", &rw) == 1)
  {
    uint64_t nonvoidside = 0; /* bit vector of which sides appear in the rel */
    fprintf (outmat, "%u", rw);
    for (unsigned int k = 0; k < rw; k++)
    {
      int w;
      uint64_t id;
      int ret = fscanf (inmat, " %" SCNu64 ":%d", &id, &w);
      ASSERT_ALWAYS (ret == 2);
      FATAL_ERROR_CHECK (id >= ncols, "Too big value of column number");
      fprintf (outmat, " %" PRIu64 ":%d", cols_data[id].new_id, w);
      nonvoidside |= ((uint64_t) 1) << cols_data[id].side;
    }
    int ret = fscanf (inmat, "\n");
    ASSERT_ALWAYS (ret == 0);
    fprintf (outmat, "\n");

    fprintf (outsideinfo, "%"PRIu64 "", i);
    for (uint64_t b = nonvoidside, s = 0; b != 0 ; b>>=1, s++)
      if (b & ((uint64_t) 1))
        fprintf (outsideinfo, " %" PRIu64 "", s);
    fprintf (outsideinfo, "\n");

    for(int side = 0, n = 0; side < poly->nb_polys; side++, nonvoidside>>=1)
    {
      for (int k = 0; k < nsm[side]; k++)
      {
        n++;
        int ret = fscanf (insm, "%s", tmp);
        ASSERT_ALWAYS (ret == 1);
        if (nonvoidside & ((uint64_t) 1))
          fprintf (outsm, "%s%s", tmp, (n == nsm_tot) ? "" : " ");
        else
          fprintf (outsm, "0%s", (n == nsm_tot) ? "" : " ");
      }
    }
    ret = fscanf (insm, "\n");
    ASSERT_ALWAYS (ret == 0);
    fprintf (outsm, "\n");

    i++;
    if (stats_test_progress (stats))
      stats_print_progress (stats, i, 0, 0, 0);
  }
  stats_print_progress (stats, i, 0, 0, 1);
  ASSERT_ALWAYS (i == nr);
  ASSERT_ALWAYS (feof(inmat));

  fclose_maybe_compressed (outsideinfo, sideinfofilename);
  fclose_maybe_compressed (inmat, matrixfilename);
  fclose_maybe_compressed (outmat, newmatrixfilename);
  fclose_maybe_compressed (insm, smfilename);
  fclose_maybe_compressed (outsm, newsmfilename);


  /* print usage of time and memory */
  print_timing_and_memory(wct0);

  free (cols_data);
  renumber_clear (tab);
  cado_poly_clear (poly);
  param_list_clear(pl);
  return EXIT_SUCCESS;
}
