/* purge --- perform singleton removal and clique removal
 * 
 * Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015
 * Cyril Bouvier, Alain Filbois, Francois Morain, Paul Zimmermann
 * 
 * This file is part of CADO-NFS.
 * 
 * CADO-NFS is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 2.1 of the License, or
 * (at your option) any later version.
 * 
 * CADO-NFS is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with CADO-NFS; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA 02110-1301, USA.
*/

/* References:
 * On the Number Field Sieve Integer Factorisation Algorithm,
 * Stefania Cavallar, PhD Thesis, University of Leiden, 2002.
 *
 * The filtering step of discrete logarithm and integer factorization
 * algorithms.
 * Cyril Bouvier, preprint, 22 pages. 2013.
 */

/*
 * Important remark:
 *   A relation corresponds to a row of the matrix and a ideal corresponds to
 *   a column of the matrix.
 *
 * This program works in two passes over the relation files:
 * - the first pass loads in memory only indexes of columns >= col_min_index
 *   and keeps a count of the weight of each column in cols_weight.
 *   Then, a first step of singleton removal is performed followed by 'nsteps'
 *   steps of singleton removal and clique removal, in order to obtained the
 *   final excess 'keep'.
 * - the second pass goes through the relations again, and dumps the remaining
 *   ones in the format needed by 'merge'.

 * This program uses the following data structures:
 * row_used[i]    - non-zero iff row i is kept (so far)
 * row_compact[i] - list of indexes of columns (greater than or equal to
 *                  col_min_index) of the row i (terminated by a -1 sentinel)
 * cols_weight [h] - weight of the column h in current rows (saturates at 256)
 */

/*
 * Exit value:
 * - 0 if enough relations
 * - 1 if an error occurred (then we should abort the factorization)
 * - 2 if not enough relations
 */

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
#include "filter_config.h"
#include "purge_matrix.h"
#include "singleton_removal.h"
#include "clique_removal.h"

/*****************************************************************************/

/* write final values to stdout */
/* This output, incl. "Final values:", is required by the script */
static void
print_final_values (purge_matrix_ptr mat, double weight)
{
  int64_t excess = purge_matrix_compute_excess (mat);
  fprintf (stdout, "Final values:\nnrows=%" PRIu64 " ncols=%" PRIu64 " "
                   "excess=%" PRId64 "\nweight=%1.0f weight*nrows=%1.2e\n",
                   mat->nrows, mat->ncols, excess, weight,
                   weight * (double) mat->nrows);
  fflush (stdout);
}

/* If nsteps is negative, then the value is chosen by the function. */
static void singletons_and_cliques_removal(purge_matrix_ptr mat, int nsteps,
                                           int64_t final_excess,
                                           double required_excess,
                                           unsigned int nthreads, int verbose)
{
  uint64_t oldnrows = 0;
  int64_t oldexcess, excess, target_excess;
  int count;

  /* First step of singletons removal */
  fprintf(stdout, "\nStep 0: only singleton removal\n");
  excess = singleton_removal (mat, nthreads, verbose);

  if (excess <= 0) /* covers case nrows = ncols = 0 */
  {
    fprintf (stdout, "number of rows < number of columns + keep\n");
    print_final_values (mat, 0);
    exit(2);
  }

  if ((double) excess < required_excess * ((double) mat->ncols))
  {
    fprintf (stdout, "(excess / ncols) = %.2f < %.2f. See -required_excess "
                     "argument.\n", ((double) excess / (double) mat->ncols),
                     required_excess);
    print_final_values (mat, 0);
    exit(2);
  }

  /* delta between the current excess and the desired final excess. */
  int64_t delta_excess = excess - final_excess;

  /* If nsteps was not given in the command line, adjust nsteps in
     [1..DEFAULT_PURGE_NSTEPS] so that each step removes at least about 1% wrt
     the number of columns */
  if (nsteps < 0)
  {
    /* If we are less than or equal to the wanted final excess, there is no need
     * for clique removal, so set nsteps to 0*/
    if (delta_excess <= 0)
      nsteps = 0;
    else if ((uint64_t) delta_excess / DEFAULT_PURGE_NSTEPS < mat->ncols / 100)
      nsteps = 1 + (100 * delta_excess) / mat->ncols;
    else
      nsteps = DEFAULT_PURGE_NSTEPS;
  }

  int64_t chunk = 0;
  if (nsteps > 0)
  {
    chunk = delta_excess / nsteps;
    fprintf(stdout, "# INFO: number of clique removal steps: %d\n", nsteps);
    fprintf(stdout, "# INFO: At each step, excess will be decreased by "
                    "%" PRId64 "\n", chunk);
    fflush (stdout);
  }
  else
    fprintf(stdout, "# INFO: No step of clique removal will be done\n");

  /* nsteps steps of clique removal + singletons removal */
  for (count = 0; count < nsteps && excess > 0; count++)
  {
    oldnrows = mat->nrows;
    oldexcess = excess;
    target_excess = excess - chunk;
    if (target_excess < final_excess)
      target_excess = final_excess;
    fprintf(stdout, "\nStep %u on %u: target excess is %" PRId64 "\n",
                    count + 1, nsteps, target_excess);
    fflush(stdout);

    /* prints some stats on columns and rows weight if verbose > 0. */
    if (verbose > 0)
    {
      purge_matrix_print_stats_columns_weight (stdout, mat, verbose);
      purge_matrix_print_stats_rows_weight (stdout, mat, verbose);
    }

    cliques_removal (mat, target_excess, nthreads, verbose);
    excess = singleton_removal (mat, nthreads, verbose);

    fprintf (stdout, "This step removed %" PRId64 " rows and decreased excess "
             "by %" PRId64 "\n", (int64_t) (oldnrows - mat->nrows),
             oldexcess - excess);
    if (oldexcess > excess)
      fprintf (stdout, "Each excess row deleted %2.2lf rows\n",
               (double) (oldnrows - mat->nrows) /
               (double) (oldexcess - excess));
  }


  /* May need an extra step of clique removal + singletons removal if excess is
     still larger than keep. It may happen due to the fact that each clique does
     not make the excess go down by one but can (rarely) left the excess
     unchanged. */
  if (excess > final_excess && nsteps > 0)
  {
    oldnrows = mat->nrows;
    oldexcess = excess;
    target_excess = final_excess;

    fprintf(stdout, "\nStep extra: target excess is %" PRId64 "\n",
                    target_excess);
    fflush(stdout);

    /* prints some stats on columns and rows weight if verbose > 0. */
    if (verbose > 0)
    {
      purge_matrix_print_stats_columns_weight (stdout, mat, verbose);
      purge_matrix_print_stats_rows_weight (stdout, mat, verbose);
    }

    cliques_removal (mat, target_excess, nthreads, verbose);
    excess = singleton_removal (mat, nthreads, verbose);

    fprintf(stdout, "This step removed %" PRId64 " rows and decreased excess "
                    "by %" PRId64 "\nEach excess row deleted %2.2lf rows\n",
                    (int64_t) (oldnrows-mat->nrows), (oldexcess-excess),
                    (double) (oldnrows-mat->nrows) / (double) (oldexcess-excess));
  }
}

/*************** Callback function called by filter_rels ********************/

/* Data struct for thread_print (called by filter_rels on pass 2) */
struct data_second_pass_s
{
  double W; /* Total weight of the matrix (counting only remaining rows) */
  purge_matrix_ptr mat;
  /* fd[0]: for printing kept relations */
  /* fd[1]: for printing deleted relations */
  FILE * fd[2];
};
typedef struct data_second_pass_s data_second_pass_t[1];
typedef struct data_second_pass_s * data_second_pass_ptr;

/* Callback function called by filter_rels on pass 2 */
void *
thread_print(data_second_pass_ptr arg, earlyparsed_relation_ptr rel)
{
  if (purge_matrix_is_row_active (arg->mat, rel->num))
  {
    arg->W += rel->nb;
    fputs(rel->line, arg->fd[0]);
  }
  else if (arg->fd[1] != NULL)
    fputs(rel->line, arg->fd[1]);
  return NULL;
}


/*********** utility functions for purge binary ****************/

  /* Build the file list (ugly). It is the concatenation of all
   *  b s p
   * where:
   *    b is the basepath (empty if not given)
   *    s ranges over all subdirs listed in the subdirlist (empty if no
   *    such list)
   *    p ranges over all paths listed in the filelist.
   */
char **filelist_from_file_with_subdirlist(const char *basepath,
					  const char *filelist,
					  const char *subdirlist)
{
    /* count the number of files in the filelist */
    int nfiles = 0;
    int nsubdirs = 0;
    char **fl = filelist_from_file(NULL, filelist, 0);
    for (char **p = fl; *p; p++, nfiles++);

    char **sl = filelist_from_file(basepath, subdirlist, 1);
    for (char **p = sl; *p; p++, nsubdirs++);

    char ** fic = (char **) malloc((nsubdirs * nfiles + 1) * sizeof(char *));
    ASSERT_ALWAYS(fic != NULL);

    char **full = fic;
    for (char **f = fl; *f; f++) {
	for (char **s = sl; *s; s++, full++) {
	    int ret = asprintf(full, "%s/%s", *s, *f);
	    ASSERT_ALWAYS(ret >= 0);
	}
    }
    *full = NULL;
    filelist_clear(fl);
    filelist_clear(sl);
    return fic;
}

static void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "filelist", "file containing a list of input files");
  param_list_decl_usage(pl, "subdirlist",
                               "file containing a list of subdirectories");
  param_list_decl_usage(pl, "basepath", "path added to all file in filelist");
  param_list_decl_usage(pl, "out", "outfile for remaining relations");
  param_list_decl_usage(pl, "nrels", "number of initial relations");
  param_list_decl_usage(pl, "col-max-index", "upper bound on the number of "
                                  "columns (must be at least the number\n"
                  "                   of prime ideals in renumber table)");
  param_list_decl_usage(pl, "col-min-index", "do not take into account columns"
                                             " with indexes <= col-min-index");
  param_list_decl_usage(pl, "keep", "wanted excess at the end of purge "
                                    "(default " STR(DEFAULT_FILTER_EXCESS) ")");
  param_list_decl_usage(pl, "nsteps", "maximal number of steps of clique "
                                      "removal (default: chosen in [1.."
                                             STR(DEFAULT_PURGE_NSTEPS) "])");
  param_list_decl_usage(pl, "required_excess", "\% of excess required at the "
                            "end of the 1st singleton removal step (default "
                            STR(DEFAULT_PURGE_REQUIRED_EXCESS) ")");
  param_list_decl_usage(pl, "outdel", "outfile for deleted relations (for DL)");
  param_list_decl_usage(pl, "t", "number of threads (default "
                                             STR(DEFAULT_PURGE_NTHREADS) ")");
  param_list_decl_usage(pl, "v", "(switch) verbose mode");
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


/*************************** main ********************************************/

int main(int argc, char **argv)
{
    char * argv0 = argv[0];
    param_list pl;
    uint64_t col_min_index_arg = UMAX(uint64_t);
    char ** input_files;
    uint64_t col_max_index_arg = 0;
    uint64_t nrows_init_arg = 0;
    purge_matrix_t mat; /* All info regarding the matrix is in this struct */
    int64_t keep = DEFAULT_FILTER_EXCESS; /* minimum final excess */
    int nsteps = -1; /* negative value means chosen by purge */
    double required_excess = DEFAULT_PURGE_REQUIRED_EXCESS;
    unsigned int nthreads = DEFAULT_PURGE_NTHREADS;
    int verbose = 0;

#ifdef HAVE_MINGW
    _fmode = _O_BINARY;		/* Binary open for all files */
#endif

    double wct0 = wct_seconds();

    param_list_init(pl);
    declare_usage(pl);
    argv++,argc--;

    param_list_configure_switch (pl, "-v", &verbose);
    param_list_configure_switch(pl, "force-posix-threads", &filter_rels_force_posix_threads);

    if (argc == 0)
      usage (pl, argv0);

    /* read all command-line parameters */
    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        /* Since we accept file names freeform, we decide to never abort
         * on unrecognized options */
        break;
    }
    /* print command-line arguments */
    verbose_interpret_parameters(pl);
    param_list_print_command_line (stdout, pl);
    fflush(stdout);

    /* read command-line parameters */
    param_list_parse_uint64(pl, "nrels", &nrows_init_arg);
    param_list_parse_uint64(pl, "col-max-index", &col_max_index_arg);
    param_list_parse_int64(pl, "keep", &keep);

    /* Only look at columns of index >= col-min-index */
    param_list_parse_uint64(pl, "col-min-index", &col_min_index_arg);


    param_list_parse_uint(pl, "t", &nthreads);
    param_list_parse_int(pl, "nsteps", &nsteps);
    param_list_parse_double(pl, "required_excess", &required_excess);

    /* These three parameters specify the set of input files, of the form
     * <base path>/<one of the possible subdirs>/<one of the possible
     * file names>
     *
     * possible subdirs are lister in the file passed as subdirlist.
     * Ditto for possible file names.
     *
     * file names need not be basenames, i.e. they may contain directory
     * components. subdirlist and basepath are optional.
     */
    const char *basepath = param_list_lookup_string(pl, "basepath");
    const char *subdirlist = param_list_lookup_string(pl, "subdirlist");
    const char *filelist = param_list_lookup_string(pl, "filelist");
    const char *purgedname = param_list_lookup_string(pl, "out");
    const char *deletedname = param_list_lookup_string(pl, "outdel");

    if (param_list_warn_unused(pl))
    {
      fprintf(stderr, "Error, unused parameters are given\n");
      usage(pl, argv0);
    }

    if (!purgedname) {
        fprintf(stderr, "Error, option -out is mandatory\n");
        exit(EXIT_FAILURE);
    }
    /* }}} */


    /*{{{ argument checking, and some statistics for things related to
     * the hash table. It needs several static parameters on the command
     * line. This is cumbersome, but while it can probably be avoided, it
     * also hard to do so efficiently */
    if ((basepath || subdirlist) && !filelist)
    {
      fprintf(stderr, "Error, -basepath / -subdirlist only valid with -filelist\n");
      usage(pl, argv0);
    }
    if ((filelist != NULL) + (argc != 0) != 1) {
      fprintf(stderr, "Error, provide either -filelist or freeform file names\n");
      usage(pl, argv0);
    }
    if (nrows_init_arg == 0)
    {
      fprintf(stderr, "Error, missing -nrels command line argument "
                      "(or nrels = 0)\n");
      usage(pl, argv0);
    }
    if (col_max_index_arg == 0)
    {
      fprintf(stderr, "Error, missing or wrong -col-max-index command line argument "
                      "(should be > 0)\n");
      usage(pl, argv0);
    }
    if (col_min_index_arg == UMAX(uint64_t))
    {
      fprintf(stderr, "Error, missing -col-min-index command line argument\n");
      usage(pl, argv0);
    }
    if (col_min_index_arg >= col_max_index_arg)
    {
      fprintf(stderr, "Error, col-min-index >= col-max-index\n");
      usage(pl, argv0);
    }
    /* If col_max_index_arg > 2^32, then we need index_t to be 64-bit */
    if (((col_max_index_arg >> 32) != 0) && sizeof(index_t) < 8)
    {
      fprintf(stderr, "Error, -col-max-index is too large for a 32-bit "
                      "program\nSee #define __SIZEOF_INDEX__ in typedefs.h\n");
      exit(EXIT_FAILURE);
    }
    if (nthreads == 0)
    {
      fprintf(stderr, "Error, -t should be non-zero\n");
      exit(EXIT_FAILURE);
    }

    /* Printing relevant information */
    comp_print_info_weight_function ();

    fprintf(stdout, "# INFO: number of rows: %" PRIu64 "\n", nrows_init_arg);
    fprintf(stdout, "# INFO: maximum possible index of a column: %" PRIu64
                    "\n", col_max_index_arg);
    fprintf(stdout, "# INFO: number of threads: %u\n", nthreads);
    fprintf(stdout, "# INFO: number of clique removal steps: ");
    if (nsteps < 0)
      fprintf(stdout, "will be chosen by the program\n");
    else
      fprintf(stdout, "%d\n", nsteps);
    ASSERT_ALWAYS(keep >= 0);
    fprintf(stdout, "# INFO: target excess: %" PRId64 "\n", keep);
    fflush (stdout);
    /*}}}*/

    /* }}} */

    purge_matrix_init (mat, nrows_init_arg, col_min_index_arg,
                       col_max_index_arg);


    set_antebuffer_path(argv0, param_list_lookup_string(pl, "path_antebuffer"));
    /* }}} */

    /*{{{ build the list of input files from the given args
     * If no filelist is given, files are on the command-line.
     * If no subdirlist is given, files are easily construct from basepath and
     * filelist.
     * If subdirlist is given, it is a little bit trickier, see
     * filelist_from_file_with_subdirlist for more details. */
    if (!filelist)
      input_files = argv;
    else if (!subdirlist)
      input_files = filelist_from_file(basepath, filelist, 0);
    else
      input_files = filelist_from_file_with_subdirlist(basepath, filelist, 
                                                       subdirlist);
    /*}}}*/

    /****************** Begin interesting stuff *************************/
    fprintf(stdout, "\nPass 1, reading and storing columns with index h >= "
                    "%" PRIu64 "\n", mat->col_min_index);

    /* first pass over relations in files */
    /* Note: Now that we no longer take a bitmap on input, all
     * relations are considered active at this point, so that we
     * do not need to pass a bitmap to filter_rels */
    filter_rels(input_files,
                (filter_rels_callback_t) &purge_matrix_set_row_from_rel,
                (void *) mat, EARLYPARSE_NEED_INDEX, NULL, NULL);

    if (mat->nrows != mat->nrows_init)
    {
      fprintf(stderr, "Error, -nrels value should match the number of scanned "
                      "relations\nexpected %" PRIu64 " relations, found "
                      "%" PRIu64 "\n", mat->nrows_init, mat->nrows);
      abort();
    }

    /* Take into account the memory allocated for all mat->row_compact[i] */
    purge_matrix_row_compact_update_mem_usage (mat);

    /* prints some stats on columns and rows weight if verbose > 0. */
    if (verbose > 0)
    {
      purge_matrix_print_stats_columns_weight (stdout, mat, verbose);
      purge_matrix_print_stats_rows_weight (stdout, mat, verbose);
    }

    /* MAIN FUNCTIONS: do singletons and cliques removal. */
    singletons_and_cliques_removal (mat, nsteps, keep, required_excess,
                                    nthreads, verbose);

    /* prints some stats on columns and rows weight if verbose > 0. */
    if (verbose > 0)
    {
      purge_matrix_print_stats_columns_weight (stdout, mat, verbose);
      purge_matrix_print_stats_rows_weight (stdout, mat, verbose);
    }

    if (mat->nrows < mat->ncols + keep)
    {
      fprintf (stdout, "number of rows < number of columns + keep\n");
      print_final_values (mat, 0);
      exit(2);
    }
    if (mat->nrows == 0 || mat->ncols == 0)
    {
      fprintf(stdout, "number of rows or number of columns is 0\n");
      print_final_values (mat, 0);
      exit(2);
    }

    /****** Pass 2: reread the relation files and write output file(s) ******/
    fprintf(stdout, "\nPass 2, reading and writing output file%s...\n",
                    deletedname == NULL ? "" : "s");
    data_second_pass_t data2;
    memset(data2, 0, sizeof(data_second_pass_t));
    data2->mat = mat;

    if (!(data2->fd[0] = fopen_maybe_compressed(purgedname, "w")))
    {
      fprintf(stderr, "Error, cannot open file %s for writing.\n", purgedname);
      exit(1);
    }

    if (deletedname != NULL)
    {
      if (!(data2->fd[1] = fopen_maybe_compressed(deletedname, "w")))
      {
        fprintf(stderr, "Error, cannot open file %s for writing.\n",
                        deletedname);
        exit(1);
      }
      /* Write the header line for the file of deleted relations. */
      fprintf(data2->fd[1], "# %" PRIu64 "\n", mat->nrows_init - mat->nrows);
    }

    /* Write the header line for the file of remaining relations:
     * compute last index i such that cols_weight[i] != 0
     */
    {
      uint64_t last_used = mat->col_max_index - 1;
      while (mat->cols_weight[last_used] == 0)
        last_used--;

      fprintf(data2->fd[0], "# %" PRIu64 " %" PRIu64 " %" PRIu64 "\n",
                            mat->nrows, last_used + 1, mat->ncols);
    }

    /* second pass over relations in files */
    filter_rels(input_files, (filter_rels_callback_t) &thread_print,
                (void *) data2, EARLYPARSE_NEED_LINE, NULL, NULL);

    /* write final values to stdout */
    /* This output, incl. "Final values:", is required by the script */
    print_final_values (mat, data2->W);

    /* Free allocated stuff */
    if (filelist)
      filelist_clear(input_files);

    fclose_maybe_compressed(data2->fd[0], purgedname);
    if (data2->fd[1])
      fclose_maybe_compressed(data2->fd[1], deletedname);

    purge_matrix_clear (mat);
    /* print usage of time and memory */
    print_timing_and_memory(wct0);

    param_list_clear(pl);

    return 0;
}
