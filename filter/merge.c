/* merge --- main program to merge relations into relation-sets (cycles)

Copyright 2008-2009 Francois Morain.
Reviewed by Paul Zimmermann, February 2009.

This file is part of CADO-NFS.

CADO-NFS is free software; you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 2.1 of the License, or (at your option)
any later version.

CADO-NFS is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
details.

You should have received a copy of the GNU Lesser General Public License
along with CADO-NFS; see the file COPYING.  If not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/* Algorithm: digest and interpolation from Cavallar02.

  Stefania Cavallar, On the Number Field Sieve Integer Factorisation Algorithm,
  PhD thesis, University of Leiden, 2002, 108 pages.
 */

#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>   /* for _O_BINARY */
#include <string.h> /* for strcmp */
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "portability.h"

#include "filter_config.h"
#include "utils_with_io.h"
#include "merge_replay_matrix.h" /* for filter_matrix_t */
#include "report.h"     /* for report_t */
#include "markowitz.h" /* for MkzInit */
#include "merge_mono.h" /* for mergeOneByOne */
#include "sparse.h"

static void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "mat", "input purged file");
  param_list_decl_usage(pl, "out", "output history file");
  param_list_decl_usage(pl, "keep", "excess to keep (default "
                                    STR(DEFAULT_FILTER_EXCESS) ")");
  param_list_decl_usage(pl, "skip", "number of heavy columns to bury (default "
                                    STR(DEFAULT_MERGE_SKIP) ")");
  param_list_decl_usage(pl, "force-bury-below-index", "columns below this index are buried early on when reading the matrix");
  param_list_decl_usage(pl, "maxlevel", "maximum number of rows in a merge "
                            "(default " STR(DEFAULT_MERGE_MAXLEVEL) ")");
  param_list_decl_usage(pl, "target_density", "stop when the average row density exceeds this value"
                            " (default " STR(DEFAULT_MERGE_TARGET_DENSITY) ")");
  param_list_decl_usage(pl, "resume", "resume from history file");
  param_list_decl_usage(pl, "mkztype", "controls how the weight of a merge is "
                            "approximated (default " STR(DEFAULT_MERGE_MKZTYPE) ")");
  param_list_decl_usage(pl, "wmstmax", "controls until when a mst is used with "
                            "-mkztype 2 (default " STR(DEFAULT_MERGE_WMSTMAX) ")");
  param_list_decl_usage(pl, "forbidden-cols", "list of columns that cannot be "
                                              "used for merges");
  param_list_decl_usage(pl, "force-posix-threads", "(switch)");
  param_list_decl_usage(pl, "path_antebuffer", "path to antebuffer program");
  param_list_decl_usage(pl, "v", "verbose level");
  param_list_decl_usage(pl, "t", "number of threads");
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
    char *argv0 = argv[0];

    filter_matrix_t mat[1];
    report_t rep[1];

    int nthreads = 1;
    int maxlevel = DEFAULT_MERGE_MAXLEVEL;
    uint32_t keep = DEFAULT_FILTER_EXCESS;
    uint32_t skip = DEFAULT_MERGE_SKIP;
    double target_density = DEFAULT_MERGE_TARGET_DENSITY;
    uint32_t mkztype = DEFAULT_MERGE_MKZTYPE;
    uint32_t wmstmax = DEFAULT_MERGE_WMSTMAX;
    unsigned int force_bury_below_index = 0;

                               /* use real MST minimum for wt[j] <= wmstmax*/

#ifdef HAVE_MINGW
    _fmode = _O_BINARY;     /* Binary open for all files */
#endif

    double tt;
    double wct0 = wct_seconds ();
    param_list pl;
    int verbose = 0;
    param_list_init (pl);
    declare_usage(pl);
    argv++,argc--;

    param_list_configure_switch(pl, "force-posix-threads", &filter_rels_force_posix_threads);
    param_list_configure_switch(pl, "v", &verbose);

#ifdef HAVE_MINGW
    _fmode = _O_BINARY;     /* Binary open for all files */
#endif

    if (argc == 0)
      usage (pl, argv0);

    for( ; argc ; ) {
      if (param_list_update_cmdline(pl, &argc, &argv)) continue;
      fprintf (stderr, "Unknown option: %s\n", argv[0]);
      usage (pl, argv0);
    }
    /* print command-line arguments */
    verbose_interpret_parameters (pl);
    param_list_print_command_line (stdout, pl);
    fflush(stdout);

    const char * purgedname = param_list_lookup_string (pl, "mat");
    const char * outname = param_list_lookup_string (pl, "out");
    /* -resume can be useful to continue a merge stopped due  */
    /* to a too small value of -maxlevel                      */
    const char * resumename = param_list_lookup_string (pl, "resume");
    const char *path_antebuffer = param_list_lookup_string(pl, "path_antebuffer");
    const char *forbidden_cols = param_list_lookup_string(pl, "forbidden-cols");

    param_list_parse_int (pl, "t", &nthreads);
#ifdef HAVE_OPENMP
    omp_set_num_threads (nthreads);
#endif

    param_list_parse_int (pl, "maxlevel", &maxlevel);
    param_list_parse_uint (pl, "keep", &keep);
    param_list_parse_uint (pl, "skip", &skip);
    param_list_parse_uint (pl, "force-bury-below-index", &force_bury_below_index);
    if (force_bury_below_index > skip) {
        fprintf(stderr, "cannot have force_bury_below_index > skip\n");
        exit(EXIT_FAILURE);
    }

    param_list_parse_double (pl, "target_density", &target_density);

    param_list_parse_uint (pl, "mkztype", &mkztype);
    param_list_parse_uint (pl, "wmstmax", &wmstmax);

    /* Some checks on command line arguments */
    if (param_list_warn_unused(pl))
    {
      fprintf(stderr, "Error, unused parameters are given\n");
      usage(pl, argv0);
    }

    if (purgedname == NULL)
    {
      fprintf(stderr, "Error, missing -mat command line argument\n");
      usage (pl, argv0);
    }
    if (outname == NULL)
    {
      fprintf(stderr, "Error, missing -out command line argument\n");
      usage (pl, argv0);
    }
    if (maxlevel <= 0 || maxlevel > MERGE_LEVEL_MAX)
    {
      fprintf (stderr, "Error: maxlevel should be positive and less than %d\n",
                       MERGE_LEVEL_MAX);
      usage (pl, argv0);
    }
    if (mkztype > 2)
    {
      fprintf (stderr, "Error: -mkztype should be 0, 1, or 2.\n");
      usage (pl, argv0);
    }

    set_antebuffer_path (argv0, path_antebuffer);

    /* Read number of rows and cols on first line of purged file */
    purgedfile_read_firstline (purgedname, &(mat->nrows), &(mat->ncols));

    /* initialize rep (i.e., mostly opens outname) and write matrix dimension */
    rep->type = 0;
    rep->outfile = fopen_maybe_compressed (outname, "w");
    ASSERT_ALWAYS(rep->outfile != NULL);

    /* Init structure containing the matrix and the heap of potential merges */
    initMat (mat, maxlevel, keep, skip, force_bury_below_index);

    /* Read all rels and fill-in the mat structure */
    tt = seconds ();
    filter_matrix_read (mat, purgedname);
    printf ("Time for filter_matrix_read: %2.2lfs\n", seconds () - tt);

    /* resume from given history file if needed */
    if (resumename != NULL)
      resume (rep, mat, resumename);

    /* Some columns can be disable so merge won't use them as pivot */
    if (forbidden_cols != NULL)
    {
      printf ("Disabling columns from %s\n", forbidden_cols);
      matR_disable_cols (mat, forbidden_cols);
    }

    mat->verbose = verbose;
    mat->wmstmax = wmstmax;
    mat->mkztype = mkztype;
    tt = seconds();
    MkzInit (mat, 1);
    printf ("Time for MkzInit: %2.2lfs\n", seconds()-tt);

    mergeOneByOne (rep, mat, maxlevel, target_density);

    fclose_maybe_compressed (rep->outfile, outname);
    printf ("Final matrix has N=%" PRIu64 " nc=%" PRIu64 " (%" PRId64 ") "
            "W=%" PRIu64 " W*N=%.2e W/N=%.2f\n",
            mat->rem_nrows, mat->rem_ncols,
            ((int64_t) mat->rem_nrows) - ((int64_t) mat->rem_ncols), mat->weight,
            compute_WN (mat), compute_WoverN (mat));
    fflush (stdout);
    MkzClear (mat, 1);
    clearMat (mat);

    param_list_clear (pl);

    printf ("Total merge time: %.2f seconds\n", seconds ());

    print_timing_and_memory (stdout, wct0);

    return 0;
}
