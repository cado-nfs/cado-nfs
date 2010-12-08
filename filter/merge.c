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
#include <string.h> /* for strcmp */

#include "utils.h" /* for gzip_open */

#include "merge_opts.h" /* for USE_MARKOWITZ */
#include "filter_matrix.h" /* for filter_matrix_t */
#include "report.h"     /* for report_t */
#ifndef USE_MARKOWITZ
# include "swar.h"      /* for initSWAR */
#else
# include "markowitz.h" /* for MkzInit */
#endif
#include "merge_mono.h" /* for mergeOneByOne */

#ifdef USE_MPI
#include "mpi.h"
#include "merge_mpi.h"
#endif

#define CWMAX_DEFAULT 100
#define RWMAX_DEFAULT 100
#define MAXLEVEL_DEFAULT 10
#define KEEP_DEFAULT 128
#define FORBW_DEFAULT 0
#define RATIO_DEFAULT 1.1
#define COVERNMAX_DEFAULT 100

static void
usage (void)
{
  fprintf (stderr, "Usage: merge [options]\n");
  fprintf (stderr, "   -v             - print some extra information\n");
  fprintf (stderr, "   -mat   xxx     - input (purged) file is xxx\n");
  fprintf (stderr, "   -out   xxx     - output (history) file is xxx\n");
  fprintf (stderr, "   -cwmax nnn     - merge columns of weight <= nnn only (default %u)\n", CWMAX_DEFAULT);
  fprintf (stderr, "   -rwmax nnn     - merge rows of weight <= nnn only (default %u)\n", RWMAX_DEFAULT);
  fprintf (stderr, "   -maxlevel nnn  - merge up to nnn rows (default %u)\n",
	   MAXLEVEL_DEFAULT);
  fprintf (stderr, "   -keep nnn      - keep an excess of nnn (default %u)\n",
	   KEEP_DEFAULT);
  fprintf (stderr, "   -forbw nnn     - controls the optimization function (default %u, see below)\n",
	   FORBW_DEFAULT);
  fprintf (stderr, "   -ratio rrr     - maximal ratio cN(final)/cN(min) with forbw=0 (default %1.1f)\n",
	   RATIO_DEFAULT);
  fprintf (stderr, "   -coverNmax nnn - with forbw=3, stop when c/N exceeds nnn (default %u)\n", COVERNMAX_DEFAULT);
  fprintf (stderr, "   -itermax nnn   - if non-zero, stop when nnn columns have been removed (cf -resume)\n");
  fprintf (stderr, "   -resume xxx    - resume from history file xxx (cf -itermax)\n");
  fprintf (stderr, "\nThe different optimization functions are, where c is the total matrix weight\n");
  fprintf (stderr, "and N the number of rows (relation-sets):\n");
  fprintf (stderr, "   -forbw 0 - optimize the matrix size N (cf -ratio)\n");
  fprintf (stderr, "   -forbw 1 - stop when the product cN is minimal\n");
  fprintf (stderr, "   -forbw 3 - stop when the ratio c/N exceeds coverNmax\n");
  exit (1);
}

int
main (int argc, char *argv[])
{
    filter_matrix_t mat[1];
    report_t rep[1];
    char *purgedname = NULL, *outname = NULL;
    char *resumename = NULL;
    int cwmax = CWMAX_DEFAULT, rwmax = RWMAX_DEFAULT;
    int maxlevel = MAXLEVEL_DEFAULT, keep = KEEP_DEFAULT;
    int verbose = 0; /* default verbose level */
    double tt;
    double ratio = RATIO_DEFAULT; /* bound on cN_new/cN to stop the merge */
    int i, forbw = FORBW_DEFAULT, coverNmax = COVERNMAX_DEFAULT;
#ifdef USE_MARKOWITZ
    int wmstmax = 7; /* use real MST minimum for wt[j] <= wmstmax */
    int mkzrnd = 0;
    int mkztype = 2;
#endif
    int itermax = 0;

    /* print comand-line arguments */
    fprintf (stderr, "%s.r%s", argv[0], CADO_REV);
    for (i = 1; i < argc; i++)
      fprintf (stderr, " %s", argv[i]);
    fprintf (stderr, "\n");
#ifdef USE_MPI
    MPI_Init(&argc, &argv);
#endif

    while(argc > 1 && argv[1][0] == '-'){
        if (argc > 2 && strcmp (argv[1], "-mat") == 0){
	    purgedname = argv[2];
	    argc -= 2;
	    argv += 2;
	}
	else if (argc > 2 && strcmp (argv[1], "-cwmax") == 0){
	    cwmax = atoi(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
	else if (argc > 2 && strcmp (argv[1], "-rwmax") == 0){
	    rwmax = atoi(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
	else if (argc > 2 && strcmp (argv[1], "-maxlevel") == 0){
	    maxlevel = atoi(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
	else if (argc > 2 && strcmp (argv[1], "-keep") == 0){
	    keep = atoi(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
	else if (argc > 1 && strcmp (argv[1], "-v") == 0){
            verbose ++;
	    argc -= 1;
	    argv += 1;
	}
	else if (argc > 2 && strcmp (argv[1], "-forbw") == 0){
	    forbw = atoi(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
	else if (argc > 2 && strcmp (argv[1], "-ratio") == 0){
	    ratio = strtod(argv[2], NULL);
	    argc -= 2;
	    argv += 2;
	}
	else if (argc > 2 && strcmp (argv[1], "-coverNmax") == 0){
	    coverNmax = atoi(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
	else if (argc > 2 && strcmp (argv[1], "-out") == 0){
	    outname = argv[2];
	    argc -= 2;
	    argv += 2;
	}
	/* -resume can be useful to continue a merge stopped due
	   to a too small value of -maxlevel */
	else if (argc > 2 && strcmp (argv[1], "-resume") == 0){
	    resumename = argv[2];
	    argc -= 2;
	    argv += 2;
	}
#ifdef USE_MARKOWITZ
	else if (argc > 2 && strcmp (argv[1], "-wmstmax") == 0){
	    wmstmax = atoi(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
	else if (argc > 2 && strcmp (argv[1], "-mkzrnd") == 0){
	    mkzrnd = atoi(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
	else if (argc > 2 && strcmp (argv[1], "-mkztype") == 0){
	    mkztype = atoi(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
#endif
	/* -itermax can be used with -resume, for example:
	   merge -itermax 1000 -out his.tmp
           merge -resume his.tmp -out his.final */
	else if (argc > 2 && strcmp (argv[1], "-itermax") == 0){
	    itermax = atoi(argv[2]);
	    argc -= 2;
	    argv += 2;
	}
	else
	  usage ();
    }

    purgedfile_stream ps;
    purgedfile_stream_init(ps);
    purgedfile_stream_openfile(ps, purgedname);

    mat->nrows = ps->nrows;
    mat->ncols = ps->ncols;
    mat->keep  = keep;
    mat->cwmax = cwmax;
    mat->rwmax = rwmax;
    mat->mergelevelmax = maxlevel;
    mat->itermax = itermax;
    
#ifdef USE_MPI
    mpi_start_proc(outname,mat,purgedfile,purgedname,forbw,ratio,coverNmax,
		   resumename);
    /* TODO: clean the mat data structure (?) */
    MPI_Finalize();
    return 0;
#endif    
    initMat (mat, 0, ps->ncols);

    tt = seconds ();
    filter_matrix_read_weights (mat, ps);
    fprintf (stderr, "Getting column weights took %2.2lf\n", seconds () - tt);
    purgedfile_stream_rewind(ps);

    /* print weight counts */
    {
      unsigned long j, *nbm;
      int w;
      nbm = (unsigned long*) malloc ((maxlevel + 1) * sizeof (unsigned long));
      memset (nbm, 0, (maxlevel + 1) * sizeof (unsigned long));
      for (j = 0; j < (unsigned long) mat->ncols; j++)
        {
          w = mat->wt[GETJ(mat, j)];
          if (w <= maxlevel)
            nbm[w] ++;
        }
      for (j = 0; j <= (unsigned long) maxlevel; j++)
        if (nbm[j] != 0)
          fprintf (stderr, "There are %lu column(s) of weight %lu\n",
                   nbm[j], j);
      free (nbm);
    }

#ifndef USE_MARKOWITZ
    fprintf (stderr, "SWAR version\n");
    initSWAR (mat);
#else
    fprintf(stderr, "Markowitz version\n");
#endif
    fillmat (mat);
    
    tt = wct_seconds ();
    filter_matrix_read (mat, ps, verbose);
    fprintf (stderr, "Time for filter_matrix_read: %2.2lf\n", wct_seconds () - tt);

    /* initialize rep, i.e., mostly opens outname */
    init_rep (rep, outname, mat, 0, MERGE_LEVEL_MAX);
    /* output the matrix dimensions in the history file */
    report2 (rep, mat->nrows, mat->ncols);

    /* resume from given history file if needed */
    if (resumename != NULL)
      resume (rep, mat, resumename);

#ifdef USE_MARKOWITZ
    mat->wmstmax = wmstmax;
    mat->mkzrnd = mkzrnd;
    mat->mkztype = mkztype;
    tt = seconds();
    MkzInit (mat);
    fprintf (stderr, "Time for MkzInit: %2.2lf\n", seconds()-tt);
#endif

#if M_STRATEGY <= 2
# ifdef USE_MARKOWITZ
    fprintf(stderr, "merge NYI for Markowitz\n");
    return 1;
# endif
    merge (mat, maxlevel, verbose, forbw);
#else
    mergeOneByOne (rep, mat, maxlevel, verbose, forbw, ratio, coverNmax);
#endif /* M_STRATEGY <= 2 */

    gzip_close (rep->outfile, outname);
    fprintf (stderr, "Final matrix has N=%d nc=%d (%d) w(M)=%lu N*w(M)=%"
	     PRIu64"\n", mat->rem_nrows, mat->rem_ncols,
	     mat->rem_nrows - mat->rem_ncols, mat->weight,
	    (uint64_t) mat->rem_nrows * (uint64_t) mat->weight);
#ifndef USE_MARKOWITZ
    closeSWAR (mat);
#else
    MkzClose (mat);
#endif
    clearMat (mat);
    return 0;
}
