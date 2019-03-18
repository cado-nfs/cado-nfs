/* merge --- new merge program

Copyright 2019 Charles Bouillaguet and Paul Zimmermann.

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

#include "cado.h"
/* the following should come after cado.h, which sets -Werror=all */
#ifdef  __GNUC__
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>  /* for _O_BINARY */
#include <gmp.h>    /* for mpn_ior_n */

#ifdef HAVE_MALLOPT
#include <malloc.h>
#endif

/* Define MARKOWITZ to use Markowitz pivoting to estimate the fill-in of
   a merge in routine merge_cost (instead of computing the real fill-in
   when adding the row of smallest weight to the other ones. */
#define MARKOWITZ

// #define BIG_BROTHER

#ifdef BIG_BROTHER
    unsigned char *touched_columns = NULL;
#endif


/* define DEBUG if printRow or copy_matrix is needed */
// #define DEBUG

/* CBOUND_INCR is the increment on the maximal cost of merges at each step.
   Setting it to 1 is optimal in terms of matrix size, but will take a very
   long time (typically 10 times more than with CBOUND_INCR=10). */
#ifndef FOR_DL
/* Experimentally on the RSA-512 matrix, CBOUND_INCR=11 gives a matrix which
   is 0.5% larger than the matrix obtained with CBOUND_INCR=1,
   and takes only 20% more time than with CBOUND_INCR=20.
   With MARKOWITZ, CBOUND_INCR=13 gives a similar number of steps than
   without MARKOWITZ and CBOUND_INCR=11. */
#ifndef MARKOWITZ
#define CBOUND_INCR 11
#else
#define CBOUND_INCR 13
#endif
#else
/* For the p180 matrix (http://caramba.loria.fr/p180.txt), CBOUND_INCR=20
   gives a matrix which is 0.2% larger than the matrix obtained with
   CBOUND_INCR=1 (with target_density=200), and takes only 17% more time than
   with CBOUND_INCR=30.
   With MARKOWITZ, CBOUND_INCR=31 gives a similar number of steps than
   without MARKOWITZ and CBOUND_INCR=20. */
#ifndef MARKOWITZ
#define CBOUND_INCR 20
#else
#define CBOUND_INCR 31
#endif
#endif

#include "portability.h"

#include "filter_config.h"
#include "utils_with_io.h"
#include "merge_replay_matrix.h" /* for filter_matrix_t */
#include "report.h"     /* for report_t */
#include "sparse.h"
#include "mst.h"
#include "transpose.h"

/* Note about variables used in the code:
 * cwmax is the (current) maximal weight of columns that will be considered
   for a merge. It starts at cwmax=2. Once we have performed *all* 2-merges,
   we increase cwmax to 3, and at each step of the algorithm, we increase it
   by 1 (not waiting for all 3-merges to be completed).
 * cbound is the maximum (current) fill-in that is allowed for a merge
   (in fact, it is a biased value to avoid negative values, one should subtract
    BIAS from cbound to get the actual value). It starts at 0, and once all
    the 2-merges have been performed (which all give a negative fill-in, thus
    they will all be allowed), we increase cbound by CBOUND_INCR at each step
    of the algorithm (where CBOUND_INCR differs for integer factorization and
    discrete logarithm).
 * j0 means that we assume that columns of index < j0 cannot have
   weight <= cwmax. It depends on cwmax (decreases when cwmax increases).
   At the first call to compute_weights(), the values j0(cwmax=2) up to
   j0(MERGE_LEVEL_MAX) are computed once for all (since the weight of a
   column usually does not decrease, the values of j0 should remain correct
   during the algorithm, but not optimal).
   In several places we use the fact that the rows are sorted by increasing
   columns: if we start from the end, we can stop at soon as j < j0.
*/

/* 0: compute_weights
   1: compute_R
   2: compute_merges
   3: apply_merges
   4: pass
   5: renumber */
double cpu_t[6] = {0};
double wct_t[6] = {0};

static int verbose = 0; /* verbosity level */

#define MARGIN 5 /* reallocate dynamic lists with increment 1/MARGIN */

/*************************** heap structures *********************************/

#ifdef USE_HEAP
#define PAGE_SIZE (1<<18) /* seems to be optimal for RSA-512 */

typedef struct {
  char** pages;          /* list of pages */
  unsigned long size;    /* number of pages allocated */
  unsigned long current; /* pages[size-1] + current is the first free cell */
} heap_struct_t;
typedef heap_struct_t heap_t[1];

heap_t global_heap;      /* global heap */
heap_t *local_heap;      /* local heap (one per thread) */

static void MAYBE_UNUSED
heap_init (heap_t h)
{
  h->pages = malloc (1 * sizeof (char*));
  h->pages[0] = (char*) malloc (PAGE_SIZE * sizeof (char));
  h->size = 1;
  h->current = 0;
}

/* s is in bytes */
static void MAYBE_UNUSED *
heap_malloc (heap_t h, size_t s)
{
  unsigned long cur = h->current;
  if (cur + s <= PAGE_SIZE)
    {
      h->current += s;
      return h->pages[h->size - 1] + cur;
    }
  /* otherwise we allocate a new page */
  h->pages = realloc (h->pages, (h->size + 1) * sizeof (char*));
  h->pages[h->size] = malloc (PAGE_SIZE * sizeof (char));
  h->current = s;
  h->size ++;
  return h->pages[h->size - 1];
}

static void MAYBE_UNUSED
heap_clear (heap_t h)
{
  for (unsigned long i = 0; i < h->size; i++)
    free (h->pages[i]);
  free (h->pages);
}
#endif

/*****************************************************************************/

static void
declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "mat", "input purged file");
  param_list_decl_usage(pl, "out", "output history file");
  param_list_decl_usage(pl, "skip", "number of heavy columns to bury (default "
				    STR(DEFAULT_MERGE_SKIP) ")");
  param_list_decl_usage(pl, "target_density", "stop when the average row density exceeds this value"
			    " (default " STR(DEFAULT_MERGE_TARGET_DENSITY) ")");
  param_list_decl_usage(pl, "force-posix-threads", "force the use of posix threads, do not rely on platform memory semantics");
  param_list_decl_usage(pl, "path_antebuffer", "path to antebuffer program");
  param_list_decl_usage(pl, "t", "number of threads");
  param_list_decl_usage(pl, "v", "verbose mode");
}

static void
usage (param_list pl, char *argv0)
{
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
}



#ifndef FOR_DL
/* sort row[0], row[1], ..., row[n-1] in non-decreasing order */
static void
sort_relation (index_t *row, unsigned int n)
{
  unsigned int i, j;

  for (i = 1; i < n; i++)
    {
      index_t t = row[i];
      if (t < row[i-1])
	{
	  row[i] = row[i-1];
	  for (j = i - 1; j > 0 && t < row[j-1]; j--)
	    row[j] = row[j-1];
	  row[j] = t;
	}
    }
}
#endif

/* callback function called by filter_rels */
static void *
insert_rel_into_table (void *context_data, earlyparsed_relation_ptr rel)
{
  filter_matrix_t *mat = (filter_matrix_t *) context_data;
  unsigned int j = 0;
  typerow_t buf[UMAX(weight_t)]; /* rel->nb is of type weight_t */

  for (unsigned int i = 0; i < rel->nb; i++)
  {
    index_t h = rel->primes[i].h;
    mat->rem_ncols += (mat->wt[h] == 0);
    mat->wt[h] += (mat->wt[h] != UMAX(col_weight_t));
    if (h < mat->skip)
	continue; /* we skip (bury) the first 'skip' indices */
#ifdef FOR_DL
    exponent_t e = rel->primes[i].e;
    /* For factorization, they should not be any multiplicity here.
       For DL we do not want to count multiplicity in mat->wt */
    buf[++j] = (ideal_merge_t) {.id = h, .e = e};
#else
    ASSERT(rel->primes[i].e == 1);
    buf[++j] = h;
#endif
  }

#ifdef FOR_DL
  buf[0].id = j;
#else
  buf[0] = j;
#endif

  /* only count the non-skipped coefficients */
  mat->tot_weight += j;

  /* sort indices to ease row merges */
#ifndef FOR_DL
  sort_relation (&(buf[1]), j);
#else
  qsort (&(buf[1]), j, sizeof(typerow_t), cmp_typerow_t);
#endif

#ifdef USE_HEAP
  mat->rows[rel->num] = heap_malloc (global_heap, (j + 1) * sizeof (typerow_t));
#else
  mat->rows[rel->num] = mallocRow (j + 1);
#endif
  compressRow (mat->rows[rel->num], buf, j);

  return NULL;
}

static void
filter_matrix_read (filter_matrix_t *mat, const char *purgedname)
{
  uint64_t nread;
  char *fic[2] = {(char *) purgedname, NULL};

  /* read all rels */
  nread = filter_rels (fic, (filter_rels_callback_t) &insert_rel_into_table,
		       mat, EARLYPARSE_NEED_INDEX, NULL, NULL);
  ASSERT_ALWAYS(nread == mat->nrows);
  mat->rem_nrows = nread;
}

static void
print_timings (char *s, double cpu, double wct)
{
  printf ("%s %.1fs (cpu), %.1fs (wct) [cpu/wct=%.1f]\n",
	  s, cpu, wct, cpu / wct);
  fflush (stdout);
}

/* renumber the columns of mat to have consecutive indices 0..ncols-1 */
static void
renumber (filter_matrix_t *mat)
{
  double cpu = seconds (), wct = wct_seconds ();

  /* compute p[j] such that ideal j is renamed to p[j] <= j */
  /* we use 64-bit words, thus need to round to an integer multiple of 64 */
  uint64_t size = 1 + (mat->ncols - 1) / 64; /* ceil(ncols/64) */
  mat->p = malloc(64 * size * sizeof (index_t));

  /* detect non-empty columns */
  int T = omp_get_max_threads();
  uint64_t * local_p[T];
  uint64_t partial[T];

  #pragma omp parallel
  {
    /* each thread observes some rows and marks the corresponding columns in local_p */
    int tid = omp_get_thread_num();
    int T = omp_get_num_threads();
    local_p[tid] = malloc(size * sizeof(uint64_t));
    memset(local_p[tid], 0, size * sizeof(uint64_t));

    /* static scheduling is OK since rows have identical weight distribution in the input matrix */
    #pragma omp for schedule(static)
    for (uint64_t i = 0; i < mat->nrows; i++)
      for (index_t l = 1; l <= matLengthRow(mat, i); l++) {
	index_t j = matCell(mat, i, l);
	uint64_t outside = j / 64;
	uint64_t inside = j & 63;
	local_p[tid][outside] |= ((uint64_t) 1u) << inside;
      }

    /* accumulate marks over all threads */
    #pragma omp for schedule(static)
    for (uint64_t j = 0; j < size; j++)
      for (int t = 1; t < T; t++)
	local_p[0][j] |= local_p[t][j];

    /************ BEGIN PREFIX SUM ********************/

    /* count non-empty columns in each chunk */
    uint64_t weight = 0;
    #pragma omp for
    for (uint64_t j = 0; j < size; j++)
      weight += __builtin_popcountll(local_p[0][j]);
    partial[tid] = weight;

    #pragma omp barrier

    /* prefix-sum the column counts (sequentially) */
    #pragma omp master
    {
      uint64_t s = 0;
      for (int t = 0; t < T; t++) {
	uint64_t tmp = partial[t];
	partial[t] = s;
	s += tmp;
      }
      mat->rem_ncols = s;
    }

    #pragma omp barrier

    /* compute the new column indices */
    uint64_t s = partial[tid];
    #pragma omp for schedule(static)
    for (uint64_t i = 0; i < size; i++) {
      uint64_t bitmask = local_p[0][i];
      for (uint64_t j = 0; j < 64; j++) {
	mat->p[i * 64 + j] = s;
	if (bitmask & (((uint64_t) 1) << j))
	  s++;
      }
    }

    /**************** END PREFIX SUM ***************/

    /* perform the renumbering in all rows */
    #pragma omp for
    for (uint64_t i = 0; i < mat->nrows; i++)
      for (index_t l = 1; l <= matLengthRow(mat, i); l++)
	matCell(mat, i, l) = mat->p[matCell(mat, i, l)];

    free(local_p[tid]);
  } /* end of parallel section */

#ifndef FOR_DL
  free(mat->p);
#else
  /* For the discrete logarithm, we keep the inverse of p, to print the
     original columns in the history file.
     Warning: for a column j of weight 0, we have p[j] = p[j'] where
     j' is the smallest column > j of positive weight, thus we only consider
     j such that p[j] < p[j+1], or j = ncols-1. */
  for (uint64_t i = 0, j = 0; j < mat->ncols; j++)
    if (mat->p[j] == i && (j + 1 == mat->ncols || mat->p[j] < mat->p[j+1]))
      mat->p[i++] = j; /* necessarily i <= j */
#endif

  /* reset ncols */
  mat->ncols = mat->rem_ncols;
  printf ("exit renumber, ncols=%" PRIu64 "\n", mat->ncols);
  fflush (stdout);

  cpu = seconds () - cpu;
  wct = wct_seconds () - wct;
  print_timings ("   renumber took", cpu, wct);
  cpu_t[5] += cpu;
  wct_t[5] += wct;
}

/* For 1 <= w <= MERGE_LEVEL_MAX, put in jmin[w] the smallest index j such that
   mat->wt[j] = w. This routine is called only once, at the first call of
   compute_weights. */

static void
compute_jmin (filter_matrix_t *mat, index_t *jmin)
{
  /* unfortunately, reduction on array sections requires OpenMP >= 4.5,
     which is not yet THAT widespread. We work around the problem */
  index_t tjmin[omp_get_max_threads()][MERGE_LEVEL_MAX + 1];

  #pragma omp parallel /* reduction(min: jmin[1:MERGE_LEVEL_MAX]) */
  {
    int T = omp_get_num_threads();
    int tid = omp_get_thread_num();

    index_t *local = tjmin[tid];

    /* first initialize to ncols */
    for (int w = 1; w <= MERGE_LEVEL_MAX; w++)
      local[w] = mat->ncols;

    #pragma omp for schedule(static)
    for (index_t j = 0; j < mat->ncols; j++) {
      col_weight_t w = mat->wt[j];
      if (0 < w && w <= MERGE_LEVEL_MAX && j < local[w])
	local[w] = j;
    }

    #pragma omp for schedule(static)
    for (int w = 1; w <= MERGE_LEVEL_MAX; w++) {
      jmin[w] = mat->ncols;
      for (int t = 0; t < T; t++)
	if (jmin[w] > tjmin[t][w])
	  jmin[w] = tjmin[t][w];
    }
  }

  jmin[0] = 1; /* to tell that jmin was initialized */

  /* make jmin[w] = min(jmin[w'], 1 <= w' <= w) */
  for (int w = 2; w <= MERGE_LEVEL_MAX; w++)
    if (jmin[w - 1] < jmin[w])
      jmin[w] = jmin[w - 1];
}

/* compute column weights (in fact, saturate to cwmax + 1 since we only need to
   know whether the weights are <= cwmax or not) */
static void
compute_weights (filter_matrix_t *mat, index_t *jmin)
{
  double cpu = seconds (), wct = wct_seconds ();
  unsigned char cwmax = mat->cwmax;

  index_t j0;
  if (jmin[0] == 0) /* jmin was not initialized */
    {
      j0 = 0;
      cwmax = MERGE_LEVEL_MAX;
    }
  else
    /* we only need to consider ideals of index >= j0, assuming the weight of
       an ideal cannot decrease (except when decreasing to zero when merged) */
    j0 = jmin[mat->cwmax];

  col_weight_t *Wt[omp_get_max_threads()];
  #pragma omp parallel
  {
    int T = omp_get_num_threads();
    int tid = omp_get_thread_num();

    /* we allocate an array of size mat->ncols, but the first j0 entries are unused */
    if (tid == 0)
	Wt[0] = mat->wt; /* trick: we use wt for Wt[0] */
    else
	Wt[tid] = malloc (mat->ncols * sizeof (col_weight_t));
    memset (Wt[tid] + j0, 0, (mat->ncols - j0) * sizeof (col_weight_t));

    /* Thread k accumulates weights in Wt[k].
     We only consider ideals of index >= j0, and put the weight of ideal j,
     j >= j0, in Wt[k][j]. */

    /* using schedule(dynamic,128) here is crucial, since during merge,
     the distribution of row lengths is no longer uniform (including
     discarded rows) */
    col_weight_t *Wtk = Wt[tid];
    #pragma omp for schedule(dynamic, 128)
    for (index_t i = 0; i < mat->nrows; i++) {
      if (mat->rows[i] == NULL) /* row was discarded */
	continue;
      for (index_t l = matLengthRow (mat, i); l >= 1; l--) {
	index_t j = matCell (mat, i, l);
	if (j < j0) /* assume ideals are sorted by increasing order */
	  break;
	else if (Wtk[j] <= cwmax)      /* (*) HERE */
	  Wtk[j]++;
      }
    }

    /* Thread k accumulates in Wt[0] the weights for the k-th block of columns,
       saturating at cwmax + 1:
       Wt[0][j] = min(cwmax+1, Wt[0][j] + Wt[1][j] + ... + Wt[nthreads-1][j]) */
    col_weight_t *Wt0 = Wt[0];
    #pragma omp for schedule(static)
    for (index_t i = j0; i < mat->ncols; i++) {
      col_weight_t val = Wt0[i];
      for (int t = 1; t < T; t++)
	if (val + Wt[t][i] <= cwmax)
	  val += Wt[t][i];
	else {
	  val = cwmax + 1;
	  break;
	}
      Wt0[i] = val;
    }

    if (tid > 0)     /* start from 1 since Wt[0] = mat->wt + j0 should be kept */
      free (Wt[tid]);
  }

  if (jmin[0] == 0) /* jmin was not initialized */
    compute_jmin (mat, jmin);

  cpu = seconds () - cpu;
  wct = wct_seconds () - wct;
  print_timings ("   compute_weights took", cpu, wct);
  cpu_t[0] += cpu;
  wct_t[0] += wct;
}

#ifdef USE_CSR
/* computes the transposed matrix for columns of weight <= cwmax
   (we only consider columns >= j0) */
static void
compute_R (filter_matrix_t *mat, index_t j0)
{
  double cpu = seconds (), wct = wct_seconds ();

  index_t *Rp = mat->Rp;
  index_t *Rq = mat->Rq;
  index_t *Rqinv = mat->Rqinv;
  uint64_t nrows = mat->nrows;
  uint64_t ncols = mat->ncols;
  int cwmax = mat->cwmax;

  /* compute the number of rows, the indices of the rowd and the row pointers */
  int T = omp_get_max_threads();
  index_t tRnz[T];
  index_t tRn[T];
  #pragma omp parallel
  {
    int T = omp_get_num_threads();
    int tid = omp_get_thread_num();
    index_t Rnz = 0;
    index_t Rn = 0;
    #pragma omp for schedule(static) nowait
    for (index_t j = j0; j < ncols; j++) {
      col_weight_t w = mat->wt[j];
      if (0 < w && w <= cwmax) {
	Rnz += w;
	Rn++;
      }
    }
    tRnz[tid] = Rnz;
    tRn[tid] = Rn;

    #pragma omp barrier

    /* prefix-sum over the T threads (sequentially) */
    #pragma omp master
    {
      index_t r = 0;
      index_t s = 0;
      for (int t = 0; t < T; t++) {
	index_t w = tRnz[t];
	index_t n = tRn[t];
	tRnz[t] = r;
	tRn[t] = s;
	r += w;
	s += n;
      }
      mat->Rn = s;
      Rp[s] = r; /* set the last row pointer */
    }

    #pragma omp barrier
    Rnz = tRnz[tid];
    Rn = tRn[tid];

    #pragma omp for schedule(static)
    for (index_t j = j0; j < ncols; j++) {
      col_weight_t w = mat->wt[j];
      if (0 < w && w <= cwmax) {
	Rq[j] = Rn;
	Rqinv[Rn] = j;
	#ifdef TRANSPOSE_EASY_WAY
	Rnz += w;
	Rp[Rn] = Rnz;
	#else
	Rp[Rn] = Rnz;
	Rnz += w;
	#endif
	Rn++;
      }
    }
  } /* end parallel section */
  index_t Rn = mat->Rn;
  index_t Rnz = Rp[Rn];


#ifdef BIG_BROTHER
  index_t n_empty = 0;
  for (index_t j = 0; j < ncols; j++)
    if (mat->wt[j] == 0)
      n_empty++;
  printf("$$$     empty: %d\n", n_empty);
  printf("$$$     light: %d\n", Rn);
#endif

  /* extract submatrix */
  index_t *Mi = aligned_alloc(64, Rnz * sizeof(index_t));
  index_t *Mj = aligned_alloc(64, Rnz * sizeof(index_t));
  index_t *Ri = aligned_alloc(64, Rnz * sizeof(index_t));
  index_t ptr = 0;

  #pragma omp parallel
  {
	index_t tptr;
	index_t slot = 0;
	#define BUFFER_SIZE 1024
	index_t row[BUFFER_SIZE] __attribute__((__aligned__(64)));
	index_t col[BUFFER_SIZE] __attribute__((__aligned__(64)));

	#pragma omp for schedule(dynamic, 128) 
	for (index_t i = 0; i < nrows; i++) {
		if (mat->rows[i] == NULL)
			continue; /* row was discarded */
		for (index_t k = matLengthRow(mat, i); k >= 1; k--) {
			index_t j = matCell (mat, i, k);
			if (j < j0)
				break;
			if (mat->wt[j] > cwmax)
				continue;
			row[slot] = i;
			col[slot] = Rq[j];
			if (slot == BUFFER_SIZE - 1) {
				#pragma omp atomic capture
				{ tptr = ptr; ptr += BUFFER_SIZE; }
				for (int j = 0; j < BUFFER_SIZE; j += 64 / sizeof(index_t)) {
					store_nontemp_64B(Mi + tptr, row + j);
					store_nontemp_64B(Mj + tptr, col + j);
					tptr += 64 / sizeof(index_t);
				}
				slot = 0;
			} else {
				slot++;
			}
		}
	}
	/* purge buffer */
	#pragma omp atomic capture
	{ tptr = ptr; ptr += slot; }
	for (index_t i = 0; i < slot; i++) {
		Mi[tptr + i] = row[i];
		Mj[tptr + i] = col[i];
	}
  }
  ASSERT(ptr == Rnz);
  
  /* finally... */
  transpose(Rnz, Mi, Mj, Rn, Rp, Ri);
  free(Mi);
  free(Mj);

  /* save */
  mat->Rn = Rn;
  mat->Ri = Ri;
  cpu = seconds () - cpu;
  wct = wct_seconds () - wct;
  print_timings ("   compute_R took", cpu, wct);
  cpu_t[1] += cpu;
  wct_t[1] += wct;
}
#else /* List of List transpose */
static void
compute_R (filter_matrix_t *mat, index_t j0)
{
  double cpu = seconds (), wct = wct_seconds ();

  /* first allocate R[j] */
  index_t Rn = 0;
  for (index_t j = j0; j < mat->ncols; j++)
    {
      col_weight_t w = mat->wt[j];
      if (0 < w && w <= mat->cwmax)
	{
	  mat->R[j] = malloc (mat->wt[j] * sizeof (index_t));
	  mat->wt[j] = 0; /* trick: we put wt[j] to 0, it will be put back
			     to its initial value in the dispatch loop below */
	  Rn ++;
	}
    }
  mat->Rn = Rn;

  /* dispatch entries */
  col_weight_t *wt = mat->wt;
  #pragma omp parallel for schedule(dynamic, 128)
  for (index_t i = 0; i < mat->nrows; i++) {
    if (mat->rows[i] == NULL)
      continue; /* row was discarded */
    for (index_t k = matLengthRow(mat, i); k >= 1; k--) {
      index_t j = matCell (mat, i, k);
      /* since ideals are sorted by increasing value, we can stop
	 when we encounter j < j0 */
      if (j < j0)
	break;
      /* we only accumulate ideals of weight <= MERGE_LEVEL_MAX */
      if (mat->wt[j] > mat->cwmax)
	continue;
      index_t s;
      #pragma omp atomic capture
      s = wt[j]++;
      mat->R[j][s] = i;
    }
  }

  cpu = seconds () - cpu;
  wct = wct_seconds () - wct;
  print_timings ("   compute_R took", cpu, wct);
  cpu_t[1] += cpu;
  wct_t[1] += wct;
}
#endif

// #define TRACE_J 1438672

static void
decrease_weight (filter_matrix_t *mat, index_t j)
{
  /* only decrease the weight if <= MERGE_LEVEL_MAX,
     since we saturate to MERGE_LEVEL_MAX+1 */
  if (mat->wt[j] <= MERGE_LEVEL_MAX) {
    /* update is enough, we do not need capture since we are not interested
       by the value of wt[j] */
    #pragma omp atomic update
    mat->wt[j]--;
#ifdef BIG_BROTHER
    touched_columns[j] = 1;
#endif
  }
}

static void
increase_weight (filter_matrix_t *mat, index_t j)
{
  /* only increase the weight if <= MERGE_LEVEL_MAX,
     since we saturate to MERGE_LEVEL_MAX+1 */
  if (mat->wt[j] <= MERGE_LEVEL_MAX) {
    #pragma omp atomic update
    mat->wt[j]++;
#ifdef BIG_BROTHER
    touched_columns[j] = 1;
#endif
  }
}

/* doit == 0: return the weight of row i1 + row i2
   doit <> 0: add row i2 to row i1 */
#ifndef FOR_DL
/* special code for factorization */
static int32_t
add_row (filter_matrix_t *mat, index_t i1, index_t i2, int doit,
	 MAYBE_UNUSED index_t j)
{
  uint32_t k1 = matLengthRow (mat, i1);
  uint32_t k2 = matLengthRow (mat, i2);
  int32_t c = 0;
  uint32_t t1 = 1, t2 = 1;
  while (t1 <= k1 && t2 <= k2)
    {
      if (mat->rows[i1][t1] == mat->rows[i2][t2])
	      t1 ++, t2 ++;
      else if (mat->rows[i1][t1] < mat->rows[i2][t2])
	      t1 ++, c ++;
      else
	      t2 ++, c ++;
    }
  c += (k1 + 1 - t1) + (k2 + 1 - t2);
  if (doit == 0)
    return c;
  /* now perform the real merge */
  index_t *t, *t0;
#ifdef USE_HEAP
  int tid = omp_get_thread_num ();
  t = heap_malloc (local_heap[tid], (c + 1) * sizeof (index_t));
#else
  t = malloc ((c + 1) * sizeof (index_t));
#endif
  t0 = t;
  *t++ = c;
  t1 = t2 = 1;
  while (t1 <= k1 && t2 <= k2)
    {
      if (mat->rows[i1][t1] == mat->rows[i2][t2])
	{
	  decrease_weight (mat, mat->rows[i1][t1]);
	  t1 ++, t2 ++;
	}
      else if (mat->rows[i1][t1] < mat->rows[i2][t2])
	*t++ = mat->rows[i1][t1++];
      else
	{
	  increase_weight (mat, mat->rows[i2][t2]);
	  *t++ = mat->rows[i2][t2++];
	}
    }
  while (t1 <= k1)
    *t++ = mat->rows[i1][t1++];
  while (t2 <= k2)
    {
      increase_weight (mat, mat->rows[i2][t2]);
      *t++ = mat->rows[i2][t2++];
    }
  ASSERT (t0 + (c + 1) == t);
#ifndef USE_HEAP
  free (mat->rows[i1]);
#endif
  mat->rows[i1] = t0;
  return c;
}
#else /* FOR_DL: j is the ideal to be merged */
#define INT32_MIN_64 (int64_t) INT32_MIN
#define INT32_MAX_64 (int64_t) INT32_MAX

static int32_t
add_row (filter_matrix_t *mat, index_t i1, index_t i2, int doit, index_t j)
{
  /* first look for the exponents of j in i1 and i2 */
  uint32_t k1 = matLengthRow (mat, i1);
  uint32_t k2 = matLengthRow (mat, i2);
  ideal_merge_t *r1 = mat->rows[i1];
  ideal_merge_t *r2 = mat->rows[i2];
  int32_t e1 = 0, e2 = 0;

  /* search by decreasing ideals as the ideal to be merged is likely large */
  for (int l = k1; l >= 1; l--)
    if (r1[l].id == j)
      {
	e1 = r1[l].e;
	break;
      }
  for (int l = k2; l >= 1; l--)
    if (r2[l].id == j)
      {
	e2 = r2[l].e;
	break;
      }

  /* we always check e1 and e2 are not zero, in order to prevent from zero
     exponents that would come from exponent overflows in previous merges */
  ASSERT_ALWAYS (e1 != 0 && e2 != 0);

  int d = (int) gcd_int64 ((int64_t) e1, (int64_t) e2);
  e1 /= -d;
  e2 /= d;
  /* we will multiply row i1 by e2, and row i2 by e1 */

  int32_t c = 0;
  uint32_t t1 = 1, t2 = 1;

  while (t1 <= k1 && t2 <= k2)
    {
      if (r1[t1].id == r2[t2].id)
	{
	  /* If exponent do not cancel, add 1 to cost.
	     Warning: we should ensure that r1[t1].e * e2 does not overflow,
	     same for r2[t2].e * e1 and the sum.
	     In fact, since the sum is computed modulo 2^32, the only bad case
	     is when the sum is a non-zero multiple of 2^32. */
	  int32_t e = r1[t1].e * e2 + r2[t2].e * e1;
	  if (e != 0)
	    c ++; /* we are sure that the sum is not zero */
	  else
	    {
	      /* We compute the sum with 64-bit integers. Since all values are
		 in [-2^31, 2^31-1], the sum is in [-2^63,2^63-2^33+2], thus
		 always fits into an int64_t. */
	      int64_t ee = (int64_t) r1[t1].e * (int64_t) e2 + (int64_t) r2[t2].e * (int64_t) e1;
	      c += ee != 0;
	    }
	  t1 ++, t2 ++;
	}
      else if (r1[t1].id < r2[t2].id)
	t1 ++, c ++;
      else
	t2 ++, c ++;
    }
  c += (k1 + 1 - t1) + (k2 + 1 - t2);

  if (doit == 0)
    return c;

  /* now perform the real merge */
  ideal_merge_t *t, *t0;
#ifdef USE_HEAP
  int tid = omp_get_thread_num ();
  t = heap_malloc (local_heap[tid], (c + 1) * sizeof (ideal_merge_t));
#else
  t = malloc ((c + 1) * sizeof (ideal_merge_t));
#endif
  t0 = t;
  (*t++).id = c; /* length of the new relation */
  t1 = t2 = 1;
  int64_t e;
  while (t1 <= k1 && t2 <= k2)
    {
      if (r1[t1].id == r2[t2].id)
	{
	  /* as above, the exponent e below cannot overflow */
	  e = (int64_t) e2 * (int64_t) r1[t1].e + (int64_t) e1 * (int64_t) r2[t2].e;
	  if (e != 0) /* exponents do not cancel */
	    {
	      (*t).id = r1[t1].id;
	      ASSERT_ALWAYS(INT32_MIN_64 <= e && e <= INT32_MAX_64);
	      (*t).e = (int32_t) e;
	      t ++;
	    }
	  else
	    decrease_weight (mat, r1[t1].id);
	  t1 ++, t2 ++;
	}
      else if (r1[t1].id < r2[t2].id)
	{
	  (*t).id = r1[t1].id;
	  e = (int64_t) e2 * (int64_t) r1[t1].e;
	  ASSERT_ALWAYS(INT32_MIN_64 <= e && e <= INT32_MAX_64);
	  (*t).e = (int32_t) e;
	  t1 ++, t ++;
	}
      else
	{
	  (*t).id = r2[t2].id;
	  e = (int64_t) e1 * (int64_t) r2[t2].e;
	  ASSERT_ALWAYS(INT32_MIN_64 <= e && e <= INT32_MAX_64);
	  (*t).e = (int32_t) e;
	  increase_weight (mat, r2[t2].id);
	  t2 ++, t ++;
	}
    }
  while (t1 <= k1)
    {
      (*t).id = r1[t1].id;
      e = (int64_t) e2 * (int64_t) r1[t1].e;
      ASSERT_ALWAYS(INT32_MIN_64 <= e && e <= INT32_MAX_64);
      (*t).e = (int32_t) e;
      t1 ++, t ++;
    }
  while (t2 <= k2)
    {
      (*t).id = r2[t2].id;
      e = (int64_t) e1 * (int64_t) r2[t2].e;
      ASSERT_ALWAYS(INT32_MIN_64 <= e && e <= INT32_MAX_64);
      (*t).e = (int32_t) e;
      increase_weight (mat, r2[t2].id);
      t2 ++, t ++;
    }
  ASSERT (t0 + (c + 1) == t);
#ifndef USE_HEAP
  free (mat->rows[i1]);
#endif
  mat->rows[i1] = t0;
  return c;
}
#endif

static void
remove_row (filter_matrix_t *mat, index_t i)
{
  int32_t w = matLengthRow (mat, i);
  for (int k = 1; k <= w; k++)
    decrease_weight (mat, rowCell(mat->rows[i], k));
#ifndef USE_HEAP
  free (mat->rows[i]);
#endif
  /* even with USE_HEAP, we need to set rows[i] to NULL to know that this
     row was discarded */
  mat->rows[i] = NULL;
}

#ifdef DEBUG
static void MAYBE_UNUSED
printRow (filter_matrix_t *mat, index_t i)
{
  if (mat->rows[i] == NULL)
    {
      printf ("row %u has been discarded\n", i);
      return;
    }
  int32_t k = matLengthRow (mat, i);
  printf ("%lu [%d]:", (unsigned long) i, k);
  for (int j = 1; j <= k; j++)
#ifndef FOR_DL
    printf (" %lu", (unsigned long) mat->rows[i][j]);
#else
    printf (" %lu^%d", (unsigned long) mat->rows[i][j].id, mat->rows[i][j].e);
#endif
  printf ("\n");
}
#endif

/* classical cost: merge the row of smaller weight with the other ones,
   and return the merge cost (taking account of cancellations).
   id is the index of a row in R. */
static int32_t
merge_cost (filter_matrix_t *mat, index_t id)
{
#ifdef USE_CSR
  index_t lo = mat->Rp[id];
  index_t hi = mat->Rp[id + 1];
  int w = hi - lo;
#else
  index_t j = id;
  int w = mat->wt[j];
  index_t lo = 0;
  index_t hi = w;
#endif

  ASSERT (1 <= w && w <= mat->cwmax);

  if (w == 1)
    return -3; /* ensure all 1-merges are processed before 2-merges with no
		  cancellation */

  /* find shortest row in the merged column */
#ifdef USE_CSR
  index_t imin = mat->Ri[lo];
#else
  index_t imin = mat->R[id][0];
#endif
  index_t cmin = matLengthRow (mat, imin);
  for (index_t k = lo + 1; k < hi; k++)
    {
#ifdef USE_CSR
      index_t i = mat->Ri[k];
#else
      index_t i = mat->R[id][k];
#endif
      index_t c = matLengthRow(mat, i);
      if (c < cmin)
	{
		imin = i;
		cmin = c;
	      }
    }

  /* we remove row imin and add it to all w-1 others: cmin * (w - 2)
     the column j disappears: -w */
  index_t c = -cmin; /* remove row imin */
  for (index_t k = lo; k < hi; k++)
    {
#ifdef USE_CSR
      index_t i = mat->Ri[k];
#else
      index_t i = mat->R[id][k];
#endif
      if (i != imin)
	      /* It is crucial here to take into account cancellations of
		 coefficients, and not to simply add the length of both
		 rows minus 2. Indeed, if row 'a' was added to two
		 relation-sets 'b' and 'c', and 'b' and 'c' are merged together,
		 all ideals from 'a' will cancel. */
	#ifndef MARKOWITZ
	  index_t j = mat->Rqinv[id];
	  c += add_row (mat, i, imin, 0, j) - matLengthRow (mat, i);
	#else /* estimation with Markowitz pivoting: might miss cancellations */
	  c += cmin - 2;
	#endif
    }
  return c;
}

/* Output a list of merges to a string.
   Assume rep->type = 0.
   size is the length of str.
   Return the number of characters written, except the final \0
   (or that would have been written if that number >= size) */
static int
#ifndef FOR_DL
sreportn (char *str, size_t size, index_signed_t *ind, int n)
#else
sreportn (char *str, size_t size, index_signed_t *ind, int n, index_t j)
#endif
{
  size_t m = 0; /* number of characters written */

  for (int i = 0; i < n; i++)
    {
      m += snprintf (str + m, size - m, "%ld", (long int) ind[i]);
      ASSERT(m < size);
      if (i < n-1)
	{
	  m += snprintf (str + m, size - m, " ");
	  ASSERT(m < size);
	}
    }
#ifdef FOR_DL
  m += snprintf (str + m, size - m, " #%lu", (unsigned long) j);
#endif
  m += snprintf (str + m, size - m, "\n");
  ASSERT(m < size);
  return m;
}

/* Perform the row additions given by the minimal spanning tree (stored in
   history[][]). */
static int
addFatherToSons (index_t history[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1],
		 filter_matrix_t *mat, int m, index_t *ind, index_t j,
		 int *father, int *sons)
{
  int i, s, t;

  for (i = m - 2; i >= 0; i--)
    {
      s = father[i];
      t = sons[i];
      if (i == 0)
	{
	  history[i][1] = ind[s];
	  ASSERT(s == 0);
	}
      else
	history[i][1] = -(ind[s] + 1);
      add_row (mat, ind[t], ind[s], 1, j);
      history[i][2] = ind[t];
      history[i][0] = 2;
    }
  return m - 2;
}

/* perform the merge described by the id-th row of R,
   computing the full spanning tree */
static int32_t
merge_do (filter_matrix_t *mat, index_t id, FILE *out)
{
  int32_t c;
#ifdef USE_CSR
  index_t j = mat->Rqinv[id];
  index_t t = mat->Rp[id];
  int w = mat->Rp[id + 1] - t;
#else
  index_t j = id;
  int w = mat->wt[j];
#endif

  ASSERT (1 <= w && w <= mat->cwmax);

  if (w == 1)
    {
      char s[MERGE_CHAR_MAX];
      int n MAYBE_UNUSED;
#ifdef USE_CSR
      index_signed_t i = mat->Ri[t]; /* only row containing j */
#else
      index_signed_t i = mat->R[id][0]; /* only row containing j */
#endif
#ifndef FOR_DL
      n = sreportn (s, MERGE_CHAR_MAX, &i, 1);
#else
      n = sreportn (s, MERGE_CHAR_MAX, &i, 1, mat->p[j]);
#endif
      ASSERT(n < MERGE_CHAR_MAX);
      fprintf (out, "%s", s);
      remove_row (mat, i);
      return -3;
    }

  /* perform the real merge and output to history file */
#ifdef USE_CSR
  index_t *ind = mat->Ri + t;
#else
  index_t *ind = mat->R[id];
#endif
  char s[MERGE_CHAR_MAX];
  int n = 0; /* number of characters written to s (except final \0) */
  int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX];
  fillRowAddMatrix (A, mat, w, ind, j);
  /* mimic MSTWithA */
  int start[MERGE_LEVEL_MAX], end[MERGE_LEVEL_MAX];
  c = minimalSpanningTree (start, end, w, A);
  /* c is the weight of the minimal spanning tree, we have to remove
     the weights of the initial relations */
  for (int k = 0; k < w; k++)
    c -= matLengthRow (mat, ind[k]);
  index_t history[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1];
  int hmax = addFatherToSons (history, mat, w, ind, j, start, end);
  for (int i = hmax; i >= 0; i--)
    {
#ifndef FOR_DL
      n += sreportn (s + n, MERGE_CHAR_MAX - n,
		     (index_signed_t*) (history[i]+1), history[i][0]);
#else
      n += sreportn (s + n, MERGE_CHAR_MAX - n,
		     (index_signed_t*) (history[i]+1), history[i][0],
		     mat->p[j]);
#endif
      ASSERT(n < MERGE_CHAR_MAX);
    }
  fprintf (out, "%s", s);
  remove_row (mat, ind[0]);
  return c;
}

/* since merge costs might be negative (for weight 2), we translate them by 3,
   so that 2-merges with no cancellation give biased cost -2+3=1, and those
   with cancellation give biased cost 0, so they will be merged first */
#define BIAS 3


/* accumulate in L all merges of (biased) cost <= cbound.
   L must be preallocated.
   L is a linear array and the merges appear by increasing cost.
   Returns the size of L. */
static int
compute_merges (index_t *L, filter_matrix_t *mat, int cbound)
{
  index_t Rn = mat->Rn;
  int * cost = malloc(Rn * sizeof(*cost));
  int T = omp_get_max_threads();
  index_t count[T][cbound + 1];
  // int Lp[cbound + 2];  cost pointers

  /* compute the cost of all candidate merges */
#ifdef USE_CSR
  #pragma omp parallel for schedule(dynamic, 64)   /* TODO: dynamic really necessary? */
  for (index_t i = 0; i < Rn; i++)
    cost[i] = merge_cost (mat, i) + BIAS;
#else
  index_t i = 0;
  index_t *Rqinv = malloc (Rn * sizeof (index_t));
  for (index_t j = 0; j < mat->ncols; j++)
    if (0 < mat->wt[j] && mat->wt[j] <= mat->cwmax)
      {
	Rqinv[i] = j;
	cost[i++] = merge_cost (mat, j) + BIAS;
      }
  ASSERT_ALWAYS(i == Rn);
#endif

  /* Yet Another Bucket Sort (sigh) : sort the candidate merges by cost. Check if worth parallelizing */
  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    index_t *tcount = &count[tid][0];

    memset(tcount, 0, (cbound + 1) * sizeof(index_t));

    #pragma omp for
    for (index_t i = 0; i < Rn; i++) {
	int c = cost[i];
	if (c <= cbound)
	  tcount[c]++;
    }
  } /* end parallel section */

  /* prefix-sum */
  int s = 0;
  for (int c = 0; c <= cbound; c++) {
     // Lp[c] = s;                     /* global row pointer in L */
     for (int t = 0; t < T; t++) {
	index_t w = count[t][c];       /* per-thread row pointer in L */
	count[t][c] = s;
	s += w;
     }
  }

 #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    index_t *tcount = &count[tid][0];

    #pragma omp for
    for (index_t i = 0; i < Rn; i++) {
      int c = cost[i];
      if (c > cbound)
	continue;
#ifdef USE_CSR
      L[tcount[c]++] = i;
#else
      L[tcount[c]++] = Rqinv[i];
#endif
    }
  } /* end parallel section */

#ifdef BIG_BROTHER
  printf("$$$     cheap: %d\n", s);
#endif

  free(cost);
#ifndef USE_CSR
  free (Rqinv);
#endif
  return s;
}


/* return the number of merges applied */
#ifdef USE_CSR
static unsigned long
apply_merges (index_t *L, index_t total_merges, filter_matrix_t *mat, FILE *out)
{
  int size = 1 + mat->nrows / 64;
  uint64_t * busy_rows = malloc(size * sizeof(*busy_rows));
  memset(busy_rows, 0, sizeof(uint64_t) * size);

  unsigned long nmerges = 0;
  unsigned long discarded = 0;
  unsigned long eventually_discarded = 0;
  int64_t fill_in = 0;
  double contention = 0;

  #pragma omp parallel reduction(+: fill_in, nmerges, discarded, eventually_discarded) reduction(max: contention)
  {
    contention = 0;
    #pragma omp for schedule(dynamic, 16)
    for (index_t it = 0; it < total_merges; it++) {
      index_t id = L[it];
      index_t lo = mat->Rp[id];
      index_t hi = mat->Rp[id + 1];

      /* merge is possible if all its rows are "available" */
      int ok = 1;
      for (index_t k = lo; k < hi; k++) {
	index_t i = mat->Ri[k];
	uint64_t x = i / 64;
	uint64_t y = i & 63;
	if (busy_rows[x] & (1ull << y)) {
	  ok = 0;
	  break;
	}
      }
      if (!ok)
	discarded++;
      if (ok) {
	double start = wct_seconds();
	#pragma omp critical
	{ /* potential merge, enter critical section */
	  /* check again, since another thread might have reserved a row */
	  contention += wct_seconds() - start;
	  for (index_t k = lo; k < hi; k++) {
	    index_t i = mat->Ri[k];
	    uint64_t x = i / 64;
	    uint64_t y = i & 63;
	    if (busy_rows[x] & (1ull << y)) {
	      ok = 0;
	      break;
	    }
	  }
	  if (!ok)
	    eventually_discarded++;
	  if (ok) /* reserve rows */
	    for (index_t k = lo; k < hi; k++) {
	      index_t i = mat->Ri[k];
	      uint64_t x = i / 64;
	      uint64_t y = i & 63;
	      busy_rows[x] |= (1ull << y);
	    }
	} /* end critical */
      }
      if (ok) {
	fill_in += merge_do(mat, id, out);
	nmerges ++;
      }
    }
  } /* end parallel section */


  mat->tot_weight += fill_in;
  /* each merge decreases the number of rows and columns by one */
  mat->rem_nrows -= nmerges;
  mat->rem_ncols -= nmerges;

  /* settings for next pass */
  if (mat->cwmax == 2) /* we first process all 2-merges */
    {
      if (nmerges == total_merges)
	mat->cwmax++;
    }
  else
    {
      if (mat->cwmax < MERGE_LEVEL_MAX)
	mat->cwmax ++;
    }

#ifdef BIG_BROTHER
  printf("$$$     discarded: %ld\n", discarded);
  printf("$$$     eventually_discarded: %ld\n", eventually_discarded);
  printf("$$$     merged: %ld\n", nmerges);
  printf("$$$     max-contention: %.2fs\n", contention);
  index_t n_rows = 0;
  for (int i = 0; i < size; i++)
    n_rows += __builtin_popcountll(busy_rows[i]);
  printf("$$$     affected-rows: %d\n", n_rows);
#endif

  free(busy_rows);
  return nmerges;
}
#else
static unsigned long
apply_merges (index_t *L, index_t total_merges, filter_matrix_t *mat, FILE *out)
{
  int size = 1 + mat->nrows / 64;
  uint64_t * busy_rows = malloc(size * sizeof(*busy_rows));
  memset(busy_rows, 0, sizeof(uint64_t) * size);

  index_t *T; /* todo list */
  T = malloc (total_merges * sizeof (index_t));
  unsigned long nmerges = 0;

  /* When the transposed matrix is in LIL-format, we need to first compute all
     independent merges before applying them. The reason is that wt[j] is used
     to know how many relations are involved in a merge: since we update wt[j]
     incrementally, it might be wrong when we consider a potential merge. */

  #pragma omp parallel
  {
    #pragma omp for schedule(dynamic, 16)
    for (index_t it = 0; it < total_merges; it++) {
      index_t id = L[it];
      index_t lo = 0;
      index_t j = id;
      index_t hi = mat->wt[j];

      /* merge is possible if all its rows are "available" */
      int ok = 1;
      for (index_t k = lo; k < hi; k++) {
	index_t i = mat->R[id][k];
	uint64_t x = i / 64;
	uint64_t y = i & 63;
	if (busy_rows[x] & (1ull << y)) {
	  ok = 0;
	  break;
	}
      }
      if (ok) {
	#pragma omp critical
	{ /* potential merge, enter critical section */
	  /* check again, since another thread might have reserved a row */
	  for (index_t k = lo; k < hi; k++) {
	    index_t i = mat->R[id][k];
	    uint64_t x = i / 64;
	    uint64_t y = i & 63;
	    if (busy_rows[x] & (1ull << y)) {
	      ok = 0;
	      break;
	    }
	  }
	  if (ok) /* reserve rows */
	    for (index_t k = lo; k < hi; k++) {
	      index_t i = mat->R[id][k];
	      uint64_t x = i / 64;
	      uint64_t y = i & 63;
	      busy_rows[x] |= (1ull << y);
	    }
	} /* end critical */
      }
      if (ok) /* put merge in T */
	{
	  index_t s;
	  #pragma omp atomic capture
	  s = nmerges ++;
	  T[s] = id;
	}
    }
  } /* end parallel section */

  /* now we apply the merges in T */
  int64_t fill_in = 0;
  #pragma omp parallel for schedule(dynamic, 16)
  for (index_t k = 0; k < nmerges; k++)
    fill_in += merge_do (mat, T[k], out);

  mat->tot_weight += fill_in;
  /* each merge decreases the number of rows and columns by one */
  mat->rem_nrows -= nmerges;
  mat->rem_ncols -= nmerges;

  /* settings for next pass */
  if (mat->cwmax == 2) /* we first process all 2-merges */
    {
      if (nmerges == total_merges)
	mat->cwmax++;
    }
  else
    {
      if (mat->cwmax < MERGE_LEVEL_MAX)
	mat->cwmax ++;
    }

  free (busy_rows);
  free (T);
  return nmerges;
}
#endif

static double
average_density (filter_matrix_t *mat)
{
  return (double) mat->tot_weight / (double) mat->rem_nrows;
}

#ifdef DEBUG
/* duplicate the matrix, where the lines of mat_copy are
   not allocated individually by malloc, but all at once */
static void MAYBE_UNUSED
copy_matrix (filter_matrix_t *mat)
{
  unsigned long weight = mat->tot_weight;
  unsigned long s = weight + mat->nrows;
  index_t *T = malloc (s * sizeof (index_t));
  index_t *p = T;
  double cpu = seconds (), wct = wct_seconds ();
  for (index_t i = 0; i < mat->nrows; i++)
    {
      if (mat->rows[i] == NULL)
	{
	  p[0] = 0;
	  p ++;
	}
      else
	{
	  memcpy (p, mat->rows[i], (mat->rows[i][0] + 1) * sizeof (index_t));
	  p += mat->rows[i][0] + 1;
	}
    }
  print_timings ("   copy_matrix took", seconds () - cpu,
		 wct_seconds () - wct);
  ASSERT_ALWAYS(p == T + s);
  free (T);
}
#endif

int
main (int argc, char *argv[])
{
    char *argv0 = argv[0];

    filter_matrix_t mat[1];
    report_t rep[1];

    int nthreads = 1;
    uint32_t skip = DEFAULT_MERGE_SKIP;
    double target_density = DEFAULT_MERGE_TARGET_DENSITY;

#ifdef HAVE_MINGW
    _fmode = _O_BINARY;     /* Binary open for all files */
#endif

    double tt;
    double cpu0 = seconds ();
    double wct0 = wct_seconds ();
    param_list pl;
    param_list_init (pl);
    declare_usage(pl);
    argv++,argc--;

    param_list_configure_switch (pl, "-v", &verbose);
    param_list_configure_switch(pl, "force-posix-threads", &filter_rels_force_posix_threads);

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

    const char *purgedname = param_list_lookup_string (pl, "mat");
    const char *outname = param_list_lookup_string (pl, "out");
    const char *path_antebuffer = param_list_lookup_string(pl, "path_antebuffer");

    param_list_parse_int (pl, "t", &nthreads);
#ifdef HAVE_OPENMP
    omp_set_num_threads (nthreads);
#endif

#if defined(HAVE_MALLOPT) && !defined(USE_HEAP) && !defined(HAVE_TCMALLOC)
#define USE_ARENAS
    /* experimentally, setting the number of arenas to twice the number of
       threads seems optimal (man mallopt says it should match the number of
       threads) */
    int arenas = 2 * nthreads;
    mallopt (M_ARENA_MAX, arenas);
#endif

    param_list_parse_uint (pl, "skip", &skip);

    param_list_parse_double (pl, "target_density", &target_density);

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

    set_antebuffer_path (argv0, path_antebuffer);

    /* Read number of rows and cols on first line of purged file */
    purgedfile_read_firstline (purgedname, &(mat->nrows), &(mat->ncols));

#if (__SIZEOF_INDEX__ == 4)
    if (mat->nrows >> 32)
      {
	fprintf (stderr, "Error, nrows = %" PRIu64 " larger than 2^32, please recompile with -D__SIZEOF_INDEX__=8\n", mat->nrows);
	exit (EXIT_FAILURE);
      }
    if (mat->ncols >> 32)
      {
	fprintf (stderr, "Error, ncols = %" PRIu64 " larger than 2^32, please recompile with -D__SIZEOF_INDEX__=8\n", mat->ncols);
	exit (EXIT_FAILURE);
      }
#endif

    /* initialize rep (i.e., mostly opens outname) and write matrix dimension */
    rep->type = 0;
    rep->outfile = fopen_maybe_compressed (outname, "w");
    ASSERT_ALWAYS(rep->outfile != NULL);

    /* some explanation about the history file */
    fprintf (rep->outfile, "# Every line starting with # is ignored.\n");
    fprintf (rep->outfile, "# A line i1 i2 ... ik means that row i1 ");
    fprintf (rep->outfile, "is added to i2, ..., ik, and row i1\n");
    fprintf (rep->outfile, "# is removed afterwards ");
    fprintf (rep->outfile, "(where row 0 is the first line in *.purged.gz).\n");
#ifdef FOR_DL
    fprintf (rep->outfile, "# A line ending with #j ");
    fprintf (rep->outfile, "means that ideal of index j should be merged.\n");
#endif

    /* initialize the matrix structure */
    initMat (mat, skip);

    /* we bury the 'skip' ideals of smallest index */
    mat->skip = skip;

#ifdef USE_HEAP
    /* allocate heaps */
    heap_init (global_heap);
    local_heap = malloc (nthreads * sizeof (heap_t));
    for (int i = 0; i < nthreads; i++)
      heap_init (local_heap[i]);
#endif

    /* Read all rels and fill-in the mat structure */
    tt = seconds ();
    filter_matrix_read (mat, purgedname);
    printf ("Time for filter_matrix_read: %2.2lfs\n", seconds () - tt);

    double cpu_after_read = seconds ();
    double wct_after_read = wct_seconds ();

    renumber (mat);

#ifdef USE_CSR
    /* Allocate the transposed matrix R in CSR format. Since Rp is of fixed
       size, we allocate it for once. However, the size of Ri will vary from
       step to step. */
    mat->Rp = malloc ((mat->ncols + 1) * sizeof (index_t));
    mat->Ri = NULL;
    mat->Rq = malloc (mat->ncols * sizeof (index_t));
    mat->Rqinv = malloc (mat->ncols * sizeof (index_t));
#else
    mat->R = malloc (mat->ncols * sizeof (index_t*));
    for (index_t j = 0; j < mat->ncols; j++)
      mat->R[j] = NULL;
#endif

#ifdef BIG_BROTHER
    touched_columns = malloc(mat->ncols * sizeof(*touched_columns));
    memset(touched_columns, 0, mat->ncols * sizeof(*touched_columns));
#endif

    printf ("Using MERGE_LEVEL_MAX=%d, CBOUND_INCR=%d",
	    MERGE_LEVEL_MAX, CBOUND_INCR);
#ifdef USE_ARENAS
    printf (", M_ARENA_MAX=%d", arenas);
#endif
#ifdef HAVE_TCMALLOC
    printf (", HAVE_TCMALLOC");
#endif
#ifdef USE_HEAP
    printf (", USE_HEAP(PAGE_SIZE=%d)", PAGE_SIZE);
#endif
    printf ("\n");

    printf ("N=%" PRIu64 " W=%" PRIu64 " W/N=%.2f cpu=%.1fs wct=%.1fs mem=%luM\n",
	    mat->rem_nrows, mat->tot_weight, average_density (mat),
	    seconds () - cpu0, wct_seconds () - wct0,
	    PeakMemusage () >> 10);
#ifdef BIG_BROTHER
    printf("$$$ N: %" PRId64 "\n", mat->nrows);
    printf("$$$ start:\n");
#endif

    fflush (stdout);

    mat->cwmax = 2;

    /* jmin[w] for 1 <= w <= MERGE_LEVEL_MAX is the smallest column of weight w
       at beginning. We set jmin[0] to 0 to tell that jmin[] was not
       initialized. */
    index_t jmin[MERGE_LEVEL_MAX + 1] = {0,};

    // copy_matrix (mat);

    unsigned long lastN, lastW;
    double lastWoverN;
    int cbound = BIAS; /* bound for the (biased) cost of merges to apply */
    int pass = 0;
    while (1)
      {

	double cpu1 = seconds (), wct1 = wct_seconds ();

	pass ++;

	/* Once cwmax >= 3, tt each pass, we increase cbound to allow more
	   merges. If one decreases CBOUND_INCR, the final matrix will be
	   smaller, but merge will take more time.
	   If one increases CBOUND_INCR, merge will be faster, but the final
	   matrix will be larger. */
	if (mat->cwmax > 2)
	  cbound += CBOUND_INCR;

	lastN = mat->rem_nrows;
	lastW = mat->tot_weight;
	lastWoverN = (double) lastW / (double) lastN;

#ifdef TRACE_J
	for (index_t i = 0; i < mat->ncols; i++)
	  {
	    if (mat->rows[i] == NULL)
	      continue;
	    for (index_t k = 1; k <= matLengthRow(mat, i); k++)
	      if (mat->rows[i][k] == TRACE_J)
		printf ("ideal %d in row %lu\n", TRACE_J, (unsigned long) i);
	  }
#endif

#ifdef BIG_BROTHER
    printf("$$$   - pass: %d\n", pass);
    printf("$$$     cwmax: %d\n", mat->cwmax);
    printf("$$$     cbound: %d\n", cbound);
#endif


	/* we only compute the weights at pass 1, afterwards they will be
	   updated at each merge */
	if (pass == 1)
	  compute_weights (mat, jmin);

	compute_R (mat, jmin[mat->cwmax]);

	double cpu2 = seconds (), wct2 = wct_seconds ();

	index_t *L = malloc(mat->Rn * sizeof(index_t));
	index_t n_possible_merges = compute_merges(L, mat, cbound);

	cpu2 = seconds () - cpu2;
	wct2 = wct_seconds () - wct2;
	if (verbose > 0)
		printf("*** compute_merges: %" PRIu64 " candidate merges\n", (uint64_t) n_possible_merges);
	print_timings ("   compute_merges took", cpu2, wct2);
	cpu_t[2] += cpu2;
	wct_t[2] += wct2;

	double cpu3 = seconds (), wct3 = wct_seconds ();

	unsigned long nmerges = apply_merges (L, n_possible_merges, mat, rep->outfile);

	cpu3 = seconds () - cpu3;
	wct3 = wct_seconds () - wct3;
	print_timings ("   apply_merges took", cpu3, wct3);
	cpu_t[3] += cpu3;
	wct_t[3] += wct3;

	free(L);
#ifdef USE_CSR
	free(mat->Ri);
#else
	for (index_t j = 0; j < mat->ncols; j++)
	  {
	    free (mat->R[j]);
	    mat->R[j] = NULL;
	  }
#endif

	cpu1 = seconds () - cpu1;
	wct1 = wct_seconds () - wct1;
	print_timings ("   pass took", cpu1, wct1);
	cpu_t[4] += cpu1;
	wct_t[4] += wct1;


#ifdef BIG_BROTHER
  int n_cols = 0;
  for (unsigned int j = 0; j < mat->ncols; j++) {
    n_cols += touched_columns[j];
    touched_columns[j] = 0;
  }
  printf("$$$     affected-columns: %d\n", n_cols);
#endif

	/* estimate current average fill-in */
	double av_fill_in = ((double) mat->tot_weight - (double) lastW)
	  / (double) (lastN - mat->rem_nrows);

	printf ("N=%" PRIu64 " W=%" PRIu64 " (%.0fMB) W/N=%.2f fill-in=%.2f cpu=%.1fs wct=%.1fs mem=%luM [pass=%d,cwmax=%d]\n",
		mat->rem_nrows, mat->tot_weight,
		9.5367431640625e-07 * (mat->rem_nrows + mat->tot_weight) * sizeof(index_t),
		(double) mat->tot_weight / (double) mat->rem_nrows, av_fill_in,
		seconds () - cpu0, wct_seconds () - wct0,
		PeakMemusage () >> 10, pass, mat->cwmax);
	fflush (stdout);

	if (average_density (mat) >= target_density)
	  break;

	if (nmerges == 0 && mat->cwmax == MERGE_LEVEL_MAX)
	  break;
      }

    fclose_maybe_compressed (rep->outfile, outname);

    if (average_density (mat) > target_density)
      {
	/* estimate N for W/N = target_density, assuming W/N = a*N + b */
	unsigned long N = mat->rem_nrows;
	double WoverN = (double) mat->tot_weight / (double) N;
	double a = (lastWoverN - WoverN) / (double) (lastN - N);
	double b = WoverN - a * (double) N;
	/* we want target_density = a*N_target + b */
	printf ("Estimated N=%" PRIu64 " for W/N=%.2f\n",
		(uint64_t) ((target_density - b) / a), target_density);
      }

#ifndef FOR_DL /* we don't do renumbering for DL */
    print_timings ("renumber       :", cpu_t[5], wct_t[5]);
#endif
    print_timings ("compute_weights:", cpu_t[0], wct_t[0]);
    print_timings ("compute_R      :", cpu_t[1], wct_t[1]);
    print_timings ("compute_merges :", cpu_t[2], wct_t[2]);
    print_timings ("apply_merges   :", cpu_t[3], wct_t[3]);
    print_timings ("pass           :", cpu_t[4], wct_t[4]);

    printf ("Final matrix has N=%" PRIu64 " nc=%" PRIu64 " (%" PRIu64
	    ") W=%" PRIu64 "\n", mat->rem_nrows, mat->rem_ncols,
	    mat->rem_nrows - mat->rem_ncols, mat->tot_weight);
    fflush (stdout);

    printf ("Before cleaning memory:\n");
    print_timing_and_memory (stdout, cpu_after_read, wct_after_read);

#ifdef FOR_DL
    free (mat->p);
#endif
#ifdef USE_CSR
    free (mat->Rp);
    free (mat->Rq);
    free (mat->Rqinv);
#else
    for (index_t j = 0; j < mat->ncols; j++)
      free (mat->R[j]);
    free (mat->R);
#endif
    clearMat (mat);

#ifdef USE_HEAP
    /* free heaps */
    heap_clear (global_heap);
    for (int i = 0; i < nthreads; i++)
      heap_clear (local_heap[i]);
    free (local_heap);
#endif

    param_list_clear (pl);

    printf ("After cleaning memory:\n");
    print_timing_and_memory (stdout, cpu_after_read, wct_after_read);

    /* print total time and memory (including reading the input matrix,
       initializing and free-ing all data) */
    print_timing_and_memory (stdout, cpu0, wct0);

    return 0;
}
