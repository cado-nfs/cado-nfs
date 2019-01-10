/* merge2 --- new merge program

Copyright 2019 Paul Zimmermann.

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
#include <fcntl.h>  /* for _O_BINARY */
#include <string.h> /* for strcmp */
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

// #define USE_MST
// #define DEBUG
// #define MEM

#include "portability.h"

#include "filter_config.h"
#include "utils_with_io.h"
#include "merge_replay_matrix.h" /* for filter_matrix_t */
#include "report.h"     /* for report_t */
#include "markowitz.h" /* for MkzInit */
#include "merge_mono.h" /* for mergeOneByOne */
#include "sparse.h"
#include "mst.h"

static void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "mat", "input purged file");
  param_list_decl_usage(pl, "out", "output history file");
  param_list_decl_usage(pl, "keep", "excess to keep (default "
                                    STR(DEFAULT_FILTER_EXCESS) ")");
  param_list_decl_usage(pl, "skip", "number of heavy columns to bury (default "
                                    STR(DEFAULT_MERGE_SKIP) ")");
  param_list_decl_usage(pl, "target_density", "stop when the average row density exceeds this value"
                            " (default " STR(DEFAULT_MERGE_TARGET_DENSITY) ")");
  param_list_decl_usage(pl, "v", "verbose level");
  param_list_decl_usage(pl, "t", "number of threads");
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
insert_rel_into_table2 (void *context_data, earlyparsed_relation_ptr rel)
{
  filter_matrix_t *mat = (filter_matrix_t *) context_data;
  unsigned int j = 0;
  typerow_t buf[UMAX(weight_t)]; /* rel->nb is of type weight_t */

  for (unsigned int i = 0; i < rel->nb; i++)
  {
    index_t h = rel->primes[i].h;
    mat->rem_ncols += (mat->wt[h] == 0);
    mat->wt[h] += (mat->wt[h] != SMAX(int32_t));
    if (h < mat->force_bury_below_index) continue;
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

  /* only count the non-buried coefficients */
  mat->tot_weight += j;

  /* sort indices to ease row merges */
#ifndef FOR_DL
  sort_relation (&(buf[1]), j);
#else
  qsort (&(buf[1]), j, sizeof(typerow_t), cmp_typerow_t);
#endif

  mat->rows[rel->num] = mallocRow (j + 1);
  compressRow (mat->rows[rel->num], buf, j);

  return NULL;
}

static void
filter_matrix_read2 (filter_matrix_t *mat, const char *purgedname)
{
  uint64_t nread;
  char *fic[2] = {(char *) purgedname, NULL};

  /* read all rels */
  nread = filter_rels (fic, (filter_rels_callback_t) &insert_rel_into_table2, mat,
		       EARLYPARSE_NEED_INDEX, NULL, NULL);
  ASSERT_ALWAYS(nread == mat->nrows);
  mat->rem_nrows = nread;
}

/********************** level-1 buckets (to compute weights) *************************/

typedef struct
{
  index_t *list;
  unsigned long alloc;
  unsigned long size;
} bucket1_t;

#define MARGIN 5 /* reallocate with increment 1/MARGIN */

/* bucket B[k][j] contains buckets for thread k, with ideals = j mod nthreads */
static bucket1_t**
init_buckets1 (int nthreads, unsigned long tot_weight)
{
  bucket1_t **B;
  unsigned long w = tot_weight / (nthreads * nthreads);

  w += w / MARGIN;
  B = malloc (nthreads * sizeof (bucket1_t*));
  for (int k = 0; k < nthreads; k++)
    {
      B[k] = malloc (nthreads * sizeof (bucket1_t));
      for (int j = 0; j < nthreads; j++)
	{
	  B[k][j].list = malloc (w * sizeof (index_t));
	  B[k][j].alloc = w;
	  B[k][j].size = 0;
	}
    }
  return B;
}

static void
clear_buckets1 (bucket1_t **B, int nthreads)
{
#ifdef MEM
  unsigned long mem = 0;
#endif
  for (int k = 0; k < nthreads; k++)
    {
      for (int j = 0; j < nthreads; j++)
	{
	  free (B[k][j].list);
#ifdef MEM
	  mem += B[k][j].alloc * sizeof (index_t);
#endif
	}
      free (B[k]);
#ifdef MEM
      mem += nthreads * sizeof (bucket1_t);
#endif
    }
  free (B);
#ifdef MEM
  mem += nthreads * sizeof (bucket1_t*);
  printf ("****** init_buckets1 allocated %luM\n", mem >> 20);
#endif
}

static void
add_bucket1 (bucket1_t *Bi, index_t j, int nthreads)
{
  int jj = j % nthreads;

  if (Bi[jj].size == Bi[jj].alloc)
    {
#ifdef DEBUG
      printf ("reallocate Bi[%d] from %lu", jj, Bi[jj].alloc);
#endif
      Bi[jj].alloc += 1 + Bi[jj].alloc / MARGIN;
#ifdef DEBUG
      printf (" to %lu\n", Bi[jj].alloc);
#endif
      Bi[jj].list = realloc (Bi[jj].list, Bi[jj].alloc * sizeof (index_t));
    }
  ASSERT(Bi[jj].size < Bi[jj].alloc);
  Bi[jj].list[Bi[jj].size] = j;
  Bi[jj].size ++;
}

/* add ideals j for row i, where Bi = B[i % nthreads] */
static void
fill_buckets1 (bucket1_t *Bi, filter_matrix_t *mat, index_t i, int nthreads)
{
  if (mat->rows[i] == NULL) /* row was discarded */
    return;
  for (uint32_t k = 1; k <= matLengthRow (mat, i); k++)
    {
      index_t j = matCell (mat, i, k);
      add_bucket1 (Bi, j, nthreads);
    }
}

static void
fill_buckets1_all (bucket1_t **B, filter_matrix_t *mat, index_t i, int nthreads)
{
  bucket1_t *Bi = B[i % nthreads];
  while (i < mat->nrows)
    {
      fill_buckets1 (Bi, mat, i, nthreads);
      i += nthreads;
    }
}

/* Apply all buckets B[i][k] to compute the column weights. Since B[i][k] contains
   ideals j such that j = k mod nthreads, and two threads use a different value of k,
   this ensures there cannot be any conflict in incrementing mat->wt[j].
   We saturate to cwmax+1, since we only need to know iff the weight is <= cwmax. */
static void
apply_buckets1_all (bucket1_t **B, filter_matrix_t *mat, index_t k, int nthreads)
{
  ASSERT(mat->cwmax < INT32_MAX);
  for (int i = 0; i < nthreads; i++)
    {
      bucket1_t Bik = B[i][k];
      for (unsigned long t = 0; t < Bik.size; t++)
	{
	  index_t j = Bik.list[t];
	  if (mat->wt[j] <= mat->cwmax)
	    mat->wt[j] ++;
	}
    }
}

/*************** level-2 buckets (to compute inverse matrix) *************************/

typedef struct
{
  index_t j; /* ideal index */
  index_t i; /* row index */
} index_pair_t;

typedef struct
{
  index_pair_t *list;
  unsigned long alloc;
  unsigned long size;
} bucket_t;

/* bucket B[k][j] contains buckets for thread k, with ideals = j mod nthreads */
static bucket_t**
init_buckets (int nthreads, unsigned long tot_weight)
{
  bucket_t **B;
  unsigned long w = tot_weight / (nthreads * nthreads);

  w += w / MARGIN;
  B = malloc (nthreads * sizeof (bucket_t*));
  for (int k = 0; k < nthreads; k++)
    {
      B[k] = malloc (nthreads * sizeof (bucket_t));
      for (int j = 0; j < nthreads; j++)
	{
	  B[k][j].list = malloc (w * sizeof (index_pair_t));
	  B[k][j].alloc = w;
	  B[k][j].size = 0;
	}
    }
  return B;
}

static void
clear_buckets (bucket_t **B, int nthreads)
{
#ifdef MEM
  unsigned long mem = 0;
#endif
  for (int k = 0; k < nthreads; k++)
    {
      for (int j = 0; j < nthreads; j++)
	{
	  free (B[k][j].list);
#ifdef MEM
	  mem += B[k][j].alloc * sizeof (index_pair_t);
#endif
	}
      free (B[k]);
#ifdef MEM
      mem += nthreads * sizeof (bucket_t);
#endif
    }
  free (B);
#ifdef MEM
  mem += nthreads * sizeof (bucket_t*);
  printf ("****** init_buckets allocated %luM\n", mem >> 20);
#endif
}

static void
add_bucket (bucket_t *Bi, index_t i, index_t j, int nthreads)
{
  int jj = j % nthreads;

  if (Bi[jj].size == Bi[jj].alloc)
    {
#ifdef DEBUG
      printf ("reallocate B[%d][%d] from %lu", (int) (i % nthreads), jj, Bi[jj].alloc);
#endif
      Bi[jj].alloc += 1 + Bi[jj].alloc / MARGIN;
#ifdef DEBUG
      printf (" to %lu\n", Bi[jj].alloc);
#endif
      Bi[jj].list = realloc (Bi[jj].list, Bi[jj].alloc * sizeof (index_pair_t));
    }
  ASSERT(Bi[jj].size < Bi[jj].alloc);
  Bi[jj].list[Bi[jj].size].i = i;
  Bi[jj].list[Bi[jj].size].j = j;
  Bi[jj].size ++;
}

/* add pairs (j,i) for row i */
static void
fill_buckets (bucket_t *Bi, filter_matrix_t *mat, index_t i, int nthreads)
{
  if (mat->rows[i] == NULL) /* row was discarded */
    return;
  for (uint32_t k = 1; k <= matLengthRow (mat, i); k++)
    {
      index_t j = matCell (mat, i, k);
      /* we only consider ideals of weight <= cwmax */
      if (mat->wt[j] <= mat->cwmax)
	add_bucket (Bi, i, j, nthreads);
    }
}

static void
fill_buckets_all (bucket_t **B, filter_matrix_t *mat, index_t i, int nthreads)
{
  bucket_t *Bi = B[i % nthreads];
  while (i < mat->nrows)
    {
      fill_buckets (Bi, mat, i, nthreads);
      i += nthreads;
    }
}

static void
add_R_entry (filter_matrix_t *mat, index_t j, index_t i)
{
  int k = mat->wt[j];

  mat->R[j][k] = i;
  mat->wt[j] = k + 1;
}

/* Apply all buckets B[i][k] to compute the inverse matrix. Since B[i][k] contains
   ideals j such that j = k mod nthreads, and two threads use a different value of k,
   this ensures there cannot be any conflict in incrementing mat->wt[j]. */
static void
apply_buckets_all (bucket_t **B, filter_matrix_t *mat, index_t k, int nthreads)
{
  ASSERT(mat->cwmax < INT32_MAX);
  for (int i = 0; i < nthreads; i++)
    {
      bucket_t Bik = B[i][k];
      for (unsigned long t = 0; t < Bik.size; t++)
	{
	  index_t j = Bik.list[t].j;
	  if (mat->wt[j] <= mat->cwmax)
	    add_R_entry (mat, j, Bik.list[t].i);
	}
    }
}

/* compute column weights (in fact, saturate to cwmax + 1 since we only need to
   know whether the weights are <= cwmax or not) */
static void
compute_weights (filter_matrix_t *mat)
{
  int nthreads;
  double cpu = seconds (), wct = wct_seconds ();

  memset (mat->wt, 0, mat->ncols * sizeof (int32_t));

#pragma omp parallel
#pragma omp master
  nthreads = omp_get_num_threads ();

  /* we first compute buckets (j,i) where thread k processes rows i = k mod nthreads,
     and bucket l corresponds to j = l mod nthreads */
  bucket1_t **B;
  B = init_buckets1 (nthreads, mat->tot_weight);
  /* thread 0 deals with rows 0, n, 2n, ...;
     thread 1 deals with rows 1, n+1, 2n+1, ...; and so on, where n = nthreads */
#pragma omp parallel for
  for (int k = 0; k < nthreads; k++)
    fill_buckets1_all (B, mat, k, nthreads);

#pragma omp parallel for schedule(static,1)
  for (int k = 0; k < nthreads; k++)
    apply_buckets1_all (B, mat, k, nthreads);

  clear_buckets1 (B, nthreads);

  printf ("   compute_weights took %.1fs (cpu), %.1fs (wct)\n",
	  seconds () - cpu, wct_seconds () - wct);
}

/* return the total weight of the matrix */
static unsigned long MAYBE_UNUSED
get_tot_weight (filter_matrix_t *mat)
{
  int nthreads;
  unsigned long tot_weight = 0;

#pragma omp parallel
#pragma omp master
  nthreads = omp_get_num_threads ();

#pragma omp parallel for schedule(static,1)
  for (int k = 0; k < nthreads; k++)
    {
      unsigned long tot_weight_thread = 0;
      for (index_t i = k; i < mat->nrows; i += nthreads)
	if (mat->rows[i] != NULL)
	  tot_weight_thread += matLengthRow (mat, i);
#pragma omp critical
      tot_weight += tot_weight_thread;
    }
  return tot_weight;
}

/* return the total weight of all columns of weight <= cwmax */
static unsigned long
get_tot_weight_columns (filter_matrix_t *mat)
{
  int nthreads;
  unsigned long tot_weight = 0;

#pragma omp parallel
#pragma omp master
  nthreads = omp_get_num_threads ();

#pragma omp parallel for schedule(static,1)
  for (int k = 0; k < nthreads; k++)
    {
      unsigned long tot_weight_thread = 0;
      for (index_t j = k; j < mat->ncols; j += nthreads)
	if (mat->wt[j] <= mat->cwmax)
	  tot_weight_thread += mat->wt[j];
#pragma omp critical
      tot_weight += tot_weight_thread;
    }
  return tot_weight;
}

static void
compute_R (filter_matrix_t *mat)
{
  int nthreads;
  double cpu = seconds (), wct = wct_seconds ();
#ifdef MEM
  unsigned long mem = 0;
#endif

  /* the inverse matrix R is already allocated, but the individual entries
     R[j] are not */
  unsigned long tot_weight = get_tot_weight_columns (mat);
  /* FIXME: we could allocate a huge chunk of memory of tot_weight * sizeof (index_t),
     to avoid small allocations, but then we cannot call free anymore on mat->R[j] */
  for (index_t j = 0; j < mat->ncols; j++)
    if (0 < mat->wt[j] && mat->wt[j] <= mat->cwmax)
      {
	mat->R[j] = malloc (mat->wt[j] * sizeof (index_t));
#ifdef MEM
	mem += mat->wt[j] * sizeof (index_t);
#endif
	mat->wt[j] = 0; /* reset to 0 */
      }
    else
      mat->R[j] = NULL; /* to please freeRj */
#ifdef MEM
  printf ("****** compute_R allocated %luM\n", mem >> 20);
#endif

#pragma omp parallel
#pragma omp master
  nthreads = omp_get_num_threads ();

  /* we first compute buckets (j,i) where thread k processes rows i = k mod nthreads,
     and bucket l corresponds to j = l mod nthreads */
  bucket_t **B;
  B = init_buckets (nthreads, tot_weight);
  /* thread 0 deals with rows 0, n, 2n, ...;
     thread 1 deals with rows 1, n+1, 2n+1, ...; and so on, where n = nthreads */
#pragma omp parallel for
  for (int k = 0; k < nthreads; k++)
    fill_buckets_all (B, mat, k, nthreads);

#pragma omp parallel for schedule(static,1)
  for (int k = 0; k < nthreads; k++)
    apply_buckets_all (B, mat, k, nthreads);

  clear_buckets (B, nthreads);

  printf ("   compute_R took %.1fs (cpu), %.1fs (wct)\n",
	  seconds () - cpu, wct_seconds () - wct);
}

typedef struct {
  index_t **list;       /* list[c] is a list of j-values with merge cost c */
  unsigned long *size;  /* size[c] is the length of list[c] */
  unsigned long *alloc; /* alloc[c] is the allocated size of list[c] */
  int cmax;             /* all costs are <= cmax */
} cost_list_t;

/* return the largest ideal of row i (excepted j) */
static index_t MAYBE_UNUSED
largest_ideal (filter_matrix_t *mat, index_t i, index_t j)
{
  uint32_t k = matLengthRow (mat, i);
  ASSERT (k > 1);
  index_t l = mat->rows[i][k];
  return (l == j) ? mat->rows[i][k-1] : l;
}

/* doit == 0: return the weight of row i1 + row i2
   doit <> 0: add row i2 to row i1 */
static int32_t
add_row (filter_matrix_t *mat, index_t i1, index_t i2, int doit)
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
  t = malloc ((c + 1) * sizeof (index_t));
  t0 = t;
  *t++ = c;
  t1 = t2 = 1;
  while (t1 <= k1 && t2 <= k2)
    {
      if (mat->rows[i1][t1] == mat->rows[i2][t2])
	t1 ++, t2 ++;
      else if (mat->rows[i1][t1] < mat->rows[i2][t2])
	*t++ = mat->rows[i1][t1++];
      else
	*t++ = mat->rows[i2][t2++];
    }
  while (t1 <= k1)
    *t++ = mat->rows[i1][t1++];
  while (t2 <= k2)
    *t++ = mat->rows[i2][t2++];
  ASSERT (t0 + (c + 1) == t);
  free (mat->rows[i1]);
  mat->rows[i1] = t0;
  return c;
}

static void
remove_row (filter_matrix_t *mat, index_t i)
{
  free (mat->rows[i]);
  mat->rows[i] = NULL;
}

static void MAYBE_UNUSED
printRow (filter_matrix_t *mat, index_t i)
{
  int32_t k = matLengthRow (mat, i);
  printf ("%lu [%d]:", (unsigned long) i, k);
  for (int j = 1; j <= k; j++)
    printf (" %lu", (unsigned long) mat->rows[i][j]);
  printf ("\n");
}

// #define TRACE_J 10637127

#ifndef USE_MST
/* classical cost: merge the row of smaller weight with the other ones:
   if out is NULL: return the merge cost
   if out <> NULL: perform the merge, output to 'out' (history file),
                   and return the merge cost */
static int32_t
merge_cost_or_do (filter_matrix_t *mat, index_t j, FILE *out)
{
  index_t imin, i;
  int32_t cmin, c;
  int32_t w = mat->wt[j];

  ASSERT (2 <= w && w <= mat->cwmax);

  imin = mat->R[j][0];
  cmin = matLengthRow (mat, imin);
#ifdef TRACE_J
  if (j == TRACE_J) printf ("TRACE_J: j=%lu i=%lu c=%d\n", j, imin, cmin);
#endif
  for (int k = 1; k < w; k++)
    {
      i = mat->R[j][k];
      c = matLengthRow (mat, i);
#ifdef TRACE_J
      if (j == TRACE_J) printf ("TRACE_J: j=%lu i=%lu c=%d\n", j, i, c);
#endif
      if (c < cmin || (c == cmin && largest_ideal (mat, i, j) <
		       largest_ideal (mat, imin, j)))
	/* in case c = cmin, we take the row with the smaller largest
	   ideal (excepted j) */
	{
	  imin = i;
	  cmin = c;
	}
    }

  if (out == NULL) /* return the merge cost, does not modify anything */
    {
      /* we remove row imin and add it to all w-1 others: cmin * (w - 2)
	 the column j disappears: -w */
      c = -cmin; /* remove row imin */
      for (int k = 0; k < w; k++)
	{
	  i = mat->R[j][k];
	  if (i != imin)
	    c += add_row (mat, i, imin, 0) - matLengthRow (mat, i);
	}
      return c;
    }
  else /* perform the real merge and output to history file */
    {
      char s0[1024], *s;
      s = s0 + sprintf (s0, "%lu", (unsigned long) imin);
      c = -cmin; /* remove row imin */
      for (int k = 0; k < w; k++)
	{
	  i = mat->R[j][k];
	  if (i != imin)
	    {
	      c -= matLengthRow (mat, i);
	      add_row (mat, i, imin, 1);
	      c += matLengthRow (mat, i);
	      s += sprintf (s, " %lu", (unsigned long) i);
	    }
	}
      /* the default sprintf function is thread-safe (see
	 https://stackoverflow.com/questions/23586682/how-to-use-printf-in-multiple-threads),
	 thus each line is written atomically in the history file
	 (the order in which lines from different threads are written is
	 non deterministic, but does not matter since the merges are
	 independent) */
      sprintf (s, "\n");
      fprintf (out, s0);
      remove_row (mat, imin);
      return c;
    }
}
#else
/* Same as reportn from report.c, but output to a string.
   Assume rep->type = 0.
   Return the number of characters written. */
static int
sreportn (char *str, index_signed_t *ind, int n, MAYBE_UNUSED index_t j)
{
  char *str0 = str;

  for (int i = 0; i < n; i++)
    {
      str += sprintf (str, "%ld", (long int) ind[i]);
      if (i < n-1)
	str += sprintf (str, " ");
    }
  str += sprintf (str, "\n");
  return str - str0;
}

/* Same as addFatherToSons() from merge_mono.c, but uses add_row instead of
   addRowsAndUpdate */
static int
addFatherToSons2 (index_t history[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1],
		  filter_matrix_t *mat, int m, index_t *ind,
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
      // addRowsAndUpdate (mat, ind[t], ind[s], ideal);
      add_row (mat, ind[t], ind[s], 1);
      history[i][2] = ind[t];
      history[i][0] = 2;
    }
  return m - 2;
}

/* compute full spanning tree */
static int32_t
merge_cost_or_do (filter_matrix_t *mat, index_t j, FILE *out)
{
  int32_t c;
  int32_t w = mat->wt[j];

  ASSERT (2 <= w && w <= mat->cwmax);

  if (out == NULL) /* return the merge cost, does not modify anything */
    return minCostUsingMST (mat, w, mat->R[j], j);
  else /* perform the real merge and output to history file */
    {
      index_t ind[MERGE_LEVEL_MAX];
      memcpy (ind, mat->R[j], w * sizeof (index_t));
      char s0[1024], *s;
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
      int hmax = addFatherToSons2 (history, mat, w, ind, start, end);
      s = s0;
      for (int i = hmax; i >= 0; i--)
	s += sreportn (s, (index_signed_t*) (history[i]+1), history[i][0], j);
      fprintf (out, s0);
      remove_row (mat, ind[0]);
      return c;
    }
}
#endif

static void
cost_list_init_aux (cost_list_t *l)
{
  l->list = NULL;
  l->size = 0;
  l->alloc = 0;
  l->cmax = -1;
}

static cost_list_t*
cost_list_init (int nthreads)
{
  cost_list_t *L;

  L = malloc (nthreads * sizeof (cost_list_t));
#pragma omp parallel for
  for (int i = 0; i < nthreads; i++)
    cost_list_init_aux (L + i);
  return L;
}

/* return the number of bytes allocated */
#ifdef MEM
static unsigned long
#else
static void
#endif
cost_list_clear_aux (cost_list_t *l)
{
#ifdef MEM
  unsigned long mem = 0;
#endif
  for (int i = 0; i <= l->cmax; i++)
    {
      free (l->list[i]);
#ifdef MEM
      mem += l->alloc[i] * sizeof (index_t);
#endif
    }
  free (l->list);
  free (l->size);
  free (l->alloc);
#ifdef MEM
  mem += l->cmax * sizeof (index_t*);
  mem += l->cmax * sizeof (unsigned long);
  mem += l->cmax * sizeof (unsigned long);
  return mem;
#endif
}

static void
cost_list_clear (cost_list_t *L, int nthreads)
{
#ifdef MEM
  unsigned long mem = 0;
#endif
#pragma omp parallel for
    for (int i = 0; i < nthreads; i++)
      {
#ifdef MEM
	unsigned long mem_thread = cost_list_clear_aux (L + i);
#pragma omp critical
	mem += mem_thread;
#else
	cost_list_clear_aux (L + i);
#endif
      }
    free (L);
#ifdef MEM
    mem += nthreads * sizeof (cost_list_t);
    printf ("****** cost_list allocated %luM\n", mem >> 20);
#endif
}

/* add pair (j,c) into l */
static void
add_cost (cost_list_t *l, index_t j, int32_t c)
{
  if (c < 0)
    c = 0; /* ensures all costs are >= 0 */
  if (c > l->cmax)
    {
      l->list = realloc (l->list, (c + 1) * sizeof (index_t*));
      l->size = realloc(l->size, (c + 1) * sizeof (unsigned long));
      l->alloc = realloc(l->alloc, (c + 1) * sizeof (unsigned long));
      for (int i = l->cmax + 1; i <= c; i++)
	{
	  l->list[i] = NULL;
	  l->size[i] = 0;
	  l->alloc[i] = 0;
	}
      l->cmax = c;
    }
  ASSERT(c <= l->cmax);
  if (l->size[c] == l->alloc[c])
    {
      unsigned long new_alloc = l->alloc[c];
      new_alloc += 1 + new_alloc / MARGIN;
      l->list[c] = realloc (l->list[c], new_alloc * sizeof (index_t));
      l->alloc[c] = new_alloc;
    }
  ASSERT(l->size[c] < l->alloc[c]);
  l->list[c][l->size[c]] = j;
  l->size[c] ++;
}

/* since merge costs might be negative (for weight 2), we translate them by 3,
   so that 2-merges with no cancellation give biased cost -2+3=1, and those
   with cancellation give biased cost 0, so they will be merged first */
#define BIAS 3

/* compute the merge costs for columns j = k mod nthreads */
static void
compute_merges_aux (cost_list_t *l, int k, int nthreads,
			 filter_matrix_t *mat)
{
#ifdef DEBUG
  int32_t cmin = INT32_MAX, cmax = INT32_MIN;
  unsigned long nmerges = 0;
#endif
  for (index_t j = k; j < mat->ncols; j += nthreads)
    if (2 <= mat->wt[j] && mat->wt[j] <= mat->cwmax)
      {
	int32_t c = merge_cost_or_do (mat, j, NULL);
	add_cost (l, j, c + BIAS);
#ifdef DEBUG
	if (c > cmax)
	  cmax = c;
	if (c < cmin)
	  cmin = c;
	nmerges ++;
#endif
      }
#ifdef DEBUG
#pragma omp critical
  printf ("l[%d] has cmin=%d cmax=%d nmerges=%lu\n", k, cmin, cmax,
	  nmerges);
#endif
}

static void
compute_merges (cost_list_t *L, int nthreads, filter_matrix_t *mat)
{
#pragma omp parallel for schedule(static,1)
  for (int i = 0; i < nthreads; i++)
    compute_merges_aux (L + i, i, nthreads, mat);

#ifdef DEBUG
  unsigned long nmerges = 0;
  for (int i = 0; i < nthreads; i++)
    for (int c = 0; c <= (L+i)->cmax; c++)
      nmerges += (L+i)->size[c];
  printf ("total number of merges: %lu\n", nmerges);
#endif
}

typedef struct {
  index_t *list;
  unsigned long size;
  unsigned long alloc;
} merge_list_t;

static void
merge_list_init (merge_list_t *l)
{
  l->list = NULL;
  l->size = 0;
  l->alloc = 0;
}

static void
merge_list_add (merge_list_t *l, index_t j)
{
  if (l->size == l->alloc)
    {
      unsigned long new_alloc = l->alloc;
      new_alloc += 1 + l->alloc / MARGIN;
      l->list = realloc (l->list, new_alloc * sizeof (index_t));
      l->alloc = new_alloc;
    }
  l->list[l->size] = j;
  l->size ++;
}

static void
merge_list_clear (merge_list_t *l)
{
  free (l->list);
}

int
largest_cmax (cost_list_t *L, int nthreads)
{
  int cmax = 0;
  for (int i = 0; i < nthreads; i++)
    if ((L + i)->cmax > cmax)
      cmax = (L + i)->cmax;
  return cmax;
}

/* apply merges of index = k mod nthreads */
static void
apply_merge_aux (index_t *l, unsigned long size, int k, int nthreads,
		 filter_matrix_t *mat, FILE *out)
{
  int32_t fill_in = 0;

  /* Note: the 2-merges are not necessarily processed in increasing fill-in,
     since those with fill-in < -2 are capped to -3. All we are sure is that
     those with fill-in < -2 (thus -4, -6, -8, ...) are processed first, then
     those with fill-in -2. */
  for (unsigned long t = k; t < size; t += nthreads)
    fill_in += merge_cost_or_do (mat, l[t], out);
#pragma omp critical
  mat->tot_weight += fill_in;
}

static void
apply_merges (cost_list_t *L, int nthreads, filter_matrix_t *mat, FILE *out)
{
  static int max_merge = 1;
  /* We first prepare a list of independent ideals of small cost to merge,
     i.e., all rows where those ideals appear are independent.
     This step is mono-thread. */
  mpz_t z; /* bit vector to detect rows not yet used */
  mpz_init (z);

  /* determine the largest cmax */
  int cmax = largest_cmax (L, nthreads);

  /* first compute the total number of possible merges */
  unsigned long total_merges = 0;
  for (int c = 0; c <= cmax; c++)
    for (int i = 0; i < nthreads; i++)
      if (c <= (L+i)->cmax)
	total_merges += (L+i)->size[c];

  merge_list_t l[1]; /* list of independent merges */
  merge_list_init (l);
  for (int c = 0; c <= cmax; c++)
    for (int i = 0; i < nthreads; i++)
      if (c <= (L+i)->cmax)
	{
	  for (unsigned long k = 0; k < (L+i)->size[c]; k++)
	    {
	      index_t j = (L+i)->list[c][k];
	      int ok = 1;
	      int w = mat->wt[j];
	      if (w > max_merge)
		{
		  printf ("First %d-merge (cost %d)\n", w, c - BIAS);
		  max_merge = w;
		}
	      for (int t = 0; t < w && ok; t++)
		ok = mpz_tstbit (z, mat->R[j][t]) == 0;
	      if (ok)
		{
		  /* mark rows used */
		  for (int t = 0; t < w; t++)
		    mpz_setbit (z, mat->R[j][t]);
		  merge_list_add (l, j);
		}
	    }
	}
  printf ("   found %lu independent merges out of %lu, with cwmax=%d\n",
	  l->size, total_merges, mat->cwmax);

  /* we increase cwmax only when the ratio of the number of independent merges
     over the number of possible merges exceeds some threshold */
  if (mat->cwmax == 2) /* we first process all 2-merges */
    {
      if (l->size == 0)
	mat->cwmax ++;
    }
  else
    {
#define THRESHOLD 0.05
      if ((double) l->size / (double) total_merges > THRESHOLD)
	{
	  if (mat->cwmax < MERGE_LEVEL_MAX)
	    mat->cwmax ++;
	}
    }

  /* We notice that the apply_merge_aux() function does not scale well when
     the number of threads is large, this is due to too many concurrent calls
     to malloc/free in add_row when performing the merges in parallel.
     We thus limit the number of threads in this part (experimentally
     MAX_THREADS=16 seems optimal to minimize the wall-clock time). */
#define MAX_THREADS 16

  /* reduce the number of threads to avoid blocking issues with malloc/free */
  if (nthreads > MAX_THREADS)
    nthreads = MAX_THREADS;

  /* now apply in parallel the independent merges */
#pragma omp parallel for schedule(static,1)
  for (int k = 0; k < nthreads; k++)
    apply_merge_aux (l->list, l->size, k, nthreads, mat, out);

  /* each merge decreases the number of rows by one */
  mat->rem_nrows -= l->size;

  merge_list_clear (l);
  mpz_clear (z);
}

static double
average_density (filter_matrix_t *mat)
{
  return (double) mat->tot_weight / (double) mat->rem_nrows;
}

int
main (int argc, char *argv[])
{
    char *argv0 = argv[0];

    filter_matrix_t mat[1];
    report_t rep[1];

    int nthreads = 1;
    uint32_t keep = DEFAULT_FILTER_EXCESS;
    uint32_t skip = DEFAULT_MERGE_SKIP;
    double target_density = DEFAULT_MERGE_TARGET_DENSITY;

#ifdef HAVE_MINGW
    _fmode = _O_BINARY;     /* Binary open for all files */
#endif

    double tt;
    double cpu0 = seconds ();
    double wct0 = wct_seconds ();
    param_list pl;
    int verbose = 0;
    param_list_init (pl);
    declare_usage(pl);
    argv++,argc--;

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

    const char *purgedname = param_list_lookup_string (pl, "mat");
    const char *outname = param_list_lookup_string (pl, "out");

    param_list_parse_int (pl, "t", &nthreads);
#ifdef HAVE_OPENMP
    omp_set_num_threads (nthreads);
#endif

    param_list_parse_uint (pl, "keep", &keep);
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

    /* Read number of rows and cols on first line of purged file */
    purgedfile_read_firstline (purgedname, &(mat->nrows), &(mat->ncols));

    /* initialize rep (i.e., mostly opens outname) and write matrix dimension */
    rep->type = 0;
    rep->outfile = fopen_maybe_compressed (outname, "w");
    ASSERT_ALWAYS(rep->outfile != NULL);

    /* Init structure containing the matrix and the heap of potential merges */
    initMat (mat, MERGE_LEVEL_MAX, keep, skip, 0);

    /* we bury the 'skip' ideals of smallest index */
    mat->force_bury_below_index = skip;

    /* Read all rels and fill-in the mat structure */
    tt = seconds ();
    filter_matrix_read2 (mat, purgedname);
    printf ("Time for filter_matrix_read: %2.2lfs\n", seconds () - tt);

    printf ("N=%lu W=%lu W/N=%.2f cpu=%.1fs wct=%.1fs mem=%luM\n",
	    mat->rem_nrows, mat->tot_weight, average_density (mat),
	    seconds () - cpu0, wct_seconds () - wct0,
	    PeakMemusage () >> 10);

    mat->cwmax = 2;

    while (average_density (mat) < target_density)
      {
	double cpu1 = seconds (), wct1 = wct_seconds ();

	compute_weights (mat);

	compute_R (mat);

	double cpu2 = seconds (), wct2 = wct_seconds ();

	cost_list_t *L = cost_list_init (nthreads);
	compute_merges (L, nthreads, mat);

	printf ("   compute_merges took %.1fs (cpu), %.1fs (wct)\n",
		seconds () - cpu2, wct_seconds () - wct2);

	double cpu3 = seconds (), wct3 = wct_seconds ();

	apply_merges (L, nthreads, mat, rep->outfile);

	printf ("   apply_merges took %.1fs (cpu), %.1fs (wct)\n",
		seconds () - cpu3, wct_seconds () - wct3);

	cost_list_clear (L, nthreads);

	printf ("   pass took %.1fs (cpu), %.1fs (wct)\n",
		seconds () - cpu1, wct_seconds () - wct1);

	printf ("N=%lu W=%lu W/N=%.2f cpu=%.1fs wct=%.1fs mem=%luM\n",
		mat->rem_nrows, mat->tot_weight,
		(double) mat->tot_weight / (double) mat->rem_nrows,
		seconds () - cpu0, wct_seconds () - wct0,
		PeakMemusage () >> 10);

      }

    fclose_maybe_compressed (rep->outfile, outname);

    clearMat (mat);

    param_list_clear (pl);

    print_timing_and_memory (stdout, wct0);

    return 0;
}
