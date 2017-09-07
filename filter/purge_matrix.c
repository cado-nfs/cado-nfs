#include "cado.h"
#include <stdio.h>
#include <stdlib.h>

#include "portability.h"

#include "utils_with_io.h"
#include "filter_config.h"
#include "purge_matrix.h"
#include "clique_removal.h" /* for computing stats on cliques */

/* If HAVE_SYNC_FETCH is not defined, we will use mutex for multithreaded
 * version of the code. May be too slow. */
#ifndef HAVE_SYNC_FETCH
static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
#endif

void
purge_matrix_init (purge_matrix_ptr mat, uint64_t nrows_init,
                   uint64_t col_min_index, uint64_t col_max_index)
{
  size_t cur_alloc;
  mat->nrows_init = nrows_init;
  mat->col_max_index = col_max_index;
  mat->nrows = mat->ncols = 0;
  mat->col_min_index = col_min_index;
  mat->tot_alloc_bytes = 0;

  /* Malloc cols_weight ans set to 0 */
  cur_alloc = mat->col_max_index * sizeof (weight_t);
  mat->cols_weight = (weight_t *) malloc(cur_alloc);
  ASSERT_ALWAYS(mat->cols_weight != NULL);
  mat->tot_alloc_bytes += cur_alloc;
  fprintf(stdout, "# MEMORY: Allocated cols_weight of %zuMB (total %zuMB "
                  "so far)\n", cur_alloc >> 20, mat->tot_alloc_bytes >> 20);
  memset(mat->cols_weight, 0, cur_alloc);

  /* Malloc sum2_row */
  cur_alloc = mat->col_max_index * sizeof (uint64_t);
  mat->sum2_row = (uint64_t *) malloc(cur_alloc);
  ASSERT_ALWAYS(mat->sum2_row != NULL);
  mat->tot_alloc_bytes += cur_alloc;
  fprintf(stdout, "# MEMORY: Allocated sum2_row of %zuMB (total %zuMB "
                  "so far)\n", cur_alloc >> 20, mat->tot_alloc_bytes >> 20);

  /* Malloc row_compact */
  cur_alloc = mat->nrows_init * sizeof (index_t *);
  mat->row_compact = (index_t **) malloc (cur_alloc);
  ASSERT_ALWAYS(mat->row_compact != NULL);
  mat->tot_alloc_bytes += cur_alloc;
  fprintf(stdout, "# MEMORY: Allocated row_compact of %zuMB (total %zuMB "
                  "so far)\n", cur_alloc >> 20, mat->tot_alloc_bytes >> 20);
}

/* Free everything and set everythin to 0. */
void
purge_matrix_clear (purge_matrix_ptr mat)
{
  free (mat->cols_weight);
  free (mat->sum2_row);
  my_malloc_free_all();
  free(mat->row_compact);

  memset (mat, 0, sizeof (purge_matrix_t));
}

void
purge_matrix_row_compact_update_mem_usage (purge_matrix_ptr mat)
{
  size_t cur_alloc = get_my_malloc_bytes();
  mat->tot_alloc_bytes += cur_alloc;
  fprintf(stdout, "# MEMORY: Allocated row_compact[i] %zuMB (total %zuMB so "
                  "far)\n", cur_alloc >> 20, mat->tot_alloc_bytes >> 20);
}


/* Set a row of a purge_matrix_t from a read relation.
 * The row that is set is decided by the relation number rel->num
 * We put in mat->row_compact only primes such that their index h is greater or
 * equal to mat->col_min_index.
 * A row in mat->row_compact is ended by a -1 (= UMAX(index_t))
 * The number of columns mat->ncols is updated is necessary.
 * The number of rows mat->ncols is increased by 1.
 *
 * The return type of this function is void * instead of void so we can use it
 * as a callback function for filter_rels.
 *
 * Not thread-safe.
 */

void *
purge_matrix_set_row_from_rel (purge_matrix_t mat, earlyparsed_relation_ptr rel)
{
  ASSERT_ALWAYS(rel->num < mat->nrows_init);

  unsigned int nb_above_min_index = 0;
  for (weight_t i = 0; i < rel->nb; i++)
    nb_above_min_index += (rel->primes[i].h >= mat->col_min_index);

  index_t *tmp_row = index_my_malloc (1 + nb_above_min_index);
  unsigned int next = 0;
  for (weight_t i = 0; i < rel->nb; i++)
  {
    index_t h = rel->primes[i].h;
    if (mat->cols_weight[h] == 0)
    {
      mat->cols_weight[h] = 1;
      (mat->ncols)++;
    }
    else if (mat->cols_weight[h] != UMAX(weight_t))
      mat->cols_weight[h]++;

    if (h >= mat->col_min_index)
      tmp_row[next++] = h;
  }

  tmp_row[next] = UMAX(index_t); /* sentinel */
  mat->row_compact[rel->num] = tmp_row;
  (mat->nrows)++;

  return NULL;
}


/* Delete a row: set mat->row_used[i] to 0, update the count of the columns
 * appearing in that row, mat->nrows and mat->ncols.
 * Warning: we only update the count of columns that we consider, i.e.,
 * columns with index >= mat->col_min_index.
 * CAREFUL: no multithread compatible with " !(--(*o))) " and
 * bit_vector_clearbit.
 */
void
purge_matrix_delete_row (purge_matrix_ptr mat, uint64_t i)
{
  index_t *row;
  weight_t *w;

  if (purge_matrix_is_row_active(mat, i))
  {
    for (row = mat->row_compact[i]; *row != UMAX(*row); row++)
    {
      w = &(mat->cols_weight[*row]);
      ASSERT(*w);
      /* Decrease only if not equal to the maximum value */
      /* If weight becomes 0, we just remove a column */
      if (*w < UMAX(*w) && !(--(*w)))
        mat->ncols--;
    }
    /* We do not free mat->row_compact[i] as it is freed later with
      my_malloc_free_all */
    purge_matrix_set_row_inactive (mat, i); /* mark as deleted */
    mat->nrows--;
  }
  else /* Row already deleted => Abort (It means there is a bug somewhere) */
  {
    fprintf (stderr, "Error, row %" PRIu64" is already deleted\n", i);
    abort();
  }
}

/***** Functions to compute mat->sum2_row array (mono and multi thread) *******/

/* sum2_row[h] = 0 if cols_weight[h] <> 2
               = i1 + i2 if cols_weight[h] = 2 and h appears in rows i1 and i2
 */

/* Code for computing mat->sum2_row --- monothread version */

static inline void
purge_matrix_compute_sum2_row_mono (purge_matrix_ptr mat)
{
  /* Reset mat->sum2_row to 0 */
  memset(mat->sum2_row, 0, mat->col_max_index * sizeof(uint64_t));
  for (uint64_t i = 0; i < mat->nrows_init; i++)
  {
    index_t h, *row_ptr;
    if (purge_matrix_is_row_active(mat, i))
      for (row_ptr = mat->row_compact[i]; (h = *row_ptr++) != UMAX(h);)
        if (mat->cols_weight[h] == 2)
          mat->sum2_row[h] += i;
  }
}

/* Code for computing mat->sum2_row --- multithread version */

typedef struct sum2_mt_data_s {
  /* Read only part */
  uint64_t begin, end;
  /* Read-write part */
  purge_matrix_ptr mat;
} sum2_mt_data_t;

void *
purge_matrix_compute_sum2_row_mt_thread (void *pt)
{
  sum2_mt_data_t *data = (sum2_mt_data_t *) pt;
  purge_matrix_ptr mat = data->mat;
  uint64_t i;

  for (i = data->begin; i < data->end; i++)
  {
    if (purge_matrix_is_row_active(mat, i))
    {
      index_t h, *row_ptr;
      for (row_ptr = mat->row_compact[i]; (h = *row_ptr++) != UMAX(h);)
      {
        if (LIKELY (mat->cols_weight[h] == 2))
        {
#ifdef HAVE_SYNC_FETCH
          __sync_add_and_fetch (mat->sum2_row + h, i);
#else /* else we use mutex to protect the addition on mat->sum2_row[h] */
          pthread_mutex_lock (&lock);
          mat->sum2_row[h] += i;
          pthread_mutex_unlock (&lock);
#endif
        }
      }
    }
  }
  return NULL;
}

static inline void
purge_matrix_compute_sum2_row_mt (purge_matrix_ptr mat, unsigned int nthreads)
{
  pthread_t *threads = NULL;
  sum2_mt_data_t *th_data = NULL;
  uint64_t nrows_per_thread, k;
  int err;

  /* Reset mat->sum2_row to 0 */
  memset(mat->sum2_row, 0, mat->col_max_index * sizeof(uint64_t));

  th_data = (sum2_mt_data_t *) malloc (nthreads * sizeof(sum2_mt_data_t));
  ASSERT_ALWAYS(th_data != NULL);
  threads = (pthread_t *) malloc (nthreads * sizeof(pthread_t));
  ASSERT_ALWAYS(threads != NULL);

  /* nrows_per_thread MUST be a multiple of BV_BITS (see how we loop on the rows
   * in purge_matrix_compute_sum2_row_mt_thread
   */
  nrows_per_thread = (mat->nrows_init / nthreads) & ((uint64_t) ~(BV_BITS - 1));
  k = 0;
  for (unsigned int i = 0; i < nthreads; i++)
  {
    th_data[i].mat = mat;
    th_data[i].begin = k;
    k += nrows_per_thread;
    th_data[i].end = k;
  }
  th_data[nthreads-1].end = mat->nrows_init;

  for (unsigned int i = 0; i < nthreads; i++)
  {
    if ((err = pthread_create(&threads[i], NULL,
                              &purge_matrix_compute_sum2_row_mt_thread,
                              (void *) &(th_data[i]))))
    {
      fprintf(stderr, "Error, pthread_create failed in purge_matrix_compute_"
                      "sum2_row_mt: %d. %s\n", err, strerror(errno));
      abort();
    }
  }
  for (unsigned int i = 0; i < nthreads; i++)
    pthread_join (threads[i], NULL);

  free(threads);
  free(th_data);
}

void
purge_matrix_compute_sum2_row (purge_matrix_ptr mat, unsigned int nthreads)
{
  if (nthreads > 1)
    purge_matrix_compute_sum2_row_mt (mat, nthreads);
  else
    purge_matrix_compute_sum2_row_mono (mat);
}

/******************* Functions to print stats ********************************/

/* These 2 functions compute and print stats on rows weight and columns weight.
 * The stats can be expensive to compute, so these functions should not be
 * called by default. */

/* Internal function */
void
print_stats_uint64 (FILE *out, uint64_t *w, uint64_t len, char name[],
                    char unit[], int verbose)
{
  uint64_t av = 0, min = UMAX(uint64_t), max = 0, std = 0, nb_nzero = 0;
  for (uint64_t i = 0; i < len; i++)
  {
    if (w[i] > 0)
    {
      nb_nzero++;
      if (w[i] < min)
        min = w[i];
      if (w[i] > max)
        max = w[i];
      av += w[i];
      std += w[i]*w[i];
    }
  }

  double av_f = ((double) av) / ((double) nb_nzero);
  double std_f = sqrt(((double) std) / ((double) nb_nzero) - av_f*av_f);

  fprintf (out, "# STATS on %s: #%s = %" PRIu64 "\n", name, name, nb_nzero);
  fprintf (out, "# STATS on %s: min %s = %" PRIu64 "\n", name, unit, min);
  fprintf (out, "# STATS on %s: max %s = %" PRIu64 "\n", name, unit, max);
  fprintf (out, "# STATS on %s: av %s = %.2f\n", name, unit, av_f);
  fprintf (out, "# STATS on %s: std %s = %.2f\n", name, unit, std_f);

  if (verbose > 1)
  {
    uint64_t *nb_w = NULL;
    nb_w = (uint64_t *) malloc ((max-min+1) * sizeof (uint64_t));
    ASSERT_ALWAYS (nb_w != NULL);
    memset (nb_w, 0, (max-min+1) * sizeof (uint64_t));
    for (uint64_t i = 0; i < len; i++)
      if (w[i] > 0)
        nb_w[w[i]-min]++;

    for (uint64_t i = 0; i < max-min+1; i++)
    {
      if (nb_w[i] > 0)
        fprintf (out, "# STATS on %s: #%s of %s %" PRIu64 " : %" PRIu64
                      "\n", name, name, unit, min+i, nb_w[i]);
    }
    free (nb_w);
  }
  fflush (out);
}

void
purge_matrix_print_stats_columns_weight (FILE *out, purge_matrix_srcptr mat,
                                         int verbose)
{
  uint64_t *w = NULL;
  index_t h, *row_ptr;
  w = (uint64_t *) malloc (mat->col_max_index * sizeof (uint64_t));
  ASSERT_ALWAYS (w != NULL);
  memset (w, 0, mat->col_max_index * sizeof (uint64_t));

  for (uint64_t i = 0; i < mat->nrows_init; i++)
    if (purge_matrix_is_row_active(mat, i))
      for (row_ptr = mat->row_compact[i]; (h = *row_ptr++) != UMAX(h);)
        w[h]++;

  print_stats_uint64 (out, w, mat->col_max_index, "cols", "weight", verbose);
  free (w);
}

void
purge_matrix_print_stats_rows_weight (FILE *out, purge_matrix_srcptr mat,
                                      int verbose)
{
  uint64_t *w = NULL;
  index_t h, *row_ptr;
  w = (uint64_t *) malloc (mat->nrows_init * sizeof (uint64_t));
  ASSERT_ALWAYS (w != NULL);
  memset (w, 0, mat->nrows_init * sizeof (uint64_t));

  for (uint64_t i = 0; i < mat->nrows_init; i++)
    if (purge_matrix_is_row_active(mat, i))
      for (row_ptr = mat->row_compact[i]; (h = *row_ptr++) != UMAX(h);)
        w[i]++;

  print_stats_uint64 (out, w, mat->nrows_init, "rows", "weight", verbose);
  free (w);
}
