#include "cado.h"

#include <string.h>
#include <stdio.h>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "portability.h"
#include "utils_with_io.h"
#include "filter_config.h"
#include "merge_replay_matrix.h"
#include "sparse.h"
#include "markowitz.h"
#include "mst.h"

/***************** memory allocation on R[j] *********************************/

static void
mallocRj (filter_matrix_t *mat, index_t j, int32_t w)
{
  if (w == 0)
    {
      /* This can happen because we do not renumber columns in purge,
	 thus the range of the j-values is larger than the number of
	 columns. We simply ignore those columns. */
      mat->R[j] = NULL;
      return;
    }
  mat->R[j] = (index_t *) malloc((w + 1) * sizeof(index_t));
  FATAL_ERROR_CHECK(mat->R[j] == NULL, "Cannot allocate memory");
  mat->R[j][0] = 0; /* last index used */
}

static void
reallocRj (filter_matrix_t *mat, index_t j, int32_t w)
{
  mat->R[j] = (index_t *) realloc (mat->R[j], (w + 1) * sizeof(index_t));
  FATAL_ERROR_CHECK(mat->R[j] == NULL, "Cannot reallocate memory");
  mat->R[j][0] = w;
}

void
freeRj (filter_matrix_t *mat, index_t j)
{
    free (mat->R[j]);
    mat->R[j] = NULL;
}

/****************************** heap management ******************************/

/* we precompute the weights for w < 256 to avoid the cost of multiple calls
   to comp_weight_function */
static float comp_weight[256];

static inline float
comp_weight_function (int32_t w MAYBE_UNUSED, int maxlevel MAYBE_UNUSED)
{
#define USE_WEIGHT_LAMBDA 0 /* LAMBDA = 0 seems to be better */
#if USE_WEIGHT_LAMBDA == 0
  /* we only count ideals of weight > maxlevel, assuming all those of weight
     <= maxlevel will be merged */
    return (w <= maxlevel) ? 0.0 : 1.0;
#elif USE_WEIGHT_LAMBDA == 1
    return powf (2.0 / 3.0, (float) (w - 2));
#elif USE_WEIGHT_LAMBDA == 2
    return powf (0.5, (float) (w - 2));
#elif USE_WEIGHT_LAMBDA == 3
    return powf (0.8, (float) (w - 2));
#elif USE_WEIGHT_LAMBDA == 4
    return 1.0 / log2f ((float) w);
#elif USE_WEIGHT_LAMBDA == 5
    return 2.0 / (float) w;
#elif USE_WEIGHT_LAMBDA == 6
    return 4.0 / (float) (w * w);
#else
#error "Invalid value of USE_WEIGHT_LAMBDA"
#endif
}

static void
heap_init (heap H, index_t nrows, int maxlevel)
{
  H->list = malloc (nrows * sizeof (index_t));
  FATAL_ERROR_CHECK(H->list == NULL, "Cannot allocate memory");
  H->size = 0;
  H->alloc = nrows;
  H->index = malloc (nrows * sizeof (index_t));
  FATAL_ERROR_CHECK(H->index == NULL, "Cannot allocate memory");
  for (unsigned long i = 0; i < nrows; i++)
    H->index[i] = UINT32_MAX;
  for (int w = 0; w < 256; w++)
    comp_weight[w] = comp_weight_function (w, maxlevel);
  printf ("Using weight function lambda=%d for clique removal\n",
          USE_WEIGHT_LAMBDA);
}

/* fill the heap of heavy rows */
void
heap_fill (filter_matrix_t *mat)
{
  for (unsigned long i = 0; i < mat->nrows; i++)
    {
      ASSERT_ALWAYS(mat->rows[i] != NULL);
      heap_push (mat->Heavy, mat, i);
    }
}

static void
heap_clear (heap H)
{
  free (H->list);
  free (H->index);
}

static float
heapRowWeight (index_t i, filter_matrix_t *mat)
{
  float W = 0.0;
  int32_t w;
  index_t j;

  ASSERT(mat->rows[i] != NULL);
  for (uint32_t k = 1; k <= matLengthRow (mat, i); k++)
    {
      j = matCell (mat, i, k);
      w = mat->wt[j];
      ASSERT(w != 0);
      if (w < 0)
        w = -w;
      if (w > 255)
        w = 255; /* saturate to 255 */
      W += comp_weight[w];
    }
  return W;
}

#if 0
static void
check_heap (heap H, filter_matrix_t *mat)
{
  for (unsigned long i = 0; i < H->size; i++)
    {
      index_t j = H->list[i];
      ASSERT_ALWAYS(H->index[j] == i);
      ASSERT_ALWAYS(mat->rows[j] != NULL);
    }
}
#endif

/* move down entry n of heap */
static void
moveDown (heap H, filter_matrix_t *mat, index_t n)
{
  index_t i = H->list[n], j;
  float w;

  w = heapRowWeight (i, mat);
  while (2 * n + 1 < H->size)
    {
      index_t left = 2 * n + 1, right = 2 * n + 2, son;
      if (right >= H->size ||
          heapRowWeight (H->list[left], mat) > heapRowWeight (H->list[right], mat))
        son = left; /* compare with left son */
      else
        son = right; /* compare with right son */
      if (heapRowWeight (H->list[son], mat) > w)
        {
          j = H->list[son];
          H->list[n] = j;
          H->index[j] = n;
          n = son;
        }
      else
        break; /* i has larger weight */
    }
  H->list[n] = i;
  H->index[i] = n;
}

/* remove relation i from heap */
void
heap_delete (heap H, filter_matrix_t *mat, index_t i)
{
  index_t n = H->index[i];

  ASSERT(H->size > 0);

  /* here, the row of cell n of the heap is invalid */

  /* put last entry in position n (if not already last) */
  H->size --;
  if (n < H->size)
    {
      H->list[n] = H->list[H->size];
      /* adjust heap */
      moveDown (H, mat, n);
    }
}

/* 1,2 -> 0, 3,4 -> 1, 5,6 -> 2, ... */
#define PARENT(i) (((i)+1)/2-1)

/* add relation i */
void
heap_push (heap H, filter_matrix_t *mat, index_t i)
{
  index_t n, j;
  float w;

  if (H->index[i] == UINT32_MAX)
    {
      ASSERT(H->size < H->alloc);

      w = heapRowWeight (i, mat);
      n = H->size;

      /* move parents down the tree as long as they have smaller weight */
      while (n > 0 && heapRowWeight (H->list[PARENT(n)], mat) < w)
        {
          j = H->list[PARENT(n)];
          H->list[n] = j;
          H->index[j] = n;
          n = PARENT(n);
        }
      /* insert new element */
      H->list[n] = i;
      H->index[i] = n;
      H->size ++;
    }
  else /* relation was already in heap */
    {
      n = H->index[i];
      w = heapRowWeight (i, mat);
      /* if new weight is smaller than old one, move down the heap */
      if ((2 * n + 1 < H->size && w < heapRowWeight (H->list[2 * n + 1], mat)) ||
          (2 * n + 2 < H->size && w < heapRowWeight (H->list[2 * n + 2], mat)))
        moveDown (H, mat, n);
      else /* move up the heap */
        {
          while (n > 0 && w > heapRowWeight (H->list[PARENT(n)], mat))
            {
              j = H->list[PARENT(n)];
              H->list[n] = j;
              H->index[j] = n;
              n = PARENT(n);
            }
          /* insert updated element */
          H->list[n] = i;
          H->index[i] = n;
        }
    }
}

/* return index i of relation with larger weight in heap H */
index_t
heap_pop (heap H, filter_matrix_t *mat MAYBE_UNUSED)
{
  ASSERT(H->size > 0);
  return H->list[0];
}

/*****************************************************************************/

int
decrS (int w)
{
  /* since MkzDecreaseColWeight is only called for mat->wt[j] >= 0
     in removeCellAndUpdate(), we always have w >= 0 here */
  ASSERT(w >= 0);
  return w - 1;
}

/* increment the absolute value of a signed number */
int
incrS (int w)
{
  return (w >= 0) ? w + 1 : w - 1;
}

/* initialize the sparse matrix mat */
void
initMat (filter_matrix_t *mat, int maxlevel, uint32_t keep,
         uint32_t nburied)
{
  mat->keep  = keep;
  mat->mergelevelmax = maxlevel;
  /* we start with cwmax = 2, and we increase it in mergeOneByOne() when
     the Markowitz queue is empty */
  mat->cwmax = 2;
  ASSERT_ALWAYS (mat->cwmax < 255); /* 255 is reserved for saturated values */
  mat->nburied = nburied;

  mat->weight = 0;
  mat->tot_weight = 0;
  mat->rem_ncols = 0;

  mat->rows = (typerow_t **) malloc (mat->nrows * sizeof (typerow_t *));
  ASSERT_ALWAYS (mat->rows != NULL);
  mat->wt = (int32_t *) malloc (mat->ncols * sizeof (int32_t));
  ASSERT_ALWAYS (mat->wt != NULL);
  memset (mat->wt, 0, mat->ncols * sizeof (int32_t));
  mat->R = (index_t **) malloc (mat->ncols * sizeof(index_t *));
  ASSERT_ALWAYS(mat->R != NULL);

  /* initialize the heap for heavy rows */
  heap_init (mat->Heavy, mat->nrows, mat->mergelevelmax);
}

void
clearMat (filter_matrix_t *mat)
{
  uint64_t i, j;

  for (i = 0; i < mat->nrows; i++)
    free (mat->rows[i]);
  free (mat->rows);
  free (mat->wt);
  for (j = 0; j < mat->ncols; j++)
    freeRj (mat, j);
  free (mat->R);

  /* free the memory for the heap of heavy rows */
  heap_clear (mat->Heavy);
}

#ifndef FOR_DL
/* Renumber the non-zero columns to contiguous values [0, 1, 2, ...]
 * We can use this function only for factorization because in this case we do
 * not need the indexes of the columns (contrary to DL where the indexes of the
 * column are printed in the history file). */
static void
renumber_columns (filter_matrix_t *mat)
{
  index_t *p = NULL;

  p = (index_t *) malloc (mat->ncols * sizeof (index_t));
  ASSERT_ALWAYS (p != NULL);

  /* first compute the mapping of column indices */
  index_t h = 0;
  for (uint64_t j = 0; j < mat->ncols; j++)
  {
    /* at this point wt[j] = 0 if ideal j never appears in relations,
     * or wt[j] > 0 if ideal j appears in relations */
    ASSERT(mat->wt[j] >= 0);
    if (mat->wt[j] > 0)
    {
      /* index j is mapped to h, with h <= j */
      p[j] = h;
#if TRACE_COL >= 0
      if (TRACE_COL == h || TRACE_COL == j)
        printf ("TRACE_COL: column %lu is renumbered into %u\n", j, h);
#endif
      mat->wt[h] = mat->wt[j];
      h++;
    }
  }
  /* h should be equal to rem_ncols, which is the number of columns with
   * non-zero weight */
  ASSERT_ALWAYS(h + mat->nburied == mat->rem_ncols);

  /* Realloc mat->wt */
  mat->wt = realloc (mat->wt, h * sizeof (int32_t));
  ASSERT_ALWAYS(mat->wt != 0);
  /* Reset mat->ncols to behave as if they were no non-zero columns. */
  mat->ncols = h;

  /* apply mapping to the rows. As p is a non decreasing function, the rows are
   * still sorted after this operation. */
  for (uint64_t i = 0; i < mat->nrows; i++)
    for (index_t j = 1; j <= matLengthRow(mat, i); j++)
      {
        index_t h = matCell (mat, i, j);
        setCell (mat->rows[i], j, p[h], 1); /* for factorization, the exponent
                                               does not matter */
      }

  free (p);
}
#endif

/* initialize Rj[j] for light columns, i.e., for those of weight <= cwmax */
static void
initMatR (filter_matrix_t *mat)
{
  index_t h;
  int32_t wmax = mat->cwmax;

#ifndef FOR_DL
  renumber_columns (mat);
#endif

#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
  for (h = 0; h < mat->ncols; h++)
    {
      int32_t w;
      w = mat->wt[h];
      if (w <= wmax)
        mallocRj (mat, h, w);
      else /* weight is larger than cwmax */
        {
          mat->wt[h] = -mat->wt[h]; // trick!!!
          mat->R[h] = NULL;
        }
    }
}

static void
recompute_weights (filter_matrix_t *mat)
{
  uint64_t i, h;

  /* reset the weights to zero */
  memset (mat->wt, 0, mat->ncols * sizeof (int32_t));

  /* recompute the weights */
  for (i = 0; i < mat->nrows; i++)
    if (mat->rows[i] != NULL) /* row is still active */
      for (unsigned int k = 1 ; k <= matLengthRow(mat, i); k++)
        {
          h = matCell(mat, i, k);
          mat->wt[h] ++;
        }
}

static void
reinitMatR (filter_matrix_t *mat)
{
  index_t h;
  int32_t wmax = mat->cwmax;

#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
  for (h = 0; h < mat->ncols; h++)
    {
      int32_t w;
      w = mat->wt[h];
      if (w <= wmax)
        {
          reallocRj (mat, h, w);
          mat->R[h][0] = 0; /* we reset the weights to 0 */
        }
      else /* weight is larger than cwmax */
        {
          mat->wt[h] = -mat->wt[h]; // trick!!!
          /* If w > wmax, since wmax is only increasing, we have already
             destroyed the column before, thus mat->R[h] should be NULL. */
          ASSERT(mat->R[h] == NULL);
        }
    }
}

static void
fillR (filter_matrix_t *mat)
{
  uint64_t i, h;

  for (i = 0; i < mat->nrows; i++)
    if (mat->rows[i] != NULL) /* row i is still active */
      for (unsigned int k = 1 ; k <= matLengthRow(mat, i); k++)
        {
          h = matCell(mat, i, k);
          if (mat->wt[h] > 0)
            {
              ASSERT(mat->R[h] != NULL);
              mat->R[h][0]++;
              mat->R[h][mat->R[h][0]] = i;
            }
        }

#if TRACE_COL >= 0
  h = TRACE_COL;
  printf ("TRACE_COL: weight of ideal %lu is %d\n", h, mat->wt[h]);
  ASSERT_ALWAYS(mat->wt[h] <= 0 || (uint32_t) mat->wt[h] == mat->R[h][0]);
#endif
}

/* Put in nbm[w] for 0 <= w < 256, the number of ideals of weight w.
   Return the number of active columns (w > 0). */
unsigned long
weight_count (filter_matrix_t *mat, uint64_t *nbm)
{
  uint64_t h, active = 0;

  for (h = 0; h < 256; h++)
    nbm[h] = 0;
  for (h = 0; h < mat->ncols; h++)
    {
      if (0 <= mat->wt[h] && mat->wt[h] < 256)
        nbm[mat->wt[h]]++;
      active += mat->wt[h] > 0;
    }
  return active;
}

void
recomputeR (filter_matrix_t *mat)
{
#ifdef TIMINGS
  trecomputeR -= seconds ();
#endif

  /* recompute the column weights */
  recompute_weights (mat);

  /* re-allocate the matrix R */
  reinitMatR (mat);

  /* fill the matrix R */
  fillR (mat);

  /* recompute the Markowitz structure */
  MkzClear (mat, 0);
  MkzInit (mat, 0);

#ifdef TIMINGS
  trecomputeR += seconds ();
#endif
}

void
matR_disable_cols (filter_matrix_t *mat, const char *infilename)
{
  FILE *file = NULL;
  char buf[256];

  file = fopen_maybe_compressed (infilename, "r");
  ASSERT_ALWAYS (file != NULL);

  int stop = 0;
  while (!stop)
  {
    if (fgets(buf, 256, file) == NULL)
      stop = 1;
    else if (buf[0] != '#')
    {
      size_t n = strnlen(buf, 256);
      ASSERT_ALWAYS(n != 256);

      index_t h;
      int ret = sscanf (buf, "%" SCNid "\n", &h);
      ASSERT_ALWAYS (ret == 1);
      if (h < mat->ncols)
      {
        if (mat->R[h] != NULL)
          freeRj (mat, h);
        if (mat->wt[h] > 0)
          mat->wt[h] = -mat->wt[h]; // trick!!!
      }
    }
  }

  fclose_maybe_compressed (file, infilename);
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
void * insert_rel_into_table (void *context_data, earlyparsed_relation_ptr rel)
{
  filter_matrix_t *mat = (filter_matrix_t *) context_data;
  unsigned int j = 0;
  typerow_t buf[UMAX(weight_t)]; /* rel->nb is of type weight_t */

  for (unsigned int i = 0; i < rel->nb; i++)
  {
    index_t h = rel->primes[i].h;
#ifdef FOR_DL
    exponent_t e = rel->primes[i].e;
    /* For factorization, they should not be any multiplicity here.
       For DL we do not want to count multiplicity in mat->wt */
    buf[++j] = (ideal_merge_t) {.id = h, .e = e};
#else
    ASSERT(rel->primes[i].e == 1);
    buf[++j] = h;
#endif
    mat->rem_ncols += (mat->wt[h] == 0);
    mat->wt[h] += (mat->wt[h] != SMAX(int32_t));
  }

#ifdef FOR_DL
  buf[0].id = rel->nb;
#else
  buf[0] = rel->nb;
#endif
  mat->tot_weight += rel->nb;

  /* sort indices to ease row merges */
#ifndef FOR_DL
  sort_relation (&(buf[1]), rel->nb);
#else
  qsort (&(buf[1]), rel->nb, sizeof(typerow_t), cmp_typerow_t);
#endif

  mat->rows[rel->num] = mallocRow (rel->nb + 1);
  compressRow (mat->rows[rel->num], buf, rel->nb);

  return NULL;
}


void
filter_matrix_read (filter_matrix_t *mat, const char *purgedname)
{
  uint64_t nread, h;
  char *fic[2] = {(char *) purgedname, NULL};

  /* read all rels */
  nread = filter_rels(fic, (filter_rels_callback_t) &insert_rel_into_table, mat,
                      EARLYPARSE_NEED_INDEX, NULL, NULL);
  ASSERT_ALWAYS(nread == mat->nrows);
  mat->rem_nrows = nread;

  /* print weight count */
  uint64_t nbm[256], total;
  total = weight_count (mat, nbm);
  for (h = 1; h <= (uint64_t) mat->mergelevelmax; h++)
    printf ("There are %" PRIu64 " column(s) of weight %" PRIu64 "\n", nbm[h], h);
  printf ("Total %" PRIu64 " columns\n", total);
  ASSERT_ALWAYS(mat->rem_ncols == mat->ncols - nbm[0]);

  int weight_buried_is_exact = 1;
  uint64_t weight_buried = 0;
  uint64_t i;
  /* Bury heavy coloumns. The 'nburied' heaviest column are buried. */
  /* Buried columns are not taken into account by merge. */
  if (mat->nburied)
  {
    uint64_t *heaviest = NULL;
    heaviest = (uint64_t *) malloc (mat->nburied * sizeof (uint64_t));
    ASSERT_ALWAYS(heaviest != NULL);
    for (i = 0; i < mat->nburied; i++)
      heaviest[i] = UMAX(uint64_t);

    /* Find the 'nburied' heaviest columns */
    for (i = 0; i < mat->ncols; i++)
    {
      int32_t w_cur = mat->wt[i];
      uint64_t j = mat->nburied;

      while (j > 0 && (heaviest[j-1] == UMAX(uint64_t) ||
                       w_cur > mat->wt[heaviest[j-1]]))
      {
        if (j < mat->nburied)
          heaviest[j] = heaviest[j-1];
        j--;
      }
      if (j < mat->nburied)
        heaviest[j] = i;
    }

    int32_t buried_max = mat->wt[heaviest[0]];
    int32_t buried_min = mat->wt[heaviest[mat->nburied-1]];

    if (buried_max == SMAX(int32_t))
      weight_buried_is_exact = 0;

    /* Compute weight of buried part of the matrix. */
    for (i = 0; i < mat->nburied; i++)
    {
      int32_t w = mat->wt[heaviest[i]];
      /* since we saturate the weights at 2^31-1, weight_buried might be less
         than the real weight of buried columns, however this can occur only
         when the number of rows exceeds 2^32-1 */
      weight_buried += w;
#if DEBUG >= 1
      fprintf(stderr, "# Burying j=%" PRIu64 " (wt = %" PRId32 ")\n",
                      heaviest[i], w);
#endif
    }
    printf("# Number of buried columns is %" PRIu64 " (min_weight=%" PRId32 ", "
           "max_weight=%" PRId32 ")\n", mat->nburied, buried_min, buried_max);
    if (weight_buried_is_exact)
      printf("# Weight of the buried part of the matrix: %" PRIu64 "\n",
             weight_buried);
    else /* weight_buried is only a lower bound */
      printf("# Weight of the buried part of the matrix is >= %" PRIu64 "\n",
             weight_buried);

    /* Remove buried columns from rows in mat structure */
    printf("# Start to remove buried columns from rels...\n");
    fflush (stdout);
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < mat->nrows; i++)
    {
      unsigned int k = 1;
      for (unsigned int j = 1; j <= matLengthRow(mat, i) ; j++)
      {
        index_t cur_id = matCell(mat, i, j);
        int32_t w = mat->wt[cur_id];
        unsigned int is_buried;
        if (w < buried_min)
          is_buried = 0;
        else if (w > buried_min)
          is_buried = 1;
        else /* w == buried_min */
        {
          unsigned int l = mat->nburied;
          is_buried = 0;
          while (l > 0 && w == mat->wt[heaviest[l-1]])
          {
            if (heaviest[l-1] == cur_id)
              is_buried = 1;
            l--;
          }
        }
        if (!is_buried) /* not a buried column */
          {
            /* we can put any exponent since we don't bury columns in DL */
            setCell(mat->rows[i], k, cur_id, 1);
            k ++;
          }
      }
      setCell(mat->rows[i], 0, k - 1, 1);
      mat->rows[i] = reallocRow (mat->rows[i], k);
    }
    printf ("# Done\n");

    /* reset to 0 the weight of the buried columns */
    for (unsigned int j = 0; j < mat->nburied; j++)
      mat->wt[heaviest[j]] = 0;

    /* compute the matrix weight */
    mat->weight = 0;
    for (i = 0; i < mat->nrows; i++)
      mat->weight += matLengthRow(mat, i);

    free (heaviest);
  }
  else

  {
    printf("# No columns were buried.\n");
    mat->weight = mat->tot_weight;
  }
  printf("# Weight of the active part of the matrix: %" PRIu64 "\n# Total "
         "weight of the matrix: %" PRIu64 "\n", mat->weight, mat->tot_weight);
  if (weight_buried_is_exact)
    ASSERT_ALWAYS (mat->weight + weight_buried == mat->tot_weight);
  else /* weight_buried is only a lower bound */
    ASSERT_ALWAYS (mat->weight + weight_buried <= mat->tot_weight);

  /* Allocate mat->R[h] if necessary */
  initMatR (mat);

  /* Re-read all rels (in memory) to fill-in mat->R */
  printf("# Start to fill-in columns of the matrix...\n");
  fflush (stdout);
  fillR (mat);
  printf ("# Done\n");
}

void
print_row(filter_matrix_t *mat, index_t i)
{
    fprintRow (stdout, mat->rows[i]);
}

void
destroyRow (filter_matrix_t *mat, index_t i)
{
    free (mat->rows[i]);
    mat->rows[i] = NULL;
}

void
remove_i_from_Rj(filter_matrix_t *mat, index_t i, index_t j)
{
  unsigned int k, n = mat->R[j][0];

  /* R[j] should not be empty */
  ASSERT(mat->R[j] != NULL);

  for (k = 1; k <= n; k++)
    if (mat->R[j][k] == i)
      {
        mat->R[j][k] = mat->R[j][n];
        mat->R[j][0] = n - 1;
        return;
      }
  ASSERT_ALWAYS(0);
}

// cell M[i, j] is incorporated in the data structure. It is used
// later on in cases where i is not already present in row[j].
void
add_i_to_Rj(filter_matrix_t *mat, index_t i, index_t j)
{
  int l;

  ASSERT(mat->R[j] != NULL);

  l = mat->R[j][0] + 1;
  reallocRj (mat, j, l);
  mat->R[j][l] = i;
}

/* return the weight of the relation obtained when adding relations i1 and i2
*/
int
weightSum(filter_matrix_t *mat, index_t i1, index_t i2, MAYBE_UNUSED index_t j)
{
  unsigned int k1, k2, w, len1, len2;

    len1 = (isRowNull(mat, i1) ? 0 : matLengthRow(mat, i1));
    len2 = (isRowNull(mat, i2) ? 0 : matLengthRow(mat, i2));
#if DEBUG >= 1
    if((len1 == 0) || (len2 == 0))
        fprintf(stderr, "i1=%d i2=%d len1=%d len2=%d\n", i1, i2, len1, len2);
#endif
#ifdef FOR_DL /* look for the exponents of j in i1 and i2*/
    int e1 = 0, e2 = 0;
    int d;
    unsigned int l;
    for (l = 1 ; l <= len1 ; l++)
        if (matCell(mat, i1, l) == j)
            e1 = mat->rows[i1][l].e;
    for (l = 1 ; l <= len2 ; l++)
        if (matCell(mat, i2, l) == j)
            e2 = mat->rows[i2][l].e;

    ASSERT (e1 != 0 && e2 != 0);

    d  = (int) gcd_int64 ((int64_t) e1, (int64_t) e2);
    e1 /= -d;
    e2 /= d;
#endif
    k1 = k2 = 1;
    w = 0;
    while((k1 <= len1) && (k2 <= len2))
    {
        if(matCell(mat, i1, k1) < matCell(mat, i2, k2))
        {
#ifdef FOR_DL
            w += (e2 * mat->rows[i1][k1].e == 0) ? 0 : 1;
#else
            w++;
#endif
            k1++;
        }
        else if(matCell(mat, i1, k1) > matCell(mat, i2, k2))
        {
#ifdef FOR_DL
            w += (e1 * mat->rows[i2][k2].e == 0) ? 0 : 1;
#else
            w++;
#endif
            k2++;
        }
        else
        {
#ifdef FOR_DL
            w += (e2*mat->rows[i1][k1].e + e1*mat->rows[i2][k2].e == 0) ? 0 : 1;
#endif
            k1++;
            k2++;
        }
    }
#ifdef FOR_DL
    // finish with k1
    for( ; k1 <= matLengthRow(mat, i1); k1++)
      w += (e2 * mat->rows[i1][k1].e == 0) ? 0 : 1;
    // finish with k2
    for( ; k2 <= matLengthRow(mat, i2); k2++)
      w += (e1 * mat->rows[i2][k2].e == 0) ? 0 : 1;
#else
    w += matLengthRow(mat, i1) + 1 - k1;
    w += matLengthRow(mat, i2) + 1 - k2;
#endif
    return w;
}

/* put in ind[0]..ind[m-1] the indices of the m (active) rows containing j */
void
fillTabWithRowsForGivenj(index_t *ind, filter_matrix_t *mat, index_t j)
{
  int ni = 0;
  unsigned int k;

  for (k = 1; k <= mat->R[j][0]; k++)
    ind[ni++] = mat->R[j][k];
}
