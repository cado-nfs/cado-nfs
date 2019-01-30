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
#include "mst.h"

/***************** memory allocation on R[j] *********************************/

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

/* Initialize the sparse matrix mat.
   If initR is 0, does not initialize R and the heap. */
void
initMat (filter_matrix_t *mat, uint32_t skip)
{
  /* we start with cwmax = 2, and we increase it in mergeOneByOne() when
     the Markowitz queue is empty */
  mat->cwmax = 2;
  ASSERT_ALWAYS (mat->cwmax < 255); /* 255 is reserved for saturated values */
  mat->skip = skip;

  mat->weight = 0;
  mat->tot_weight = 0;
  mat->rem_ncols = 0;

  mat->rows = (typerow_t **) malloc (mat->nrows * sizeof (typerow_t *));
  ASSERT_ALWAYS (mat->rows != NULL);
  mat->wt = (int32_t *) malloc (mat->ncols * sizeof (int32_t));
  ASSERT_ALWAYS (mat->wt != NULL);
  memset (mat->wt, 0, mat->ncols * sizeof (int32_t));
}

void
clearMat (filter_matrix_t *mat)
{
  uint64_t i, j;

#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
  for (i = 0; i < mat->nrows; i++)
    free (mat->rows[i]);
  free (mat->rows);
  free (mat->wt);
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
  for (j = 0; j < mat->ncols; j++)
    free (mat->R[j]);
  free (mat->R);
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
    mat->rem_ncols += (mat->wt[h] == 0);
    mat->wt[h] += (mat->wt[h] != SMAX(int32_t));
    if (h < mat->skip) continue;
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

  /* do as if the coefficients we're discarding early on are still here
   */
  mat->tot_weight += rel->nb;

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

int cmp_u64(const uint64_t * a, const uint64_t * b)
{
    return (*a > *b) - (*b > *a);
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
