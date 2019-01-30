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

/*****************************************************************************/

/* Initialize the sparse matrix mat. */
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

int cmp_u64(const uint64_t * a, const uint64_t * b)
{
    return (*a > *b) - (*b > *a);
}

void
print_row(filter_matrix_t *mat, index_t i)
{
    fprintRow (stdout, mat->rows[i]);
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
