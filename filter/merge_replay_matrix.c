#include "cado.h"

#include <string.h>
#include <stdio.h>

#include "portability.h"
#include "utils_with_io.h"
#include "filter_config.h"
#include "merge_replay_matrix.h"
#include "sparse.h"

int
decrS (int w)
{
  return (w >= 0 ? w - 1 : w + 1);
}

/* increment the absolute value of a signed number */
int
incrS (int w)
{
  return (w >= 0 ? w + 1 : w - 1);
}

/* initialize the sparse matrix mat */
void
initMat (filter_matrix_t *mat, int maxlevel, uint32_t keep,
         uint32_t nburied)
{
  mat->keep  = keep;
  mat->mergelevelmax = maxlevel;
  mat->cwmax = 2 * maxlevel;
  ASSERT_ALWAYS (mat->cwmax < 255);
  mat->nburied = nburied;

  mat->weight = 0;
  mat->tot_weight = 0;
  mat->rem_ncols = 0;

  mat->rows = (typerow_t **) malloc (mat->nrows * sizeof (typerow_t *));
  ASSERT_ALWAYS (mat->rows != NULL);
  mat->wt = (int32_t *) malloc (mat->ncols * sizeof (int32_t));
  memset (mat->wt, 0, mat->ncols * sizeof (int32_t));
  mat->R = (index_t **) malloc (mat->ncols * sizeof(index_t *));
  ASSERT_ALWAYS(mat->R != NULL);
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
    free (mat->R[j]);
  free (mat->R);
}

/* initialize Rj[j] for light columns, i.e., for those of weight <= cwmax */
void
InitMatR (filter_matrix_t *mat)
{
  index_t h;
  int32_t w;

  for(h = 0; h < mat->ncols; h++)
  {
    w = mat->wt[h];
    if (w <= mat->cwmax)
    {
      mat->R[h] = (index_t *) malloc((w + 1) * sizeof(index_t));
      FATAL_ERROR_CHECK(mat->R[h] == NULL, "Cannot allocate memory");
      mat->R[h][0] = 0; /* last index used */
    }
    else /* weight is larger than cwmax */
    {
      mat->wt[h] = -mat->wt[h]; // trick!!!
      mat->R[h] = NULL;
    }
  }
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
        {
          free (mat->R[h]);
          mat->R[h] = NULL;
        }
        if (mat->wt[h] > 0)
          mat->wt[h] = -mat->wt[h]; // trick!!!
      }
    }
  }

  fclose_maybe_compressed (file, infilename);
}

/* callback function called by filter_rels */
void * insert_rel_into_table (void *context_data, earlyparsed_relation_ptr rel)
{
  filter_matrix_t *mat = (filter_matrix_t *) context_data;

  /* XXX For now can't use my_malloc, because rows are realloc later */
  mat->rows[rel->num] = (typerow_t*) malloc ((rel->nb + 1) * sizeof (typerow_t));
  FATAL_ERROR_CHECK(mat->rows[rel->num] == NULL, "Cannot allocate memory");
  matLengthRow(mat, rel->num) = rel->nb;

  for (unsigned int i = 0; i < rel->nb; i++)
  {
    index_t h = rel->primes[i].h;
    exponent_t e = rel->primes[i].e;
    mat->tot_weight++;
    /* For factorization, they should not be any multiplicity here.
       For DL we do not want to count multiplicity in mat->wt */
#ifndef FOR_DL
    ASSERT_ALWAYS (e == 1);
#endif
    if (mat->wt[h] == 0)
    {
      mat->wt[h] = 1;
      mat->rem_ncols++;
    }
    else if (mat->wt[h] != SMAX(int32_t))
      mat->wt[h]++;

    setCell(mat->rows[rel->num][i+1], h, e);
  }

  /* sort indices to ease row merges */
  qsort(&(mat->rows[rel->num][1]), rel->nb, sizeof(typerow_t), cmp_typerow_t);

  return NULL;
}


void
filter_matrix_read (filter_matrix_t *mat, const char *purgedname)
{
  uint64_t nread, i, h;
  char *fic[2] = {(char *) purgedname, NULL};

  /* read all rels */
  nread = filter_rels(fic, (filter_rels_callback_t) &insert_rel_into_table, mat,
                      EARLYPARSE_NEED_INDEX, NULL, NULL);
  ASSERT_ALWAYS(nread == mat->nrows);
  mat->rem_nrows = nread;

  /* print weight count */
  uint64_t *nbm;
  nbm = (uint64_t *) malloc ((mat->mergelevelmax + 1) * sizeof (uint64_t));
  memset (nbm, 0, (mat->mergelevelmax + 1) * sizeof (uint64_t));
  for (h = 0; h < mat->ncols; h++)
    if (mat->wt[h] <= mat->mergelevelmax)
      nbm[mat->wt[h]]++;
  for (h = 1; h <= (uint64_t) mat->mergelevelmax; h++)
    printf ("There are %" PRIu64 " column(s) of weight %" PRIu64 "\n", nbm[h], h);
  ASSERT_ALWAYS(mat->rem_ncols == mat->ncols - nbm[0]);
  free (nbm);

  /* Bury heavy coloumns. The 'nburied' heaviest column are buried. */
  /* Buried columns are not took into account by merge. */
  uint64_t weight_buried = 0;
  if (mat->nburied)
  {
    uint64_t *heaviest = NULL;
    heaviest = (uint64_t *) malloc (mat->nburied * sizeof (uint64_t));
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

    /* Compute weight of buried part of the matrix. */
    for (i = 0; i < mat->nburied; i++)
    {
      int32_t w = mat->wt[heaviest[i]];
      weight_buried += w;
#if DEBUG >= 1
      fprintf(stderr, "# Burying j=%" PRIu64 " (wt = %" PRId32 ")\n",
                      heaviest[i], w);
#endif
    }
    printf("# Number of buried columns is %" PRIu64 " (min_weight=%" PRId32 ", "
           "max_weight=%" PRId32 ")\n", mat->nburied, buried_min, buried_max);
    printf("# Weight of the buried part of the matrix: %" PRIu64 "\n",
           weight_buried);

    /* Remove buried column from rows in mat structure */
    printf("# Start to remove buried columns from rels...\n");
    fflush (stdout);
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
          mat->weight++;
          mat->rows[i][k++] = mat->rows[i][j];
        }
      }

      matLengthRow(mat, i) = k-1;
      mat->rows[i] = (typerow_t*) realloc (mat->rows[i], k * sizeof (typerow_t));
    }
    printf ("# Done\n");

    free (heaviest);
  }
  else
  {
    printf("# No columns were buried.\n");
    mat->weight = mat->tot_weight;
  }
  printf("# Weight of the active part of the matrix: %" PRIu64 "\n# Total "
         "weight of the matrix: %" PRIu64 "\n", mat->weight, mat->tot_weight);
  ASSERT_ALWAYS (mat->weight + weight_buried == mat->tot_weight);

  /* Allocate mat->R[h] if necessary */
  InitMatR (mat);

  /* Re-read all rels (in memory) to fill-in mat->R */
  printf("# Start to fill-in columns of the matrix...\n");
  fflush (stdout);
  for (i = 0; i < mat->nrows; i++)
  {
    for(unsigned int k = 1 ; k <= matLengthRow(mat, i); k++)
    {
      h = matCell(mat, i, k);
      if(mat->wt[h] > 0)
      {
        mat->R[h][0]++;
        mat->R[h][mat->R[h][0]] = i;
      }
    }
  }
    printf ("# Done\n");
}

void
print_row(filter_matrix_t *mat, int i)
{
    fprintRow (stdout, mat->rows[i]);
}

void
destroyRow(filter_matrix_t *mat, int i)
{
    free(mat->rows[i]);
    mat->rows[i] = NULL;
}

void
freeRj(filter_matrix_t *mat, int j)
{
    free(mat->R[j]);
    mat->R[j] = NULL;
}

// Don't touch to R[j][0]!!!!
void
remove_i_from_Rj(filter_matrix_t *mat, int i, int j)
{
  // be dumb for a while
  unsigned int k;

  if(mat->R[j] == NULL)
  {
	  fprintf(stderr, "Row %d already empty\n", j);
	  return;
  }
  for(k = 1; k <= mat->R[j][0]; k++)
	  if((int) mat->R[j][k] == i)
    {
#if DEBUG >= 2
      fprintf(stderr, "Removing row %d from R[%d]\n", i, j);
#endif
      mat->R[j][k] = UMAX(index_t);
      break;
	  }
}

// cell M[i, j] is incorporated in the data structure. It is used
// later on in cases where i is not already present in row[j].
void
add_i_to_Rj(filter_matrix_t *mat, int i, int j)
{
  // be semi-dumb for a while
  unsigned int k;

#if DEBUG >= 2
  fprintf(stderr, "Adding row %d to R[%d]\n", i, j);
#endif
  for(k = 1; k <= mat->R[j][0]; k++)
    if(mat->R[j][k] == UMAX(index_t))
      break;
  if(k <= mat->R[j][0])
  {
  // we have found a place where it is -1
  mat->R[j][k] = i;
  }
  else
  {
#if DEBUG >= 2
    fprintf(stderr, "WARNING: reallocing things in add_i_to_Rj for R[%d]\n", j);
#endif
    int l = mat->R[j][0]+2;
    mat->R[j] = (index_t *)realloc(mat->R[j], l * sizeof(index_t));
    mat->R[j][l-1] = i;
    mat->R[j][0] = l-1;
  }
}

// Remove j from R[i] and crunch R[i].
void
remove_j_from_row(filter_matrix_t *mat, int i, int j)
{
  unsigned int k;

#if TRACE_ROW >= 0
    if(i == TRACE_ROW){
	fprintf(stderr, "TRACE_ROW: remove_j_from_row i=%d j=%d\n", i, j);
    }
#endif
#if DEBUG >= 2
    fprintf(stderr, "row[%d]_b=", i);
    print_row(mat, i);
    fprintf(stderr, "\n");
#endif
    for(k = 1; k <= matLengthRow(mat, i); k++)
	if((int) matCell(mat, i, k) == j)
	    break;
    ASSERT(k <= matLengthRow(mat, i));
    // crunch
    for(++k; k <= matLengthRow(mat, i); k++)
        matCell(mat, i, k-1) = matCell(mat, i, k);
    matLengthRow(mat, i) -= 1;
    mat->weight--;
#if DEBUG >= 2
    fprintf(stderr, "row[%d]_a=", i);
    print_row(mat, i);
    fprintf(stderr, "\n");
#endif
}

/* return the weight of the relation obtained when adding relations i1 and i2
*/
int
weightSum(filter_matrix_t *mat, int i1, int i2, MAYBE_UNUSED int32_t j)
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
        if ((int) matCell(mat, i1, l) == j)
            e1 = mat->rows[i1][l].e;
    for (l = 1 ; l <= len2 ; l++)
        if ((int) matCell(mat, i2, l) == j)
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
    // finish with k1
    for( ; k1 <= matLengthRow(mat, i1); k1++)
#ifdef FOR_DL
            w += (e2 * mat->rows[i1][k1].e == 0) ? 0 : 1;
#else
            w++;
#endif
    // finish with k2
    for( ; k2 <= matLengthRow(mat, i2); k2++)
#ifdef FOR_DL
            w += (e1 * mat->rows[i2][k2].e == 0) ? 0 : 1;
#else
            w++;
#endif
    return w;
}

/* put in ind[0]..ind[m-1] the indices of the m (active) rows containing j,
   and return the total weight of all those rows */
int
fillTabWithRowsForGivenj(int32_t *ind, filter_matrix_t *mat, int32_t j)
{
  int ni = 0, w = 0;
  unsigned int k;
  index_t i;

  for (k = 1; k <= mat->R[j][0]; k++)
    if ((i = mat->R[j][k]) != UMAX(index_t))
      {
        ind[ni++] = i;
        w += matLengthRow(mat, i);
        if (ni == mat->wt[j])
          break;
      }
  return w;
}
