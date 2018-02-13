// TODO: use unified compact lists...!
// TODO: reintroduce mat->weight

/* try to align C and C++ code */
#define MIMIC

/* If one wants to change the R data structure, please check the diff of
   revision 1568, which clearly identifies the places where R is used. */

#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "portability.h"
#include "utils_with_io.h"

#include "filter_config.h"
#include "merge_replay_matrix.h"
#include "sparse.h"
#include "mst.h"
#include "report.h"

#include "markowitz.h"
#include "merge_mono.h"

#define DEBUG 0

#if DEBUG >= 1
/* total number of row additions in merge */
static unsigned long row_additions = 0;
#endif

/* Given an ideal of weight m, returns the best index i, 0 <= i < m,
   such that merging the i-th relation containing the ideal with all
   others (m-1) gives the smallest total weight.
   Note: for m >= 4 this might not be the best possible merge. For example
   for m=4 the best possible merge might be: 0-1, 1-2, 2-3 which is not
   of this shape (but for m <= 3 all possible merges are of this shape). */
static int
findBestIndex(filter_matrix_t *mat, int m, index_t *ind, index_t ideal)
{
    /* not mallocing A[][] to speed up things */
    int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX], i, j, imin, wmin, w;

    ASSERT(m <= MERGE_LEVEL_MAX);
    if (m == 2)
      return 0;
    fillRowAddMatrix(A, mat, m, ind, ideal);
    // iterate over all vertices
    imin = -1;
    wmin = SMAX(int);
    for(i = 0; i < m; i++){
	// compute the new total weight if i is used as pivot
	w = 0;
	for(j = 0; j < m; j++)
	    // note A[i][i] = 0
	    w += A[i][j];
#if DEBUG >= 1
	printf ("W[%d]=%d\n", i, w);
#endif
	if (w < wmin)
          {
	    imin = i;
	    wmin = w;
          }
    }
    return imin;
}

#if DEBUG >= 2
// We should have mat->wt[j] == m == nchk
static void
checkCoherence(filter_matrix_t *mat, int m, index_t j)
{
    int nchk = 0, k;

    nchk = mat->R[j][0];
    ASSERT(nchk == (mat->wt[j] >= 0 ? mat->wt[j] : -mat->wt[j]));
    if(m != -1){
	if(nchk != m){
	    printf ("HYPERCHECK: mat->R[%d][0]=%ld, m=%d\n", j,
                     (long int) mat->R[j][0], m);
	    printf ("Gasp: nchk=%d\n", nchk);
	}
	ASSERT(nchk == m);
    }
}
#endif

/* Check that the rows are sorted, and contain only once each ideal.
   This should normally always be the case, but we have to prevent from
   users who build their own corrupted purge file,
   see https://gforge.inria.fr/tracker/index.php?func=detail&aid=21330. */
static void
checkMat (filter_matrix_t *mat)
{
  for (index_t i = 0; i < mat->nrows; i++)
    for (index_t k = 2; k <= matLengthRow(mat, i); k++)
      if (!(matCell(mat, i, k - 1) < matCell(mat, i, k)))
        {
          fprintf (stderr, "Error, invalid row %lu:\n",
                   (unsigned long) i);
          fprintRow (stderr, mat->rows[i]);
          fprintf (stderr, "\nPlease fix your purge file\n");
          exit (1);
        }
}

//////////////////////////////////////////////////////////////////////
// making things independent of the real data structure used
//////////////////////////////////////////////////////////////////////

/* remove column j and update matrix */
static void
removeColumnAndUpdate(filter_matrix_t *mat, index_t j)
{
    MkzRemoveJ (mat, j);
    freeRj (mat, j);
}

// The cell [i, j] may be incorporated to the data structure, at least
// if j is not too heavy, etc.
static void
addCellAndUpdate(filter_matrix_t *mat, index_t i, index_t j)
{
#if TRACE_COL >= 0
    int oldw = mat->wt[j];
#endif
    int w = MkzIncrCol(mat, j);

    /* FIXME: if an ideal has weight exactly cwmax before a merge,
       and keeps the same weight afterwards, depending on the order in which
       we apply the row[i1] += row[i2] operations, it might be that its weight
       exceeds cwmax in the middle of the merge, and thus it will be disabled,
       and we will lose a possible merge. */

    if (w >= 0)
      {
	if(w > mat->cwmax)
          {
	    // the weight of column w exceeds cwmax, thus we remove it
	    removeColumnAndUpdate(mat, j);
	    mat->wt[j] = -w;
          }
	else // update R[j] by adding i
          add_i_to_Rj (mat, i, j);
      }
#if TRACE_COL >= 0
    if (j == TRACE_COL && mat->wt[j] != oldw)
      printf ("TRACE_COL: addCellAndUpdate changed weight of %d from %d to %d\n",
              j, oldw, mat->wt[j]);
#endif
}

// remove the cell (i,j), and updates matrix correspondingly.
static void
removeCellAndUpdate(filter_matrix_t *mat, index_t i, index_t j)
{
#if TRACE_ROW >= 0
    if(i == TRACE_ROW){
	printf ("TRACE_ROW: removeCellAndUpdate i=%d j=%d\n", i, j);
    }
#endif
    if(mat->wt[j] < 0){
	// if mat->wt[j] is already < 0, we don't care about
	// decreasing, updating, etc. except when > 2
	return;
    }
#if TRACE_COL >= 0
    int oldw = mat->wt[j];
#endif
    MkzDecreaseColWeight(mat, j);
#if TRACE_COL >= 0
    if (j == TRACE_COL && mat->wt[j] != oldw)
      printf ("TRACE_COL: removeCellAndUpdate changed weight of %d from %d to %d\n",
              j, oldw, mat->wt[j]);
#endif
    // update R[j] by removing i
    remove_i_from_Rj(mat, i, j);
}

//////////////////////////////////////////////////////////////////////
// now, these guys are generic...!

// for all j in row[i], removes j and update data
// if final, also update the Markowitz counts, free memory, and remove from
// heap of heavy rows
static void
removeRowAndUpdate(filter_matrix_t *mat, index_t i, int final)
{
    unsigned int k;

#if TRACE_ROW >= 0
    if(i == TRACE_ROW)
	printf ("TRACE_ROW: removeRowAndUpdate i=%d\n", i);
#endif
    mat->weight -= matLengthRow(mat, i);
    for(k = 1; k <= matLengthRow(mat, i); k++){
#if TRACE_COL >= 0
      unsigned int w = mat->wt[matCell(mat, i, k)];
#endif
	removeCellAndUpdate(mat, i, matCell(mat, i, k));
#if TRACE_COL >= 0
	if(matCell(mat, i, k) == TRACE_COL)
          printf ("TRACE_COL: removeRowAndUpdate removes %d from R_%d (weight %d -> %d)\n",
                  TRACE_COL, i, w, mat->wt[TRACE_COL]);
#endif
    }

    if (final)
      {
        /* free the memory occupied by row i */
        destroyRow (mat, i);
        mat->rem_nrows--;
        /* remove row i from the heap of heavy rows */
        heap_delete (mat->Heavy, mat, i);
      }
}

#ifdef FOR_DL
// All entries M[i, j] are potentially added to the structure.
static void
addOneRowAndUpdate(filter_matrix_t *mat, index_t i)
{
  unsigned int k;

  mat->weight += matLengthRow(mat, i);
#if TRACE_ROW >= 0
  if (i == TRACE_ROW)
    printf ("TRACE_ROW: addOneRowAndUpdate i=%d\n", i);
#endif
  for(k = 1; k <= matLengthRow(mat, i); k++)
    {
#if TRACE_COL >= 0
      if (matCell(mat, i, k) == TRACE_COL)
        printf ("TRACE_ROW: addOneRowAndUpdate i=%d j=%d (weight %d -> %d)\n",
                i, matCell(mat, i, k), mat->wt[matCell(mat, i, k)],
                mat->wt[matCell(mat, i, k)] + 1);
#endif
      addCellAndUpdate(mat, i, matCell(mat, i, k));
    }
}

// realize mat[i1] += mat[i2] and update the data structure.
// len could be the real length of row[i1]+row[i2] or -1.
static void
addRowsAndUpdate(filter_matrix_t *mat, index_t i1, index_t i2, index_t j)
{
    // cleaner one, that shares addRowsData() to prepare the next move...!
    // i1 is to disappear, replaced by a new one
    removeRowAndUpdate(mat, i1, 0); /* fake removal */
    // we know the length of row[i1]+row[i2]
    addRows(mat->rows, i1, i2, j);
    addOneRowAndUpdate(mat, i1);
}
#else
/* special version for factorization */
static void
addRowsAndUpdate (filter_matrix_t *mat, index_t i1, index_t i2,
                  index_t j MAYBE_UNUSED)
{
    uint32_t k, k1, k2, l1, l2, len;
    typerow_t **rows = mat->rows, *tmp, *r1 = rows[i1], *r2 = rows[i2];

    l1 = rowLength (rows, i1);
    l2 = rowLength (rows, i2);

    /* removeRowAndUpdate(mat, i1, 0) */
    mat->weight -= l1;
    for (k1 = 1; k1 <= l1; k1++)
      removeCellAndUpdate (mat, i1, rowCell(r1, k1));

    /* addRows(mat->rows, i1, i2, j), i.e.,
       addRowsUpdateIndex(mat->rows, NULL, i1, i2, j) */
    len = 1 + l1 + l2;
    tmp = (typerow_t *) malloc (len * sizeof(typerow_t));
    k = k1 = k2 = 1;

    while ((k1 <= l1) && (k2 <= l2))
      if (rowCell(r1, k1) < rowCell(r2, k2))
        tmp[k++] = rowCell(r1, k1), k1++;
      else if (rowCell(r1, k1) > rowCell(r2, k2))
        tmp[k++] = rowCell(r2, k2), k2++;
      else
        k1++, k2++;
    for(; k1 <= l1;)
      tmp[k++] = rowCell(r1, k1), k1++;
    for(; k2 <= l2;)
      tmp[k++] = rowCell(r2, k2), k2++;
    tmp[0] = k - 1;
#if defined(FOR_DLP) || __SIZEOF_INDEX__ == 4 || __SIZEOF_INDEX__ == 8
    free (r1);
    rows[i1] = (k == len) ? tmp : realloc (tmp, k * sizeof(typerow_t));
#else
    if (k-1 != l1)
      rows[i1] = realloc (rows[i1], k * sizeof(typerow_t));
    compressRow (rows[i1], tmp, k - 1);
    free (tmp);
#endif

    /* addOneRowAndUpdate(mat, i1) */
    mat->weight += k - 1;
    for (r1 = rows[i1], k1 = 1; k1 < k; k1++)
      addCellAndUpdate (mat, i1, rowCell(r1, k1));
}
#endif

static void
removeRowDefinitely(report_t *rep, filter_matrix_t *mat, index_t i)
{
    report1(rep, i, -1);
    removeRowAndUpdate(mat, i, 1);
}

/* Remove column j, and update the matrix.
   This function is only called for singleton ideals. */
static void
removeColDefinitely(report_t *rep, filter_matrix_t *mat, index_t j)
{
  ASSERT(mat->R[j][0] == 1);

#if TRACE_COL >= 0 || TRACE_ROW >= 0
  if (j == TRACE_COL || mat->R[j][1] == (index_t) TRACE_ROW)
    printf ("TRACE: removeColDefinitely j=%d: remove row %d, weight %d\n",
            j, mat->R[j][1], mat->wt[j]);
#endif
  removeRowDefinitely (rep, mat, mat->R[j][1]);
  removeColumnAndUpdate (mat, j);
  mat->rem_ncols --;
}

static int
removeSingletons(report_t *rep, filter_matrix_t *mat)
{
  index_t j;
  int nr=0;

  for (j = 0; j < mat->ncols; j++)
    if (mat->wt[j] == 1) {
      removeColDefinitely (rep, mat, j);
      nr++;
    }
  return nr;
}

/* Try all combinations of merging one row with all others (m-1) ones
   to find the smaller one; resists to m==1.
 */
static void
tryAllCombinations(report_t *rep, filter_matrix_t *mat, int m, index_t *ind,
                   index_t j)
{
    int i, k;

    if(m == 1)
      printf ("Warning: m=1 in tryAllCombinations\n");
    i = findBestIndex(mat, m, ind, j);
#if DEBUG >= 1
    printf ("Minimal is i=%d (%d %d)\n", i, ind[0], ind[1]);
#endif
    for(k = 0; k < m; k++)
	if(k != i){
	    addRowsAndUpdate(mat, ind[k], ind[i], j);
#if DEBUG >= 1
	    printf ("new row[%d]=", ind[k]);
	    print_row (mat, ind[k]);
	    printf ("\n");
#endif
	}

    /* swap ind[0] and ind[i] */
    index_t itmp = ind[i];
    ind[i] = ind[0];
    ind[0] = itmp;

#if DEBUG >= 1
    printf ("=> new_ind: %d %d\n", ind[0], ind[1]);
#endif
    reportn(rep, (index_signed_t*) ind, m, j);
    removeRowAndUpdate(mat, ind[0], 1);
}

/*
  Input: (father[i], sons[i]) for 0 <= i < m-1 are the edges of the
  minimal spanning tree.
  Output: perform the corresponding merges, and stores them in the history.
  father[i], sons[i] are integers in [0..m-1]
*/
int
addFatherToSons(index_t history[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1],
		filter_matrix_t *mat, int m, index_t *ind, index_t ideal,
		int *father, int *sons)
{
  int i, s, t;

  /* we process the entries in backward order, to avoid modifying a parent
     before adding it to one of its sons, but we put the entries in history[]
     in increasing order: history[] is processed backwards in MSTWithA() */
  for (i = m - 2; i >= 0; i--)
    {
      s = father[i];
      t = sons[i];
      if (i == 0)
        {
          history[i][1] = ind[s]; /* to destroy root */
          ASSERT(s == 0);
        }
      else
        history[i][1] = -(ind[s] + 1); /* not to destroy row ind[s] */
      addRowsAndUpdate (mat, ind[t], ind[s], ideal);
      history[i][2] = ind[t];
      history[i][0] = 2;
    }
  return m - 2;
}

static void
MSTWithA (report_t *rep, filter_matrix_t *mat, int m, index_t *ind, index_t j,
          int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX])
{
    index_t history[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1];
    int start[MERGE_LEVEL_MAX], end[MERGE_LEVEL_MAX];
    int hmax, i;

    minimalSpanningTree (start, end, m, A);
#if DEBUG >= 1
    printMST (start, end, m, ind);
#endif
    hmax = addFatherToSons (history, mat, m, ind, j, start, end);
    for(i = hmax; i >= 0; i--)
      reportn (rep, (index_signed_t*) (history[i]+1), history[i][0], j);
    removeRowAndUpdate(mat, ind[0], 1);
}

static void
useMinimalSpanningTree (report_t *rep, filter_matrix_t *mat, int m,
                        index_t *ind, index_t j)
{
    int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX];

    fillRowAddMatrix (A, mat, m, ind, j);
    MSTWithA (rep, mat, m, ind, j, A);
}

/* return the list of all alive ideals involved in merge of column j, and
   put its length in 'len' */
static index_t*
get_alive_ideals (filter_matrix_t *mat, int m, index_t *ind, index_t j,
                  int *len)
{
  index_t *l, jj;
  int k, s, t, n;

  *len = -m; /* we don't consider ideal j */
  for (k = 0; k < m; k++)
    *len += matLengthRow (mat, ind[k]);
  l = malloc (*len * sizeof (index_t));
  for (s = k = 0; k < m; k++)
    {
      index_t i = ind[k];
      n = matLengthRow (mat, i);
      for (t = 1; t <= n; t++)
        {
          jj = matCell(mat, i, t);
          if (jj != j && MkzIsAlive (mat->MKZA, jj))
            l[s++] = jj;
        }
    }
  ASSERT_ALWAYS(s <= *len);

  /* sort list elements */
  qsort (l, s, sizeof(index_t), cmp_index);

  /* remove duplicates */
  for (k = t = 0; t < s; t++)
    /* invariant: l[0] < l[1] < ... < l[k-1] */
    if (k == 0 || l[k-1] < l[t])
      l[k++] = l[t];

  *len = k;
  return l;
}

static void
findOptimalCombination (report_t *rep, filter_matrix_t *mat, int m,
                        index_t *ind, index_t j)
{
  index_t *l;
  int len, i, k;

  /* get the list of all alive ideals involved in this merge */
  l = get_alive_ideals (mat, m, ind, j, &len);

  /* we can use here two algorithms:
     (a) tryAllCombinations tries to merge row i with all other (m-1) rows,
         for each i, 0 <= i < m, and keeps the best i
     (b) useMinimalSpanningTree computes a minimal spanning tree
     Both have complexity O(m^2), and useMinimalSpanningTree is always better
     or equal. Moreover useMinimalSpanningTree is better to find cancellations
     (if any) between rows. */

  if (m <= 2)
    tryAllCombinations (rep, mat, m, ind, j);
  else
    useMinimalSpanningTree (rep, mat, m, ind, j);

  /* compute all still alive ideals in l */
  for (i = k = 0; k < len; k++)
    if (MkzIsAlive (mat->MKZA, l[k]))
      l[i++] = l[k];
  /* update the Markowitz costs of the remaining ideals */
  MkzUpdateN (mat, l, i);

  free (l);
}

#if DEBUG >= 1
static void
checkWeights (filter_matrix_t *mat)
{
  int *W, k, minw = INT_MAX;
  index_t i, j;

  W = (int*) malloc (mat->ncols * sizeof(int));
  for (j = 0; j < mat->ncols; j++)
    W[j] = 0;
  for(i = 0; i < mat->nrows; i++)
    if(!isRowNull(mat, i))
      {
        for(k = 1; k <= lengthRow(mat, i); k++)
          {
            j = cell(mat, i, k);
            ASSERT_ALWAYS(0 <= j && j < mat->ncols);
            W[j] ++;
          }
      }
  for (j = 0; j < mat->ncols; j++)
    {
      static int count = 0;
      if (W[j] != mat->wt[j])
        {
          if (mat->wt[j] >= 0)
            {
              fprintf (stderr, "Wrong weight for column %d: expected %d, got %d\n",
                       j, mat->wt[j], W[j]);
              exit (1);
            }
          else if (W[j] <= mat->mergelevelmax)
            {
              if (count++ <= 10)
                fprintf (stderr, "Warning, column %d: wt=%d W=%d\n",
                         j, mat->wt[j], W[j]);
            }
        }
      if (W[j] == 1 && mat->wt[j] == 1)
        {
          fprintf (stderr, "Error, weight=1 for column %d\n", j);
          exit (1);
        }
      if (W[j] != 0 && mat->wt[j] >= 0 && W[j] < minw)
        minw = W[j];
      if (W[j] == 2 && mat->wt[j] == 2)
        {
          static int count = 0;
          if (count ++ <= 10)
            fprintf (stderr, "Warning, column %d: wt=%d W=%d\n",
                     j, mat->wt[j], W[j]);
        }
    }
  printf ("Minimum non-zero column weight: %d\n", minw);
  free (W);
}
#endif


// j has weight m, which should coherent with mat->wt[j] == m
static void
mergeForColumn (report_t *rep, filter_matrix_t *mat, int m, index_t j)
{
    index_t ind[MERGE_LEVEL_MAX];
    int ni;

#if TRACE_COL >= 0
    if (j == TRACE_COL)
      printf ("TRACE_COL: merge column %lu\n", (unsigned long) j);
#endif

    ASSERT(mat->wt[j] == m);
    ASSERT(sizeof(index_t) == sizeof(mat->R[j][0]));

    /* each m-merge leads to m-1 additions of rows */
#if DEBUG >= 1
    row_additions += m - 1;
    if (m > MERGE_LEVEL_MAX)
      {
	fprintf (stderr, "Error: m=%d > MERGE_LEVEL_MAX=%d\n",
                 m, MERGE_LEVEL_MAX);
	exit(1);
      }
    printf ("Treating column %d of weight %d\n", j, mat->wt[j]);
#if DEBUG >= 2
    printf ("Status before next j=%d to start\n", j);
    // the corresponding rows are in R[j], skipping 1st cell
    checkCoherence(mat, m, j);
#endif
#endif

    memcpy (ind, mat->R[j] + 1, mat->R[j][0] * sizeof (index_t));

    /* now ind[0], ..., ind[m-1] are the m rows containing j */

#if DEBUG >= 1
    printf (" %d", j);
    printf (" => the %d rows are:\n", m);
    for(int k = 0; k < m; k++){
	printf ("row[%d]=", ind[k]);
	print_row (mat, ind[k]);
	printf ("\n");
    }
    printf ("\n");
#endif

    findOptimalCombination (rep, mat, m, ind, j);
    mat->rem_ncols--;
    removeColumnAndUpdate (mat, j);

    /* update the status of rows in heap of heavy rows */
    for (ni = 0; ni < m; ni++)
      if (mat->rows[ind[ni]] != NULL) /* row is still active */
        heap_push (mat->Heavy, mat, ind[ni]);
}

// Delete at most niremmax superfluous rows such that nrows-ncols >= keep.
// Let's say we are at merge-level m.
// Return the number of removed rows.
static int
deleteSuperfluousRows (report_t *rep, filter_matrix_t *mat,
                       int niremmax, int m)
{
  int keep = mat->keep;
  int nirem = 0;

  if (m <= 2)
    return 0; /* it is not worth removing rows during level-2 merges */

  if ((int) (mat->rem_nrows - mat->rem_ncols) <= keep)
    return 0;

  if (niremmax > (int) (mat->rem_nrows - mat->rem_ncols) - keep)
    niremmax = (mat->rem_nrows - mat->rem_ncols) - keep;

  printf ("Removing %d excess rows; current nrows=%" PRIu64 "\n",
          niremmax, mat->rem_nrows);
  for (nirem = 0; nirem < niremmax && mat->Heavy->size > 0; nirem++)
    {
      index_t i = heap_pop (mat->Heavy, mat);
      removeRowDefinitely (rep, mat, i);
    }
  return nirem;
}

#define REPORT_INCR 3.0

static int
number_of_superfluous_rows (filter_matrix_t *mat)
{
    int excess = mat->rem_nrows - mat->rem_ncols;

    ASSERT(excess >= 0);

    /* the following formula is close to optimal experimentally for RSA-100,
       with an initial excess of abut 100%, maxlevel=20 and a target density
       of 170 (with USE_WEIGHT_LAMBDA = 0) */
    return (int) (0.145 * (double) excess * REPORT_INCR / 3.0);
}

static void
print_memory_usage (filter_matrix_t *mat)
{
  unsigned long mem_rows, mem_wt, mem_cols, mem_MKZ;

  mem_rows = mat->nrows * sizeof (typerow_t*);
  for (unsigned long i = 0; i < mat->nrows; i++)
    if (mat->rows[i] != NULL)
#if defined(FOR_DL) || __SIZEOF_INDEX__ == 4 || __SIZEOF_INDEX__ == 8
      mem_rows += (rowLength(mat->rows, i) + 1) * sizeof (typerow_t);
#else
      mem_rows += (rowLength(mat->rows, i) + 1) * __SIZEOF_INDEX__;
#endif
  mem_wt = mat->ncols * sizeof (int32_t);
  mem_cols = mat->ncols * sizeof (index_t*);
  for (unsigned long j = 0; j < mat->ncols; j++)
    if (mat->R[j] != NULL)
      mem_cols += (mat->R[j][0] + 1) * sizeof (index_t);
  mem_MKZ = (mat->MKZQ[0] + 1) * 2 * sizeof(index_t); /* memory for MKZQ */
  mem_MKZ += (mat->ncols + 1) * sizeof(index_t); /* memory for MKZA */
  printf ("# memory usage: rows %luMB, cols %luMB, wt %luMB, MKZ %luMB, "
          "tot %lu MB (VmPeak %ld)\n",
          mem_rows >> 20, mem_cols >> 20, mem_wt >> 20, mem_MKZ >> 20,
          (mem_rows + mem_cols + mem_wt + mem_MKZ) >> 20,
          PeakMemusage () >> 10);
}

static inline void
print_report (filter_matrix_t *mat)
{
  long vmsize = Memusage();
  long vmrss = Memusage2();
  static double t0 = -1;
  char sbuf0[16];
  char sbuf1[16];
  if (t0 < 0) t0 = wct_seconds();
  printf ("%.1f N=%" PRIu64 " (%" PRId64 ") W=%" PRIu64 " W*N=%.2e "
          "W/N=%.2f #Q=%lu [%s / %s]\n",
          wct_seconds() - t0,
          mat->rem_nrows,
          ((int64_t) mat->rem_nrows) - ((int64_t) mat->rem_ncols),
          mat->weight, (double) compute_WN (mat), compute_WoverN (mat),
          (unsigned long) MkzQueueCardinality (mat),
          size_disp(vmrss<<10, sbuf0),
          size_disp(vmsize<<10, sbuf1));
  if (mat->verbose)
    print_memory_usage (mat);
  fflush (stdout);
}

struct merge_stats_s
{
  uint64_t nb;
  int32_t w_min;
  int32_t w_max;
  float w_av;
};

typedef struct merge_stats_s merge_stats_t;

merge_stats_t *
merge_stats_malloc (int maxlevel)
{
  merge_stats_t *data = NULL;
  data = (merge_stats_t *) malloc ((maxlevel + 1) * sizeof (merge_stats_t));
  ASSERT_ALWAYS (data != NULL);
  for (int i = 0; i <= maxlevel; i++)
  {
    data[i].nb = 0;
    data[i].w_min = INT32_MAX;
    data[i].w_max = INT32_MIN;
    data[i].w_av = 0.0;
  }

  return data;
}

/* This function assumes that 0 <= level <= maxlevel */
static inline int
merge_stats_is_first_merge (merge_stats_t *data, int level)
{
  return (data[level].nb == 0);
}

/* This function assumes that 0 <= level <= maxlevel */
static inline void
merge_stats_update (merge_stats_t *data, int level, int32_t weight)
{
  data[level].nb++;
  data[level].w_max = MAX (data[level].w_max, weight);
  data[level].w_min = MIN (data[level].w_min, weight);
  data[level].w_av += (float) weight;
}

void
merge_stats_print (merge_stats_t *data, int maxlevel)
{
  uint64_t nbtot = 0;
  for (int level = 1; level <= maxlevel; level++)
  {
    nbtot += data[level].nb;
    if (data[level].nb > 0)
      printf ("Number of %d-merges: %" PRIu64 " (weight min %" PRId32 ", "
              "max %" PRId32 ", avg %.2f)\n", level,
              data[level].nb, data[level].w_min, data[level].w_max,
              data[level].w_av / ((float) data[level].nb));
    else
      printf ("Number of %d-merges: %" PRIu64 "\n", level, data[level].nb);
  }
  printf ("Total number of merges: %" PRIu64 "\n", nbtot);
}

void
mergeOneByOne (report_t *rep, filter_matrix_t *mat, int maxlevel,
               double target_density)
{
  index_t j;
  index_signed_t mkz;
  uint64_t WN_prev, WN_cur, WN_min;
  double WoverN;
  unsigned int ncost = 0;
  int m;
  merge_stats_t *merge_data;

  /* first check the input matrix is valid */
  checkMat (mat);

  printf ("# Using %s to compute the merges\n", __func__);

  /* fill the heap of heavy rows */
  heap_fill (mat);

  // clean things
#ifdef MIMIC
  for( ; removeSingletons (rep, mat); )
      recomputeR(mat);
#else
  removeSingletons (rep, mat);
#endif

  merge_data = merge_stats_malloc (maxlevel);

  WN_min = WN_cur = compute_WN(mat);
  WoverN = compute_WoverN (mat);

  /* report lines are printed for each multiple of report_incr */
  double report_incr = REPORT_INCR;
  double report_next = ceil (WoverN / report_incr) * report_incr;

  int64_t excess;
  /* stop when excess = keep and (WoverN >= target or queue is empty) */
  while ((excess = mat->rem_nrows - mat->rem_ncols) > mat->keep ||
         (WoverN < target_density || MkzQueueCardinality (mat) > 0))
  {
    if (excess > mat->keep &&
        (WoverN >= target_density || MkzQueueCardinality (mat) == 0))
      { /* we hope removing the remaining excess will decrease the average
           row density and thus enable a few more merges */
        printf ("Removing final excess, nrows=%" PRIu64 "\n",
                mat->rem_nrows);
        deleteSuperfluousRows (rep, mat, excess - mat->keep, INT_MAX);
        printf ("Removing singletons, nrows=%" PRIu64 "\n", mat->rem_nrows);
        removeSingletons (rep, mat);
        recomputeR (mat);
        WoverN = compute_WoverN (mat);
        continue; /* this cannot loop forever since the excess decreases */
      }

    /* now excess = keep or (WoverN < target_density and #Q > 0) */

    if (excess == mat->keep &&
        (WoverN >= target_density || MkzQueueCardinality (mat) == 0))
      break; /* we cannot do any more merge */

    /* now WoverN < target_density and #Q > 0 */

    /* Do one merge */
    int ret = MkzPopQueue (&j, &mkz, mat);
    ASSERT_ALWAYS(ret != 0);

    m = mat->wt[j];
    ASSERT_ALWAYS(m >= 0);

    /* m=0 can happen if two ideals have weight 2, and are in the same 2
       relations: when we merge one of them, the other one will have weight
       0 after the merge */
    if (m == 0)
      goto next_merge;

    uint64_t weight0 = mat->weight;
    if (m == 1) /* singleton ideal */
        removeColDefinitely(rep, mat, j);
    else /* m >= 2 */
      {
        fprintf(rep->outfile, "#\n");
        mergeForColumn (rep, mat, m, j);
      }
    index_signed_t real_mkz = mat->weight - weight0;
    if (m > 1 && merge_stats_is_first_merge (merge_data, m))
      {
        fprintf (rep->outfile, "## First %d-merge\n", m);
        print_report (mat);
        printf ("First %d-merge, estimated cost %ld, real cost %ld\n",
                m, (long) mkz, (long) real_mkz);
        fflush (stdout);
      }

    /* Update values and report if necessary */
    merge_stats_update (merge_data, m, real_mkz);
    WoverN = compute_WoverN (mat);
    WN_prev = WN_cur;
    WN_cur = compute_WN(mat);
    if (WN_cur > WN_prev)
      ncost++;
    else
      ncost = 0;
    if (WN_cur < WN_min)
      WN_min = WN_cur;

    if (WoverN >= report_next)
    {
      int ni2rem;

      print_report (mat);
      report_next += report_incr;
      ni2rem = number_of_superfluous_rows (mat);
      deleteSuperfluousRows (rep, mat, ni2rem, m);
      removeSingletons (rep, mat);
      /* even if we had cwmax = mergelevelmax before, or at the very beginning,
         it is important to recompute the R matrix (and the Markowitz queue),
         since some ideals that had weight > mergelevelmax before could now
         have weight <= mergelevelmax */
      if (mat->cwmax == mat->mergelevelmax)
        recomputeR (mat);
    }

  next_merge:
    /* if we have not yet reached maxlevel:
       (a) if we have reached the target excess, we increase cwmax to maxlevel
       (b) if we have not yet reached keep, and the queue is empty, we increase
           cwmax by 1 only */
    if (mat->cwmax < mat->mergelevelmax &&
        (MkzQueueCardinality (mat) == 0 || excess == mat->keep))
      {
        if (excess == mat->keep)
          mat->cwmax = mat->mergelevelmax;
        else /* MkzQueueCardinality (mat) == 0 */
          mat->cwmax ++;
        recomputeR (mat);
      }
  }

  /* now excess = keep and (WoverN >= target or queue is empty) */

  if (WoverN >= target_density)
    printf ("Reached target density W/N=%.2f.\n", WoverN);
  else
    printf ("Heap is empty, stopping. Rerun with larger maxlevel if "
            "more merges are needed\n");

#if DEBUG >= 1
  checkWeights (mat);
#endif

  merge_stats_print (merge_data, maxlevel);
  free (merge_data);

  printf ("Total number of removed columns: %" PRId64 "\n",
          mat->ncols - mat->rem_ncols);

#ifdef TIMINGS
  for (int m = 2; m < MERGE_LEVEL_MAX; m++)
    {
      if (tfill[m] != 0 || tmst[m] != 0)
        printf ("m=%d: ncalls=%.0f tfill=%.2f tmst=%.2f\n",
                m, nfill[m], tfill[m], tmst[m]);
      tfill[0] += tfill[m];
      tmst[0] += tmst[m];
      nfill[0] += nfill[m];
    }
  printf ("Total: ncalls=%.0f tfill=%.2f tmst=%.2f\n", nfill[0], tfill[0],
          tmst[0]);
  printf ("Time for recomputeR: %.2f\n", trecomputeR);
#endif
}

//////////////////////////////////////////////////////////////////////
//
// Resume section: very much inspired by replay.c, of course...!
//
//////////////////////////////////////////////////////////////////////

// A line is "i i1 ... ik".
// If i >= 0 then
//     row[i] is to be added to rows i1...ik and destroyed at the end of
//     the process.
//     Works also is i is alone (hence: destroyed row).
// If i < 0 then
//     row[-i-1] is to be added to rows i1...ik and NOT destroyed.
//
static void
doAllAdds(report_t *rep, filter_matrix_t *mat, char *str)
{
  index_t j;
  index_signed_t ind[MERGE_LEVEL_MAX], i0;
  int ni, sg, k;

  ni = parse_hisfile_line (ind, str, &j);  
  
  if (ind[0] < 0)
    {
      sg = -1;
      i0 = -ind[0]-1;
    }
  else
    {
      sg = 1;
      i0 = ind[0];
    }

  for (k = 1; k < ni; k++)
      addRowsAndUpdate(mat, ind[k], i0, j);

  reportn(rep, ind, ni, j);
  
  if (sg > 0)
    removeRowAndUpdate(mat, i0, 1);
}

// resumename is a file of the type mergehis.
// TODO: Compiles, but not really tested with Markowitz...!
void
resume(report_t *rep, filter_matrix_t *mat, const char *resumename)
{
    FILE *resumefile = fopen(resumename, "r");
    char str[RELATION_MAX_BYTES];
    unsigned long addread = 0;
    index_t nactivej;
    index_t j;
    char * rp;
    int REPORT = 10000;

    printf ("Resuming computations from %s\n", resumename);
    // skip first line containing nrows ncols
    rp = fgets(str, RELATION_MAX_BYTES, resumefile);
    ASSERT_ALWAYS(rp);

    printf ("Reading row additions\n");
    while(fgets(str, RELATION_MAX_BYTES, resumefile)){
	addread++;
	if((addread % (10 * REPORT)) == 0)
	    printf ("%lu lines read at %2.2lf\n", addread, seconds());
	if(str[strlen(str)-1] != '\n'){
	    printf ("Gasp: not a complete a line!");
	    printf (" I stop reading and go to the next phase\n");
	    break;
	}
        doAllAdds(rep, mat, str);
    }
    fclose(resumefile);
    for(j = 0; j < mat->ncols; j++)
	if((mat->wt[j] == 0) && MkzIsAlive(mat->MKZA, j))
	    // be sure j was removed...
	    removeColumnAndUpdate(mat, j);
    nactivej = 0;
    for(j = 0; j < mat->ncols; j++)
	if(mat->wt[j] != 0)
	    nactivej++;
    mat->rem_ncols = nactivej;
    printf ("At the end of resume, we have");
    printf (" nrows=%" PRIu64 " ncols=%" PRIu64 " (%" PRIu64 ")\n",
            mat->rem_nrows, mat->rem_ncols, mat->rem_nrows - mat->rem_ncols);
}
