#include "cado.h" // IWYU pragma: keep
#include <limits.h>

#include "filter_config.h"
#include "sparse.h"
#include "mst.h"
#include "macros.h"

//////////////////////////////////////////////////////////////////////
// heap data structure to store the remaining edges to be considered
// in Prim's algorithm

// #define DEBUG 1

/* naive implementation of Prim's algorithm:
   we put in start[i] and end[i] the values s and t of each edge (s, t) */
int
minimalSpanningTree (int *start, int *end, int m,
		     int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX])
{
  int n, k, i, j, w = 0, imin, jmin, wmin;
  int S[MERGE_LEVEL_MAX], T[MERGE_LEVEL_MAX];

  ASSERT_ALWAYS(m - 1 < MERGE_LEVEL_MAX);
  ASSERT_ALWAYS(m >= 2);

  /* S is the set of vertices in the current tree, T is the set of remaining
     vertices */
  S[0] = 0; /* S = {0} */
  for (i = 1; i < m; i++)
    T[i-1] = i; /* T = {1, 2, ..., m-1} */
  n = 1;     /* number of vertices in S */
  k = m - 1; /* number of vertices in T */
  while (k)
    {
      int s, t;
      /* find the edge with minimal weight from S to T */
      wmin = INT_MAX;
      imin = jmin = -1;
      for (i = 0; i < n; i++)
        for (j = 0; j < k; j++)
          if (A[S[i]][T[j]] < wmin)
            {
              imin = i;
              jmin = j;
              wmin = A[S[i]][T[j]];
            }
      ASSERT(imin != -1 && jmin != -1);
      s = S[imin];
      t = T[jmin];
      w += wmin;
      S[n] = t;
      T[jmin] = T[k - 1];
      start[n-1] = s;
      end[n-1] = t;
      n++;
      k--;
    }
  return w;
}

/* given an ideal of weight m, fills the m x m matrix A so that
   A[i][j] is the weight of the sum of the i-th and j-th rows
   containing the ideal, for 0 <= i, j < m */
void
fillRowAddMatrix(int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX], filter_matrix_t *mat,
                 int m, index_t *ind, index_t ideal)
{
    int i, j;

    /* A[i][i] is not used, thus we don't initialize it. */
    for(i = 0; i < m; i++)
	for(j = i+1; j < m; j++){
	    A[i][j] = weightSum (mat, ind[i], ind[j], ideal);
	    A[j][i] = A[i][j];
	}
}

int
minCostUsingMST (filter_matrix_t *mat, int m, index_t *ind, index_t j)
{
    int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX], w;
    int sons[MERGE_LEVEL_MAX];
    int father[MERGE_LEVEL_MAX];

    fillRowAddMatrix (A, mat, m, ind, j);
    w = minimalSpanningTree (father, sons, m, A);
    /* w is the cost of all merges, we should subtract the cost of the
       initial relations */
    for (int i = 0; i < m; i++)
      w -= matLengthRow(mat, ind[i]);
    return w;
}

#ifdef DEBUG
void
printMST (int *father, int *sons, int m, index_t *ind)
{
  for (int i = 0; i < m; i++)
    printf ("father=%d(%lu) son=%d(%lu)\n", father[i],
            (unsigned long) ind[father[i]], sons[i],
            (unsigned long) ind[sons[i]]);
}
#endif
