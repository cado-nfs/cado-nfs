#include "cado.h"
#include <stdio.h>
#include <string.h>

#include "portability.h"
#include "filter_config.h"
#include "utils.h"
#include "merge_replay_matrix.h"
#include "sparse.h"
#include "mst.h"

//////////////////////////////////////////////////////////////////////
// heap data structure to store the remaining edges to be considered
// in Prim's algorithm

// #define DEBUG 1

/* The queue is a heap of 32-bit entries, with the upper 16 bits being the
   weight of the edge, the next 8 bits being the end vertex, and the last
   8 bits the start vertex. In such a way by comparing two values we compare
   the weight. The heap starts at Q[1], Q[0] being the number of elements of
   the heap, thus the heap values are Q[1]...Q[Q[0]]. */

#define isQueueEmpty(Q) (Q[0] == 0)
#define W(x) ((x) >> 16)
#define S(x) ((x) & 255)
#define T(x) (((x) >> 8) & 255)

#define MAX_QUEUE_SIZE (MERGE_LEVEL_MAX*MERGE_LEVEL_MAX/2)

#if DEBUG >= 1
static void
printQueue (index_t *Q)
{
  for (unsigned int i = 1; i <= Q[0]; i++)
    fprintf (stderr, "w=%d s=%d t=%d\n", W(Q[i]), S(Q[i]), T(Q[i]));
}
#endif

static void
upQueue (index_t *Q, int k)
{
  index_t x = Q[k];

  while ((k > 1) && (x < Q[k/2]))
    {
      // we are at level > 0 and the father is > son
      // the father replaces the son
      Q[k] = Q[k/2];
      k /= 2;
    }
  // we found the place of x
  Q[k] = x;
}

static void
addEdge (index_t *Q, int s, int t, int Auv)
{
  Q[0] ++;
  ASSERT(Q[0] < MAX_QUEUE_SIZE);
  Q[Q[0]] = (Auv << 16) | (t << 8) | s;
  upQueue (Q, Q[0]);
}

// Move Q[k] down.
static void
downQueue (index_t *Q, unsigned int k)
{
  index_t x = Q[k];
  unsigned int j;

  /* the left son of k is 2k, the right son is 2k+1 */
  while (2*k <= Q[0])
    {
      // k has at least a left son
      j = 2*k;
      if (j < Q[0])
        // k has a right son
        if (Q[j] > Q[j+1])
          j++;
      // at this point, Q[j] is the smallest son
      if (x <= Q[j])
        break; /* x is smaller than smallest son */
      else
        {
          // the father takes the place of the son
          Q[k] = Q[j];
          k = j;
        }
    }
  // we found the place of x
  Q[k] = x;
}

static void
popQueue(int *s, int *t, index_t *Q)
{
  *s = S(Q[1]);
  *t = T(Q[1]);
  Q[1] = Q[Q[0]];
  Q[0]--;
  downQueue (Q, 1);
}

// Add all edges (u,v) with v in V
static void
addAllEdges (index_t *Q, int u, int *V, int nV,
             int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX])
{
  for (int i = 0; i < nV; i++)
    addEdge (Q, u, V[i], A[u][V[i]]);
}

/* naive implementation of Prim's algorithm:
   we put in start[i] and end[i] the values s and t of each edge (s, t) */
static int
minimalSpanningTreePrimNaive (int *start, int *end,
                              int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX], int m)
{
  int n, k, i, j, w = 0, imin, jmin, wmin;
  int S[MERGE_LEVEL_MAX], T[MERGE_LEVEL_MAX];

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

/* Return the weight of the minimal spanning tree.
   For each (s,t) which is part of the tree, we have
   start[i] = s and end[i] = t for 0 <= i < m-1. */
static int
minimalSpanningTreeWithPrim(int *start, int *end,
			    int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX], int m)
{
  /* the queue contains at most m-1 edges when nU=1 (those connected to the
     root node), then at most 2*(m-2) when nU=2 (we remove one node and add
     m-2), ... More generally when nU=k it is at most k*(m-k) + (k-1)*(k-2)/2.
     The maximum is for k=m-1 or m-2 when it equals m^2/2-3/2*m+2 < m^2/2. */
    index_t Q[MAX_QUEUE_SIZE];
    int u, s, t, i, nU, nV, w;
    int V[MERGE_LEVEL_MAX]; /* edges remaining to be dealt with */
    int index[MERGE_LEVEL_MAX];

    // nodes of T
    for(i = 0; i < m; i++){
        V[i] = i;
        index[i] = i; /* V[index[i]] = i if i is in V, -1 otherwise */
    }
    u = 0;
    index[u] = -1;
    index[V[m-1]] = u;
    V[u] = V[m-1];
    nU = 1;   /* number of edges already in the MST */
    nV = m-1; /* number of edges remaining to be dealt with */
    Q[0] = 0; /* heap is empty initially */
    addAllEdges (Q, u, V, nV, A);
#if DEBUG >= 1
    printQueue(Q);
#endif
    w = 0;
    while (!isQueueEmpty(Q)) /* while queue is non empty */
      {
        popQueue (&s, &t, Q); // pop queue
#if DEBUG >= 1
	fprintf(stderr, "Popping a = (%d, %d) of weight %d\n", s, t, A[s][t]);
#endif
	if (index[t] != -1){
	    // t does not close a cycle
	    // T[u] <- T[u] union (s, t)
#if DEBUG >= 1
	    fprintf(stderr, "new edge: (%d, %d) of weight %d\n",s,t,A[s][t]);
#endif
	    w += A[s][t];
            start[nU - 1] = s;
            end[nU - 1] = t;
            ASSERT(V[index[t]] == t);
            index[V[nV-1]] = index[t];
            V[index[t]] = V[nV-1];
            index[t] = -1;
	    nV--;
            nU++;
	    if (nV == 0)
		break;
            addAllEdges (Q, t, V, nV, A);
#if DEBUG >= 1
	    printQueue(Q);
#endif
	}
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
minimalSpanningTree(int *start, int *end,
		    int m, int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX])
{
  int ret;

  if (m <= 25)
    ret = minimalSpanningTreePrimNaive (start, end, A, m);
  else
    ret = minimalSpanningTreeWithPrim (start, end, A, m);
  return ret;
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
