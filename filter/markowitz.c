#include "cado.h"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "portability.h"
#include "filter_config.h"
#include "utils_with_io.h"
#include "merge_replay_matrix.h"
#include "sparse.h"
#include "mst.h"
#include "report.h"
#include "markowitz.h"

#define MKZ_DEBUG 0
// #define MKZ_TIMINGS 1

#if MKZ_TIMINGS
double tmkzup, tmkzdown, tmkzupdown, tmkzcount;
#endif

// Again, a priority queue as a heap...!
// Q[0] contains the number of items in Q[], so that useful part of Q
// is Q[1..Q[0]]

// Q[2*i] contains j-jmin = dj
// Q[2*i+1] contains the Markowitz count for j

// A[j] gives u s.t. Q[2*u] = j

// get r-th component of Q[i]
#define MkzGet(Q, i, r) (Q[((i)<<1)+(r)])
#define MkzSet(Q, i, r, val) (Q[((i)<<1)+(r)] = (val))

inline int
MkzIsAlive(index_t *A, index_t dj)
{
  return A[dj] != MKZ_INF;
}

// (Q, A)[k1] <- (Q, A)[k2]
static void
MkzAssign(index_t *Q, index_t *A, index_t k1, index_t k2)
{
    index_t dj = MkzGet(Q, k2, 0);

    MkzSet(Q, k1, 0, dj);
    MkzSet(Q, k1, 1, MkzGet(Q, k2, 1)); // could be simplified...!
    A[dj] = k1;
}

#if MKZ_DEBUG
static void MAYBE_UNUSED
MkzPrintQueue(index_t *Q)
{
    int level = 0;
    index_t i, imax = 1;

    fprintf(stderr, "L0:");
    for(i = 1; i <= Q[0]; i++){
	fprintf(stderr, " [%d, %d]", MkzGet(Q, i, 1), MkzGet(Q, i, 0));
	if(i == imax){
	    imax = (imax<<1)+1;
	    fprintf(stderr, "\nL%d:", ++level);
	}
    }
    fprintf(stderr, "\n");
}
#endif

static void
MkzUpQueue(index_t *Q, index_t *A, index_t k)
{
    index_t dj = MkzGet(Q, k, 0), count = MkzGet(Q, k, 1);
#if MKZ_TIMINGS
    double tt = seconds();
#endif

    while((k > 1) && (MkzGet(Q, k/2, 1) >= count)){
	// we are at level > 0 and the father is >= son
	// the father replaces the son
      MkzAssign(Q, A, k, k/2);
	k /= 2;
    }
    // we found the place of (dj, count)
    MkzSet(Q, k, 0, dj);
    MkzSet(Q, k, 1, count);
    A[dj] = k;
#if MKZ_TIMINGS
    tmkzup += (seconds()-tt);
#endif
}

// Move Q[k] down, by keeping the structure of Q as a heap, i.e.,
// each node has a smaller cost than its two left and right nodes
static void
MkzDownQueue(index_t *Q, index_t *A, index_t k)
{
    index_t dj = MkzGet(Q, k, 0), count = MkzGet(Q, k, 1), j;
#if MKZ_TIMINGS
    double tt = seconds();
#endif

    while ((j = 2*k) <= Q[0]){ /* node k has at least one left son 2k */
	if(j < Q[0])
	    // node k has also a right son
	    if(MkzGet(Q, j, 1) > MkzGet(Q, j+1, 1))
		j++;
	// at this point, Q[j] is the son with the smallest "count"
	if(count <= MkzGet(Q, j, 1)) /* Q[k] has smaller cost than both sons */
	    break;
	else{
	    // the father takes the place of the "smaller" son
            MkzAssign(Q, A, k, j);
	    k = j;
	}
    }
    // we found the place of (dj, count)
    MkzSet(Q, k, 0, dj);
    MkzSet(Q, k, 1, count);
    A[dj] = k;
#if MKZ_TIMINGS
    tmkzdown += (seconds()-tt);
#endif
}

// (Q, A)[k] has just arrived, but we have to move it in the heap, so that
// it finds its place.
static void
MkzMoveUpOrDown(index_t *Q, index_t *A, index_t k)
{
#if MKZ_TIMINGS
    double tt = seconds();
#endif

    // move new node up or down
    if(k == 1)
	// rare event!
	MkzDownQueue(Q, A, k);
    else{
	// k has a father
	if(MkzGet(Q, k/2, 1) > MkzGet(Q, k, 1))
	    // we have to move up
	    MkzUpQueue(Q, A, k);
	else
	    MkzDownQueue(Q, A, k);
    }
#if MKZ_TIMINGS
    tmkzupdown += (seconds()-tt);
#endif
}

#if MKZ_DEBUG >= 1
static int MAYBE_UNUSED
MkzIsHeap(index_t *Q)
{
    index_t k;

    for(k = 1; k <= Q[0]/2; k++){
	// k has a left son
	if(MkzGet(Q, k, 1) > MkzGet(Q, 2*k, 1)){
	    fprintf(stderr, "Pb: father=%d > lson=%d\n",
		    MkzGet(Q, k, 1), MkzGet(Q, 2*k, 1));
	    return 0;
	}
	if(k < Q[0]/2){
	    // k has a right son
	    if(MkzGet(Q, k, 1) > MkzGet(Q, 2*k+1, 1)){
		fprintf(stderr, "Pb: father=%d > rson=%d\n",
			MkzGet(Q, k, 1), MkzGet(Q, 2*k+1, 1));
		return 0;
	    }
	}
    }
    return 1;
}
#endif

/* here we count a cost k for an ideal of weight k */
static int
Cavallar (filter_matrix_t *mat, index_t j)
{
  return mat->wt[j];
}

/* This functions returns the difference in matrix element when we add the
   lightest row with ideal j to all other rows. If ideal j has weight w,
   and the lightest row has weight w0:
   * we remove w elements corresponding to ideal j
   * we remove w0-1 other elements for the lightest row
   * we add w0-1 elements to the w-1 other rows
   Thus the result is (w0-1)*(w-2) - w = (w0-2)*(w-2) - 2.
   Note: we could also take into account the "cancelled ideals", i.e., the
   ideals apart the pivot j that got cancelled "by chance". However some
   experiments show that this does not improve (if any) the result of merge.
   A possible explanation is that the contribution of cancelled ideals follows
   a normal distribution, and that on the long run we mainly see the average.
*/
static int
pureMkz(filter_matrix_t *mat, index_t j)
{
    int w = mat->wt[j];

    if (w <= 1)
      return -4; /* ensures that empty columns and singletons are removed earlier */
    else if (w == 2)
      return -2;
    else
      {
        index_t i, w0;

        /* approximate traditional Markowitz count: we assume we add the
	   lightest row to all other rows */
        i = mat->R[j][1];
        w0 = matLengthRow(mat, i);
        for (unsigned int k = 2; k <= mat->R[j][0]; k++)
          {
            i = mat->R[j][k];
            if (matLengthRow(mat, i) < w0)
              w0 = matLengthRow(mat, i);
          }
	return (w0 - 2) * (w - 2) - 2;
      }
}

/* this function takes into account cancelled ideals "by chance" for
   w <= mat->wmstmax, and is identical to pureMkz() for larger weights.
   Thus for mat->wmstmax = 1, it should be identical to pureMkz(). */
static int
lightColAndMkz (filter_matrix_t *mat, index_t j)
{
    int wj = mat->wt[j];

    if (wj <= 1)
      return -4; /* like pureMkz */
    else if (wj <= mat->wmstmax)
      {
        index_t *ind = (index_t*) mat->R[j] + 1;
        if (wj == 2)
          return weightSum (mat, ind[0], ind[1], j)
            - matLengthRow (mat, ind[0]) - matLengthRow (mat, ind[1]);
        else
          return minCostUsingMST (mat, wj, ind, j);
      }
    // real traditional Markowitz count
    return pureMkz (mat, j);
}

/* return the cost of merging column j (the smaller, the better) */
static int32_t
MkzCount(filter_matrix_t *mat, index_t j)
{
  int32_t cost;
  switch(mat->mkztype){
  case MKZTYPE_LIGHT:
    cost = lightColAndMkz(mat, j);
    break;
  case MKZTYPE_PURE:
    cost = pureMkz(mat, j);
    break;
  case MKZTYPE_CAVALLAR:
  default:
      /* for the double-matrix trick, we count k for an ideal of weight k */
    cost = Cavallar(mat, j);
    /* fall through */
  }
  /* ensure cost >= 0, since the MKZ data structure is now unsigned */
  if (cost < 0)
    cost = 0;
  return cost;
}

void
MkzClear (filter_matrix_t *mat, int verbose)
{
  if (verbose)
    printf ("Max Markowitz count: %lu\n",
	    (unsigned long) MkzGet(mat->MKZQ, mat->MKZQ[0], 1));
#if MKZ_TIMINGS
  printf ("MKZT: up=%d down=%d updown=%d count=%d\n",
          (int)tmkzup, (int)tmkzdown, (int)tmkzupdown, (int)tmkzcount);
#endif
    free (mat->MKZQ);
    free (mat->MKZA);
}

/* increment the weight of column j (in absolute value) */
int
MkzIncrCol(filter_matrix_t *mat, index_t j)
{
    int ind;

#if MKZ_DEBUG >= 1
    fprintf(stderr, "Incr: wt(%d) was %d\n", j, mat->wt[j]);
#endif
    ind = mat->wt[j] = incrS(mat->wt[j]);
    return ind;
}

/* update the Markowitz cost of column j */
void
MkzUpdate (filter_matrix_t *mat, index_t j)
{
    index_t adr = mat->MKZA[j];
    index_t mkz;

    ASSERT(adr != MKZ_INF);
    /* compute the new Markowitz cost */
    mkz = MkzCount (mat, j);
    /* update new cost */
    MkzSet (mat->MKZQ, adr, 1, mkz);
    /* move it up or down in the heap */
    MkzMoveUpOrDown (mat->MKZQ, mat->MKZA, adr);
}

#if 0 /* parallel version */
/* update in parallel j[0], j[1], ..., j[n-1] */
void
MkzUpdateN (filter_matrix_t *mat, index_t *j, int n)
{
  index_t *mkz;
  int i;

  mkz = malloc (n * sizeof (index_t));
  ASSERT_ALWAYS(mkz != NULL);
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
  for (i = 0; i < n; i++)
    mkz[i] = MkzCount (mat, j[i]);
  /* the following is sequential */
  for (i = 0; i < n; i++)
    {
      uint32_t adr = mat->MKZA[j[i]];
      MkzSet (mat->MKZQ, adr, 1, mkz[i]);
      MkzMoveUpOrDown (mat->MKZQ, mat->MKZA, adr);
    }
  free (mkz);
}
#else
/* update in parallel j[0], j[1], ..., j[n-1] */
void
MkzUpdateN (filter_matrix_t *mat, index_t *j, int n)
{
  for (int i = 0; i < n; i++)
    MkzUpdate (mat, j[i]);
}
#endif

/*
   Updates:
   - mat->wt[j] (weight of column j)

   We arrive here when mat->wt[j] > 0.

*/
void
MkzDecreaseColWeight(filter_matrix_t *mat, index_t j)
{
    index_t dj = j;

#if MKZ_DEBUG >= 1
    fprintf(stderr, "Decreasing col %d; was %d\n", j, mat->wt[dj]);
#endif
    mat->wt[dj] = decrS(mat->wt[dj]);
}
