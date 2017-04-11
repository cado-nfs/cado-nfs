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

static void
MkzInsert(index_t *Q, index_t *A, index_t dj, index_t count)
{
    Q[0]++;
    MkzSet(Q, Q[0], 0, dj);
    MkzSet(Q, Q[0], 1, count);
    A[dj] = Q[0];
    MkzUpQueue(Q, A, Q[0]);
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

// Remove (Q, A)[k].
static void
MkzDelete(index_t *Q, index_t *A, index_t k)
{
#if MKZ_DEBUG >= 1
    fprintf(stderr, "MKZ: deleting (Q, A)[%d]=[%d, %d]\n", k,
	    MkzGet(Q, k, 0), MkzGet(Q, k, 1));
#endif
    // we put Q[Q[0]] in Q[k]
    MkzAssign (Q, A, k, Q[0]);
    Q[0]--;
    MkzMoveUpOrDown(Q, A, k);
}

#if 0
// Remove (Q, A)[k].
void
MkzRemove(index_t *dj, index_t *mkz, index_t *Q, uint32_t *A, index_t k)
{
    *dj = MkzGet(Q, k, 0);
    *mkz = MkzGet(Q, k, 1);
    A[*dj] = MKZ_INF;
    MkzAssign(Q, A, k, Q[0]);
    Q[0]--;
    MkzMoveUpOrDown(Q, A, k);
}
#endif

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

#if MKZ_DEBUG >= 1
static void
MkzCheck(filter_matrix_t *mat)
{
  uint64_t dj;
  int maxlevel = mat->mergelevelmax;

  for(dj = 0; dj < mat->ncols; dj++)
  {
    if (0 < mat->wt[dj] && mat->wt[dj] <= maxlevel)
	  {
      if(MkzGet(mat->MKZQ, mat->MKZA[dj], 0) != (index_t) dj)
      {
        fprintf(stderr, "GASP: %" PRId32 " <> %" PRIu64 " in MkzCheck\n",
                        MkzGet(mat->MKZQ, mat->MKZA[dj], 0), dj);
        exit (1);
      }
    }
  }
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
    switch(mat->mkztype){
    case MKZTYPE_LIGHT:
	return lightColAndMkz(mat, j);
    case MKZTYPE_PURE:
	return pureMkz(mat, j);
    case MKZTYPE_CAVALLAR:
    default:
      /* for the double-matrix trick, we count k for an ideal of weight k */
	return Cavallar(mat, j);
    }
}

/* pop the top element from the heap, return 0 iff heap is empty */
int
MkzPopQueue(index_t *dj, index_signed_t *mkz, filter_matrix_t *mat)
{
  index_t *Q = mat->MKZQ;
  index_t *A = mat->MKZA;

  if (Q[0] == 0)
    return 0;

  /* Q[0] contains the number of items in Q[], thus the first element is
     stored in Q[2..3] */
  *dj = MkzGet(Q, 1, 0);
  *mkz = MkzGet(Q, 1, 1);
  while (mat->wt[*dj] > mat->mergelevelmax)
    {
      /* remove heavy column */
      MkzDelete (Q, A, 1);
      A[*dj] = MKZ_INF;

      if (MkzQueueCardinality (mat) == 0)
        return 0;

      *dj = MkzGet(Q, 1, 0);
      *mkz = MkzGet(Q, 1, 1);
    }
  A[*dj] = MKZ_INF; /* already done in MkzRemoveJ, but if we don't do it,
                       we get A[j1]=A[j2] for some j1 <> j2 */
  if (Q[0] > 1)
    {
      MkzAssign(Q, A, 1, Q[0]); /* move entry of index Q[0] to index 1 */
      MkzDownQueue(Q, A, 1);    /* reorder heap structure */
    }
  Q[0]--;                   /* decrease number of entries in Q,A */
  return 1;
}

void
MkzInit (filter_matrix_t *mat, int verbose)
{
    int32_t *mkz;
    index_t j;
    uint64_t sz = 0;
    int maxlevel = mat->mergelevelmax;

#if MKZ_TIMINGS
    tmkzup = tmkzdown = tmkzupdown = tmkzcount = 0.0;
#endif
    if (verbose)
      {
        fprintf (stderr, "Entering initMarkowitz (type=%d", mat->mkztype);
        if (mat->mkztype == MKZTYPE_LIGHT)
          fprintf (stderr, ", wmstmax=%d", mat->wmstmax);
        fprintf (stderr, ")\n");
      }

    // compute number of eligible columns in the heap
    for (j = 0; j < mat->ncols; j++)
      if (0 < mat->wt[j] && mat->wt[j] <= maxlevel)
        sz++;

    // Allocating heap MKZQ
    size_t tmp_alloc = (sz + 1) * 2 * sizeof(index_t);
    if (verbose)
      fprintf (stderr, "Allocating heap for %" PRIu64 " columns (%zuMB)\n", sz,
               tmp_alloc >> 20);
    mat->MKZQ = (index_t *) malloc (tmp_alloc);
    ASSERT_ALWAYS(mat->MKZQ != NULL);
    mat->MKZQ[0] = 0;
    mat->MKZQ[1] = (index_t) sz; // why not?
    
    // every j needs a pointer (MKZA)
    tmp_alloc = (mat->ncols + 1) * sizeof(index_t);
    if (verbose)
      fprintf (stderr, "Allocating pointers to heap: %zuMB\n",
               tmp_alloc >> 20);
    mat->MKZA = (index_t *) malloc (tmp_alloc);
    ASSERT_ALWAYS(mat->MKZA != NULL);

    mkz = malloc (mat->ncols * sizeof (int32_t));
    ASSERT_ALWAYS(mkz != NULL);

    /* since the computation of the Markowitz cost is read-only,
       we can perform it in parallel */
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for(j = 0; j < mat->ncols; j++)
      if (0 < mat->wt[j] && mat->wt[j] <= maxlevel)
        mkz[j] = MkzCount(mat, j);
      else
        mkz[j] = INT32_MIN;

    /* insertion in the heap cannot be done in parallel */
    for(j = 0; j < mat->ncols; j++)
      if (mkz[j] != INT32_MIN)
        MkzInsert(mat->MKZQ, mat->MKZA, j, mkz[j]);
      else
        mat->MKZA[j] = MKZ_INF;

    free (mkz);
    
#if MKZ_DEBUG >= 1
    MkzCheck(mat);
    fprintf(stderr, "Initial queue is\n");
    MkzPrintQueue(mat->MKZQ);
#endif
}

void
MkzClear (filter_matrix_t *mat, int verbose)
{
  if (verbose)
    fprintf(stderr, "Max Markowitz count: %lu\n",
	    (unsigned long) MkzGet(mat->MKZQ, mat->MKZQ[0], 1));
#if MKZ_TIMINGS
    fprintf(stderr, "MKZT: up=%d down=%d updown=%d count=%d\n",
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

/* remove column j from and update matrix */
void
MkzRemoveJ(filter_matrix_t *mat, index_t j)
{
    mat->wt[j] = 0;

    /* This can happen when maxlevel < weight[j] <= cwmax initially, thus
       A[j] was initialized to MKZ_INF, but because of a merge it becomes
       larger than cwmax. */
    if (mat->MKZA[j] == MKZ_INF)
      return;

#if MKZ_DEBUG >= 1
    fprintf(stderr, "Removing col %d of weight %d\n", j, mat->wt[j]);
#endif
    // remove j from the QA structure
    MkzDelete (mat->MKZQ, mat->MKZA, mat->MKZA[j]);
    mat->MKZA[j] = MKZ_INF;
#if MKZ_DEBUG >= 1
    MkzIsHeap(mat->MKZQ);
#endif
}
