#include "cado.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include "cado_poly.h"
#include "polyselect_qroots.h"

void
polyselect_qroots_init (polyselect_qroots_ptr R)
{
  R->alloc = 0;
  R->size = 0;
  R->q = NULL;
  R->nr = NULL;
  R->roots = NULL;
}

void
polyselect_qroots_realloc (polyselect_qroots_ptr R, unsigned long newalloc)
{
  ASSERT (newalloc >= R->size);
  R->alloc = newalloc;
  R->q = realloc (R->q, newalloc * sizeof (unsigned int));
  if (R->q == NULL)
    {
      fprintf(stderr, "Error, cannot reallocate memory in %s\n", __func__);
      exit(1);
    }
  R->nr = realloc (R->nr, newalloc * sizeof (unsigned int));
  if (R->nr == NULL)
    {
      fprintf(stderr, "Error, cannot reallocate memory in %s\n", __func__);
      exit(1);
    }
  R->roots = realloc (R->roots, newalloc * sizeof (uint64_t*));
  if (R->roots == NULL)
    {
      fprintf(stderr, "Error, cannot reallocate memory in %s\n", __func__);
      exit(1);
    }
}

/* reorder by decreasing number of roots (nr) 
 */

struct qr_flat {
    unsigned int q;
    unsigned int nr;
    uint64_t * roots;
};

int qr_flat_compare_nroots(const struct qr_flat * a, const struct qr_flat * b)
{
    return (a->nr < b->nr) - (b->nr < a->nr);
}

void
polyselect_qroots_rearrange (polyselect_qroots_ptr R)
{
    if (R->size <= 1) return;
#if 1
    /* honestly, it's embarrassing to have an insertion sort. */
    /* We'll keep it on just for debugging, since change in the ordering
     * leads to a a change in the polynomials that are examined!
     */
    unsigned int i, j, k, max, tmpq, tmpnr;
    uint64_t *tmpr = malloc (MAX_DEGREE * sizeof (uint64_t));

    for (i = 0; i < R->size; i ++) {
        max = i;
        for (j = i+1; j < R->size; j++) {
            if (R->nr[j] > R->nr[max]) {
                max = j;
            }
        }

        tmpq = R->q[i];
        tmpnr = R->nr[i];
        for (k = 0; k < MAX_DEGREE; k ++)
            tmpr[k] = R->roots[i][k];

        R->q[i] = R->q[max];
        R->nr[i] = R->nr[max];
        for (k = 0; k < MAX_DEGREE; k ++)
            R->roots[i][k] = R->roots[max][k];

        R->q[max] = tmpq;
        R->nr[max] = tmpnr;
        for (k = 0; k < MAX_DEGREE; k ++)
            R->roots[max][k] = tmpr[k];
    }
    free (tmpr);
#else
       struct qr_flat * tmp = malloc(R->size * sizeof(struct qr_flat));
       for(unsigned int i = 0 ; i < R->size ; i++) {
           tmp[i].q = R->q[i];
           tmp[i].nr = R->nr[i];
           tmp[i].roots = R->roots[i];
       }
       typedef int (*sortfunc_t) (const void *, const void *);
       qsort(tmp, R->size, sizeof(struct qr_flat), (sortfunc_t) qr_flat_compare_nroots);
       for(unsigned int i = 0 ; i < R->size ; i++) {
           R->q[i] = tmp[i].q;
           R->nr[i] = tmp[i].nr;
           R->roots[i] = tmp[i].roots;
       }
       free(tmp);
#endif
}

void
polyselect_qroots_add (polyselect_qroots_ptr R, unsigned int q, unsigned int nr, uint64_t *roots)
{
  unsigned int i;

  if (nr == 0)
    return;
  if (R->size == R->alloc)
    polyselect_qroots_realloc (R, R->alloc + R->alloc / 2 + 1);
  R->q[R->size] = q;
  R->nr[R->size] = nr;
  R->roots[R->size] = malloc (MAX_DEGREE * sizeof (uint64_t));
  if (R->roots[R->size] == NULL)
    {
      fprintf(stderr, "Error, cannot allocate memory in %s\n", __func__);
      exit(1);
    }
  for (i = 0; i < nr; i++)
    R->roots[R->size][i] = roots[i];
  R->size ++;
}

void
polyselect_qroots_print (polyselect_qroots_srcptr R)
{
  unsigned int i, j;
  for (i = 0; i < R->size; i++) {
    fprintf (stderr, "q: %u, r: ", R->q[i]);
    for (j = 0; j < R->nr[i]; j ++)
      fprintf (stderr, "%" PRIu64 " ", R->roots[i][j]);
    fprintf (stderr, "\n");
  }
}

void
polyselect_qroots_clear (polyselect_qroots_ptr R)
{
  unsigned int i;

  free (R->q);
  free (R->nr);
  for (i = 0; i < R->size; i++)
    free (R->roots[i]);
  free (R->roots);
}

