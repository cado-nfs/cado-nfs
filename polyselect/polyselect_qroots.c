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
    fprintf (stderr, "Error, cannot reallocate memory in roots_realloc\n");
    exit (1);
  }
  R->nr = realloc (R->nr, newalloc * sizeof (unsigned int));
  if (R->nr == NULL)
  {
    fprintf (stderr, "Error, cannot reallocate memory in roots_realloc\n");
    exit (1);
  }
  R->roots = realloc (R->roots, newalloc * sizeof (uint64_t*));
  if (R->roots == NULL)
  {
    fprintf (stderr, "Error, cannot reallocate memory in roots_realloc\n");
    exit (1);
  }
}

/* reorder by decreasing number of roots (nr) */
void
polyselect_qroots_rearrange (polyselect_qroots_ptr R)
{
  if (R->size > 1) {
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
  }
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
    fprintf (stderr, "Error, cannot allocate memory in roots_add\n");
    exit (1);
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

