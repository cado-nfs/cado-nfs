#include "cado.h"
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>

#include "polyselect_proots.h"

/* init polyselect_proots_t */
void
polyselect_proots_init (polyselect_proots_ptr R,
              unsigned long size )
{
  R->size = size;

  /* length of nr&roots is known now. lengths of roots[i] are TBD. */
  /* +1 for R->nr for end guard in collision_on_each_sq */
  R->nr = (uint8_t *) malloc ((size + 1) * sizeof (*(R->nr)));
  R->roots = (uint64_t **) malloc (size * sizeof (*(R->roots)));

  if (R->nr == NULL || R->roots == NULL) {
    fprintf (stderr, "Error, cannot allocate memory in polyselect_proots_init().\n");
    exit (1);
  }
}


/* add a root to polyselect_proots_t */
void
polyselect_proots_add ( polyselect_proots_ptr R,
             unsigned long nr,
             uint64_t *roots,
             unsigned long index )
{
  unsigned int i;
  R->nr[index] = nr;

  if (nr != 0) {
    R->roots[index] = (uint64_t *) malloc (nr * sizeof (uint64_t));
    if (R->roots[index] == NULL) {
      fprintf (stderr, "Error, cannot allocate memory in polyselect_proots_add\n");
      exit (1);
    }

    for (i = 0; i < nr; i++)
      R->roots[index][i] = roots[i];
  }
  else
    R->roots[index] = NULL;
}


/* print roots */
void
polyselect_proots_print (polyselect_proots_srcptr R)
{
  unsigned int i, j;
  for (i = 0; i < R->size; i++) {
    if (R->nr[i] == 0) {
      fprintf (stderr, "NULL\n");
    }
    else {
      for (j = 0; j < R->nr[i]; j ++)
        fprintf (stderr, "%" PRIu64 " ", R->roots[i][j]);
      fprintf (stderr, "\n");
    }
  }
}


/* clear roots */
void
polyselect_proots_clear (polyselect_proots_ptr R)
{
  unsigned int i;

  free (R->nr);
  for (i = 0; i < R->size; i++)
    free (R->roots[i]);
  free (R->roots);
}


