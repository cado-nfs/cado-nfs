#include "cado.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "polyselect_priority_queue.h"

void polyselect_priority_queue_init(polyselect_priority_queue_ptr Q, size_t n)
{
    Q->size = 0;
    Q->data = NULL;
    polyselect_priority_queue_resize(Q, n);
}

void polyselect_priority_queue_clear(polyselect_priority_queue_ptr Q)
{
    free(Q->data);
}

void polyselect_priority_queue_reset(polyselect_priority_queue_ptr Q)
{
    for(size_t i = 0 ; i < Q->size ; i++) {
        Q->data[i] = NAN;
    }
}
void polyselect_priority_queue_resize(polyselect_priority_queue_ptr Q, size_t n)
{
    Q->data = (double *) realloc(Q->data, n * sizeof(double));
    for(size_t i = Q->size ; i < n ; i++) {
        Q->data[i] = NAN;
    }
    Q->size = n;
}

/* Insert a value into a sorted array of length len.
   Returns 1 if element was inserted, 0 if it was too big */
int polyselect_priority_queue_push(polyselect_priority_queue_ptr Q, double x)
{
    /* XXX write the O(log(n)) push. Presently it's silly O(n)
     */
  int result = 0;
  if (Q->size == 0)
    return 0;
  double * w = Q->data + Q->size - 1;
  double * w0 = Q->data;
  if (isnan(*w) || x < *w)
    {
      for (; w > w0 && (isnan(w[-1]) || x < w[-1]); w--)
          *w = w[-1];
      *w = x;
      result = 1;
    }
  return result;
}

/* Insert a value into a sorted array of length len.
   Returns 1 if element was inserted, 0 if it was too big */
void polyselect_priority_queue_merge(polyselect_priority_queue_ptr Q, polyselect_priority_queue_srcptr Q0)
{
    /* XXX this is ugly, complexity wise... */
    for(size_t i = 0 ; i < Q0->size && !isnan(Q0->data[i]) ; i++)
        polyselect_priority_queue_push(Q, Q0->data[i]);
}

double polyselect_priority_queue_top(polyselect_priority_queue_srcptr Q)
{
    return Q->data[0];
}


void polyselect_priority_queue_snprintf(polyselect_priority_queue_srcptr Q, char * s, size_t n, const char * format, const char * delim)
{
    size_t c = 0;
    if (n) *s = '\0';
    for (size_t i = 0; i < Q->size && !isnan(Q->data[i]); i++) {
        if (i) c += snprintf(s + c, n - c, "%s", delim);
        snprintf(s + c, n - c, format, Q->data[i]);
    }
}
