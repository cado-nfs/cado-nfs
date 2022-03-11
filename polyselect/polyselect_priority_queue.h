#ifndef POLYSELECT_PRIORITY_QUEUE_H_
#define POLYSELECT_PRIORITY_QUEUE_H_

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

struct polyselect_priority_queue_s {
  unsigned long size;
  double * data;
};
typedef struct polyselect_priority_queue_s polyselect_priority_queue_t[1];
typedef struct polyselect_priority_queue_s * polyselect_priority_queue_ptr;
typedef const struct polyselect_priority_queue_s * polyselect_priority_queue_srcptr;

extern void polyselect_priority_queue_init(polyselect_priority_queue_ptr, size_t);
extern void polyselect_priority_queue_reset(polyselect_priority_queue_ptr Q);
extern void polyselect_priority_queue_resize(polyselect_priority_queue_ptr Q, size_t n);
extern void polyselect_priority_queue_clear(polyselect_priority_queue_ptr);
extern int polyselect_priority_queue_push(polyselect_priority_queue_ptr, double x);
extern double polyselect_priority_queue_top(polyselect_priority_queue_srcptr);
extern void polyselect_priority_queue_merge(polyselect_priority_queue_ptr Q, polyselect_priority_queue_srcptr Q0);

extern void polyselect_priority_queue_snprintf(polyselect_priority_queue_srcptr Q, char * s, size_t n, const char * format, const char * delim);
#ifdef __cplusplus
}
#endif

#endif	/* POLYSELECT_PRIORITY_QUEUE_H_ */
