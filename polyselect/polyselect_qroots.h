#ifndef POLYSELECT_QROOTS_H_
#define POLYSELECT_QROOTS_H_

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* structure to store q roots */
struct polyselect_qroots_s {
  unsigned int alloc;   /* allocated size */
  unsigned int size;    /* used size */
  unsigned int *q;
  unsigned int *nr;     /* number of roots of x^d = N (mod q) */
  uint64_t **roots;     /* roots of (m0+x)^d = N (mod q^2) */
};
typedef struct polyselect_qroots_s polyselect_qroots_t[1];
typedef struct polyselect_qroots_s * polyselect_qroots_ptr;
typedef const struct polyselect_qroots_s * polyselect_qroots_srcptr;

void polyselect_qroots_init (polyselect_qroots_ptr);
void polyselect_qroots_realloc (polyselect_qroots_ptr, unsigned long);
void polyselect_qroots_add (polyselect_qroots_ptr, unsigned int, unsigned int, uint64_t*);
void polyselect_qroots_print (polyselect_qroots_srcptr);
void polyselect_qroots_rearrange (polyselect_qroots_ptr R);
void polyselect_qroots_clear (polyselect_qroots_ptr);


#ifdef __cplusplus
}
#endif

#endif	/* POLYSELECT_QROOTS_H_ */
