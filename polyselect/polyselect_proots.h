#ifndef POLYSELECT_PROOTS_H_
#define POLYSELECT_PROOTS_H_

#include <stdint.h>
#include <gmp.h>

#include "polyselect_poly_header.h"
#include "gmp_aux.h"
#include "polyselect_main_data.h"

#ifdef __cplusplus
extern "C" {
#endif

/* structure to store P roots */
struct polyselect_proots_s {
  unsigned long size;    /* used size -- this is the same as main->lenPrimes */
  uint8_t *nr;     /* number of roots of x^d = N (mod p) */
  uint64_t **roots; /* roots of (m0+x)^d = N (mod p^2) */
};
typedef struct polyselect_proots_s polyselect_proots_t[1];
typedef struct polyselect_proots_s * polyselect_proots_ptr;
typedef const struct polyselect_proots_s * polyselect_proots_srcptr;

extern void polyselect_proots_init (polyselect_proots_ptr, int, unsigned long);
extern void polyselect_proots_add (polyselect_proots_ptr, unsigned long, uint64_t*, unsigned long);
extern void polyselect_proots_print (polyselect_proots_srcptr);
extern void polyselect_proots_clear (polyselect_proots_ptr);

extern unsigned long polyselect_proots_compute(polyselect_proots_ptr R, polyselect_poly_header_srcptr header, polyselect_main_data_srcptr main, gmp_randstate_ptr rstate);

#ifdef __cplusplus
}
#endif

#endif	/* POLYSELECT_PROOTS_H_ */
