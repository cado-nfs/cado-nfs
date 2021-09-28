#ifndef POLYSELECT_COLLISIONS_H_
#define POLYSELECT_COLLISIONS_H_

#include <gmp.h>
#include "polyselect_locals.h"
#include "polyselect_shash.h"
#include "gmp_aux.h"

#ifdef __cplusplus
extern "C" {
#endif

/* find collisions between "P" primes, return number of loops */
extern unsigned long
collision_on_p(polyselect_shash_ptr H, polyselect_thread_locals_ptr loc);


/* collision on special-q, call collision_on_batch_sq */
extern void
collision_on_sq(unsigned long c, polyselect_shash_ptr H, polyselect_thread_locals_ptr loc);


#ifdef __cplusplus
}
#endif

#endif	/* POLYSELECT_COLLISIONS_H_ */
