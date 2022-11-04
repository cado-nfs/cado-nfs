#ifndef POLYSELECT_COLLISIONS_H_
#define POLYSELECT_COLLISIONS_H_

#include <gmp.h>
#include "gmp_aux.h"
#include "polyselect_thread.h"
#include "polyselect_shash.h"

#ifdef __cplusplus
extern "C" {
#endif


/* find collisions between "P" primes, return number of loops */

/* This _conductor version divide the work into pieces, and arranges so
 * that the different threads of the same team as the calling thread
 * participate in these work pieces collectively.
 */
extern unsigned long
collision_on_p_conductor(polyselect_thread_ptr thread);


/* collision on special-q, call collision_on_batch_sq */
extern void
collision_on_sq_conductor(unsigned long c, polyselect_thread_ptr thread);

/* TODO: reinstate the "not-conductor" version, for single-threaded work
 * (might be a good thing to have for debugging).
 */

#ifdef __cplusplus
}
#endif

#endif	/* POLYSELECT_COLLISIONS_H_ */
