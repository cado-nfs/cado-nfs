#ifndef CADO_POLYSELECT_COLLISIONS_HPP
#define CADO_POLYSELECT_COLLISIONS_HPP

#include <gmp.h>
#include "gmp_aux.h"
#include "polyselect_thread.hpp"
#include "polyselect_shash.h"

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

#endif	/* CADO_POLYSELECT_COLLISIONS_HPP */
