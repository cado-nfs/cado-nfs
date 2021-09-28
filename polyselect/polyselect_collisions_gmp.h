#ifndef POLYSELECT_COLLISIONS_GMP_H_
#define POLYSELECT_COLLISIONS_GMP_H_

/* XXX I think that there is no real need for this interface, it should
 * go eventually.
 */
#ifdef __cplusplus
extern "C" {
#endif

#include <gmp.h>
#include "polyselect_locals.h"
#include "polyselect_shash.h"
#include "gmp_aux.h"

#ifdef __cplusplus
extern "C" {
#endif

/* find collisions between "P" primes, return number of loops */
extern unsigned long
gmp_collision_on_p(polyselect_thread_locals_ptr loc);


/* collision on special-q, call collision_on_batch_sq */
extern void
gmp_collision_on_sq(unsigned long c, polyselect_thread_locals_ptr loc);

#ifdef __cplusplus
}
#endif


#ifdef __cplusplus
}
#endif

#endif	/* POLYSELECT_COLLISIONS_GMP_H_ */
