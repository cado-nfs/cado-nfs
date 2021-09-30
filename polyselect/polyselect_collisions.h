#ifndef POLYSELECT_COLLISIONS_H_
#define POLYSELECT_COLLISIONS_H_

#include <gmp.h>
#include "gmp_aux.h"
#include "polyselect_locals.h"
#include "polyselect_hash.h"
#include "polyselect_shash.h"

#ifdef __cplusplus
extern "C" {
#endif

/* find collisions between "P" primes, return number of loops */
extern unsigned long
collision_on_p(polyselect_shash_ptr H,
        polyselect_hash_match_t match,
        polyselect_thread_locals_ptr loc);


/* collision on special-q, call collision_on_batch_sq */
extern void
collision_on_sq(unsigned long c, polyselect_shash_ptr H,
        polyselect_hash_match_t match,
        polyselect_thread_locals_ptr loc);

extern void polyselect_proots_dispatch_to_shash_flat(
        polyselect_shash_ptr H,
        const uint32_t * Primes,
        size_t lenPrimes,
        const unsigned long * roots_per_prime,
        const uint8_t * number_of_roots_per_prime,
        int64_t umax);

extern void polyselect_proots_dispatch_to_hash_flat(
        polyselect_hash_ptr H,
        const uint32_t * Primes,
        size_t lenPrimes,
        const unsigned long * roots_per_prime,
        const uint8_t * number_of_roots_per_prime,
        int64_t umax,
        unsigned long q,
        mpz_srcptr rq,
        polyselect_thread_locals_ptr loc);

extern void polyselect_proots_dispatch_to_shash_notflat(
        polyselect_shash_ptr H,
        const uint32_t * Primes,
        size_t lenPrimes,
        uint64_t * const * roots_per_prime,
        const uint8_t * number_of_roots_per_prime,
        int64_t umax);

extern void polyselect_proots_dispatch_to_hash_notflat(
        polyselect_hash_ptr H,
        const uint32_t * Primes,
        size_t lenPrimes,
        uint64_t * const * roots_per_prime,
        const uint8_t * number_of_roots_per_prime,
        int64_t umax,
        unsigned long q,
        mpz_srcptr rq,
        polyselect_thread_locals_ptr loc
        );

#ifdef __cplusplus
}
#endif

#endif	/* POLYSELECT_COLLISIONS_H_ */
