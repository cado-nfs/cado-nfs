#ifndef POLYSELECT_STR_H
#define POLYSELECT_STR_H

#include <stdint.h>     // int64_t
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <gmp.h>
#include <limits.h> /* for ULONG_MAX */
#include "macros.h"     // LIKELY

/* This is used in the collisions calls */
#define POLYSELECT_HASH_ALLOC_RATIO 8UL

#if ULONG_MAX == 4294967295UL
#define LEN_SPECIAL_Q 57
#else
#define LEN_SPECIAL_Q 59
#endif
//#define DEBUG_HASH_TABLE
extern const unsigned int SPECIAL_Q[LEN_SPECIAL_Q];

struct polyselect_thread_s;

/* hash table slots */
struct polyselect_hash_slot_s
{
  int64_t i;               /* contains the values of r such that p^2
                              divides N - (m0 + r)^2 */
  uint32_t p;              /* contains the primes */
};


/* We seem to have different codes that react in subtly different ways on
 * the collisions, which is why we resort to using function pointers
 * here. Not absolutely sure that it's _really_ needed.
 */
typedef void (*polyselect_hash_match_t)(
              unsigned long p1, unsigned long p2, const int64_t i,
              uint64_t q, mpz_srcptr rq,
              struct polyselect_thread_s *);

/* hash table structure */
struct polyselect_hash_s
{
  struct polyselect_hash_slot_s *slot;
  unsigned int alloc;      /* total allocated size */
  unsigned int size;       /* number of entries in hash table */
#ifdef DEBUG_HASH_TABLE
  unsigned long coll;
  unsigned long coll_all;
#endif
};
typedef struct polyselect_hash_s polyselect_hash_t[1];
typedef struct polyselect_hash_s * polyselect_hash_ptr;
typedef const struct polyselect_hash_s * polyselect_hash_srcptr;

/* inline functions */

/* declarations */

extern const unsigned int SPECIAL_Q[];

#ifdef __cplusplus
extern "C" {
#endif

extern size_t expected_memory_usage_for_primes(unsigned long P);

extern void polyselect_hash_init (polyselect_hash_ptr, unsigned int);

extern void polyselect_hash_add (polyselect_hash_ptr, unsigned long, int64_t,
              unsigned long, mpz_srcptr,
              struct polyselect_thread_s * loc);

extern void polyselect_hash_clear (polyselect_hash_ptr);


#ifdef __cplusplus
}
#endif

#endif
