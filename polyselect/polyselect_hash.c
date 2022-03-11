#define TOKENPASTE(x, y) x ## y
#define TOKENPASTE2(x, y) TOKENPASTE(x, y)
#define LABEL_UNIQUE TOKENPASTE2(Label, __LINE__)

/* Data struct used for polyselect */
#include "cado.h" // IWYU pragma: keep
#include <float.h> // DBL_MAX
#include <math.h> // log
#include <string.h> // memset
#include <pthread.h> // pthread_mutex_lock
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "polyselect_hash.h"
#include "polyselect_thread.h"
#include "cado_poly.h"
#include "getprime.h"   // getprime
#include "misc.h"   // nprimes_interval
#include "macros.h"
#include "polyselect_match.h"

/* XXX this does not belong here! */
/* The following are primes used as factors in the special-q.
   Warning: if you add larger primes, you should ensure that any product
   still fits in an "unsigned long" (cf routine collision_on_sq).
   Here on a 32-bit machine the largest product of 4 primes is
   241*251*257*263 < 2^32, and on a 64-bit machine the largest product of 8
   primes is 233*239*241*251*257*263*269*271 < 2^64.
   LEN_SPECIAL_Q is defined in the header.
*/
#if ULONG_MAX == 4294967295UL
const unsigned int SPECIAL_Q[LEN_SPECIAL_Q] = {
  2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
  31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
  73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
  127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
  179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
  233, 239, 241, 251, 257, 263, 0 };
#else /* 64-bit */
const unsigned int SPECIAL_Q[LEN_SPECIAL_Q] = {
  2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
  31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
  73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
  127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
  179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
  233, 239, 241, 251, 257, 263, 269, 271, 0 };
#endif

//#define LESS_P

size_t expected_memory_usage_for_primes(unsigned long P)
{
  unsigned long Pmax = 2*P;
#ifdef LESS_P // if impatient for root finding
  Pmax = P + P/2;
#endif
  unsigned long maxprimes = nprimes_interval(P, Pmax);
  return maxprimes * sizeof (uint32_t);
}
    

/* init hash table */
void
polyselect_hash_init (polyselect_hash_ptr H, unsigned int init_size)
{
  memset(H, 0, sizeof(*H));
  H->alloc = init_size;
  H->slot = (struct polyselect_hash_slot_s *) malloc (H->alloc * sizeof (struct polyselect_hash_slot_s));
  if (H->slot == NULL) {
    fprintf (stderr, "Error, cannot allocate memory in polyselect_hash_init\n");
    exit (1);
  }

  for (unsigned int j = 0; j < H->alloc; j ++) {
    H->slot[j].i = 0;
    H->slot[j].p = 0;
  }
  H->size = 0;
#ifdef DEBUG_HASH_TABLE
  H->coll = 0;
  H->coll_all = 0;
#endif
}

/* rq is a root of N = (m0 + rq)^d mod (q^2) */
void
polyselect_hash_add (polyselect_hash_ptr H, unsigned long p, int64_t i,
          unsigned long q, mpz_srcptr rq,
          polyselect_thread_ptr loc)
{
  uint32_t h;

  ASSERT(H->size < H->alloc);

  h = (uint32_t) i % H->alloc;

#ifdef DEBUG_HASH_TABLE
  if (H->slot[h].i != 0)
    H->coll ++;
#endif
  while (H->slot[h].i != 0)
  {
      if (H->slot[h].i == i) {
          polyselect_match_info_ptr job;
          if (dllist_is_empty(&loc->empty_job_slots)) {
              job = malloc(sizeof(polyselect_match_info_t));
              polyselect_match_info_init(job, H->slot[h].p, p, i, q, rq, loc);
          } else {
              /* recycle an old one ! */
              struct dllist_head * ptr = dllist_get_first_node(&loc->empty_job_slots);
              dllist_pop(ptr);
              job = dllist_entry(ptr, struct polyselect_match_info_s, queue);
              polyselect_match_info_set(job, H->slot[h].p, p, i, q, rq, loc);
          }

          dllist_push_back(&loc->async_jobs, &job->queue);
      }
    if (UNLIKELY(++h == H->alloc))
      h = 0;
#ifdef DEBUG_HASH_TABLE
    H->coll_all ++;
#endif
  }
  H->slot[h].p = p;
  H->slot[h].i = i;
  H->size ++;
}

#if 0
void
polyselect_hash_find_collision_multi (const polyselect_hash_t * H, unsigned long p, int64_t i,
          unsigned long q, mpz_srcptr rq,
          polyselect_thread_ptr loc)
{
    /* XXX This sounds really really hard to code indeed. Not sure we'll
     * find a way that doesn't reimplement everything.
     *
     * Piggy-backing on top of the shash code is likely more robust than
     * to try to invent a reduction step for open hash tables (which
     * doesn't really seem to work, I'm afraid...)
     */
}
#endif

#if 0
/* no longer needed. But there seem to be some subtle differences. Do we
 * want them unified? */
/* rq is a root of N = (m0 + rq)^d mod (q^2) */
void
polyselect_hash_add_gmp (polyselect_hash_ptr H, uint32_t p, int64_t i, mpz_srcptr m0, mpz_srcptr ad,
              unsigned long d, mpz_srcptr N, uint64_t q, mpz_srcptr rq,
              polyselect_thread_locals_ptr loc)
{
  unsigned long h;

  if (H->size >= H->alloc)
    polyselect_hash_grow (H);
  if (i >= 0)
    h = ((int)i) % H->alloc;
  else
  {
    h = H->alloc - ( ((int)(-i)) % H->alloc);
    if (h == H->alloc)
      h = 0;
  }

  while (H->slot[h].p != 0)
  {
    if (m0 != NULL && H->slot[h].i== i && H->slot[h].p != p) {
      (*H->match) (H->slot[h].p, p, i, m0, ad, d, N, q, rq, loc);
    }
    if (++h == H->alloc)
      h = 0;
  }
  H->slot[h].p = p;
  H->slot[h].i = i;
  H->size ++;
}
#endif

void
polyselect_hash_clear (polyselect_hash_ptr H)
{
  free (H->slot);
}

#if 0
/* This is used only by the old polyselect_hash_add_gmp ; weird.
 */
static void
polyselect_hash_grow (polyselect_hash_ptr H)
{
  H->alloc = 2 * H->alloc + 16;
  H->slot = (struct polyselect_hash_slot_s *) realloc (H->slot, H->alloc * sizeof (struct polyselect_hash_slot_s));
  if (H->slot == NULL) {
    fprintf (stderr, "Error, cannot allocate memory in polyselect_hash_init\n");
    exit (1);
  }
  memset (H->slot + H->size, 0, sizeof(struct polyselect_hash_slot_s) * (H->alloc - H->size));
}
#endif
