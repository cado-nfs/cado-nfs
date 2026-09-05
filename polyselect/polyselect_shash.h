#ifndef CADO_POLYSELECT_SHASH_H
#define CADO_POLYSELECT_SHASH_H

#define xxxDEBUG_SHASH_OVERRUNS

#ifdef DEBUG_SHASH_OVERRUNS
#include <stdio.h>
#endif

#include <stdint.h>     // int64_t
#include <stdlib.h>
#include <inttypes.h>
#include <gmp.h>
#include <limits.h> /* for ULONG_MAX */
#include "macros.h"     // LIKELY

#define POLYSELECT_SHASH_ALLOC_RATIO 4UL

#ifdef __cplusplus
extern "C" {
#endif

struct polyselect_thread_s;

/* For the moment, this value is static. It's highly critical for
   performance. A priori, the right choice is guided by the number of TLB
   of the CPU being used. The values below are probably misled by
   other aspects ; quoting lingering comment:

   * 10 (or 9) seems the best for an Intel nehalem (?).
   * 6 seems the best for Intel Core2 (?).
   * 7 seems the best for AMD (?).
   So, the next optimization will include real time tests to
   evaluate the best value.
   NOTA: The good range is between 6 and 10. Don't use values <= 4!
   Values >= 12 are not interesting.
*/
#define LN2SHASH_NBUCKETS 8
#define polyselect_SHASH_NBUCKETS (1<<LN2SHASH_NBUCKETS)

struct polyselect_shash_s
{
  uint64_t *current[polyselect_SHASH_NBUCKETS+1]; /* +1 for guard */
  uint64_t *base[polyselect_SHASH_NBUCKETS+1];    /* +1 for guard */
  uint64_t *mem;
  uint32_t *pmem;
  uint32_t alloc;      /* total allocated size */
  uint32_t balloc;     /* allocated size for each bucket */
#ifdef DEBUG_SHASH_OVERRUNS
  /* this structure is attached to one thread doing the work. Let's add
   * some debug info.
   */
  unsigned int it, nt;
  double average_hits;
  double my_expected_hits;
#endif
};
typedef struct polyselect_shash_s polyselect_shash_t[1];
typedef struct polyselect_shash_s * polyselect_shash_ptr;
typedef const struct polyselect_shash_s * polyselect_shash_srcptr;

#ifndef EMIT_ADDRESSABLE_shash_add
# if __GNUC__ && !__GNUC_STDC_INLINE__
extern inline
# else
inline
# endif
#endif
void
polyselect_shash_add (polyselect_shash_t H, uint64_t i, int64_t ppl MAYBE_UNUSED)
{
  const unsigned int ib = i & (polyselect_SHASH_NBUCKETS - 1);
  *(H->current[ib])++ = i >> LN2SHASH_NBUCKETS;
#ifdef DEBUG_SHASH_OVERRUNS
  const uint64_t * cur = H->current[ib];
  const uint64_t * b0 = H->base[ib];
  const uint64_t * b1 = H->base[ib+1];
  if (UNLIKELY(cur >= b1)) {
      fprintf (stderr, "polyselect_Shash bucket %u is full (%zd/%zd)"
              " on i=0x%" PRIx64 " p2=0x%" PRIx64 "."
              " Called from thread %u/%u. Average hit rate %.1f, actual %.1f\n",
              ib, cur-b0, b1-b0,
              i, ppl,
              H->it, H->nt, H->average_hits, H->my_expected_hits);
      // H->current[i & (polyselect_SHASH_NBUCKETS - 1)]--;
      ASSERT_ALWAYS((cur-b0) < 2*(b1-b0));
      // exit (1);
  }
#endif
}

#ifndef EMIT_ADDRESSABLE_shash_add
# if __GNUC__ && !__GNUC_STDC_INLINE__
extern inline
# else
inline
# endif
#endif
void
polyselect_shash2_add (polyselect_shash_t H, uint64_t i, uint32_t p, int64_t ppl MAYBE_UNUSED)
{
  const unsigned int ib = i & (polyselect_SHASH_NBUCKETS - 1);
  H->pmem[H->current[ib] - H->mem] = p;
  *H->current[ib]++ = i >> LN2SHASH_NBUCKETS;
#ifdef DEBUG_SHASH_OVERRUNS
  const uint64_t * cur = H->current[ib];
  const uint64_t * b0 = H->base[ib];
  const uint64_t * b1 = H->base[ib+1];
  if (UNLIKELY(cur >= b1)) {
      fprintf (stderr, "polyselect_Shash2 bucket %u is full (%zd/%zd)"
              " on i=0x%" PRIx64 " p2=0x%" PRIx64 "."
              " Called from thread %u/%u. Average hit rate %.1f, actual %.1f\n",
              ib, cur-b0, b1-b0,
              i, ppl,
              H->it, H->nt, H->average_hits, H->my_expected_hits);
      ASSERT_ALWAYS((cur-b0) < 2*(b1-b0));
      // exit (1);
  }
#endif
}

extern void polyselect_shash_init (polyselect_shash_ptr, unsigned int);
extern void polyselect_shash_init_multi (polyselect_shash_t *, unsigned int, unsigned int);
extern void polyselect_shash_reset (polyselect_shash_ptr);
extern size_t polyselect_shash_size(polyselect_shash_srcptr);
/* polyselect_shash_find_collision is deprecated because it should rather
 * be replaced with calls such as
 *
 * polyselect_thread_team_post_work(thread->team, thread, polyselect_DCS_notflat
_subtask, &found);
 *
 */
extern int polyselect_shash_find_collision (polyselect_shash_srcptr H)
#ifndef EXPOSE_DEPRECATED_polyselect_shash_find_collision
    ATTRIBUTE_DEPRECATED
#endif
    ;
extern int polyselect_shash_find_collision_multi(const polyselect_shash_t * H, unsigned int multi, uint32_t k0, uint32_t k1);
extern void polyselect_shash_clear (polyselect_shash_ptr);
extern void polyselect_shash_clear_multi (polyselect_shash_t * H, unsigned int multi);
extern int polyselect_shash2_find_collision_multi(const polyselect_shash_t * H, unsigned int multi, uint32_t k0, uint32_t k1,
        unsigned long q, mpz_srcptr rq, struct polyselect_thread_s * thread);
extern void
polyselect_shash_reset_multi (polyselect_shash_t * H, unsigned int multi);

#ifdef __cplusplus
}
#endif

#endif	/* CADO_POLYSELECT_SHASH_H */
