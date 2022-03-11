#ifndef POLYSELECT_SHASH_H_
#define POLYSELECT_SHASH_H_

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
   of the CPU being used. The values below are probably misleaded by
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
polyselect_shash_add (polyselect_shash_t H, uint64_t i)
{
  *(H->current[i & (polyselect_SHASH_NBUCKETS - 1)])++ = i;
#if 1
  if (UNLIKELY(H->current[i & (polyselect_SHASH_NBUCKETS - 1)] >= H->base[(i & (polyselect_SHASH_NBUCKETS - 1)) + 1]))
    {
      fprintf (stderr, "polyselect_Shash bucket %" PRIu64 " is full.\n",
               i & (polyselect_SHASH_NBUCKETS - 1));
      H->current[i & (polyselect_SHASH_NBUCKETS - 1)]--;
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
polyselect_shash2_add (polyselect_shash_t H, uint64_t i, uint32_t p)
{
  unsigned int ib = i & (polyselect_SHASH_NBUCKETS - 1);
  H->pmem[H->current[ib] - H->mem] = p;
  *H->current[ib]++ = i;
#if 1
  if (UNLIKELY(H->current[i & (polyselect_SHASH_NBUCKETS - 1)] >= H->base[(i & (polyselect_SHASH_NBUCKETS - 1)) + 1]))
    {
      fprintf (stderr, "polyselect_Shash2 bucket %" PRIu64 " is full.\n",
               i & (polyselect_SHASH_NBUCKETS - 1));
      exit (1);
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

#endif	/* POLYSELECT_SHASH_H_ */
