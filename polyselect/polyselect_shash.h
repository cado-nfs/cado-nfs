#ifndef POLYSELECT_SHASH_H_
#define POLYSELECT_SHASH_H_

#include <stdint.h>     // int64_t
#include <stdlib.h>
#include <inttypes.h>
#include <gmp.h>
#include <limits.h> /* for ULONG_MAX */
#include "macros.h"     // LIKELY

#ifdef __cplusplus
extern "C" {
#endif


/* For the moment, this value is static. But it's highly critical for
   performance in the polyselect_shash table cribble:
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
  uint32_t alloc;      /* total allocated size */
  uint32_t balloc;     /* allocated size for each bucket */
};
typedef struct polyselect_shash_s polyselect_shash_t[1];
typedef struct polyselect_shash_s * polyselect_shash_ptr;
typedef const struct polyselect_shash_s * polyselect_shash_srcptr;

#ifndef INLINE
# if __GNUC__ && !__GNUC_STDC_INLINE__
#  define INLINE extern inline
# else
#  define INLINE inline
# endif
#endif

#ifndef EMIT_ADDRESSABLE_shash_add
INLINE
#endif
void
polyselect_shash_add (polyselect_shash_t H, uint64_t i)
{
  *(H->current[i & (polyselect_SHASH_NBUCKETS - 1)])++ = i;
  if (UNLIKELY(H->current[i & (polyselect_SHASH_NBUCKETS - 1)] >= H->base[(i & (polyselect_SHASH_NBUCKETS - 1)) + 1]))
    {
      fprintf (stderr, "polyselect_Shash bucket %" PRIu64 " is full.\n",
               i & (polyselect_SHASH_NBUCKETS - 1));
      exit (1);
    }
}

void polyselect_shash_init (polyselect_shash_ptr, unsigned int);
void polyselect_shash_reset (polyselect_shash_ptr);
size_t polyselect_shash_size(polyselect_shash_srcptr);
int polyselect_shash_find_collision (polyselect_shash_srcptr);
void polyselect_shash_clear (polyselect_shash_ptr);

#ifdef __cplusplus
}
#endif

#endif	/* POLYSELECT_SHASH_H_ */
