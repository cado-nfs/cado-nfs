#include "cado.h"
#include <pthread.h>
#include "polyselect_shash.h"
#include "misc.h"

static pthread_mutex_t lock=PTHREAD_MUTEX_INITIALIZER;
/* init_size is an approximation of the number of entries */
void
polyselect_shash_init (polyselect_shash_ptr H, unsigned int init_size)
{
  unsigned int init_size0 = init_size;

  /* round up to multiple of polyselect_SHASH_NBUCKETS */
  init_size = 1 + (init_size - 1) / polyselect_SHASH_NBUCKETS;
  init_size += init_size / 8 + 128; /* use 12.5% margin */
  if (init_size > init_size0)
    init_size = init_size0;
  H->alloc = init_size * (polyselect_SHASH_NBUCKETS + 1) + 8;
  /* + init_size for guard for the last buckets to avoid seg fault */
  /* + 8 for extreme guard (ASM X86 needs 8, C needs 5 when init_size is too small */
  H->mem = (uint64_t*) malloc (H->alloc * sizeof (uint64_t));
  if (!H->mem)
    {
      fprintf (stderr, "Error, cannot allocate memory in polyselect_shash_init\n");
      exit (1);
    }
  H->balloc = init_size;
}

size_t polyselect_shash_size(polyselect_shash_srcptr H)
{
    size_t r = 0;
    for(size_t i = 0 ; i < polyselect_SHASH_NBUCKETS ; i++)
        r += H->current[i] - H->base[i];
    return r;
}

void
polyselect_shash_reset (polyselect_shash_ptr H)
{
  H->base[0] = H->current[0] = H->mem;
  for (int j = 1; j <= polyselect_SHASH_NBUCKETS; j++)
    H->base[j] = H->current[j] = H->base[j-1] + H->balloc;
  /* Trick for prefetch T in polyselect_shash_find_collision after the end
     of the last bucket. Each H->base[j] has balloc entries of type uint64_t,
     where balloc >= 128.

     XXX several things are odd here
     What is "the last bucket" ?
     Is it [polyselect_SHASH_NBUCKETS-1] ?
     Is it [polyselect_SHASH_NBUCKETS] ?
     If the latter, then "the end of the last bucket" would be at position
     H->base[polyselect_SHASH_NBUCKETS] + balloc, and the reason why this
     doesn't overflow is that we have a +8 in H->alloc in the function
     above.
     If the former, then the place where we're doing the memset agrees
     with the description "after the end of the last bucket". There are
     several ways to argue that this memset doesn't overrun the buffer,
     including the one above, or the aforementioned +8. But then, the
     fact of allocating H->alloc with (polyselect_SHASH_NBUCKETS + 1)
     times the init_size would probably be a bug.
   */
  memset (H->base[polyselect_SHASH_NBUCKETS], 0, sizeof(**H->base) * 8);
}


/* return non-zero iff there is a collision */
int
polyselect_shash_find_collision (polyselect_shash_srcptr H)
{
  static uint32_t size = 0, mask;
  uint64_t *Hj, *Hjm;
  uint32_t *T;
  uint32_t k;

  /* XXX what does this maro do? Documentation needed!
   */
#define polyselect_SHASH_RESEARCH(TH,I)				\
  do {							\
    key = ((I) >> 32) + (I);				\
    if (UNLIKELY(*TH)) do {				\
      if (UNLIKELY(*TH == key)) { free (T); return 1; }	\
    } while (*(++TH));					\
    *TH = key;						\
  } while (0)

  /* XXX what does this macro do? Documentation needed!
   */
#define polyselect_SHASH_TH_I(TH,I,IND)			\
  do {						\
    I = Hj[IND];				\
    TH = T + ((I >> LN2SHASH_NBUCKETS) & mask); \
    __builtin_prefetch(TH, 1, 3);		\
  } while (0)					\

  pthread_mutex_lock (&lock);
  if (!size) {
    size = H->balloc << 1;
    /* round up to power of 2 */
    size --;
    while (size & (size - 1))
      size &= size - 1;
    size <<= 2;
    ASSERT_ALWAYS((size & (size - 1)) == 0);
    mask = size - 1;
    size += 16; /* Guard to avoid to test the end of polyselect_hash_table when ++TH */
  }

  pthread_mutex_unlock (&lock);
  T = (uint32_t*) malloc (size * sizeof(*T));
  for (k = 0; k < polyselect_SHASH_NBUCKETS; k++) {
    Hj = H->base[k];
    Hjm = H->current[k];
    if (Hj == Hjm) continue;
    memset (T, 0, size * sizeof(*T));
    /* Here, a special guard at the end of polyselect_shash_init allows
       until Hjm[polyselect_SHASH_BUCKETS-1] + 5.
       So, it's not needed to test if Hj + 4 < Hjm to avoid prefetch problem. */
    uint32_t *Th0, *Th1, *Th2, *Th3, *Th4;
    uint64_t i0, i1, i2, i3, i4;
    unsigned int key;

    polyselect_SHASH_TH_I(Th0, i0, 0);
    polyselect_SHASH_TH_I(Th1, i1, 1);
    polyselect_SHASH_TH_I(Th2, i2, 2);
    polyselect_SHASH_TH_I(Th3, i3, 3);
    polyselect_SHASH_TH_I(Th4, i4, 4);
    Hj += 5;
    while (LIKELY(Hj < Hjm)) {
      __builtin_prefetch(((void *) Hj) + 0x280, 0, 3);
      polyselect_SHASH_RESEARCH(Th0, i0); polyselect_SHASH_TH_I(Th0, i0, 0);
      polyselect_SHASH_RESEARCH(Th1, i1); polyselect_SHASH_TH_I(Th1, i1, 1);
      polyselect_SHASH_RESEARCH(Th2, i2); polyselect_SHASH_TH_I(Th2, i2, 2);
      polyselect_SHASH_RESEARCH(Th3, i3); polyselect_SHASH_TH_I(Th3, i3, 3);
      polyselect_SHASH_RESEARCH(Th4, i4); polyselect_SHASH_TH_I(Th4, i4, 4);
      Hj += 5;
    }
    switch (Hj - Hjm) { /* no break: it's NOT an error! */
    // coverity[unterminated_case]
    case 0: polyselect_SHASH_RESEARCH(Th4, i4); no_break();
    // coverity[unterminated_case]
    case 1: polyselect_SHASH_RESEARCH(Th3, i3); no_break();
    // coverity[unterminated_case]
    case 2: polyselect_SHASH_RESEARCH(Th2, i2); no_break();
    // coverity[unterminated_case]
    case 3: polyselect_SHASH_RESEARCH(Th1, i1); no_break();
    case 4: polyselect_SHASH_RESEARCH(Th0, i0); // no_break();
    }
  }
  free (T);
  return 0;
}
#undef polyselect_SHASH_TH_I
#undef polyselect_SHASH_RESEARCH

/* we know that i is in the bucket that we're currently examining, so it
 * has all its LN2SHASH_NBUCKETS low bits in common with others in the
 * same bucket. We're going to do a secondary in-memory dispatch of this
 * i into an open hash table, with plentiful storage that should be
 * sufficient to keep long runs to a minimum.  The starting point in that
 * open hash table is given by the bits that follow the LN2SHASH_NBUCKETS
 * low bits (sufficiently many of these bits so that on average, we
 * expect only very few contiguous entries).
 *
 * The key that gets actually stored in the open hash table is used as a
 * tie-breaker for match that end up being inserted nearby. nearby values
 * share their low LN2SHASH_NBUCKETS bits as well as the high bits of
 * their insert location in the open hash table (not the low bits, since
 * overruns _may_ happen). This leaves a few different bits (typically
 * something like 32-(LN2SHASH_NBUCKETS=8)-(log2(open hash table
 * size)=approx 16) = less than 8 in the low part, plus something like
 * log2(2*4*P^2/2^32) = up to approx 20 bits in the high part.
 *
 * A "sure" tie breaker would arrange these bits so that there's no
 * overlap. But on the other hand (and because false positives are
 * actually _not_ a problem), we can go along with the simpler addition
 * of the low and high part of i, and hope for the best.
 *
 * Note that this algorithm does not correctly deal with the situation
 * where the integer 0 is added to the hash table. Fortunately, for our
 * application I think that this can't happen.
 */
static inline uint64_t transform(uint64_t i, uint32_t mask)
{
    uint32_t bnum = (i >> LN2SHASH_NBUCKETS) & mask;
    uint32_t key = (uint32_t) ((i >> 32) + i);
    // uint32_t key = (i >> LN2SHASH_NBUCKETS) / (mask + 1);
    return ((uint64_t) bnum << 32) + key;
}

/* return non-zero iff there is a collision */
#define PREFETCH 256
int
MAYBE_UNUSED polyselect_shash_find_collision_old (polyselect_shash_t H)
{
  static uint32_t size = 0, mask;
  uint64_t data[PREFETCH], *pdata, *edata, *ldata;
  uint32_t *Th;
  uint32_t key;
  uint64_t *Hj, *Hjm;
  uint32_t *T;
  uint64_t i;
  unsigned int j, k, l;

  if (!size) {
    size = H->balloc << 1;
    /* round up to power of 2 */
    size --;
    while (size & (size - 1))
      size &= size - 1;
    size <<= 2;
    ASSERT_ALWAYS((size & (size - 1)) == 0);
    mask = (size - 1);
    size += 16;
  }
  T = (uint32_t*) malloc (size * sizeof(*T));
  edata = data + PREFETCH;
  for (k = 0; k < polyselect_SHASH_NBUCKETS; k++) {
    Hj = H->base[k];
    Hjm = H->current[k];
    if (Hj == Hjm) continue;
    memset (T, 0, size * sizeof(*T));
    pdata = data;
    j = Hjm - Hj;
    for (l = 0; l < PREFETCH * 8; l += 64)
      __builtin_prefetch(((void *) data) + l, 1, 3);
    if (j > PREFETCH) j = PREFETCH;
    for (l = 0; l < j; l++) {
      i = Hj[l];
      key = (uint32_t) ((i >> 32) + i);
      i = (i >> LN2SHASH_NBUCKETS) & mask;
      __builtin_prefetch (T + i, 1, 3);
      i = (i << 32) + key;
      data[l] = i;
    }
    Hj += j;
    if (LIKELY(j == PREFETCH)) {
      while (LIKELY(Hj != Hjm)) {
	__builtin_prefetch(((void *) Hj) + 0x280, 0, 3);
	i = *pdata;
	key = (uint32_t) i;
	Th = T + (i >> 32);
	if (UNLIKELY(*Th)) do {
	  if (UNLIKELY(*Th == key)) { free (T); return 1; }
	} while (*(++Th));
	*Th = key;
	i = *Hj++;
	key = (uint32_t) ((i >> 32) + i);
	i = (i >> LN2SHASH_NBUCKETS) & mask;
	__builtin_prefetch (T + i, 1, 3);
	i = (i << 32) + key;
	*pdata++ = i;
	if (UNLIKELY (pdata == edata)) pdata = data;
      }
      ldata = pdata;
    }
    else
      ldata = data + l;
    do {
      i = *pdata++;
      key = (uint32_t) i;
      Th = T + (i >> 32);
      if (UNLIKELY (pdata == edata)) pdata = data;
      if (UNLIKELY(*Th)) do {
	if (UNLIKELY(*Th == key)) { free (T); return 1; }
      } while (*(++Th));
      *Th = key;
    } while (LIKELY(ldata != pdata));
  }
  free (T);
  return 0;
}

void
polyselect_shash_clear (polyselect_shash_ptr H)
{
  free (H->mem);
}

