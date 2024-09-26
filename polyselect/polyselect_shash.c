#include "cado.h"
#include <pthread.h>
#include "polyselect_shash.h"
#include "polyselect_match.h"
#include "misc.h"

/*
 * This is an implementation of a quick hash table with integer entries
 * that is meant to look for collisions among its integer entries. The
 * expected usage range is for integers up to about 50 bits, with at most
 * about 2^25 integers pushed to the hash table.
 *
 */

/* Entries are first pushed to buckets based on their low bits.
 *
 * Then buckets are examined one by one to look for collisions in them.
 *
 * For each integer i that went to a bucket, we know that i is in the
 * bucket that we're currently examining, so it has all its
 * LN2SHASH_NBUCKETS low bits in common with others in the same bucket.
 * We're going to do a secondary in-memory dispatch of this i into an
 * open hash table, with plentiful storage that should be sufficient to
 * keep long runs to a minimum.  The starting point in that open hash
 * table is given by the bits that follow the LN2SHASH_NBUCKETS low bits
 * (sufficiently many of these bits so that on average, we expect only
 * very few contiguous entries).
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

/* init_size is an approximation of the number of entries */
void
polyselect_shash_init (polyselect_shash_ptr H, unsigned int init_size)
{
    polyselect_shash_init_multi((polyselect_shash_t *) H, init_size, 1);
}


/* Same, but allocate several tables together. Such tables **MUST** be
 * cleared with polyselect_shash_clear_multi ; intermediary resizing in
 * order to reduce the "multi" parameter is possible, provided it never
 * exceeds the origin
 */

/* We want an alloc size S so that for any integer k such that 1 <= k <=
 * multi, we have (with N=polyselect_SHASH_NBUCKETS)
 *
 *  - S can be divided in k areas of size floor(S/k), aligned at an
 *  address that is a multiple of [large_alignment_constraint]
 */
#define SHASH_LARGE_ALIGNMENT_CONSTRAINT        0x40
/*
 *  - Each of these areas has room for N+1 buckets, + 8 bytes
 *  [ The +8 dates from 4778875a2 and is not explained ]
 */
#define SHASH_ALLOC_K_OFFSET    8
/*
 *  - each of the N+1 buckets has room for
 *      o(init_size/k/N) entries, with o(x) >= 1.25*x+128.
 *  [ The +1 buckets dates from ef98432867 (and 6b2a7998cf) and is not explained ]
 */
#define SHASH_BALLOC_CONSTANT_MARGIN    128
#define SHASH_BALLOC_INVMUL_MARGIN      4 /* x means 1+1/x */
/*
 *  - (optional) Each bucket start is aligned to a certain value (which
 *  can be 1.
 */
#define SHASH_SMALL_ALIGNMENT_CONSTRAINT        1
/*
 * The last condition is arranged by setting o(x) =
 * iceildiv(1.25*x+128,0x40)*small_alignment_constraint, which is upper
 * bounded by 1.25*x+8+small_alignment_constraint, which is sufficient to
 * get an upper bound afterwards.  However it's easy enough to loop over
 * the possible values of k as well.
 */
void
polyselect_shash_init_multi (polyselect_shash_t * H, unsigned int init_size, unsigned int multi)
{
    size_t alloc = 0;
    for(unsigned int k = 1 ; k <= multi ; k++) {
        size_t balloc = iceildiv(init_size, polyselect_SHASH_NBUCKETS * k);
        balloc += balloc / SHASH_BALLOC_INVMUL_MARGIN;
        balloc += SHASH_BALLOC_CONSTANT_MARGIN;
        balloc = next_multiple_of(balloc, SHASH_SMALL_ALIGNMENT_CONSTRAINT);
        size_t alloc_k = balloc * (polyselect_SHASH_NBUCKETS+1);
        alloc_k += SHASH_ALLOC_K_OFFSET;
        alloc_k = next_multiple_of(alloc_k, SHASH_LARGE_ALIGNMENT_CONSTRAINT);
        if (alloc_k * k >= alloc)
            alloc = alloc_k * k;
    }
    H[0]->alloc = alloc;
    H[0]->mem = (uint64_t*) malloc (H[0]->alloc * sizeof (uint64_t));
    H[0]->pmem = (uint32_t*) malloc (H[0]->alloc * sizeof (uint32_t));
    if (!H[0]->mem || !H[0]->pmem)
    {
        fprintf(stderr, "Error, cannot allocate memory in %s\n", __func__);
        exit (1);
    }

    polyselect_shash_reset_multi(H, multi);
}

void
polyselect_shash_reset_multi (polyselect_shash_t * H, unsigned int k)
{
    /* The alloc_k estimate that we had computed before is by
     * construction less than the result of the computation below.
     */
    size_t alloc_k = H[0]->alloc / k;
    alloc_k -= alloc_k & (SHASH_LARGE_ALIGNMENT_CONSTRAINT - 1);
    size_t balloc = (alloc_k - SHASH_ALLOC_K_OFFSET) / (polyselect_SHASH_NBUCKETS+1);
    balloc -= balloc & (SHASH_SMALL_ALIGNMENT_CONSTRAINT - 1);
    for(unsigned int i = 0 ; i < k ; i++) {
        H[i]->mem = H[0]->mem + i * alloc_k;
        H[i]->pmem = H[0]->pmem + i * alloc_k;
        H[i]->base[0] = H[i]->current[0] = H[i]->mem;
        /* This is only used for the spacing of the pointers */
        H[i]->balloc = balloc;
        for (int j = 1; j <= polyselect_SHASH_NBUCKETS; j++)
            H[i]->base[j] = H[i]->current[j] = H[i]->base[j-1] + H[i]->balloc;
        /*
         * Trick for prefetch T in polyselect_shash_find_collision after
         * the end of the last bucket. Each H->base[j] has balloc entries
         * of type uint64_t, where balloc >= 128.
         *
         * XXX several things are odd here
         * What is "the last bucket" ?
         * Is it [polyselect_SHASH_NBUCKETS-1] ?
         * Is it [polyselect_SHASH_NBUCKETS] ?
         * If the latter, then "the end of the last bucket" would be at position
         * H->base[polyselect_SHASH_NBUCKETS] + balloc, and the reason why this
         * doesn't overflow is that we have a +8 in H->alloc in the function
         * above.
         * If the former, then the place where we're doing the memset agrees
         * with the description "after the end of the last bucket". There are
         * several ways to argue that this memset doesn't overrun the buffer,
         * including the one above, or the aforementioned +8. But then, the
         * fact of allocating H->alloc with (polyselect_SHASH_NBUCKETS + 1)
         * times the init_size would probably be a bug.
         *
         */
        memset (H[i]->base[polyselect_SHASH_NBUCKETS], 0, sizeof(**H[0]->base) * 8);
    }
}


size_t polyselect_shash_size(polyselect_shash_srcptr H)
{
    size_t r = 0;
    for(size_t i = 0 ; i < polyselect_SHASH_NBUCKETS ; i++)
        r += H->current[i] - H->base[i];
    return r;
}

size_t polyselect_shash_size_fragment(polyselect_shash_srcptr H, size_t i0, size_t i1)
{
    size_t r = 0;
    for(size_t i = i0 ; i < i1 ; i++)
        r += H->current[i] - H->base[i];
    return r;
}

void
polyselect_shash_reset (polyselect_shash_ptr H)
{
    ASSERT_ALWAYS(H->alloc);
    polyselect_shash_reset_multi((polyselect_shash_t *) H, 1);
}

static inline size_t polyselect_shash_secondary_table_size(polyselect_shash_srcptr H, unsigned int multi)
{
    size_t size = next_power_of_2(H->balloc * multi) << 2;
    ASSERT_ALWAYS((size & (size - 1)) == 0);
    return size;
}

/* This is the most versatile version.
 *
 * Search for collision in the data that is stored in H[0]..H[multi-1],
 * looking *only* at buckets number k for k0 <= k < k1
 */
int
polyselect_shash_find_collision_multi(const polyselect_shash_t * H, unsigned int multi, uint32_t k0, uint32_t k1)
{

  /* XXX what does this macro do? Documentation needed!
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

  /* We implicitly assume that all H[0] .. to H[multi-1] share identical
   * allocation parameters.
   */
  uint32_t size = polyselect_shash_secondary_table_size(H[0], multi);
  uint32_t mask = size - 1;
  size += 16; /* Guard to avoid to test the end of polyselect_hash_table when ++TH */

  uint32_t * T = (uint32_t*) malloc (size * sizeof(*T));
  for (uint32_t k = k0; k < k1; k++) {
      memset (T, 0, size * sizeof(*T));
      for(unsigned int j = 0 ; j < multi ; j++) {
          const uint64_t * Hj = H[j]->base[k];
          const uint64_t * Hjm = H[j]->current[k];
          if (Hj == Hjm) continue;
          /* Here, a special guard at the end of polyselect_shash_init
           * allows
             until Hjm[polyselect_SHASH_BUCKETS-1] + 5.  So, it's not
             needed to test if Hj + 4 < Hjm to avoid prefetch problem. */
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
  }
  free (T);
  return 0;
}
#undef polyselect_SHASH_TH_I
#undef polyselect_SHASH_RESEARCH

/* {{{ polyselect_shash_find_collision (old version)
 * return non-zero iff there is a collision given the data that was
 * stored in the shash buckets 
 */
int
polyselect_shash_find_collision (polyselect_shash_srcptr H)
{
  return polyselect_shash_find_collision_multi((const polyselect_shash_t *) H, 1, 0, polyselect_SHASH_NBUCKETS);
}
/* }}} */


static inline uint64_t transform(uint64_t i, uint32_t mask)
{
    uint32_t bnum = (i >> LN2SHASH_NBUCKETS) & mask;
    uint32_t key = (uint32_t) ((i >> 32) + i);
    // uint32_t key = (i >> LN2SHASH_NBUCKETS) / (mask + 1);
    return ((uint64_t) bnum << 32) + key;
}

/* return non-zero iff there is a collision */
int
polyselect_shash_find_collision_simple_multi(const polyselect_shash_t * H, unsigned int multi, uint32_t k0, uint32_t k1)
{
    // uint64_t data0;
    uint32_t *Th;

    uint32_t size = polyselect_shash_secondary_table_size(H[0], multi);
    uint32_t mask = size - 1;
    size += 16; /* Guard to avoid to test the end of polyselect_hash_table when ++TH */

    uint32_t * T = (uint32_t*) malloc (size * sizeof(*T));
    for (uint32_t k = k0; k < k1; k++) {
        memset (T, 0, size * sizeof(*T));
        for(unsigned int j = 0 ; j < multi ; j++) {
            const uint64_t * Hj = H[j]->base[k];
            const uint64_t * Hjm = H[j]->current[k];
            if (Hj == Hjm) continue;

            for( ; Hj != Hjm ; Hj++) {
                uint64_t ix = transform(*Hj, mask);

                Th = T + (ix >> 32);
                for( ; *Th ; Th++) {
                    if (*Th == (uint32_t) ix) {
                        free (T);
                        return 1;
                    }
                }
                *Th = ix;
            }
        }
    }
    free (T);
    return 0;
}

int
polyselect_shash2_find_collision_multi(const polyselect_shash_t * H, unsigned int multi, uint32_t k0, uint32_t k1,
        unsigned long q, mpz_srcptr rq, polyselect_thread_ptr thread)
{
    // uint64_t data0;

    uint32_t size = polyselect_shash_secondary_table_size(H[0], multi);
    uint32_t mask = size - 1;
    size += 16; /* Guard to avoid to test the end of polyselect_hash_table when ++TH */

    struct slot {
        int64_t i;
        uint32_t p;
    };
    struct slot * T = (struct slot*) malloc (size * sizeof(struct slot));

    for (uint32_t k = k0; k < k1; k++) {
        memset (T, 0, size * sizeof(struct slot));
        for(unsigned int j = 0 ; j < multi ; j++) {
            const uint64_t * Hj = H[j]->base[k];
            const uint64_t * Hjm = H[j]->current[k];
            if (Hj == Hjm) continue;

            for( ; Hj != Hjm ; Hj++) {
                int64_t i = *Hj;
                uint32_t p2 = H[j]->pmem[Hj - H[j]->mem];
                uint64_t ix = transform(i, mask);

                struct slot * Th = T + (ix >> 32);

                for( ; Th->i ; Th ++) {
                    if (Th->i == i) {
                        /* This is a collision */
                        polyselect_match_info_ptr job;
                        uint32_t p1 = Th->p;
                        if (dllist_is_empty(&thread->empty_job_slots)) {
                            job = malloc(sizeof(polyselect_match_info_t));
                            polyselect_match_info_init(job, p1, p2, i, q, rq, thread);
                        } else {
                            /* recycle an old one ! */
                            struct dllist_head * ptr = dllist_get_first_node(&thread->empty_job_slots);
                            dllist_pop(ptr);
                            job = dllist_entry(ptr, struct polyselect_match_info_s, queue);
                            polyselect_match_info_set(job, p1, p2, i, q, rq, thread);
                        }
                        dllist_push_back(&thread->async_jobs, &job->queue);
                    }
                }
                Th->i = i;
                Th->p = p2;
            }
        }
    }
    free (T);
    return 0;
}

int
polyselect_shash_find_collision_simple (polyselect_shash_srcptr H)
{
  return polyselect_shash_find_collision_simple_multi((const polyselect_shash_t *) H, 1, 0, polyselect_SHASH_NBUCKETS);
}

/* return non-zero iff there is a collision */

/* not entirely clear that this was correctly adapted to the new
 * semantics */
#define PREFETCH 256
int
MAYBE_UNUSED polyselect_shash_find_collision_old_multi (const polyselect_shash_t * H, unsigned int multi, uint32_t k0, uint32_t k1)
{
  uint64_t data[PREFETCH], *pdata, *edata, *ldata;
  uint32_t *Th;
  uint32_t key;
  uint64_t *Hj, *Hjm;
  uint64_t i;
  unsigned int l;

  uint32_t size = polyselect_shash_secondary_table_size(H[0], multi);
  uint32_t mask = size - 1;
  size += 16; /* Guard to avoid to test the end of polyselect_hash_table when ++TH */


  uint32_t * T = (uint32_t*) malloc (size * sizeof(*T));
  edata = data + PREFETCH;
  for (uint32_t k = k0; k < k1; k++) {
    memset (T, 0, size * sizeof(*T));
    for(unsigned int j = 0 ; j < multi ; j++) {
        Hj = H[j]->base[k];
        Hjm = H[j]->current[k];
        if (Hj == Hjm) continue;
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
  }
  free (T);
  return 0;
}

int
polyselect_shash_find_collision_old (polyselect_shash_srcptr H)
{
  return polyselect_shash_find_collision_old_multi((const polyselect_shash_t *) H, 1, 0, polyselect_SHASH_NBUCKETS);
}


void
polyselect_shash_clear_multi (polyselect_shash_t * H, unsigned int multi MAYBE_UNUSED)
{
  free (H[0]->mem);
  free (H[0]->pmem);
}
void
polyselect_shash_clear (polyselect_shash_ptr H)
{
    ASSERT_ALWAYS(H->alloc);
    polyselect_shash_clear_multi((polyselect_shash_t *) H, 1);
}

void polyselect_shash2_print_fragment(polyselect_shash_srcptr H, size_t i0, size_t i1)
{
    for(size_t i = i0 ; i < i1 ; i++) {
        fprintf(stderr, "%zu:", i);
        for(const uint64_t * p = H->base[i] ; p != H->current[i] ; p++) {
            fprintf(stderr, " (%" PRId64 ",%" PRIu32 ")",
                    (int64_t) *p,
                    H->pmem[p-H->mem]);
        }
        fprintf(stderr, "\n");
    }
}

void polyselect_shash_print_fragment(polyselect_shash_srcptr H, size_t i0, size_t i1)
{
    for(size_t i = i0 ; i < i1 ; i++) {
        fprintf(stderr, "%zu:", i);
        for(const uint64_t * p = H->base[i] ; p != H->current[i] ; p++) {
            fprintf(stderr, " (%" PRId64 ",*)",
                    (int64_t) *p);
        }
        fprintf(stderr, "\n");
    }
}

