#include "cado.h" // IWYU pragma: keep

/* This compilation units reacts to TRACK_CODE_PATH and uses macros
 * such as WHERE_AM_I_UPDATE.
 * This compilation unit _must_ produce different object files depending
 * on the value of TRACK_CODE_PATH.
 * The WHERE_AM_I_UPDATE macro itself is defined in las-where-am-i.hpp
 */

#include <cmath>                   // for sqrt
#include <cstdint>                 // for uint32_t, uint64_t, UINT32_C
#include <cinttypes>               // for PRIu32
#include <cstdlib>                 // for size_t, NULL, free, realloc
#include <cstring>                 // for memset
#include <new>                     // for bad_alloc
#ifdef TRACE_K
#include <type_traits>             // for is_same
#endif

#include "bucket.hpp"              // for bucket_array_t, bucket_update_t
#include "memory.h"             // free_aligned
#include "verbose.h"             // verbose_output_print

#include "bucket-push-update.hpp"  // for bucket_single::push_update
#include "fb-types.h"              // for slice_index_t, FBPRIME_FORMAT, fbp...
#include "fb.hpp"                  // for fb_factorbase, fb_factorbase::slicing
#include "iqsort.h"                // for QSORT
#include "las-config.h"            // for BUCKET_REGIONS, LOG_BUCKET_REGIONS
#include "las-where-am-i.hpp"      // for where_am_I, WHERE_AM_I_UPDATE
#include "las-memory.hpp"          // for las_memory_accessor
#include "las-where-am-i-proxy.hpp" // IWYU pragma: keep    // for where_am_I
#include "macros.h"                // for MAYBE_UNUSED, ASSERT_ALWAYS, UNLIKELY
#ifdef TRACE_K
#include "las-output.hpp"          // for TRACE_CHANNEL
#endif

template <int LEVEL, typename HINT> struct bucket_update_t;


static size_t
bucket_misalignment(const size_t sz, const size_t sr MAYBE_UNUSED) {
  size_t size = sz; 
#ifdef HAVE_SSE2
  /* Round up to a multiple of CACHELINESIZE to make SSE2 aligned accesses
     work */
  size = iceildiv(size, 16U) * 16U;
#endif
  return size;
}

/* Set the read and write pointers of the buckets back to the respective bucket
   start, and set nr_slices back to 0. */
template <int LEVEL, typename HINT>
void
bucket_array_t<LEVEL, HINT>::reset_pointers()
{
  aligned_medium_memcpy (bucket_write, bucket_start, size_b_align);
  aligned_medium_memcpy (bucket_read,  bucket_start, size_b_align);
  nr_slices = 0;
}

template <int LEVEL, typename HINT>
bucket_array_t<LEVEL, HINT>::~bucket_array_t()
{
  if (used_accessor) used_accessor->physical_free (big_data, big_size);
  free (slice_index);
  free_aligned(slice_start);
  free_aligned(bucket_read);
  free_aligned(bucket_start);
  free_pagealigned(bucket_write);
}

template <int LEVEL, typename HINT>
void
bucket_array_t<LEVEL, HINT>::move(bucket_array_t<LEVEL, HINT> &other)
{
#define MOVE_ENTRY(x, zero) do {x = other.x; other.x = zero;} while(0)
  MOVE_ENTRY(big_data, NULL);
  MOVE_ENTRY(big_size, 0);
  MOVE_ENTRY(bucket_write, NULL);
  MOVE_ENTRY(bucket_start, NULL);
  MOVE_ENTRY(bucket_read, NULL);
  MOVE_ENTRY(slice_index, NULL);
  MOVE_ENTRY(slice_start, NULL);
  MOVE_ENTRY(n_bucket, 0);
  MOVE_ENTRY(size_b_align, 0);
  MOVE_ENTRY(nr_slices, 0);
  MOVE_ENTRY(alloc_slices, 0);
#undef MOVE_ENTRY
}

/* Allocate enough memory to be able to store _n_bucket buckets, each of at
   least min_bucket_size entries. If enough (or more) memory was already
   allocated, does not shrink the allocation. */
template <int LEVEL, typename HINT>
void
bucket_array_t<LEVEL, HINT>::allocate_memory(
        las_memory_accessor & memory,
        const uint32_t new_n_bucket,
        const double fill_ratio,
        int logI,
        const slice_index_t prealloc_slices)
{
    used_accessor = &memory;
  /* Don't try to allocate anything, nor print a message, for sieving levels
     where the corresponding factor base part is empty. */
  if (fill_ratio == 0.)
    return;

  /* We'll allocate bucket regions of size 2^LOG_BUCKET_REGIONS[LEVEL],
   * and those will be used to cover lines of size 2^logI.
   *
   * If LOG_BUCKET_REGIONS[LEVEL] > logI, then we'll have both
   * even and odd lines in there, so 75% of the locations will receive
   * updates.
   *
   * If LOG_BUCKET_REGIONS[LEVEL] <= logI, it's different. The
   * line ordinate for bucket region number N will be N, maybe shifted
   * right by a few bits (logI - LOG_BUCKET_REGIONS[LEVEL]).
   * Bucket regions for which this line ordinate is even will receive
   * updates for 50% of the locations only, in contrast to 100% when the
   * ordinate is even.
   */

  const size_t Q = 0.25 * fill_ratio * BUCKET_REGIONS[LEVEL];

  size_t bs_even, bs_odd;
  size_t new_big_size;
  /* The number of updates has a variance which is very often close to
   * its mean. By aiming at a number of updates that is at most
   * N+ndev*sqrt(N), the probability to overrun is at most the probability
   * that this random variable exceeds its mean by ndev standard deviations
   * or more. Now it's a sum of random variables (number of hits per line
   * for p) that are independent but not identically distributed.
   * Probably there are some tail bounds which can tell me how rare this
   * is, depending on ndev.
   */
  const int ndev = NB_DEVIATIONS_BUCKET_REGIONS;

  uint32_t bitmask_line_ordinate = 0;
  if (LOG_BUCKET_REGIONS[LEVEL] <= logI) {
      bitmask_line_ordinate = UINT32_C(1) << (logI - LOG_BUCKET_REGIONS[LEVEL]);
      ASSERT_ALWAYS(new_n_bucket % 2 == 0);
      bs_even = 2 * Q + ndev * sqrt(2 * Q);
      bs_odd  = 4 * Q + ndev * sqrt(4 * Q);
      bs_even = bucket_misalignment(bs_even, sizeof(update_t));
      bs_odd  = bucket_misalignment(bs_odd, sizeof(update_t));
      new_big_size = (bs_even + bs_odd) * (new_n_bucket / 2) * sizeof(update_t);
  } else {
      bs_even = bs_odd = 3 * Q + ndev * sqrt(3 * Q);
      bs_even = bucket_misalignment(bs_even, sizeof(update_t));
      bs_odd  = bucket_misalignment(bs_odd, sizeof(update_t));
      new_big_size = bs_odd * new_n_bucket * sizeof(update_t);
  }

  /* add 1 megabyte to each bucket array, so as to deal with overflowing
   * buckets */
  new_big_size += 1 << 20;

  const size_t new_size_b_align = ((sizeof(void *) * new_n_bucket + 0x3F) & ~((size_t) 0x3F));

  if (new_big_size > big_size) {
    if (big_data != NULL)
      memory.physical_free (big_data, big_size);
    if (bitmask_line_ordinate) {
        verbose_output_print(0, 3, "# [%d%c] Allocating %zu bytes for %" PRIu32 " buckets of %zu to %zu update entries of %zu bytes each\n",
                             LEVEL, HINT::rtti[0],
                             new_big_size, new_n_bucket,
                             bs_even, bs_odd,
                             sizeof(update_t));
    } else {
        verbose_output_print(0, 3, "# [%d%c] Allocating %zu bytes for %" PRIu32 " buckets of %zu update entries of %zu bytes each\n",
                             LEVEL, HINT::rtti[0],
                             new_big_size, new_n_bucket, 
                             bs_even,
                             sizeof(update_t));
    }
    big_size = new_big_size;
    big_data = (update_t *) memory.physical_alloc (big_size, 1);
    void * internet_of_things MAYBE_UNUSED = NULL;
  }

  if (!big_data)
      throw std::bad_alloc();

  // bucket_size = new_bucket_size;
  n_bucket = new_n_bucket;

  if (new_size_b_align > size_b_align) {
    int const all_null =  (bucket_write == NULL) &&
                    (bucket_start == NULL) &&
                    (bucket_read == NULL);
    int const none_null = bucket_write &&
                    bucket_start &&
                    bucket_read;
    ASSERT_ALWAYS(all_null || none_null);
    if (none_null) {
        verbose_output_print(0, 1, "# [%d%c] Changing bucket allocation from %zu bytes to %zu bytes\n", LEVEL, HINT::rtti[0], size_b_align, new_size_b_align);
        free_pagealigned(bucket_write);
        free_aligned(bucket_start);
        free_aligned(bucket_read);
    }
    size_b_align = new_size_b_align;
    bucket_write = (update_t **) malloc_pagealigned (size_b_align);
    /* bucket_start is allocated as an array of n_bucket+1 pointers */
    size_t const alloc_bstart = MAX(size_b_align, (n_bucket+1) * sizeof(void*));
    bucket_start = (update_t **) malloc_aligned (alloc_bstart, 0x40);
    bucket_read = (update_t **) malloc_aligned (size_b_align, 0x40);
    memset(bucket_write, 0, size_b_align);
    memset(bucket_start, 0, alloc_bstart);
    memset(bucket_read, 0, size_b_align);
    if (none_null) {
        /* must refresh this pointer too */
        free_slice_start();
    }
  }

  /* This requires size_b_align to have been set to the new value */
  if (prealloc_slices > alloc_slices)
    realloc_slice_start(prealloc_slices - alloc_slices);

  /* Spread bucket_start pointers equidistantly over the big_data array */
  update_t * cur = big_data;
  bucket_start[0] = cur;
  for (uint32_t i = 0; bucket_start[i]=cur, i < n_bucket; i++) {
      cur += ((i & bitmask_line_ordinate) != 0) ? bs_odd : bs_even;
  }
  reset_pointers();
#ifdef SAFE_BUCKET_ARRAYS
  verbose_output_print(0, 0, "# WARNING: SAFE_BUCKET_ARRAYS is on !\n");
#endif
}

template <int LEVEL, typename HINT>
void
bucket_array_t<LEVEL, HINT>::free_slice_start()
{
  free_aligned(slice_start); slice_start = NULL;
  free(slice_index); slice_index = NULL;
  alloc_slices = 0;
}

template <int LEVEL, typename HINT>
void
bucket_array_t<LEVEL, HINT>::realloc_slice_start(const size_t extra_space)
{
  const size_t new_alloc_slices = alloc_slices + extra_space;
  if (alloc_slices)
  verbose_output_print(0, 3, "# [%d%c] Reallocating BA->slice_start from %zu entries to %zu entries\n",
                       LEVEL, HINT::rtti[0],
                       alloc_slices, new_alloc_slices);

  const size_t old_size = size_b_align * alloc_slices;
  const size_t new_size = size_b_align * new_alloc_slices;
  slice_start = (update_t **) realloc_aligned(slice_start, old_size, new_size, 0x40);
  ASSERT_ALWAYS(slice_start != NULL);
  slice_index = (slice_index_t *) realloc(slice_index, new_alloc_slices * sizeof(slice_index_t));
  ASSERT_ALWAYS(slice_index != NULL);
  memset(slice_index + alloc_slices, 0, extra_space * sizeof(slice_index_t));
  alloc_slices = new_alloc_slices;
}

/* Returns how full the fullest bucket is, as a fraction of its size */
template <int LEVEL, typename HINT>
double
bucket_array_t<LEVEL, HINT>::max_full (unsigned int * fullest_index) const
{
  double max = 0;
  for (unsigned int i = 0; i < n_bucket; ++i)
    {
      double const j = (double) nb_of_updates (i) / room_allocated_for_updates(i);
      if (max < j) {
          max = j;
          if (fullest_index) *fullest_index = i;
      }
    }
  return max;
}
template <int LEVEL, typename HINT>
double
bucket_array_t<LEVEL, HINT>::average_full() const
{
    size_t a = 0;
    for (unsigned int i = 0; i < n_bucket; ++i)
        a += nb_of_updates (i);
    return (double) a / (bucket_start[n_bucket] - bucket_start[0]);
}

template <int LEVEL, typename HINT>
void
bucket_array_t<LEVEL, HINT>::log_this_update (
    const update_t update MAYBE_UNUSED,
    const uint64_t offset MAYBE_UNUSED,
    const uint64_t bucket_number MAYBE_UNUSED,
    where_am_I & w MAYBE_UNUSED) const
{
#if defined(TRACE_K)
    size_t  const(&BRS)[FB_MAX_PARTS] = BUCKET_REGIONS;
    unsigned int const saveN = w->N;
    /* flatten the (N,x) coordinate as if relative to a unique array of
     * level-1 bucket regions */
    unsigned int const x = update.x % BRS[1];
    unsigned int const N = w->N +
        bucket_number*BRS[LEVEL]/BRS[1] + (update.x / BRS[1]);

    WHERE_AM_I_UPDATE(w, x, x);
    WHERE_AM_I_UPDATE(w, N, N);

    if (trace_on_spot_Nx(w->N, w->x)) {
        verbose_output_print (TRACE_CHANNEL, 0,
            "# Pushed hit at location (x=%u, side %d), from factor base entry "
            "(slice_index=%u, slice_offset=%u, p=%" FBPRIME_FORMAT "), "
            "to BA<%d>[%u]\n",
            (unsigned int) w->x, w->side, (unsigned int) w->i,
            (unsigned int) w->h, w->p, LEVEL, (unsigned int) w->N);
        if (std::is_same<HINT,longhint_t>::value) {
          verbose_output_print (TRACE_CHANNEL, 0,
             "# Warning: did not check divisibility during downsorting p=%"
             FBPRIME_FORMAT "\n", w->p);
        } else {
          ASSERT_ALWAYS(test_divisible(w));
        }
    }
    WHERE_AM_I_UPDATE(w, N, saveN);
#endif
}




#if 0
/* no longer used */
/* A compare function suitable for sorting updates in order of ascending x
   with qsort() */
static inline int
bucket_cmp_x (const bucket_update_t<1, longhint_t>  *a, const bucket_update_t<1, longhint_t>  *b)
{
  if (a->x < b->x)
    return -1;
  if (a->x == b->x)
    return 0;
  return 1;
}
#endif

template <int LEVEL, typename HINT>
void
bucket_single<LEVEL, HINT>::sort()
{
//  qsort (start, write - start, sizeof (bucket_update_t<1, longhint_t> ),
//	 (int(*)(const void *, const void *)) &bucket_cmp_x);
#define islt(a,b) ((a)->x < (b)->x)
  QSORT(update_t, start, write - start, islt);
#undef islt  
}

void
bucket_primes_t::purge (const bucket_array_t<1, shorthint_t> &BA,
              const int i, fb_factorbase::slicing const & fb, const unsigned char *S)
{
    for (slice_index_t i_slice = 0; i_slice < BA.get_nr_slices(); i_slice++) {
        const slice_index_t slice_index = BA.get_slice_index(i_slice);
        for(auto const & it : BA.slice_range(i, i_slice)) {
            if (UNLIKELY(S[it.x] != 255)) {
                fbprime_t const p = fb[slice_index].get_prime(it.hint);
                push_update(bucket_update_t<1, primehint_t>(it.x, p, 0, 0));
            }
        }
    }
}

/*
template <typename HINT>
static inline bucket_update_t<1, longhint_t>
to_longhint(const bucket_update_t<1, HINT> &update, slice_index_t slice_index);
*/

static inline
bucket_update_t<1, longhint_t>
to_longhint(const bucket_update_t<1, shorthint_t> &update,
                         const slice_index_t slice_index)
{
  return bucket_update_t<1, longhint_t> (update.x, 0, update.hint, slice_index);
}

static inline
bucket_update_t<1, longhint_t>
to_longhint(const bucket_update_t<1, longhint_t> &update,
                        const slice_index_t slice_index MAYBE_UNUSED)
{
  return update;
}


template <typename HINT>
void
bucket_array_complete::purge(
    const bucket_array_t<1, HINT> &BA,
    const int i, const unsigned char *S)
{
    for (slice_index_t i_slice = 0; i_slice < BA.get_nr_slices(); i_slice++) {
        const slice_index_t slice_index = BA.get_slice_index(i_slice);
        for(auto const & it : BA.slice_range(i, i_slice)) {
            if (UNLIKELY(S[it.x] != 255)) {
                push_update(to_longhint(it, slice_index));
            }
        }
    }
}

template void
bucket_array_complete::purge<shorthint_t>(
    const bucket_array_t<1, shorthint_t> &BA,
    const int i, const unsigned char *S);

template void
bucket_array_complete::purge<longhint_t>(
    const bucket_array_t<1, longhint_t> &BA,
    const int i, const unsigned char *S);


#if defined(HAVE_SSE2) && defined(SMALLSET_PURGE)

template <typename HINT, int SIZE>
void
bucket_array_complete::purge_1 (
    const bucket_array_t<1, HINT> &BA, const int i,
    const std::vector<typename bucket_update_t<1, HINT>::br_index_t> &survivors)
{
    smallset<SIZE, typename bucket_update_t<1, HINT>::br_index_t> surv_set(survivors);

    for (slice_index_t i_slice = 0; i_slice < BA.get_nr_slices(); i_slice++) {
        const slice_index_t slice_index = BA.get_slice_index(i_slice);
        for(auto const & it : BA.slice_range(i, i_slice)) {
            if (UNLIKELY(surv_set.contains(it.x))) {
                push_update(to_longhint(it, slice_index));
            }
        }
    }
}

template <typename HINT>
void
bucket_array_complete::purge (
    const bucket_array_t<1, HINT> &BA, const int i,
    const unsigned char *S,
    const std::vector<typename bucket_update_t<1, HINT>::br_index_t> &survivors)
{
  const size_t items_per_size = smallset<1, typename bucket_update_t<1, HINT>::br_index_t>::nr_items;
  if (survivors.size() == 0)
    return;
  switch ((survivors.size() + items_per_size - 1) / items_per_size) {
    case 1: purge_1<HINT, 1>(BA, i, survivors); break;
    case 2: purge_1<HINT, 2>(BA, i, survivors); break;
    case 3: purge_1<HINT, 3>(BA, i, survivors); break;
    case 4: purge_1<HINT, 4>(BA, i, survivors); break;
    case 5: purge_1<HINT, 5>(BA, i, survivors); break;
    case 6: purge_1<HINT, 6>(BA, i, survivors); break;
    case 7: purge_1<HINT, 7>(BA, i, survivors); break;
    case 8: purge_1<HINT, 8>(BA, i, survivors); break;
    default: purge(BA, i, S);
  }
}

template
void
bucket_array_complete::purge<shorthint_t> (
    const bucket_array_t<1, shorthint_t> &, int,
    const unsigned char *,
    const std::vector<typename bucket_update_t<1, shorthint_t>::br_index_t> &);
#endif

template<int INPUT_LEVEL>
void
downsort(fb_factorbase::slicing const & fbs MAYBE_UNUSED,
        bucket_array_t<INPUT_LEVEL - 1, longhint_t> &BA_out,
        const bucket_array_t<INPUT_LEVEL, shorthint_t> &BA_in,
        uint32_t bucket_number, where_am_I & w)
{
    /* Time recording for this function is done by the caller
     * (downsort_wrapper)
     */
    /* Rather similar to purging, except it doesn't purge */
    for (slice_index_t i_slice = 0; i_slice < BA_in.get_nr_slices(); i_slice++) {
        const slice_index_t slice_index = BA_in.get_slice_index(i_slice);
        WHERE_AM_I_UPDATE(w, i, slice_index);
        for(auto const & it : BA_in.slice_range(bucket_number, i_slice)) {
            WHERE_AM_I_UPDATE(w, p, fbs[slice_index].get_prime(it.hint));
            WHERE_AM_I_UPDATE(w, h, it.hint);
            BA_out.push_update(it.x, 0, it.hint, slice_index, w);
        }
    }
}

template<int INPUT_LEVEL>
void
downsort(fb_factorbase::slicing const & /* unused */,
        bucket_array_t<INPUT_LEVEL - 1, longhint_t> &BA_out,
         const bucket_array_t<INPUT_LEVEL, longhint_t> &BA_in,
         uint32_t bucket_number, where_am_I & w) 
{
    /* longhint updates don't write slice end pointers, so there must be
       exactly 1 slice per bucket */
    ASSERT_ALWAYS(BA_in.get_nr_slices() == 1);

    for(auto const & it : BA_in.slice_range(bucket_number, 0)) {
        BA_out.push_update(it.x, 0, it.hint, it.index, w);
    }
}

template<int INPUT_LEVEL>
void
downsort(fb_factorbase::slicing const & fbs,
        bucket_array_t<INPUT_LEVEL - 1, logphint_t> &BA_out,
        const bucket_array_t<INPUT_LEVEL, emptyhint_t> &BA_in,
        uint32_t bucket_number, where_am_I & w)
{
  /* Time recording for this function is done by the caller
   * (downsort_wrapper)
   */
  /* Rather similar to purging, except it doesn't purge */
  for (slice_index_t i_slice = 0; i_slice < BA_in.get_nr_slices(); i_slice++) {
    const slice_index_t slice_index = BA_in.get_slice_index(i_slice);
    // WHERE_AM_I_UPDATE(w, i, slice_index);
    for (auto const & it : BA_in.slice_range(bucket_number, i_slice)) {
        logphint_t const h = fbs[slice_index].get_logp();
        BA_out.push_update(it.x, h, w);
    }
  }
}

template<int INPUT_LEVEL>
void
downsort(fb_factorbase::slicing const & /* unused */,
        bucket_array_t<INPUT_LEVEL - 1, logphint_t> &BA_out,
        const bucket_array_t<INPUT_LEVEL, logphint_t> &BA_in,
        uint32_t bucket_number, where_am_I & w) 
{
    /* longhint updates don't write slice end pointers, so there must be
       exactly 1 slice per bucket */
    ASSERT_ALWAYS(BA_in.get_nr_slices() == 1);
    for (auto const & it : BA_in.slice_range(bucket_number, 0)) {
        BA_out.push_update(it.x, it, w);
    }
}

/* Explicitly instantiate the versions of downsort() that we'll need:
   downsorting shorthint from level 3 and level 2, and downsorting
   longhint from level 2. */
template
void
downsort<2>(fb_factorbase::slicing const &,
        bucket_array_t<1, longhint_t> &BA_out,
        const bucket_array_t<2, shorthint_t> &BA_in,
        uint32_t bucket_number, where_am_I & w);

template
void
downsort<3>(fb_factorbase::slicing const &,
        bucket_array_t<2, longhint_t> &BA_out,
        const bucket_array_t<3, shorthint_t> &BA_in,
        uint32_t bucket_number, where_am_I & wr);

template
void
downsort<2>(fb_factorbase::slicing const &,
        bucket_array_t<1, longhint_t> &BA_out,
        const bucket_array_t<2, longhint_t> &BA_in,
        uint32_t bucket_number, where_am_I & w);

template
void
downsort<2>(fb_factorbase::slicing const &,
        bucket_array_t<1, logphint_t> &BA_out,
        const bucket_array_t<2, emptyhint_t> &BA_in,
        uint32_t bucket_number, where_am_I & w);

template
void
downsort<3>(fb_factorbase::slicing const &,
        bucket_array_t<2, logphint_t> &BA_out,
        const bucket_array_t<3, emptyhint_t> &BA_in,
        uint32_t bucket_number, where_am_I & wr);

template
void
downsort<2>(fb_factorbase::slicing const &,
        bucket_array_t<1, logphint_t> &BA_out,
        const bucket_array_t<2, logphint_t> &BA_in,
        uint32_t bucket_number, where_am_I & w);

/* Instantiate concrete classes that we need or some methods do not get
   compiled and cause "undefined reference" errors during linking. */
template class bucket_single<1, primehint_t>;
template class bucket_single<1, longhint_t>;
template class bucket_array_t<1, shorthint_t>;
template class bucket_array_t<2, shorthint_t>;
template class bucket_array_t<3, shorthint_t>;
template class bucket_array_t<1, longhint_t>;
template class bucket_array_t<2, longhint_t>;
template class bucket_array_t<1, emptyhint_t>;
template class bucket_array_t<2, emptyhint_t>;
template class bucket_array_t<3, emptyhint_t>;
template class bucket_array_t<1, logphint_t>;
template class bucket_array_t<2, logphint_t>;

