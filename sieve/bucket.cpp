#include "cado.h"
#include <inttypes.h>
#include <stdlib.h>   // for malloc and friends
#include <string.h>   // for memcpy
#include <new>        // for std::bad_alloc
#include <gmp.h>
#if defined(HAVE_SSE2)
#include <emmintrin.h>
#endif
#include "las-types.hpp"
#include "bucket.hpp"
#include "portability.h"
#include "memory.h"
#include "las-config.h"
#include "las-debug.hpp"
#include "iqsort.h"
#include "verbose.h"
#include "ularith.h"
#include "smallset.hpp"

/* sz is the size of a bucket for an array of buckets. In bytes, a bucket
   size is sz * sr, with sr = sizeof of one element of the bucket (a record).
   If sz * sr is a multiple of the page, and if the buckets array is filled
   *UNIFORMLY & REGULARLY*, all the buckets in the array have the possibility
   to reach their end of the memory pages at the same time.
   In this case, the system must manage multiple TLB misses & must provide
   multiple physical pages of memory in same time = really a bad idea!
   So, in this case, I increase lightly sz.
   NB: I use in fact sz*sr*4 and not sz*sr, because an half or a quart
   of the buckets which reach at the same time the end of their page is
   sufficient to slow down the code.
*/
static size_t
bucket_misalignment(const size_t sz, const size_t sr) {
  size_t size; 
  if ((sz * sr * 4) & (pagesize() - 1))
    size = sz;
  else 
    size = sz + 8 / sr + ((8 % sr) ? 1 : 0);
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
bucket_array_t<LEVEL, HINT>::bucket_array_t()
  : big_data(NULL), big_size(0), bucket_write(NULL), bucket_start(NULL),
  bucket_read(NULL), slice_index(NULL), slice_start(NULL), n_bucket(0),
  bucket_size(0), size_b_align(0), nr_slices(0), alloc_slices(0)
{
}

template <int LEVEL, typename HINT>
bucket_array_t<LEVEL, HINT>::~bucket_array_t()
{
  physical_free (big_data, big_size);
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
  MOVE_ENTRY(bucket_size, 0);
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
bucket_array_t<LEVEL, HINT>::allocate_memory(const uint32_t new_n_bucket,
                                const double fill_ratio,
                                const slice_index_t prealloc_slices)
{
  /* Don't try to allocate anything, nor print a message, for sieving levels
     where the corresponding factor base part is empty. */
  if (fill_ratio == 0.)
    return;

  const size_t min_bucket_size = fill_ratio * BUCKET_REGIONS[LEVEL];
  const size_t new_bucket_size = bucket_misalignment(min_bucket_size, sizeof(update_t));
  size_t new_big_size = new_bucket_size * new_n_bucket * sizeof(update_t);
  /* add 1 megabyte to each bucket array, so as to deal with overflowing
   * buckets */
  new_big_size += 1 << 20;
  const size_t new_size_b_align = ((sizeof(void *) * new_n_bucket + 0x3F) & ~((size_t) 0x3F));

  if (new_big_size > big_size) {
    if (big_data != NULL)
      physical_free (big_data, big_size);
    verbose_output_print(0, 3, "# Allocating %zu bytes for %" PRIu32 " buckets of %zu update entries of %zu bytes each\n",
                         new_big_size, new_n_bucket, new_bucket_size, sizeof(update_t));
    big_size = new_big_size;
    big_data = (update_t *) physical_malloc (big_size, 1);
    void * internet_of_things MAYBE_UNUSED = NULL;
  }

  if (!big_data)
      throw std::bad_alloc();

  bucket_size = new_bucket_size;
  n_bucket = new_n_bucket;

  if (new_size_b_align > size_b_align) {
    size_b_align = new_size_b_align;
    bucket_write = (update_t **) malloc_pagealigned (size_b_align);
    bucket_start = (update_t **) malloc_aligned (size_b_align, 0x40);
    bucket_read = (update_t **) malloc_aligned (size_b_align, 0x40);
  }

  /* This requires size_b_align to have been set to the new value */
  if (prealloc_slices > alloc_slices)
    realloc_slice_start(prealloc_slices - alloc_slices);

  /* Spread bucket_start pointers equidistantly over the big_data array */
  for (uint32_t i_bucket = 0; i_bucket < n_bucket; i_bucket++) {
    bucket_start[i_bucket] = big_data + i_bucket * bucket_size;
  }
  reset_pointers();
#ifdef SAFE_BUCKETS
  verbose_output_print(0, 0, "# WARNING: SAFE_BUCKETS is on !\n");
#endif
}

template <int LEVEL, typename HINT>
void
bucket_array_t<LEVEL, HINT>::realloc_slice_start(const size_t extra_space)
{
  const size_t new_alloc_slices = alloc_slices + extra_space;
  verbose_output_print(0, 3, "# Reallocating BA->slice_start from %zu entries to %zu entries\n",
                       alloc_slices, new_alloc_slices);

  const size_t old_size = size_b_align * alloc_slices;
  const size_t new_size = size_b_align * new_alloc_slices;
  slice_start = (update_t **) realloc_aligned(slice_start, old_size, new_size, 0x40);
  ASSERT_ALWAYS(slice_start != NULL);
  slice_index = (slice_index_t *) realloc(slice_index, new_alloc_slices * sizeof(slice_index_t));
  ASSERT_ALWAYS(slice_index != NULL);
  alloc_slices = new_alloc_slices;
}

/* Returns how full the fullest bucket is, as a fraction of its size */
template <int LEVEL, typename HINT>
double
bucket_array_t<LEVEL, HINT>::max_full (unsigned int * fullest_index) const
{
  unsigned int max = 0;
  for (unsigned int i = 0; i < n_bucket; ++i)
    {
      unsigned int j = nb_of_updates (i);
      if (max < j) {
          max = j;
          if (fullest_index) *fullest_index = i;
      }
    }
  return (double) max / (double) bucket_size;
}

// Replace std::is_same()::value that is not yet available.
template <typename T, typename U>
struct is_same
{
    static const bool value = false;
};

template <typename T>
struct is_same<T, T> { static const bool value = true; };


template <int LEVEL, typename HINT>
void
bucket_array_t<LEVEL, HINT>::log_this_update (
    const update_t update MAYBE_UNUSED,
    const uint64_t offset MAYBE_UNUSED,
    const uint64_t bucket_number MAYBE_UNUSED,
    where_am_I & w MAYBE_UNUSED) const
{
#if defined(TRACE_K)
    size_t (&BRS)[FB_MAX_PARTS] = BUCKET_REGIONS;
    unsigned int saveN = w.N;
    /* flatten the (N,x) coordinate as if relative to a unique array of
     * level-1 bucket regions */
    unsigned int x = update.x % BRS[1];
    unsigned int N = w.N +
        bucket_number*BRS[LEVEL]/BRS[1] + (update.x / BRS[1]);

    WHERE_AM_I_UPDATE(w, x, x);
    WHERE_AM_I_UPDATE(w, N, N);

    if (trace_on_spot_Nx(w.N, w.x)) {
        verbose_output_print (TRACE_CHANNEL, 0,
            "# Pushed hit at location (x=%u, side %d), from factor base entry "
            "(slice_index=%u, slice_offset=%u, p=%" FBPRIME_FORMAT "), "
            "to BA<%d>[%u]\n",
            (unsigned int) w.x, w.side, (unsigned int) w.i,
            (unsigned int) w.h, w.p, LEVEL, (unsigned int) w.N);
        if (is_same<HINT,longhint_t>::value) {
          verbose_output_print (TRACE_CHANNEL, 0,
             "# Warning: did not check divisibility during downsorting p=%"
             FBPRIME_FORMAT "\n", w.p);
        } else {
          ASSERT_ALWAYS(test_divisible(w));
        }
    }
    WHERE_AM_I_UPDATE(w, N, saveN);
#endif
}



/* Instantiate concrete classes that we need or some methods do not get
   compiled and cause "undefined reference" errors during linking. */
template class bucket_array_t<1, shorthint_t>;
template class bucket_array_t<2, shorthint_t>;
template class bucket_array_t<3, shorthint_t>;
template class bucket_array_t<1, longhint_t>;
template class bucket_array_t<2, longhint_t>;


/* A compare function suitable for sorting updates in order of ascending x
   with qsort() */
int
bucket_cmp_x (const bucket_update_t<1, longhint_t>  *a, const bucket_update_t<1, longhint_t>  *b)
{
  if (a->x < b->x)
    return -1;
  if (a->x == b->x)
    return 0;
  return 1;
}

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
              const int i, const fb_part *fb, const unsigned char *S)
{
  for (slice_index_t i_slice = 0; i_slice < BA.get_nr_slices(); i_slice++) {
    const slice_index_t slice_index = BA.get_slice_index(i_slice);
    const bucket_update_t<1, shorthint_t> *it = BA.begin(i, i_slice);
    const bucket_update_t<1, shorthint_t> * const end_it = BA.end(i, i_slice);

    for ( ; it != end_it ; it++) {
      if (UNLIKELY(S[it->x] != 255)) {
        const fb_slice_interface *slice = fb->get_slice(slice_index);
        ASSERT_ALWAYS(slice != NULL);
        fbprime_t p = slice->get_prime(it->hint);
        push_update(bucket_update_t<1, primehint_t>(it->x, p, 0, 0));
      }
    }
  }
}

template <typename HINT>
static inline bucket_update_t<1, longhint_t>
to_longhint(const bucket_update_t<1, HINT> &update, slice_index_t slice_index);

template<>
bucket_update_t<1, longhint_t>
to_longhint<shorthint_t>(const bucket_update_t<1, shorthint_t> &update,
                         const slice_index_t slice_index)
{
  return bucket_update_t<1, longhint_t> (update.x, 0, update.hint, slice_index);
}

template<>
bucket_update_t<1, longhint_t>
to_longhint<longhint_t>(const bucket_update_t<1, longhint_t> &update,
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
    const bucket_update_t<1, HINT> *it = BA.begin(i, i_slice);
    const bucket_update_t<1, HINT> * const end_it = BA.end(i, i_slice);

    for ( ; it != end_it ; it++) {
      if (UNLIKELY(S[it->x] != 255)) {
        push_update(to_longhint(*it, slice_index));
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
    const bucket_update_t<1, HINT> *it = BA.begin(i, i_slice);
    const bucket_update_t<1, HINT> * const end_it = BA.end(i, i_slice);

    for ( ; it != end_it ; it++) {
      if (UNLIKELY(surv_set.contains(it->x))) {
        push_update(to_longhint(*it, slice_index));
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

template class bucket_single<1, primehint_t>;
template class bucket_single<1, longhint_t>;

template<int INPUT_LEVEL>
void
downsort(bucket_array_t<INPUT_LEVEL - 1, longhint_t> &BA_out,
         const bucket_array_t<INPUT_LEVEL, shorthint_t> &BA_in,
         uint32_t bucket_number, where_am_I & w)
{
  /* Rather similar to purging, except it doesn't purge */
  for (slice_index_t i_slice = 0; i_slice < BA_in.get_nr_slices(); i_slice++) {
    const slice_index_t slice_index = BA_in.get_slice_index(i_slice);
    WHERE_AM_I_UPDATE(w, i, slice_index);
    const bucket_update_t<INPUT_LEVEL, shorthint_t> *it = BA_in.begin(bucket_number, i_slice);
    const bucket_update_t<INPUT_LEVEL, shorthint_t> * const end_it = BA_in.end(bucket_number, i_slice);

    for ( ; it != end_it ; it++) {
      WHERE_AM_I_UPDATE(w, p,
          w.psi->sides[w.side].fb->get_slice(slice_index)->get_prime(it->hint));
      BA_out.push_update(it->x, 0, it->hint, slice_index, w);
    }
  }
}

template<int INPUT_LEVEL>
void
downsort(bucket_array_t<INPUT_LEVEL - 1, longhint_t> &BA_out,
         const bucket_array_t<INPUT_LEVEL, longhint_t> &BA_in,
         uint32_t bucket_number, where_am_I & w) 
{
  /* longhint updates don't write slice end pointers, so there must be
     exactly 1 slice per bucket */
  ASSERT_ALWAYS(BA_in.get_nr_slices() == 1);
  const bucket_update_t<INPUT_LEVEL, longhint_t> *it = BA_in.begin(bucket_number, 0);
  const bucket_update_t<INPUT_LEVEL, longhint_t> * const end_it = BA_in.end(bucket_number, 0);

  for ( ; it != end_it ; it++) {
    BA_out.push_update(it->x, 0, it->hint, it->index, w);
  }
}

/* Explicitly instantiate the versions of downsort() that we'll need:
   downsorting shorthint from level 3 and level 2, and downsorting
   longhint from level 2. */
template
void
downsort<2>(bucket_array_t<1, longhint_t> &BA_out,
            const bucket_array_t<2, shorthint_t> &BA_in,
            uint32_t bucket_number, where_am_I & w);

template
void
downsort<3>(bucket_array_t<2, longhint_t> &BA_out,
            const bucket_array_t<3, shorthint_t> &BA_in,
            uint32_t bucket_number, where_am_I & wr);

template
void
downsort<2>(bucket_array_t<1, longhint_t> &BA_out,
            const bucket_array_t<2, longhint_t> &BA_in,
            uint32_t bucket_number, where_am_I & w);

void
sieve_checksum::update(const unsigned int other)
{
    unsigned long r;
    ularith_addmod_ul_ul(&r, checksum, other, checksum_prime);
    checksum = r;
}

void
sieve_checksum::update(const unsigned char *data, const size_t len)
{
    mpz_t mb;
    unsigned int new_checksum;

    mpz_init(mb);
    mpz_import(mb, len, -1, sizeof(unsigned char), -1, 0, data);
    new_checksum = mpz_tdiv_ui(mb, checksum_prime);
    mpz_clear(mb);
    this->update(new_checksum);
}
