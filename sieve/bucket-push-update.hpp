#include <stdio.h>
#ifndef BUCKET_PUSH_UPDATE_HPP_
#define BUCKET_PUSH_UPDATE_HPP_

#include "bucket.hpp"
#include "las-where-am-i.hpp"       // WHERE_AM_I_UPDATE
#include "macros.h"

template<int LEVEL, typename HINT>
inline void
bucket_array_t<LEVEL, HINT>::push_update(const int i, const update_t& update)
{
#ifdef SAFE_BUCKET_ARRAYS
    if (bucket_write[i] >= bucket_start[i + 1]) {
        fprintf(stderr, "# Warning: hit end of bucket nb %d\n", i);
        ASSERT_ALWAYS(0);
        return;
    }
#endif
    *bucket_write[i]++ = update;
}

template<int LEVEL, typename HINT>
inline void
bucket_single<LEVEL, HINT>::push_update(const update_t& update)
{
#ifdef SAFE_BUCKETS_SINGLE
    if (start + _size <= write) {
        fprintf(stderr, "# Warning: hit end of bucket\n");
        ASSERT_ALWAYS(0);
        write--;
    }
#endif
    *(write++) = update;
}

template<int LEVEL, typename HINT>
inline const typename bucket_single<LEVEL, HINT>::update_t&
bucket_single<LEVEL, HINT>::get_next_update()
{
#ifdef SAFE_BUCKETS_SINGLE
    ASSERT_ALWAYS(read < write);
#endif
    return *read++;
}

template<int LEVEL, typename HINT>
inline void
bucket_array_t<LEVEL, HINT>::push_update(
       const uint64_t offset, const fbprime_t p,
       const slice_offset_t slice_offset, const slice_index_t slice_index,
      where_am_I& w MAYBE_UNUSED)
  {
      int logB = LOG_BUCKET_REGIONS[LEVEL];
      const uint64_t bucket_number = offset >> logB;
      ASSERT_EXPENSIVE(bucket_number < n_bucket);
      update_t update(offset & ((UINT64_C(1) << logB) - 1), p, slice_offset, slice_index);
      WHERE_AM_I_UPDATE(w, i, slice_index);
#if defined(TRACE_K)
      log_this_update(update, offset, bucket_number, w);
#endif
      push_update(bucket_number, update);
  }


#endif	/* BUCKET_PUSH_UPDATE_HPP_ */
