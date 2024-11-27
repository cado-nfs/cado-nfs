#ifndef BUCKET_PUSH_UPDATE_HPP_
#define BUCKET_PUSH_UPDATE_HPP_

#ifdef SAFE_BUCKET_ARRAYS
#include <cstdio>
#endif

#include "bucket.hpp"
#include "las-where-am-i.hpp"       // WHERE_AM_I_UPDATE
#include "macros.h"

template<int LEVEL, typename HINT>
inline void
bucket_array_t<LEVEL, HINT>::push_update(const int i, const update_t& update, where_am_I& w MAYBE_UNUSED)
{
#ifdef SAFE_BUCKET_ARRAYS
    if (bucket_write[i] >= bucket_start[i + 1]) {
        fprintf(stderr, "# Warning: hit end of bucket nb %d\n", i);
        ASSERT_ALWAYS(0);
        return;
    }
#endif
#if defined(TRACE_K)
    log_this_update(update, i, w);
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

#endif	/* BUCKET_PUSH_UPDATE_HPP_ */
