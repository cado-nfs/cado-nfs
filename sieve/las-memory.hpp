#ifndef LAS_MEMORY_HPP_
#define LAS_MEMORY_HPP_

#include <stddef.h>                    // for size_t
#include <list>                        // for list
#include <map>                         // for map
#include <set>                         // for set
#include <stack>                       // for stack
#include "macros.h"                    // for ATTR_ASSUME_ALIGNED
#include "las-config.h"                // for BUCKET_REGION, MEMSET_MIN
#include "lock_guarded_container.hpp"  // for lock_guarded_container


/* This structure is shared by threads that have the same memory binding.
 * It is in charge of providing momory-bound zones, e.g. for buckets, or
 * bucket regions.
 *
 * alloc()/free() calls are not meant to be fast, here.
 */

class las_memory_accessor {
    lock_guarded_container<std::map<size_t, std::stack<void *>>> frequent_regions_pool;
    std::list<void*> large_pages_for_pool;

    /* large memory chunks follow the same logic as in utils/memory.c,
     * but we reimplement it here so as to stick to one memory binding
     * only.
     *
     * A large memory area may be returned as follows, in decreasing
     * priority order:
     *
     *  - if support is available, via mmap(...,MAP_HUGETLB). If it
     *  succeeds, we get a memory are in multiples of 2G, and call that
     *  an "mmapped region".
     *
     *  - if support is available, via malloc_aligned() +
     *  madvise(MADV_HUGEPAGE). If it succeeds, we call that a "malloced
     *  region"
     *
     *  - otherwise, via malloc_aligned(), and it is still called a
     *  "default malloced region".
     *
     * How the requested size is rounded up depends on the mode.
     */
    lock_guarded_container<std::set<void*>> was_mmapped; // used on free()

    void touch(void *, size_t);
    public:

    static inline size_t bucket_region_size() {
        /* round to next multiple of 128 */
        return (((BUCKET_REGION + MEMSET_MIN) - 1) | 127) + 1;
    }
    void * alloc_frequent_size(size_t);
    void free_frequent_size(void *, size_t);
    inline unsigned char * alloc_bucket_region() { return (unsigned char *) alloc_frequent_size(bucket_region_size()); }
    inline void free_bucket_region(unsigned char * p) { free_frequent_size((void *) p, bucket_region_size()); }

    void * physical_alloc(size_t, bool = false) ATTR_ASSUME_ALIGNED(256);
    void physical_free(void*, size_t);

    las_memory_accessor() = default;
    las_memory_accessor(las_memory_accessor const &) = delete;
    las_memory_accessor(las_memory_accessor&&) = default;
    ~las_memory_accessor();
};

#endif	/* LAS_MEMORY_HPP_ */
