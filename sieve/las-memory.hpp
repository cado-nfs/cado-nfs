#ifndef CADO_LAS_MEMORY_HPP
#define CADO_LAS_MEMORY_HPP

#include <cstddef>

#include <list>
#include <map>
#include <set>
#include <stack>
#include <memory>

#include "macros.h"
#include "las-config.hpp"
#include "lock_guarded_container.hpp"
#include "memory.h"


/* This structure is shared by threads that have the same memory binding.
 * It is in charge of providing memory-bound zones, e.g. for buckets, or
 * bucket regions.
 *
 * alloc()/free() calls are not meant to be fast, here.
 */

class las_memory_accessor {
    lock_guarded_container<std::map<size_t, std::stack<void *>>> frequent_regions_pool;
    std::list<unique_aligned_array<char>> large_pages_for_pool;

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

    static size_t bucket_region_size() {
        /* round to next multiple of 128 */
        return (((BUCKET_REGION + MEMSET_MIN) - 1) | 127) + 1;
    }
    void * alloc_frequent_size(size_t);
    void free_frequent_size(void *, size_t);
    unsigned char * alloc_bucket_region() { return (unsigned char *) alloc_frequent_size(bucket_region_size()); }
    void free_bucket_region(unsigned char * p) { free_frequent_size((void *) p, bucket_region_size()); }

    void * physical_alloc(size_t, bool = false) ATTR_ASSUME_ALIGNED(256);
    void physical_free(void*, size_t);

    struct unique_frequent_array_deleter {
        size_t size = 0;
        las_memory_accessor * a = nullptr;
        template<typename T>
        void operator()(T * p) const {
            a->free_frequent_size(p, size * sizeof(T));
        }
    };
    template<typename T>
        using unique_frequent_array = std::unique_ptr<T, unique_frequent_array_deleter>;

    template<typename T>
    unique_frequent_array<T> make_unique_frequent_array(size_t size)
    {
        return { static_cast<T *>(alloc_frequent_size(size * sizeof(T))),
                 { size, this } };
    }

    struct unique_physical_array_deleter {
        size_t size = 0;
        las_memory_accessor * a = nullptr;
        template<typename T>
        void operator()(T * p) const {
            a->physical_free(p, size * sizeof(T));
        }
    };
    template<typename T>
        using unique_physical_array = std::unique_ptr<T, unique_physical_array_deleter>;

    template<typename T>
    unique_physical_array<T> make_unique_physical_array(size_t size, bool affect = false)
    {
        return { static_cast<T *>(physical_alloc(size * sizeof(T), affect)),
                 { size, this } };
    }


    las_memory_accessor() = default;
    las_memory_accessor(las_memory_accessor const &) = delete;
    las_memory_accessor(las_memory_accessor&&) = default;
    las_memory_accessor& operator=(las_memory_accessor const &) = delete;
    las_memory_accessor& operator=(las_memory_accessor&&) = default;
    ~las_memory_accessor() = default;
};

#endif	/* CADO_LAS_MEMORY_HPP */
