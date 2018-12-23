#include "cado.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <errno.h>
#include <limits.h>
#include <unistd.h>

#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif

#ifdef HAVE_SYS_MMAN_H
#include <sys/mman.h>
#endif

#include "portability.h"

#ifndef LARGE_PAGE_SIZE
#define LARGE_PAGE_SIZE (2UL*1024*1024)
#endif

#include "las-memory.hpp"

// memory.h provides some back-ends that we use here.
#include "memory.h"

las_memory_accessor::~las_memory_accessor()
{
    std::lock_guard<std::mutex> dummy(bucket_regions_pool.mutex());
    // ASSERT_ALWAYS(bucket_regions_pool.empty());
    for(auto x : large_pages_for_pool) free_aligned(x);
    large_pages_for_pool.clear();
}

unsigned char * las_memory_accessor::alloc_bucket_region()
{
    std::lock_guard<std::mutex> dummy(bucket_regions_pool.mutex());
    if (bucket_regions_pool.empty()) {
        /* allocate some more */
        unsigned char * w = static_cast<unsigned char*>(malloc_aligned(LARGE_PAGE_SIZE, LARGE_PAGE_SIZE));
        ASSERT_ALWAYS(w != 0);
        for(size_t s = 0 ; s + bucket_region_size() <= LARGE_PAGE_SIZE ; s += bucket_region_size())
            bucket_regions_pool.push(w + s);
        large_pages_for_pool.push_back(w);
    }
    unsigned char * v = bucket_regions_pool.top();
    bucket_regions_pool.pop();
    return v;
}

void las_memory_accessor::free_bucket_region(unsigned char * v)
{
    if (!v) return;
    std::lock_guard<std::mutex> dummy(bucket_regions_pool.mutex());
    bucket_regions_pool.push(v);
}

void * las_memory_accessor::physical_alloc(size_t size, bool affect)
{
#if defined(HAVE_MMAP) && defined(MAP_HUGETLB)
    {
        size_t nr_pages = iceildiv(size, LARGE_PAGE_SIZE);
        size_t rounded_up_size = nr_pages * LARGE_PAGE_SIZE; 
        /* Start by trying mmap() */
        void *m = mmap (NULL, rounded_up_size, PROT_READ|PROT_WRITE, MAP_ANONYMOUS | MAP_PRIVATE | MAP_HUGETLB, -1, 0);
        if (m == MAP_FAILED) {
            // Commented out because it's spammy
            // perror("mmap failed");
        } else {
            {
                std::lock_guard<std::mutex> dummy(was_mmapped.mutex());
                was_mmapped.insert(m);
            }
            if (affect) touch(m, size);
            return m;
        }
    }
#endif


#ifdef MADV_HUGEPAGE
    size_t nr_pages = iceildiv(size, LARGE_PAGE_SIZE);
    size_t rounded_up_size = nr_pages * LARGE_PAGE_SIZE; 
    /* If mmap() didn't work, try aligned malloc() with madvise() */
    void *m = malloc_aligned(rounded_up_size, LARGE_PAGE_SIZE);
    int r;
    static int printed_error = 0;
    do {
        r = madvise(m, rounded_up_size, MADV_HUGEPAGE);
    } while (r == EAGAIN);
    if (r != 0 && !printed_error) {
        perror("madvise failed");
        printed_error = 1;
    }
#else
    /* If all else fails, return regular page-aligned memory */
    void * m = malloc_pagealigned(size);
#endif

    if (affect) touch(m, size);
    return m;

}

void las_memory_accessor::touch(void * p, size_t x)
{
    size_t i, m;
#ifdef HAVE_SSE2
    const __m128i a = (__m128i) {0, 0};
#endif    
    i = ((size_t) p + 15) & (~15ULL);
    m = ((size_t) p + x - 1) & (~15ULL);
    while (i < m) {
#ifdef HAVE_SSE2
        _mm_stream_si128((__m128i *)i, a);
#else
        *(unsigned char *) i = 0;
#endif
        i += pagesize ();
    }
}

void las_memory_accessor::physical_free(void * p, size_t size MAYBE_UNUSED)
{
#if defined(HAVE_MMAP) && defined(MAP_HUGETLB)
    {
        size_t nr_pages = iceildiv(size, LARGE_PAGE_SIZE);
        size_t rounded_up_size = nr_pages * LARGE_PAGE_SIZE; 
        std::lock_guard<std::mutex> dummy(was_mmapped.mutex());
        auto it = was_mmapped.find(p);
        if (it != was_mmapped.end()) {
            munmap((void *) p, rounded_up_size);
            was_mmapped.erase(it);
            return;
        }
    }
#endif
    /* otherwise it was simply obtained with malloc (aligned to large
     * page size or simply page size, depending on compile-time support.
     */
#ifdef MADV_HUGEPAGE
    free_aligned(p);
#else
    free_pagealigned(p);
#endif
}
