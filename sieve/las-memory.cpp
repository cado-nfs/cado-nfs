#include "cado.h" // IWYU pragma: keep

#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif
#include <cerrno>
#include <cstdio>
#ifdef HAVE_SYS_MMAN_H
// IWYU pragma: no_include <bits/mman-map-flags-generic.h>
#include <sys/mman.h>
#endif
#include <mutex>
#include <utility>

#include "las-memory.hpp"
#include "memory.h"
#include "misc.h"
#include "verbose.h"
#include "portability.h"
#include "macros.h"

#ifndef LARGE_PAGE_SIZE
#define LARGE_PAGE_SIZE (2UL*1024*1024)
#endif

const size_t small_size_cutoff = 4096;

void * las_memory_accessor::alloc_frequent_size(size_t size)
{
    if (size <= small_size_cutoff)
        return malloc_aligned(size, 128);
    if (size == 0)
        return NULL;
    size_t const rsize = next_power_of_2(size);
    std::lock_guard<std::mutex> const dummy(frequent_regions_pool.mutex());
    auto & pool(frequent_regions_pool[rsize]);
    if (rsize > LARGE_PAGE_SIZE) {
        /* This should hardly ever occur, most probably never, but in
         * ultra weird cases like #30058, this does happen. Well, in such
         * a case let's use plain malloc...
         */
        verbose_fmt_print(1, 1, "# Extraordinarily large (but temporary) allocation"
                " for an area of size {}\n", rsize);
        return malloc_aligned(rsize, LARGE_PAGE_SIZE);
    }
    if (pool.empty()) {
        /* allocate some more */
        verbose_fmt_print(1, 2,
                "# Allocating new large page dedicated to returning"
                " memory areas of size {}\n", rsize);
        auto w = make_unique_aligned_array<char>(LARGE_PAGE_SIZE, LARGE_PAGE_SIZE);
        ASSERT_ALWAYS(w);
        for(size_t s = 0 ; s + rsize <= LARGE_PAGE_SIZE ; s += rsize)
            pool.push(static_cast<void *>(w.get() + s));
        large_pages_for_pool.push_back(std::move(w));
    }
    void * v = pool.top();
    pool.pop();
    return v;
}

void las_memory_accessor::free_frequent_size(void * v, size_t size)
{
    if (size <= small_size_cutoff) {
        free_aligned(v);
        return;
    }
    if (!v) return;
    size_t const rsize = next_power_of_2(size);
    if (rsize > LARGE_PAGE_SIZE) {
        free_aligned(v);
        return;
    }
    std::lock_guard<std::mutex> const dummy(frequent_regions_pool.mutex());
    auto & pool(frequent_regions_pool[rsize]);
    pool.push(v);
}

void * las_memory_accessor::physical_alloc(size_t size, bool affect)
{
#if defined(HAVE_MMAP) && defined(MAP_HUGETLB)
    {
        size_t const nr_pages = iceildiv(size, LARGE_PAGE_SIZE);
        size_t const rounded_up_size = nr_pages * LARGE_PAGE_SIZE; 
        /* Start by trying mmap() */
        void *m = mmap (NULL, rounded_up_size, PROT_READ|PROT_WRITE, MAP_ANONYMOUS | MAP_PRIVATE | MAP_HUGETLB, -1, 0);
        if (m == MAP_FAILED) {
            // Commented out because it's spammy
            // perror("mmap failed");
        } else {
            {
                std::lock_guard<std::mutex> const dummy(was_mmapped.mutex());
                was_mmapped.insert(m);
            }
            if (affect) touch(m, size);
            return m;
        }
    }
#endif


#ifdef MADV_HUGEPAGE
    size_t const nr_pages = iceildiv(size, LARGE_PAGE_SIZE);
    size_t const rounded_up_size = nr_pages * LARGE_PAGE_SIZE; 
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
        size_t const nr_pages = iceildiv(size, LARGE_PAGE_SIZE);
        size_t const rounded_up_size = nr_pages * LARGE_PAGE_SIZE; 
        std::lock_guard<std::mutex> const dummy(was_mmapped.mutex());
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
