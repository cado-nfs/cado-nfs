#ifndef CADO_UTILS_MEMALLOC_HPP
#define CADO_UTILS_MEMALLOC_HPP

#include <cstddef>
#include "typedefs.h"

#include <memory>
#include <mutex>
#include <list>

/********************** own memory allocation routines ***********************/

/* Rationale: calling one malloc() for each read relation is expensive, since 
   malloc() allocates some extra information to keep track of every memory 
   blocks. Instead, we allocate memory in big blocks of 1MB. */

template<typename T, size_t block_size = (1U << 20)>
struct simple_minded_chunk_allocator {
    T * alloc(size_t n) {
        const std::scoped_lock dummy(m);
        if (used_in_last_block + n > block_size) {
            blocks.emplace_back(new T[block_size]);
            used_in_last_block = 0;
            total_allocated_bytes += block_size * sizeof(T);
        }
        T * p = blocks.back().get() + used_in_last_block;
        used_in_last_block += n;
        return p;
    }
    size_t get_allocated_bytes() const {
        std::scoped_lock dummy(m);
        return total_allocated_bytes;
    }
    private:
    std::list<std::unique_ptr<T[]>> blocks;
    size_t used_in_last_block = block_size;
    size_t total_allocated_bytes = 0;
    mutable std::mutex m;
};

/* these are just proxies to static instances of
 * simple_minded_chunk_allocator
 */
index_t * index_my_malloc (size_t);
ideal_merge_t * ideal_merge_my_malloc (size_t);
size_t get_my_malloc_bytes ();

#endif /* CADO_UTILS_MEMALLOC_HPP */
