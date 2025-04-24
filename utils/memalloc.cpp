#include "cado.h" // IWYU pragma: keep

#include <cstddef>

#include <list>
#include <memory>
#include <mutex>

#include "memalloc.h" // for BLOCK_SIZE
#include "typedefs.h" // for index_t ideal_merge_t

/* memory blocks are allocated of that # of index_t's */
#define BLOCK_SIZE (1U<<20U)

template<typename T, size_t block_size = BLOCK_SIZE>
struct simple_minded_chunk_allocator {
    std::list<std::unique_ptr<T[]>> blocks;
    size_t used_in_last_block = block_size;
    size_t total_allocated_bytes = 0;
    std::mutex m;
    T * alloc(size_t n) {
        const std::lock_guard<std::mutex> dummy(m);
        if (used_in_last_block + n > block_size) {
            blocks.emplace_back(new T[block_size]);
            used_in_last_block = 0;
            total_allocated_bytes += block_size * sizeof(T);
        }
        T * p = blocks.back().get() + used_in_last_block;
        used_in_last_block += n;
        return p;
    }
};

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static simple_minded_chunk_allocator<index_t, BLOCK_SIZE> index_pool;
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static simple_minded_chunk_allocator<ideal_merge_t, BLOCK_SIZE> ideal_merge_pool;

index_t * index_my_malloc(size_t n)
{
    return index_pool.alloc(n);
}

ideal_merge_t * ideal_merge_my_malloc(size_t n)
{
    return ideal_merge_pool.alloc(n);
}

size_t get_my_malloc_bytes()
{
    size_t s = 0;
    {
        const std::lock_guard<std::mutex> dummy(index_pool.m);
        s += index_pool.total_allocated_bytes;
    }
    {
        const std::lock_guard<std::mutex> dummy(ideal_merge_pool.m);
        s += ideal_merge_pool.total_allocated_bytes;
    }
    return s;
}
