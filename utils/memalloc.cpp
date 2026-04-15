#include "cado.h" // IWYU pragma: keep

#include <cstddef>

#include "memalloc.hpp" // for BLOCK_SIZE
#include "typedefs.h" // for index_t ideal_merge_t

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static simple_minded_chunk_allocator<index_t> index_pool;
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static simple_minded_chunk_allocator<ideal_merge_t> ideal_merge_pool;

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
    return index_pool.get_allocated_bytes() + ideal_merge_pool.get_allocated_bytes();
}
