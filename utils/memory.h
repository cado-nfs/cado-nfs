#ifndef CADO_UTILS_MEMORY_H
#define CADO_UTILS_MEMORY_H

#include <stdlib.h>
#ifdef __cplusplus
#include <memory>
#endif

#include "macros.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void * malloc_check(size_t x);
ATTRIBUTE((malloc)) extern void * malloc_aligned(size_t size, size_t alignment);
ATTRIBUTE((warn_unused_result))
void * realloc_aligned(void * p, size_t old_size, size_t new_size,
                       size_t alignment);
extern void free_aligned(void * ptr);

extern void * aligned_alloc(size_t alignment, size_t size);

extern void * malloc_pagealigned(size_t sz) ATTR_ASSUME_ALIGNED(32);
extern void free_pagealigned(void * ptr);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
/* std::allocator<T>::pointer has always been the same as T *, and the
 * former is now deprecated in favor of the latter
 */
template <typename T, int align = sizeof(T)>
class aligned_allocator : public std::allocator<T>
{
    using super = std::allocator<T>;

  public:
    template <typename U> struct rebind {
        using other = aligned_allocator<U, align>;
    };
    T * allocate(size_t n) const
    {
        return static_cast<T *>(malloc_aligned(n * sizeof(T), align));
    }
    void deallocate(T * p, size_t) const { return free_aligned(p); }
    template <typename X> T * allocate(size_t const n, X const *) const
    {
        return allocate(n);
    }
};

struct unique_aligned_array_deleter {
    template <typename T> void operator()(T * p) { free_aligned(p); }
};
template <typename T>
using unique_aligned_array = std::unique_ptr<T[], unique_aligned_array_deleter>;
template <typename T>
inline unique_aligned_array<T> make_unique_aligned_array(size_t size, size_t align)
{
    return unique_aligned_array<T> { static_cast<T *>(malloc_aligned(size * sizeof(T), align)) };
}

#endif

#endif
