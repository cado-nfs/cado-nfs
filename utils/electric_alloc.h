#ifndef CADO_UTILS_ELECTRIC_ALLOC_H
#define CADO_UTILS_ELECTRIC_ALLOC_H
#include "macros.h"

/* This file is a debugging aid. It carries good chances of working on
 * POSIX system, but I wouldn't bet much on it, since mmap is kind of
 * strongly tied to the OS.
 *
 * To use, simply include this .h file, and use electric_alloc and
 * electric_free for for allocation/free routines.
 *
 * The vanilla electric_free needs the size of the allocated area. If
 * this is an inconvenient, try the _nosize versions below. Never tested.
 */

/* By default we protect overruns. Undefine this macro to protect
 * underruns instead */
#define PROTECT_OVERRUN

#include <fcntl.h>
#include <sys/mman.h>

#ifndef __APPLE__
#ifndef MAP_ANONYMOUS
#error "Please define _GNU_SOURCE or _BSD_SOURCE on top of the translation unit"
#endif
#endif

#ifdef __cplusplus
#include <cstdlib>

#include <new>  /* for std::bad_alloc */
#endif

CADO_INLINE
void * electric_alloc(size_t s)
{
    /* Use the method of the good old days from electric fence. */
    char *p;
    size_t r = 8192;        /* Any multiple of the page size will do. */
    unsigned int multip = (s+r-1)/r;
    p = (char *)
#ifndef __APPLE__
    mmap(0, (multip + 1) * r, PROT_READ | PROT_WRITE,
                MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
#else
    mmap(0, (multip + 1) * r, PROT_READ | PROT_WRITE,
                MAP_PRIVATE | MAP_ANON, -1, 0);
#endif
    // could please valgrind ?
    // memset(p, 0, (multip + 1) * r);
#ifdef PROTECT_OVERRUN
    p += (multip + 1) * r;
    mprotect((void*) (p-r), r, PROT_NONE);
    p -= r+s;
#else
    /* protect underrun */
    mprotect(p, r, PROT_NONE);
    p += r;
#endif
    return (void *) p;
}

CADO_INLINE
void electric_free(void * p0, size_t s)
{
    char * p = (char *) p0;
    size_t r = 8192;
    unsigned int multip = (s+r-1)/r;
#ifdef PROTECT_OVERRUN
    p += s;
    mprotect((void*) p, r, PROT_READ | PROT_WRITE);
    p -= multip * r;
    // XXX new!
    mprotect((void*) p, multip * r, PROT_READ | PROT_WRITE);
#else
    // XXX new!
    mprotect((void*) p, multip * r, PROT_READ | PROT_WRITE);
    p -= r;
    mprotect(p, r, PROT_READ | PROT_WRITE);
#endif
    munmap(p, (multip + 1) * r);
}

CADO_INLINE
void electric_mark_readonly(void * p0, size_t s)
{
    char * p = (char *) p0; // W: C-style casts are discouraged; use static_cast (fix available)
    size_t r = 8192; // W: variable 'r' of type 'size_t' (aka 'unsigned long') can be declared 'co…
    unsigned int multip = (s+r-1)/r; // W: variable 'multip' of type 'unsigned int' can be declare…
#ifdef PROTECT_OVERRUN
    p += s;
    mprotect((void*) (p - multip * r), multip * r, PROT_READ); // W: C-style casts are discouraged; use stati…
#else
    mprotect(p, multip * r, PROT_READ);
#endif
}

CADO_INLINE
void electric_mark_readwrite(void * p0, size_t s)
{
    char * p = (char *) p0; // W: C-style casts are discouraged; use static_cast (fix available)
    size_t r = 8192; // W: variable 'r' of type 'size_t' (aka 'unsigned long') can be declared 'co…
    unsigned int multip = (s+r-1)/r; // W: variable 'multip' of type 'unsigned int' can be declare…
#ifdef PROTECT_OVERRUN
    p += s;
    mprotect((void*) (p - multip * r), multip * r, PROT_READ | PROT_WRITE); // W: C-style casts are discouraged; use stati…
#else
    mprotect(p, multip * r, PROT_READ | PROT_WRITE);
#endif
}





CADO_INLINE
void * electric_alloc_nosize(size_t s)
{
    void * ptr = electric_alloc(s + sizeof(s));
    *(size_t *)ptr = s;
    return (void *) (((size_t *) ptr) + 1);
}

CADO_INLINE
void electric_free_nosize(void * p0)
{
    p0 = (void*) (((size_t *)p0)-1);
    size_t s = * (size_t *) p0;
    electric_free(p0, s + sizeof(s));
}

#ifdef  __cplusplus
#include <memory>
#include <new>
#include <type_traits>

namespace cado {
template<typename T>
struct ElectricAllocator {
    using value_type = T;

    [[nodiscard]] T* allocate(std::size_t n) // W: no header providing "std::size_t" is…
    {
        auto * p = static_cast<T*>(electric_alloc(n * sizeof(T)));
        if (p == nullptr)
            throw std::bad_alloc();
        return p;
    }

    void deallocate(T* p, std::size_t n) noexcept
    {
        electric_free(p, n * sizeof(T));
    }

};

template<typename T>
struct electric_deleter {
    void operator()(T* ptr) const {
        if (ptr) {
            ptr->~T();
            electric_free(ptr, sizeof(T));
        }
    }
};

template<typename T>
struct electric_deleter<T[]> {
    size_t size;

    explicit electric_deleter(size_t s = 0) : size(s) {}

    void operator()(T* ptr) const {
        if (ptr) {
            for (size_t i = size; i > 0; --i)
                ptr[i - 1].~T();
            electric_free(ptr, size * sizeof(T));
        }
    }
};

template<typename T>
using electric_unique_ptr = std::unique_ptr<T, electric_deleter<T>>;


template<typename T, typename... Args>
T* electric_new(Args&&... args) {
    void* p = electric_alloc(sizeof(T));
    if (!p) throw std::bad_alloc();
    try {
        return new (p) T(std::forward<Args>(args)...);
    } catch (...) {
        electric_free(p, sizeof(T));
        throw;
    }
}

template<typename T>
T* electric_new_array(size_t size) {
    if (size == 0) return nullptr;
    
    void* p = electric_alloc(sizeof(T) * size);
    if (!p) throw std::bad_alloc();
    
    T* arr = static_cast<T*>(p);
    size_t i = 0;
    try {
        for (; i < size; ++i) {
            new (&arr[i]) T();
        }
    } catch (...) {
        // Rollback on throw
        while (i > 0)
            arr[--i].~T();
        electric_free(p, sizeof(T) * size);
        throw;
    }
    return arr;
}

template<typename T, typename... Args>
    requires (!std::is_array_v<T>)
electric_unique_ptr<T> electric_make_unique(Args&&... args) {
    return { electric_new<T>(std::forward<Args>(args)...) };
}

template<typename T>
    requires std::is_unbounded_array_v<T>
electric_unique_ptr<T> electric_make_unique(size_t size) {
    using Elem = std::remove_extent_t<T>;
    return { electric_new_array<Elem>(size), electric_deleter<T>(size) };
}

template<typename T, typename... Args>
    requires std::is_bounded_array_v<T>
void electric_make_unique(Args&&...) = delete;

template<typename T>
    requires std::is_unbounded_array_v<T>
void electric_mark_readonly(electric_unique_ptr<T> & u) {
    auto size = u.get_deleter().size;
    using Elem = std::remove_extent_t<T>;
    ::electric_mark_readonly(u.get(), size * sizeof(Elem));
}
} // namespace cado


#endif

#endif  /* CADO_UTILS_ELECTRIC_ALLOC_H */
