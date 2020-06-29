#ifndef CADO_UTILS_ELECTRIC_ALLOC_H_
#define CADO_UTILS_ELECTRIC_ALLOC_H_

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
#include <new>  /* for std::bad_alloc */
#include <stdlib.h>
#endif

static inline
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

static inline
void electric_free(void * p0, size_t s)
{
    char * p = (char *) p0;
    size_t r = 8192;
    unsigned int multip = (s+r-1)/r;
#ifdef PROTECT_OVERRUN
    p += s;
    mprotect((void*) p, r, PROT_READ | PROT_WRITE);
    p -= multip * r;
#else
    p -= r;
    mprotect(p, r, PROT_READ | PROT_WRITE);
#endif
    munmap(p, (multip + 1) * r);
}

static inline
void * electric_alloc_nosize(size_t s)
{
    void * ptr = electric_alloc(s + sizeof(s));
    *(size_t *)ptr = s;
    return (void *) (((size_t *) ptr) + 1);
}

static inline
void electric_free_nosize(void * p0)
{
    p0 = (void*) (((size_t *)p0)-1);
    size_t s = * (size_t *) p0;
    electric_free(p0, s + sizeof(s));
}

#ifdef  __cplusplus
template<typename T> inline T * electric_new(size_t s) {
    T * res = (T*) electric_alloc(sizeof(T)*s);
    if (!res)
        throw std::bad_alloc();
    return res;
}
template<typename T> inline void electric_delete(T * p, size_t s) {
    electric_free(p, sizeof(T)*s);
}
#endif

#endif  /* CADO_UTILS_ELECTRIC_ALLOC_H_ */
