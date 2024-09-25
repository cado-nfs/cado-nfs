#include "cado.h" // IWYU pragma: keep
// IWYU pragma: no_include <mm_malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#ifdef HAVE_SYS_MMAN_H
#include <sys/mman.h>
#endif
/* For MinGW Build */
#if defined(_WIN32) || defined(_WIN64)
#include <windows.h>
#endif
#include "macros.h"
#include "memory.h"
#include "portability.h"

#ifndef LARGE_PAGE_SIZE
#define LARGE_PAGE_SIZE (2UL*1024*1024)
#endif

void
*malloc_check (const size_t x)
{
    void *p;
    p = malloc (x);
    if (p == NULL)
      {
        fprintf (stderr, "Error, malloc of %zu bytes failed\n", x);
        fflush (stderr);
        abort ();
      }
    return p;
}

#ifndef MAP_ANONYMOUS
#define MAP_ANONYMOUS MAP_ANON
#endif

/* Not everybody has posix_memalign. In order to provide a viable
 * alternative, we need an ``aligned free'' matching the ``aligned
 * malloc''. We rely on posix_memalign if it is available, or else fall
 * back on ugly pointer arithmetic so as to guarantee alignment. Note
 * that not providing the requested alignment can have some troublesome
 * consequences. At best, a performance hit, at worst a segv (sse-2
 * movdqa on a pentium4 causes a GPE if improperly aligned).
 */

void *malloc_aligned(size_t size, size_t alignment)
{
#ifdef HAVE_POSIX_MEMALIGN
    void *res = NULL;
    int rc = posix_memalign(&res, alignment, size);
    // ASSERT_ALWAYS(rc == 0);
    DIE_ERRNO_DIAG(rc != 0, "malloc_aligned(%s)", "");
    return res;
#else
    char * res;
    res = malloc(size + sizeof(size_t) + alignment);
    res += sizeof(size_t);
    size_t displ = alignment - ((uintptr_t) res) % alignment;
    res += displ;
    memcpy(res - sizeof(size_t), &displ, sizeof(size_t));
    ASSERT_ALWAYS((((uintptr_t) res) % alignment) == 0);
    return (void*) res;
#endif
}

/* Reallocate aligned memory.
   p must have been allocated via malloc_aligned() or realloc_aligned().
   old_size must be equal to the size parameter of malloc_aligned(), or
   to the new_size parameter of realloc_aligned(), of the function that
   allocated p. */

void *
realloc_aligned(void * p, const size_t old_size, const size_t new_size,
                const size_t alignment)
{
  if (p == NULL)
    return malloc_aligned(new_size, alignment);
#ifdef HAVE_POSIX_MEMALIGN
  /*  Alas, there is no posix_realloc_aligned(). Try to realloc(); if it
      happens to result in the desired alignment, there is nothing left
      to do. If it does not result in the desired alignment, then we
      actually do two data copies: one as part of realloc(), and another
      below. Let's hope this happens kinda rarely. */
  p = realloc(p, new_size);
  if (((uintptr_t) p) % alignment == 0) {
    return p;
  }
#else
  /* Without posix_memalign(), we always alloc/copy/free */
#endif
  /* We did not get the desired alignment, or we don't have posix_memalign().
     Allocate new memory with the desired alignment and copy the data */
  void * const alloc_p = malloc_aligned(new_size, alignment);
  memcpy(alloc_p, p, MIN(old_size, new_size));
  /* If we have posix_memalign(), then p was allocated by realloc() and can be
     freed with free(), which is what free_aligned() does. If we don't have
     posix_memalign(), then p was allocated by malloc_aligned() or
     realloc_aligned(), so using free_aligned() is correct again. */
  free_aligned(p);
  return alloc_p;
}


void free_aligned(void * p)
{
#ifdef HAVE_POSIX_MEMALIGN
    free((void *) p);
#else
    if (p == NULL)
      return;
    const char * res = (const char *) p;
    size_t displ;
    memcpy(&displ, res - sizeof(size_t), sizeof(size_t));
    res -= displ;
    res -= sizeof(size_t);
    free((void *)res);
#endif
}

/* if aligned_alloc is not provided, define it */
#ifndef HAVE_ALIGNED_ALLOC
void*
aligned_alloc (size_t alignment, size_t size)
{
  return malloc_aligned (size, alignment);
}
#endif

void *malloc_pagealigned(size_t sz)
{
    long ps = pagesize();
    ASSERT_ALWAYS(ps > 0);
    void *p = malloc_aligned (sz, ps);
    ASSERT_ALWAYS(p != NULL);
    return p;
}

void free_pagealigned(void * p)
{
    free_aligned(p);
}
