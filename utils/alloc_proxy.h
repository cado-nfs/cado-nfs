#ifndef CADO_ALLOC_PROXY_H
#define CADO_ALLOC_PROXY_H

#ifdef __cplusplus
#include <cstdlib>
#endif

#include "macros.h"

// scan-headers: stop here

/* This header file defines mynew/mynew/mymalloc/myfree, where deletions
 * require the size to be provided -- this makes it easier to plug in
 * memory debuggers.
 */
#define xxxUSE_ELECTRIC_ALLOC

/* announce C prototypes with appropriate linkage. */
#ifdef __cplusplus
extern "C" {
#endif
CADO_INLINE void * mymalloc(size_t s);
CADO_INLINE void myfree(void * p, size_t s);
#ifdef __cplusplus
}
#endif

#ifdef USE_ELECTRIC_ALLOC
#include "electric_alloc.h"
CADO_INLINE void * mymalloc(size_t s) { return electric_alloc(s); }
CADO_INLINE void myfree(void * p, size_t s) { electric_free(p, s); }
#else
CADO_INLINE void * mymalloc(size_t s) { return malloc(s); }
CADO_INLINE void myfree(void * p, size_t s MAYBE_UNUSED) { free(p); }
#endif

#ifdef __cplusplus
#ifdef USE_ELECTRIC_ALLOC
template<typename T>
inline T * mynew(size_t s) { return electric_new<T>(s); }
template<typename T>
inline void mydelete(T * & p, size_t s) { electric_delete(p,s); p = NULL; }
#else
template<typename T>
inline T * mynew(size_t s) { return new T[s]; }
template<typename T>
inline void mydelete(T * & p, size_t s MAYBE_UNUSED) { delete[] p; p = NULL; }
#endif
#endif

#endif	/* CADO_ALLOC_PROXY_H */
