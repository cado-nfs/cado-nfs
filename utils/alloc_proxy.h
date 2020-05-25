#include <stdlib.h>
#ifndef ALLOC_PROXY_H_
#define ALLOC_PROXY_H_

#include "macros.h"

/* This header file defines mynew/mynew/mymalloc/myfree, where deletions
 * require the size to be provided -- this makes it easier to plug in
 * memory debuggers.
 */
#define xxxUSE_ELETRIC_ALLOC

/* announce C prototypes with appropriate linkage. */
#ifdef __cplusplus
extern "C" {
#endif
static inline void * mymalloc(size_t s);
static inline void myfree(void * p, size_t s);
#ifdef __cplusplus
}
#endif

#ifdef USE_ELETRIC_ALLOC
#include "electric_alloc.h"
static inline void * mymalloc(size_t s) { return electric_alloc(s); }
static inline void myfree(void * p, size_t s) { electric_free(p, s); }
#else
static inline void * mymalloc(size_t s) { return malloc(s); }
static inline void myfree(void * p, size_t s MAYBE_UNUSED) { free(p); }
#endif

#ifdef __cplusplus
#ifdef USE_ELETRIC_ALLOC
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

#endif	/* ALLOC_PROXY_H_ */
