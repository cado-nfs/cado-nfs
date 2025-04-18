#ifndef CADO_UTILS_MEMALLOC_H
#define CADO_UTILS_MEMALLOC_H

#include <stddef.h>     // size_t
#include "typedefs.h"

/********************** own memory allocation routines ***********************/

/* Rationale: calling one malloc() for each read relation is expensive, since 
   malloc() allocates some extra information to keep track of every memory 
   blocks. Instead, we allocate memory in big blocks of size BLOCK_SIZE. */


#ifdef __cplusplus
extern "C" {
#endif

index_t * index_my_malloc (size_t);
ideal_merge_t * ideal_merge_my_malloc (size_t);
size_t get_my_malloc_bytes ();

#ifdef __cplusplus
}
#endif

#endif /* CADO_UTILS_MEMALLOC_H */
