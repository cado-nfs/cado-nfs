#ifndef FILTER_MEMALLOC_H_
#define FILTER_MEMALLOC_H_

#include <stddef.h>     // size_t
#include "typedefs.h"

/********************** own memory allocation routines ***********************/

/* Rationale: calling one malloc() for each read relation is expensive, since 
   malloc() allocates some extra information to keep track of every memory 
   blocks. Instead, we allocate memory in big blocks of size BLOCK_SIZE. */


/* memory blocks are allocated of that # of index_t's */
#define BLOCK_SIZE (1<<20)
#define BLOCK_INDEX_BYTES (BLOCK_SIZE * sizeof(index_t))
#define BLOCK_IDEALMERGE_BYTES (BLOCK_SIZE * sizeof(ideal_merge_t))
#define INIT_NB_BLOCK (1<<16)

#ifdef __cplusplus
extern "C" {
#endif

index_t * index_my_malloc (unsigned int);
ideal_merge_t * idealmerge_my_malloc (unsigned int);
void my_malloc_free_all ();
size_t get_my_malloc_bytes ();

#ifdef __cplusplus
}
#endif

#endif /* FILTER_MEMALLOC_H_ */
