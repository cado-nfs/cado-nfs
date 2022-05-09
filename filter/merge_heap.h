#ifndef MERGE_HEAP_H_
#define MERGE_HEAP_H_

#include "merge_replay_matrix.h"  // typerow_t, index_t

/* 
 * Custom memory management for arrays of typerow_t (i.e. sparse matrices)
 * where rows are frequently allocated/deleted.
 */


#ifdef __cplusplus
extern "C" {
#endif

/* Setup data structures. Must be called outside of a parallel region */
void heap_setup();

int heap_config_get_PAGE_DATA_SIZE();

/* Allocate space for row i, holding s coefficients in row[1:s+1] (row[0] == s).
   This writes s in row[0]. The size of the row must not change afterwards. 
   Thread-safe. Usually very fast. */
typerow_t * heap_alloc_row (index_t i, size_t s);

/* Shrinks the last allocated row. It must be the last one allocated by this thread. 
   Lock-free and constant-time. */
void heap_resize_last_row (typerow_t *row, index_t new_size);

/* given the pointer provided by heap_alloc_row, mark the row as deleted.
   The memory is not released. Lock-free and constant-time. */
void heap_destroy_row (typerow_t *row);

/* Reclaim all the memory occupied by destroyed rows. 
   Must be called outside of a parallel region.  May move rows in memory and update row pointers. */
void heap_garbage_collection(typerow_t **rows);

/* release all memory. This is technically not necessary, because the "malloc"
   allocations are internal to the process, and all space allocated to the
   process is reclaimed by the OS on termination. However, doing this enables
   valgrind to check the absence of leaks.
*/
void heap_clear ();

#ifdef __cplusplus
}
#endif

#endif	/* MERGE_HEAP_H_ */
