#ifndef MERGE_HEAP_H_
#define MERGE_HEAP_H_

#include "merge_replay_matrix.h"

#ifdef __cplusplus
extern "C" {
#endif

typerow_t * heap_alloc_row (index_t i, size_t s);
void heap_resize_last_row (typerow_t *row, index_t new_size);
void heap_destroy_row (typerow_t *row);
void heap_clear ();
void heap_setup();
int heap_config_get_PAGE_DATA_SIZE();
void full_garbage_collection(filter_matrix_t *mat);

#ifdef __cplusplus
}
#endif

#endif	/* MERGE_HEAP_H_ */
