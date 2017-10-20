#ifndef LIST_INT64_VECTOR_INDEX_H
#define LIST_INT64_VECTOR_INDEX_H 

#include <stdio.h>
#include <stdint.h>
#include "utils.h"
#include "int64_vector.h"
#include "mat_int64.h"

#ifdef __cplusplus
extern "C" {
#endif


/*
 * int64_vector with its index.
 */
typedef struct
{
  int64_vector_t vec;
  uint64_t index;
}  s_int64_vector_index_t;

typedef s_int64_vector_index_t int64_vector_index_t[1];
typedef s_int64_vector_index_t * int64_vector_index_ptr;
typedef const s_int64_vector_index_t * int64_vector_index_srcptr;

/*
 * List of int64_vector_index.
 */
typedef struct
{
  int64_vector_index_t * v;
  unsigned int length;
  unsigned int alloc;
  unsigned int vector_dim;
}  s_list_int64_vector_index_t;

typedef s_list_int64_vector_index_t list_int64_vector_index_t[1];
typedef s_list_int64_vector_index_t * list_int64_vector_index_ptr;
typedef const s_list_int64_vector_index_t * list_int64_vector_index_srcptr;

/*
 * Initialise a list of int64_vector_index.
 * 
 * list: the list.
 */
void list_int64_vector_index_init(list_int64_vector_index_ptr list,
    unsigned int vector_dim);

/*
 * Add a vector_index in the list.
 *
 * list: the list of vector.
 * v: the vector we add.
 */
void list_int64_vector_index_add_int64_vector_index(list_int64_vector_index_ptr list, int64_vector_srcptr v, uint64_t index);

/*
 * Delete the list.
 *
 * list: the list.
 */
void list_int64_vector_index_clear(list_int64_vector_index_ptr list);

/*
 * Write a list in a file.
 *
 * file: the file in which we write.
 * lest: the list.
 */
void list_int64_vector_index_fprintf(FILE * file, list_int64_vector_index_srcptr list);

/*
 * Write a list in a file, with # before each lines.
 *
 * file: the file in which we write.
 * lest: the list.
 */
void list_int64_vector_index_fprintf_comment(FILE * file,
    list_int64_vector_index_srcptr list);

/*
 * Sort the list by the last coordinate of each vectors.
 */
void list_int64_vector_index_sort_last(list_int64_vector_index_ptr list);

/*
 * Delete the element at pos.
 */
void list_int64_vector_index_delete_int64_vector_index(
    list_int64_vector_index_ptr list, unsigned int pos);

/*
 * Remove duplicate vectors.
 */
void list_int64_vector_index_remove_duplicate(
    list_int64_vector_index_ptr list);

/*
 * Remove duplicate vectors and sort the list.
 */
void list_int64_vector_index_remove_duplicate_sort(
    list_int64_vector_index_ptr list);

void list_int64_vector_index_set(
    list_int64_vector_index_ptr list_new,
    list_int64_vector_index_srcptr list_old);

#ifdef __cplusplus
}
#endif
#endif /* LIST_INT64_VECTOR_INDEX_H */
