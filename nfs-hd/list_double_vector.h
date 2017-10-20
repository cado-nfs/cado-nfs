#ifndef LIST_DOUBLE_VECTOR_H
#define LIST_DOUBLE_VECTOR_H 

#include <stdio.h>
#include "utils.h"
#include "list.h"
#include "mat_double.h"

#define DEFAULT_LENGTH_LIST_DOUBLE_VECTOR 3

/* TODO: STL */
#ifdef __cplusplus
extern "C" {
#endif

/*
 * Initialise a list of int64_vector.
 * 
 * list: the list.
 */
void list_double_vector_init(list_double_vector_ptr list);

/*
 * Add a vector in the list.
 *
 * list: the list of vector.
 * v: the vector we add.
 */
void list_double_vector_add_double_vector(list_double_vector_ptr list, double_vector_srcptr v);

/*
 * Delete the list.
 *
 * list: the list.
 */
void list_double_vector_clear(list_double_vector_ptr list);

/*
 * Write a list in a file.
 *
 * file: the file in which we write.
 * lest: the list.
 */
void list_double_vector_fprintf(FILE * file, list_double_vector_srcptr list);

/*
 * Extract all the vector (in column) of the matrix and put them in a list.
 *
 * list: the list of vector.
 * matrix: the matrix.
 */
void list_double_vector_extract_mat_int64(list_double_vector_ptr list,
    mat_int64_srcptr matrix);


#ifdef __cplusplus
}
#endif


#endif /* LIST_DOUBLE_VECTOR_H */
