#ifndef LIST_INT64_VECTOR_H
#define LIST_INT64_VECTOR_H 

#include <stdio.h>
#include <stdint.h>
#include "utils.h"
#include "int64_vector.h"
#include "mat_int64.h"
#include "list.h"

/* question/TODO: should this be a list<uint64_t> ? what for by the way ?
 * */
#ifdef __cplusplus
extern "C" {
#endif


/*
 * Initialise a list of int64_vector.
 * 
 * list: the list.
 */
void list_int64_vector_init(list_int64_vector_ptr list,
    unsigned int vector_dim);

/*
 * Add a vector in the list.
 *
 * list: the list of vector.
 * v: the vector we add.
 */
unsigned int list_int64_vector_add_int64_vector(list_int64_vector_ptr list,
    int64_vector_srcptr v);

/*
 * Delete the list.
 *
 * list: the list.
 */
void list_int64_vector_clear(list_int64_vector_ptr list);

/*
 * Set length of the list to zero, do not delete the elements.
 */
void list_int64_vector_delete_elements(list_int64_vector_ptr list);
/*
 * Write a list in a file.
 *
 * file: the file in which we write.
 * lest: the list.
 */
void list_int64_vector_fprintf(FILE * file, list_int64_vector_srcptr list);

/*
 * Write a list in a file with # before each lines.
 */
void list_int64_vector_fprintf_comment(FILE * file,
    list_int64_vector_srcptr list);
/*
 * Extract all the vector (in column) of the matrix and put them in a list.
 *
 * list: the list of vector.
 * matrix: the matrix.
 */
void list_int64_vector_extract_mat_int64(list_int64_vector_ptr list,
    mat_int64_srcptr matrix);

/*
 * TODO: continue here.
 * Return 1 a vector is inside the polytop formed by all the vector in the
 * list. Works only with vector of dimension 2.
 * This code is inspired from PNPOLY
 * (http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html).
 */
int int64_vector_in_polytop_list_int64_vector(int64_vector_srcptr vec,
    list_int64_vector_srcptr list);

/*
 * Sort vectors of the list by last coordinate.
 */
void list_int64_vector_sort_last(list_int64_vector_ptr list);

/*
 * Assume that list is sorted by last coordinate. Delete all the elements with
 *  a last coordinate strictly less than val.
 */
void list_int64_vector_delete_last_coordinate(list_int64_vector_ptr list,
    int64_t val);

void list_int64_vector_remove_duplicate(list_int64_vector_ptr list);

void list_int64_vector_delete_int64_vector(list_int64_vector_ptr list,
    unsigned int pos);

#ifdef __cplusplus
}
#endif
#endif /* LIST_INT64_VECTOR_H */
