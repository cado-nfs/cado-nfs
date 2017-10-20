#ifndef LIST_H
#define LIST_H 

#include "int64_vector.h"
#include "double_vector.h"

/* TODO: eeek ? Are there really "lists", a.k.a. data structures with
 * guaranteed constant-time insertion ??
 */
/*
 * List of int64_vector.
 */
typedef struct
{
  unsigned int alloc;
  unsigned int length;
  unsigned int vector_dim;
  int64_vector_t * v;
}  s_list_int64_vector_t;

typedef s_list_int64_vector_t list_int64_vector_t[1];
typedef s_list_int64_vector_t * list_int64_vector_ptr;
typedef const s_list_int64_vector_t * list_int64_vector_srcptr;

/*
 * List of double vector.
 */
typedef struct
{
  double_vector_t * v;
  unsigned int length;
}  s_list_double_vector_t;

typedef s_list_double_vector_t list_double_vector_t[1];
typedef s_list_double_vector_t * list_double_vector_ptr;
typedef const s_list_double_vector_t * list_double_vector_srcptr;

#endif /* LIST_H */
