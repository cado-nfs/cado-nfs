#ifndef VECTOR_H
#define VECTOR_H 

// pragma no prototypes
#include <stdint.h>

/*
 * int64_vector.
 */
typedef struct
{
  unsigned int dim;
  int64_t * c;
} int64_vector_struct_t;

typedef int64_vector_struct_t int64_vector_t[1];
typedef int64_vector_struct_t * int64_vector_ptr;
typedef const int64_vector_struct_t * int64_vector_srcptr;

/*
 * double_vector.
 */
typedef struct
{
  unsigned int dim;
  double * c;
} double_vector_struct_t;

typedef double_vector_struct_t double_vector_t[1];
typedef double_vector_struct_t * double_vector_ptr;
typedef const double_vector_struct_t * double_vector_srcptr;


#endif /* VECTOR_H */
