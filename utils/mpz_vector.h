#ifndef MPZ_VECTOR_H_
#define MPZ_VECTOR_H_

#include <stdio.h> // FILE
#include <stdint.h>
#include <gmp.h>
#include "mpz_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

/* TODO: This should be an std::vector<cxx_mpz>.
 */
typedef struct 
{
  unsigned int dim; /* dimension of the vector */
  mpz_t *c;         /* its coordinates */
} mpz_vector_struct_t;

typedef mpz_vector_struct_t mpz_vector_t[1];
typedef mpz_vector_struct_t * mpz_vector_ptr;
typedef const mpz_vector_struct_t * mpz_vector_srcptr;

/* Management of the structure: init, clear, set and swap. */
void mpz_vector_init (mpz_vector_ptr, unsigned int);
void mpz_vector_clear(mpz_vector_ptr);
void mpz_vector_swap (mpz_vector_ptr, mpz_vector_ptr);
void mpz_vector_set (mpz_vector_ptr, mpz_vector_srcptr);
void mpz_vector_setcoordinate (mpz_vector_ptr, unsigned int, mpz_srcptr);
void mpz_vector_setcoordinate_ui (mpz_vector_ptr, unsigned int, unsigned int);
void mpz_vector_setcoordinate_si (mpz_vector_ptr, unsigned int, int);
void mpz_vector_setcoordinate_uint64 (mpz_vector_ptr, unsigned int, uint64_t);
void mpz_vector_setcoordinate_int64 (mpz_vector_ptr, unsigned int, int64_t);
int mpz_vector_is_coordinate_zero (mpz_vector_srcptr, unsigned int);
/*
  Return 0 if a and b are equal,
  -1 if a is smaller and 1 if a is bigger.
*/
int mpz_vector_cmp (mpz_vector_srcptr a, mpz_vector_srcptr b);
void mpz_vector_fprintf(FILE * file, mpz_vector_srcptr v);

/* Implementation of dot product and norm (skew and non-skew version) */
void mpz_vector_dot_product (mpz_ptr, mpz_vector_srcptr, mpz_vector_srcptr);
void mpz_vector_skew_dot_product (mpz_ptr, mpz_vector_srcptr, mpz_vector_srcptr, mpz_srcptr);
void mpz_vector_norm (mpz_ptr, mpz_vector_srcptr);
void mpz_vector_skew_norm (mpz_ptr, mpz_vector_srcptr, mpz_srcptr);

/* Convert from mpz_vector_t to mpz_poly */
void mpz_vector_get_mpz_poly (mpz_poly_ptr, mpz_vector_srcptr);

/* Operations on vectors */
void mpz_vector_submul (mpz_vector_ptr r, mpz_srcptr q, mpz_vector_srcptr v);

/* Lagrange algo to reduce lattice of rank 2 */
void mpz_vector_Lagrange (mpz_vector_ptr, mpz_vector_ptr,
                          mpz_vector_srcptr, mpz_vector_srcptr, mpz_srcptr);
void mpz_vector_reduce_with_max_skew (mpz_vector_ptr, mpz_vector_ptr, mpz_ptr,
                                      mpz_vector_srcptr, mpz_vector_srcptr, mpz_srcptr, int);

#ifdef __cplusplus
}
#endif

#endif	/* MPZ_VECTOR_H_ */
