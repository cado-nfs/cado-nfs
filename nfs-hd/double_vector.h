#ifndef DOUBLE_VECTOR_H
#define DOUBLE_VECTOR_H 

#include <stdio.h>
#include <math.h>
#include "utils.h"
#include "sieving_bound.h"
#include "vector.h"

// TODO: STL

#ifdef __cplusplus
extern "C" {
#endif


/*
 * Initialise a vector.
 *
 * v: vector.
 * d: dimension of the vector.
 */
void double_vector_init(double_vector_ptr v, unsigned int d);

/*
 * Delete a vector.
 *
 * v: vector.
 */
void double_vector_clear(double_vector_ptr v);

/*
 * Swap v0 and v1, assume that v0 and v1 have the same dimension.
 *
 * v0: first vector.
 * v1: second vector.
 */
void double_vector_swap(double_vector_ptr v0, double_vector_ptr v1);

/*
 * Copy s in v. Assume that s and v have the same dimension.
 *
 * v: the root vector.
 * s: the new vector.
 */ 
void double_vector_set(double_vector_ptr v, double_vector_srcptr s);

/*
 * Set all the coordinate of v to 0.
 *
 * v: the vector.
 */
void double_vector_set_zero(double_vector_ptr v);
/*
 * Get the ith coefficient of v.
 *
 * v: vector.
 * i: index of the coefficient.
 */
static inline double double_vector_getcoordinate(double_vector_srcptr v,
    unsigned int i)
{
  ASSERT(i < v->dim);
  return v->c[i];
}

/*
 * Set the ith coefficient of v.
 *
 * v: vector.
 * i: index of the coefficient.
 * d: value of the coefficient.
 */
static inline void double_vector_setcoordinate(double_vector_ptr v,
    unsigned int i, double d)
{
  ASSERT(i < v->dim);
  v->c[i] = d;
}

/*
 * Return 1 if the vector are equal, 0 otherwise.
 *
 * a: first vector.
 * b: second vector.
 */
int double_vector_equal(double_vector_srcptr a, double_vector_srcptr b);

/*
 * a = b + c
 *
 * a: the result.
 * b: the first vector.
 * c: the second vector.
 */
void double_vector_add(double_vector_ptr a, double_vector_srcptr b,
    double_vector_srcptr c);

/*
 * a = b - c
 *
 * a: the result.
 * b: the first vector.
 * c: the second vector.
 */
void double_vector_sub(double_vector_ptr a, double_vector_srcptr b,
    double_vector_srcptr c);

/*
 * Write a vector in a file.
 *
 * file: the file.
 * v: the vector.
 */
void double_vector_fprintf(FILE * file, double_vector_srcptr v);

/*
 * Compute the square of the norm of the vector a.
 *
 * a: vector.
 */
double double_vector_norml2sqr(double_vector_srcptr a);

/*
 * Compute the norm of a vector a.
 *
 * a: vector.
 */
double double_vector_norml2(double_vector_srcptr a);

/*
 * Compute the dot product between v0 and v1.
 *
 * v0: first vector.
 * v1: second vector.
 */
double double_vector_dot_product(double_vector_srcptr v0,
    double_vector_srcptr v1);

/*
 * Return 1 if a vector is in the sieving region, 0 otherwise.
 *
 * v: vector.
 * H: sieving bound that gives a sieving region.
 */
int double_vector_in_sieving_region(double_vector_srcptr v,
    sieving_bound_srcptr H);

/*
 * Transform an int64_vector in a double_vector.
 *
 * v_d: the double vector.
 * v_i: the int64 vector.
 */
void int64_vector_to_double_vector(double_vector_ptr v_d,
    int64_vector_srcptr v_i);

/*
 * res is the orthogonal projection of v on u.
 *
 * res: an orthogonal vector.
 * u: the vector on which we project v.
 * v: the projected vector.
 */
double double_vector_orthogonal_projection(double_vector_ptr res,
    double_vector_srcptr u, double_vector_srcptr v);

#ifdef __cplusplus
}
#endif

#endif /* DOUBLE_VECTOR_H */
