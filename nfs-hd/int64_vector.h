#ifndef INT64_VECTOR_H
#define INT64_VECTOR_H

#include <stdio.h>
#include <math.h>
#include "utils.h"
#include "sieving_bound.h"
#include "vector.h"

/* TODO: should be std::vector<int64_t> */

#ifdef __cplusplus
extern "C" {
#endif


/*
 * Initialise a vector.
 *
 * v: vector.
 * d: dimension of the vector.
 */
void int64_vector_init(int64_vector_ptr v, unsigned int d);

/*
 * Delete a vector.
 *
 * v: vector.
 */
void int64_vector_clear(int64_vector_ptr v);

/*
 * Swap v0 and v1, assume that v0 and v1 have the same dimension.
 *
 * v0: first vector.
 * v1: second vector.
 */
void int64_vector_swap(int64_vector_ptr v0, int64_vector_ptr v1);

/*
 * Copy s in v. Assume that s and v have the same dimension.
 *
 * v: the root vector.
 * s: the new vector.
 */ 
void int64_vector_set(int64_vector_ptr v, int64_vector_srcptr s);

/*
 * Set all the coordinate of v to 0.
 *
 * v: the vector.
 */

void int64_vector_set_zero(int64_vector_ptr v);
/*
 * Get the ith coefficient of v.
 *
 * v: vector.
 * i: index of the coefficient.
 */
static inline int64_t int64_vector_getcoordinate(int64_vector_srcptr v,
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
static inline void int64_vector_setcoordinate(int64_vector_ptr v,
    unsigned int i, int64_t d)
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
int int64_vector_equal(int64_vector_srcptr a, int64_vector_srcptr b);

/*
 * a = b + c
 *
 * a: the result.
 * b: the first vector.
 * c: the second vector.
 */
void int64_vector_add(int64_vector_ptr a, int64_vector_srcptr b,
    int64_vector_srcptr c);

/*
 * a = b - c
 *
 * a: the result.
 * b: the first vector.
 * c: the second vector.
 */
void int64_vector_sub(int64_vector_ptr a, int64_vector_srcptr b,
    int64_vector_srcptr c);

/*
 * a = c * b
 *
 * a: the result.
 * b: the first vector.
 * c: the factor.
 */
void int64_vector_mul(int64_vector_ptr a, int64_vector_srcptr b,
                      int64_t c);

/*
 * a = b / c
 */
void int64_vector_div(int64_vector_ptr a, int64_vector_srcptr b,
    int64_t c);

/*
 * a = b + d * c
 */
void int64_vector_addmul(int64_vector_ptr a, int64_vector_srcptr b,
    int64_vector_srcptr c, int64_t d);
/*
 * Write a vector in a file.
 *
 * file: the file.
 * v: the vector.
 */
void int64_vector_fprintf(FILE * file, int64_vector_srcptr v);

/*
 * Add one to the first coordinate of v while v is in the sieving region.
 *
 * v: vector.
 * H: sieving bound that describes the sieving region.
 */
void int64_vector_add_one(int64_vector_ptr v, sieving_bound_srcptr H);

/*
 * Add one the the ith coordinate of v while v is in the sieving region.
 *
 * v: vector.
 * H: sieving bound that describes the sieving region.
 */
unsigned int int64_vector_add_one_i(int64_vector_ptr v, unsigned int i,
                                    sieving_bound_srcptr H);

/*
 * Build a mpz_vector from a int64_vector.
 *
 * a: the new mpz_vector.
 * b: the old int64_vector.
 */
void int64_vector_to_mpz_vector(mpz_vector_ptr a, int64_vector_srcptr b);

/*
 * Compute the square of the norm of the vector a.
 *
 * a: vector.
 */
uint64_t int64_vector_norml2sqr(int64_vector_srcptr a);

/*
 * Compute the norm of a vector a.
 *
 * a: vector.
 */
double int64_vector_norml2(int64_vector_srcptr a);

/*
 * Compute the dot product between v0 and v1.
 *
 * v0: first vector.
 * v1: second vector.
 */
int64_t int64_vector_dot_product(int64_vector_srcptr v0,
    int64_vector_srcptr v1);

/*
 * Return 1 if a vector is in the sieving region, 0 otherwise.
 *
 * v: vector.
 * H: sieving bound that gives a sieving region.
 */
int int64_vector_in_sieving_region(int64_vector_srcptr v,
    sieving_bound_srcptr H);

/*
 * Transform a double_vector in an int64_vector.
 *
 * v_i: the int64 vector.
 * v_d: the double vector.
 */
void double_vector_in_int64_vector(int64_vector_ptr v_i,
    double_vector_srcptr v_d);

/*
 * Compute the gcd between all the coefficients of v.
 */
int64_t int64_vector_gcd(int64_vector_srcptr v);

/*
 * u = v / gcd(v)
 */
void int64_vector_reduce(int64_vector_ptr u, int64_vector_srcptr v);

#ifdef __cplusplus
}
#endif
#endif // INT64_VECTOR_H
