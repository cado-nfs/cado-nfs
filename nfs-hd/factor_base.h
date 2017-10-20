#ifndef FACTOR_BASE_H
#define FACTOR_BASE_H

#include <stdio.h>
#include <stdint.h>
#include "ideal.h"

#ifdef __cplusplus
extern "C" {
#endif


/*
 * Representation of a factor base.
 */
typedef struct {
  //For ideal of degree 1.
  ideal_1_t * factor_base_1;
  uint64_t number_element_1;
  //For ideal of degree > 1.
  ideal_u_t * factor_base_u;
  uint64_t number_element_u;
  //For ideal with a projective root.
  ideal_pr_t * factor_base_pr;
  uint64_t number_element_pr;
} s_factor_base_t;

typedef s_factor_base_t factor_base_t[1];
typedef s_factor_base_t * factor_base_ptr;
typedef const s_factor_base_t * factor_base_srcptr;

/*
 * Initialise a factor base.
 *
 * factor_base: the factor base.
 * number_element: number of element (generally, an upper bound).
 * t: dimension of the lattice.
 */
void factor_base_init(factor_base_ptr factor_base, uint64_t number_element_1,
    uint64_t number_element_u, uint64_t number_element_pr);

/*
 * Realloc the factor base. The new number of element must be less than the old
 *  number of element stored in the factor base.
 *
 *  factor_base: the factor base.
 *  new_number_element: the number of element in the factor base.
 */
void factor_base_realloc(factor_base_ptr factor_base,
    uint64_t number_element_1, uint64_t number_element_u,
    uint64_t number_element_pr);

/*
 * Delete the factor base bound.
 *
 * factor_base: the factor base.
 */
void factor_base_clear(factor_base_ptr factor_base, unsigned int t);

/*
 * Write the factor base bound.
 *
 * file: the file in which we want to write.
 * factor_base: the factor base we want to write.
 */
void factor_base_fprintf(FILE * file, factor_base_srcptr factor_base,
    unsigned int t);

/*
 * Set an ideal in part at an index.
 *
 * factor_base: the factor base.
 * index: index in the factor base.
 * r: the r of the ideal (r, h).
 * h: the h of the ideal (r, h).
 */
static inline void factor_base_set_ideal_1_part(factor_base_ptr factor_base,
    uint64_t index, uint64_t r, mpz_poly_srcptr h, unsigned int t)
{
  ASSERT(index < factor_base->number_element_1);

  ideal_1_init(factor_base->factor_base_1[index]);
  ideal_1_set_part(factor_base->factor_base_1[index], r, h, t);
}

/*
 * Set an ideal at an index.
 *
 * factor_base: the factor base.
 * index: index in the factor base.
 * ideal: the ideal.
 * t: dimension of the lattice.
 */
static inline void factor_base_set_ideal_1(factor_base_ptr factor_base,
    uint64_t index, ideal_1_srcptr ideal, unsigned int t)
{
  ASSERT(index < factor_base->number_element_1);

  ideal_1_init(factor_base->factor_base_1[index]);
  ideal_1_set(factor_base->factor_base_1[index], ideal, t);
}

/*
 * Set an ideal in part at an index.
 *
 * factor_base: the factor base.
 * index: index in the factor base.
 * r: the r of the ideal (r, h).
 * h: the h of the ideal (r, h).
 */
static inline void factor_base_set_ideal_u_part(factor_base_ptr factor_base,
    uint64_t index, uint64_t r, mpz_poly_srcptr h, unsigned int t)
{
  ASSERT(index < factor_base->number_element_u);

  ideal_u_init(factor_base->factor_base_u[index]);
  ideal_u_set_part(factor_base->factor_base_u[index], r, h, t);
}

/*
 * Set an ideal at an index.
 *
 * factor_base: the factor base.
 * index: index in the factor base.
 * ideal: the ideal.
 * t: dimension of the lattice.
 */
static inline void factor_base_set_ideal_u(factor_base_ptr factor_base,
    uint64_t index, ideal_u_srcptr ideal, unsigned int t)
{
  ASSERT(index < factor_base->number_element_u);

  ideal_u_init(factor_base->factor_base_u[index]);
  ideal_u_set(factor_base->factor_base_u[index], ideal, t);
}

/*
 * Set an ideal in part at an index.
 *
 * factor_base: the factor base.
 * index: index in the factor base.
 * r: the r of the ideal (r, h).
 */
static inline void factor_base_set_ideal_pr(factor_base_ptr factor_base,
    uint64_t index, uint64_t r, unsigned int t)
{
  ASSERT(index < factor_base->number_element_pr);

  ideal_pr_init(factor_base->factor_base_pr[index]);
  ideal_pr_set_part(factor_base->factor_base_pr[index], r, t);
}


#ifdef __cplusplus
}
#endif
#endif  /* FACTOR_BASE_H */
