#ifndef SIEVING_INTERVAL_H
#define SIEVING_INTERVAL_H

#include <stdio.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif


typedef struct
{
  unsigned int * h;
  unsigned int t;
}  s_sieving_bound_t;

typedef s_sieving_bound_t sieving_bound_t[1];
typedef s_sieving_bound_t * sieving_bound_ptr;
typedef const s_sieving_bound_t * sieving_bound_srcptr;

/*
 * Initialise a sieving bound.
 *
 * H: the sieving bound.
 * t: dimension of the lattice.
 */
void sieving_bound_init(sieving_bound_ptr H, unsigned int t);

/*
 * Set the ith value of a sieving bound.
 *
 * H: the sieving bound.
 * i: index of the value.
 * value: the value of h[i].
 */
void sieving_bound_set_hi(sieving_bound_ptr H, unsigned int i, unsigned
    int value);

/*
 * Return the number of elements contained in the sieving region defined by H
 *  and set the result in nb.
 *
 * H: the sieving bound.
 */
uint64_t sieving_bound_number_element(sieving_bound_srcptr H);

/*
 * Delete a sieving bound.
 *
 * H: the sieving bound.
 */
void sieving_bound_clear(sieving_bound_ptr H);

/*
 * To write a sieving bound in a file.
 *
 * file: the file.
 * H: the sieving bound.
 */
void sieving_bound_fprintf(FILE * file, sieving_bound_srcptr H);

/*
 * To write a sieving bound with bounds in a file.
 *
 * file: the file.
 * H: the sieving bound.
 */
void sieving_bound_fprintf_detailed(FILE * filew, sieving_bound_srcptr H);

/*
 * To write, with # before, a sieving bound with bounds in a file.
 *
 * file: the file.
 * H: the sieving bound.
 */
void sieving_bound_fprintf_detailed_comment(FILE * filew,
    sieving_bound_srcptr H);

#ifdef __cplusplus
}
#endif
#endif /* SIEVING_INTERVAL_H */
