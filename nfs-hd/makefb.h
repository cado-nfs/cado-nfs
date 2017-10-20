#ifndef MAKEFB_H
#define MAKEFB_H

#include "utils.h"
#include <stdint.h>
#include "factor_base.h"

#ifdef __cplusplus
extern "C" {
#endif


/*
 * Make the factor bases for the V number fiels.
 *
 * fb: the V factor bases.
 * f: the V polynomials.
 * fbb: the V factor base bound.
 * t: dimension of the lattice.
 * lpb: the V large prime bound.
 * V: number of number fields.
 */
void makefb(factor_base_t * fb, cado_poly_srcptr f, uint64_t * fbb,
    unsigned int t, unsigned int * lpb);

/*
 * Read the factor base from a file.
 *
 * file: the file.
 * fb: the factor base.
 * fbb: the factor base bound.
 * lpb: the large prime bound.
 * f: polynomial that define the number field.
 */
void read_factor_base(FILE * file, factor_base_t * fb, uint64_t * fbb,
    unsigned int * lpb, cado_poly_srcptr f, double * log2_base, unsigned int t);


#ifdef __cplusplus
}
#endif
#endif /* MAKEFB_H */
