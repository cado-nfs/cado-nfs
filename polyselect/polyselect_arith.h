#ifndef POLYSELECT_ARITH_H
#define POLYSELECT_ARITH_H

#include <gmp.h>
#include <stdint.h> // uint64_t
#include "polyselect_qroots.h"
#include "polyselect_poly_header.h"
#include "gmp_aux.h"

/* declarations */

#ifdef __cplusplus
extern "C" {
#endif

extern unsigned long invert (unsigned long, unsigned long);

extern unsigned long roots_lift (uint64_t*, mpz_srcptr, unsigned long, mpz_srcptr,
                          unsigned long, unsigned long int);


/* These functions don't really belong here. These are basic
 * combinatorics computations, we should expose them as such.
 */
extern void first_comb (unsigned long, unsigned long *);

extern unsigned long next_comb (unsigned long, unsigned long, unsigned long *);

extern void print_comb (unsigned long, unsigned long *);

extern unsigned long number_comb (polyselect_qroots_srcptr SQ_R, unsigned long k, unsigned long lq);

extern unsigned long binomial (unsigned long, unsigned long);

extern void comp_sq_roots (polyselect_poly_header_srcptr, polyselect_qroots_ptr, gmp_randstate_ptr);

extern void crt_sq (mpz_ptr, mpz_ptr, unsigned long *, unsigned long *, unsigned long);

extern uint64_t return_q_norq (polyselect_qroots_srcptr, unsigned long *, unsigned long);

#ifdef __cplusplus
}
#endif

#endif
