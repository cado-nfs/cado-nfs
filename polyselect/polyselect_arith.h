#ifndef POLYSELECT_ARITH_H
#define POLYSELECT_ARITH_H

#include <gmp.h>
#include <stdint.h> // uint64_t
#include "polyselect_str.h"
#include "gmp_aux.h"

/* declarations */

#ifdef __cplusplus
extern "C" {
#endif

unsigned long invert (unsigned long, unsigned long);

unsigned long roots_lift (uint64_t*, mpz_srcptr, unsigned long, mpz_srcptr,
                          unsigned long, unsigned long int);

void first_comb (unsigned long, unsigned long *);

unsigned long next_comb (unsigned long, unsigned long, unsigned long *);

void print_comb (unsigned long, unsigned long *);

unsigned long number_comb (polyselect_qroots_srcptr SQ_R, unsigned long k, unsigned long lq);

unsigned long binom (unsigned long, unsigned long);

void comp_sq_roots (polyselect_poly_header_srcptr, polyselect_qroots_ptr, gmp_randstate_ptr);

void crt_sq (mpz_ptr, mpz_ptr, unsigned long *, unsigned long *, unsigned long);

uint64_t return_q_rq (polyselect_qroots_srcptr, unsigned long *, unsigned long,
                      mpz_ptr, mpz_ptr);

uint64_t return_q_norq (polyselect_qroots_srcptr, unsigned long *, unsigned long);

#ifdef __cplusplus
}
#endif

#endif
