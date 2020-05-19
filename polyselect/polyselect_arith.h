#ifndef POLYSELECT_ARITH_H
#define POLYSELECT_ARITH_H

#include "polyselect_str.h"

/* declarations */

#ifdef __cplusplus
extern "C" {
#endif

unsigned long invert (unsigned long, unsigned long);

unsigned long roots_lift (uint64_t*, mpz_t, unsigned long, mpz_t,
                          unsigned long, unsigned long int);

void first_comb (unsigned long, unsigned long *);

unsigned long next_comb (unsigned long, unsigned long, unsigned long *);

void print_comb (unsigned long, unsigned long *);

unsigned long number_comb (qroots_t SQ_R, unsigned long k, unsigned long lq);

unsigned long binom (unsigned long, unsigned long);

void comp_sq_roots (header_t, qroots_t);

void crt_sq (mpz_t, mpz_t, unsigned long *, unsigned long *, unsigned long);

uint64_t return_q_rq (qroots_t, unsigned long *, unsigned long,
                      mpz_t, mpz_t);

uint64_t return_q_norq (qroots_t, unsigned long *, unsigned long);

#ifdef __cplusplus
}
#endif

#endif
