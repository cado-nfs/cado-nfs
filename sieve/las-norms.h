#ifndef LAS_NORMS_H_
#define LAS_NORMS_H_

#include <stdint.h>
#include "las-types.h"

#ifdef __cplusplus
extern "C" {
#endif

/*  initializing norms */
/* Knowing the norm on the rational side is bounded by 2^(2^k), compute
   lognorms approximations for k bits of exponent + NORM_BITS-k bits
   of mantissa */
void
init_norms (sieve_info_ptr si);


/*  initialize norms for bucket regions */
/* Initialize lognorms on the rational side for the bucket_region
 * number N.
 * For the moment, nothing clever, wrt discarding (a,b) pairs that are
 * not coprime, except for the line j=0.
 */
void
init_rat_norms_bucket_region (unsigned char *S, const int N, sieve_info_ptr si);

/* Initialize lognorms on the algebraic side for the bucket
 * number N.
 * Only the survivors of the rational sieve will be initialized, the
 * others are set to 255. Case GCD(i,j)!=1 also gets 255.
 * return the number of reports (= number of norm initialisations)
 */
int
init_alg_norms_bucket_region (unsigned char *alg_S, 
                              const unsigned char *rat_S, const int N, 
                              sieve_info_ptr si);

/* This prepares the auxiliary data which is used by
 * init_rat_norms_bucket_region and init_alg_norms_bucket_region
 */
void sieve_info_init_norm_data(sieve_info_ptr si, mpz_srcptr q0);

void sieve_info_clear_norm_data(sieve_info_ptr si);

void sieve_info_update_norm_data(sieve_info_ptr si);

/* Determine whether a sieve entry with sieve residue S1 on sieving side 1
   and sieve residue S2 on sieving side 2 is likely smooth. 
   The array entry C1[S1] is initialized by sieve_info_init_lognorm() 
   to something similar to 
   -log(Pr[norm on side 1 with sieve residue S1 is smooth]),
   similar for C2, S2. Assuming the two probabilities are independent enough,
   we can estimate the neg log of the probability that both sides are smooth 
   by C1[S1] + C2[S2]. 
   If that sum does not exceed a theshold, the corresponding sieve entry is 
   a sieve survivor. 
   Alternative: have a bit array telling whether (S1,S2) is likely smooth */
static inline int 
sieve_info_test_lognorm (const unsigned char *C1, 
                         const unsigned char *C2, 
                         const unsigned char S1,
                         const unsigned char S2,
                         const unsigned char threshold)
{
  return C1[S1] + C2[S2] <= threshold;
}


#ifdef __cplusplus
}
#endif

#endif	/* LAS_NORMS_H_ */
