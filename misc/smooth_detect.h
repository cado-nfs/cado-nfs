#ifndef __SMOOTH_DETECT_H__
#define __SMOOTH_DETECT_H__

#include <gmp.h>

/*
 * Smoothness detector: take a list of pairs of integers and a smoothness
 * bound, test them for smoothness and returns one which is smooth.
 * It does not claim exhaustivity: if a candidate does not look promising
 * after a few ECMs, it is just skipped.
 *
 * The list of candidates is not given as a table, but a function must be
 * passed that returns the next candidate (computed or read from a
 * stream, whatever). A void * parameter is passed to this function, so
 * that it can have a context, and update it after each call.
 *
 * The next_cand() function must create candidates with something like:
 *   cand_set_original_values(C, u0, v0);
 * where u0 and v0 are the two integers that must be checked for
 * simultaneous smoothness.
 * Note that the candidate C is initialized and cleared by the caller.
 *
 * An example of usage is given in descent_init_Fp.c .
 */


// Type and basic functions for candidates
typedef struct {                                                                
  mpz_t u0, v0;      // original integers to be tested for smoothness
  mpz_t u, v;        // unfactored parts of u, v, must be composite or 1
  unsigned int lpu, lpv; // bitsize of largest primes found so far on both sides
  double effort;     // sum of the B1 already tried on these numbers
  unsigned long id;  // for the caller to remember who is who
} cand_s;                                                                       
typedef cand_s cand_t[1];

// Type for tuning parameters for smooth_detect. 
//   min_effort: the effort at the start (effort = sum of the B1 already tried)
//   max_effort: when this value is reached, don't increase effort anymore
//   max_pool_size: number of candidates to keep in mind at the same time
//   verbose: 0 = no verbose, 1 verbose
//   minB1: first B1 to try, in the ECM chain.
// Default values: {2000, +inf, 10, 1, 100.0}.
typedef struct {
  double min_effort;
  double max_effort;
  unsigned int max_pool_size;
  int verbose;
  double minB1;
} smooth_detect_param_s;

#ifdef __cplusplus
extern "C" {
#endif
void cand_init(cand_t c);
void cand_clear(cand_t c);
void cand_set(cand_t c, const cand_t d);
void cand_set_original_values(cand_t c, const mpz_t u0, const mpz_t v0,
        unsigned long id);
void cand_set_presieved_values(cand_t c, const mpz_t u0, const mpz_t v0,
        const mpz_t u, const mpz_t v,
        unsigned int lpu, unsigned int lpv,
        unsigned long id);


// The main exported function. Smooth candidate is put in C.
// Last argument is for changing default strategy. NULL can be passed.
void smooth_detect(cand_t C, int (*next_cand)(cand_t, void *),
        void *param_next_cand, unsigned long bound,
        const smooth_detect_param_s* param);
#ifdef __cplusplus
}
#endif



#endif   /* __SMOOTH_DETECT_H__ */
