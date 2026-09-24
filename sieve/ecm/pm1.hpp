#ifndef CADO_SIEVE_ECM_PM1_HPP
#define CADO_SIEVE_ECM_PM1_HPP

#include "arith/modredc_ul.h"
#include "arith/modredc_15ul.h"
#include "arith/modredc_2ul2.h"
#include "arith/mod_mpz.h"
#include "cxx_mpz.hpp"
#include "stage2.h"

struct pm1_plan_t {
    cxx_mpz E;              /* The exponent for stage 1 */
    unsigned int exp2 = 0;  /* Exponent of 2 in stage 1 primes */
    unsigned int B1 = 0;
    stage2_plan_t stage2 {};
};

/* P-1 with the given plan. The arithmetic is done with the arithxx layer
 * that suits the size of the modulus. The factor found, or 1, is written
 * to f. Returns 1 if backtracking was used, 0 otherwise.
 */
int pm1(modintredcul_t f, const modulusredcul_t m, const pm1_plan_t * plan);
int pm1(modintredc15ul_t f, const modulusredc15ul_t m, const pm1_plan_t * plan);
int pm1(modintredc2ul2_t f, const modulusredc2ul2_t m, const pm1_plan_t * plan);
int pm1(modintmpz_t f, const modulusmpz_t m, const pm1_plan_t * plan);

void pm1_make_plan (pm1_plan_t *, unsigned int, unsigned int, int);
void pm1_clear_plan (pm1_plan_t *);

/* Notes:
  If s_1 * s_2 = eulerphi(d), we need s_1 precomputed values, s_2 passes and
  ~(B2-B1)/d steps in each pass. Assuming we can compute V_{k_1}(x+1/x)
  and V_{k_2}(x+1/x) with about 1 multiply per value by using a common
  addition chain for S_1 \cup S_2, we can reduce the precomputation cost to
  only s_1 + s_2 instead of eulerphi(P) if we do only 1 pass. This allows for
  a larger d: with only 1 pass, sqrt(B2-B1) is optimal, with a flexible number
  of passes, (B2-B1)^{2/3} is optimal.

  In one pass we process one k_2 \in S_2 and thus primes p = i*d +- j with
  j = k_1 + k_2, k_1 \in S_1. We'd like to be able to bound these j by d/2
  so that we get compact "block" with no overlap between consecutive blocks.
  This prevents having to start at lower i/continue to larger i than
  necessary, to be able to write all p, B_1 < p <= B_2, as
  p = i*d +- (k_1 + k_2). However, there seems to be no way to write the
  smallest positive representatives <d/2 of units mod d as a sum of two sets:
  e.g. for d=30, we'd like {1,7,11,13}.
*/

#endif	/* CADO_SIEVE_ECM_PM1_HPP */
