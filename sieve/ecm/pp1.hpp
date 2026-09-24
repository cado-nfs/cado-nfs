#ifndef CADO_SIEVE_ECM_PP1_HPP
#define CADO_SIEVE_ECM_PP1_HPP

#include "arith/modredc_ul.h"
#include "arith/modredc_15ul.h"
#include "arith/modredc_2ul2.h"
#include "arith/mod_mpz.h"
#include "bytecode.h"
#include "stage2.h"

struct pp1_plan_t {
  bytecode bc;          /* Bytecode for stage 1 */
  unsigned int exp2;    /* Exponent of 2 in stage 1 primes */
  unsigned int B1;
  stage2_plan_t stage2;
};

/* P+1 with the given plan, starting from 2/7 (pp1_27) or 6/5 (pp1_65).
 * The arithmetic is done with the arithxx layer that suits the size of
 * the modulus. The factor found, or 1, is written to f. Returns 1 if
 * backtracking was used, 0 otherwise.
 */
int pp1_27(modintredcul_t f, const modulusredcul_t m, const pp1_plan_t * plan);
int pp1_27(modintredc15ul_t f, const modulusredc15ul_t m, const pp1_plan_t * plan);
int pp1_27(modintredc2ul2_t f, const modulusredc2ul2_t m, const pp1_plan_t * plan);
int pp1_27(modintmpz_t f, const modulusmpz_t m, const pp1_plan_t * plan);
int pp1_65(modintredcul_t f, const modulusredcul_t m, const pp1_plan_t * plan);
int pp1_65(modintredc15ul_t f, const modulusredc15ul_t m, const pp1_plan_t * plan);
int pp1_65(modintredc2ul2_t f, const modulusredc2ul2_t m, const pp1_plan_t * plan);
int pp1_65(modintmpz_t f, const modulusmpz_t m, const pp1_plan_t * plan);

void pp1_make_plan (pp1_plan_t *, unsigned int, unsigned int, int);
void pp1_clear_plan (pp1_plan_t *);

#endif	/* CADO_SIEVE_ECM_PP1_HPP */
