#ifndef CADO_SIEVE_ECM_PP1_HPP
#define CADO_SIEVE_ECM_PP1_HPP

#include "bytecode.h"
#include "stage2.h"

struct pp1_plan_t {
  bytecode bc;          /* Bytecode for stage 1 */
  unsigned int exp2;    /* Exponent of 2 in stage 1 primes */
  unsigned int B1;
  stage2_plan_t stage2;
};

/* P+1 with the given plan, modulo m, with an arithxx layer, starting from
 * 2/7 (pp1_27) or 6/5 (pp1_65). The factor found, or 1, is written to f.
 * Returns 1 if backtracking was used, 0 otherwise. pp1.cpp instantiates
 * these for modredc64, modredc96, modredc126 and mod_mpz_new.
 */
template<typename layer>
int pp1_27(typename layer::Integer & f, typename layer::Modulus const & m,
        pp1_plan_t const & plan);
template<typename layer>
int pp1_65(typename layer::Integer & f, typename layer::Modulus const & m,
        pp1_plan_t const & plan);

void pp1_make_plan (pp1_plan_t *, unsigned int, unsigned int, int);
void pp1_clear_plan (pp1_plan_t *);

#endif	/* CADO_SIEVE_ECM_PP1_HPP */
