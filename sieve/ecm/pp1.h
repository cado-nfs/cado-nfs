#ifndef SIEVE_ECM_PP1_H_
#define SIEVE_ECM_PP1_H_

#include "modredc_ul.h"
#include "modredc_15ul.h"
#include "modredc_2ul2.h"
#include "mod_mpz.h"
#include "bytecode.h"
#include "stage2.h"

typedef struct {
  bytecode bc;          /* Bytecode for stage 1 */
  unsigned int exp2;    /* Exponent of 2 in stage 1 primes */
  unsigned int B1;
  stage2_plan_t stage2;
} pp1_plan_t;

#ifdef __cplusplus
extern "C" {
#endif


int pp1_27_ul (modintredcul_t, const modulusredcul_t, const pp1_plan_t *);
int pp1_65_ul (modintredcul_t, const modulusredcul_t, const pp1_plan_t *);
int pp1_stage2_ul (residueredcul_t, const residueredcul_t,
                   const stage2_plan_t *, const residueredcul_t,
                   const modulusredcul_t);

int pp1_27_15ul (modintredc15ul_t, const modulusredc15ul_t, const pp1_plan_t *);
int pp1_65_15ul (modintredc15ul_t, const modulusredc15ul_t, const pp1_plan_t *);
int pp1_stage2_15ul (residueredc15ul_t, const residueredc15ul_t,
                     const stage2_plan_t *, const residueredc15ul_t,
                     const modulusredc15ul_t);

int pp1_27_2ul2 (modintredc2ul2_t, const modulusredc2ul2_t, const pp1_plan_t *);
int pp1_65_2ul2 (modintredc2ul2_t, const modulusredc2ul2_t, const pp1_plan_t *);
int pp1_stage2_2ul2 (residueredc2ul2_t, const residueredc2ul2_t,
                     const stage2_plan_t *, const residueredc2ul2_t,
                     const modulusredc2ul2_t);

int pp1_27_mpz (modintmpz_t, const modulusmpz_t, const pp1_plan_t *);
int pp1_65_mpz (modintmpz_t, const modulusmpz_t, const pp1_plan_t *);
int pp1_stage2_mpz (residuempz_t, const residuempz_t,
                    const stage2_plan_t *, const residuempz_t,
                    const modulusmpz_t);
void pp1_make_plan (pp1_plan_t *, const unsigned int, const unsigned int, int);

void pp1_clear_plan (pp1_plan_t *);

#ifdef __cplusplus
}
#endif

#endif	/* SIEVE_ECM_PP1_H_ */
