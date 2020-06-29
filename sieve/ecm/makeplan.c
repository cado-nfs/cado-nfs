#include "cado.h" // IWYU pragma: keep
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "bytecode.h"
#include "ec_arith_cost.h"
#include "facul_ecm.h"
#include "getprime.h"   // getprime
#include "pm1.h"
#include "pp1.h"
#include "stage2.h"         // for stage2_clear_plan, stage2_make_plan, stag...
#include "verbose.h"             // verbose_output_print
#include "macros.h"

/* By default we compress the bytecode chain.
 * You can set it to 0 to disabled it.
 */
#define BYTECODE_COMPRESS 1

/* Costs of operations for P+1 (for Lucas chain computed with PRAC algorithm) */
prac_cost_t pp1_opcost = { .DBL = 1, .dADD = 1. };

/* Costs of operations for interpreting PRAC chains on Montgomery curves
 * The cost are the number of modular multiplications and squarings.
 */
prac_cost_t ec_montgomery_opcost = { .DBL = MONTGOMERY_DBL,
                                     .dADD = MONTGOMERY_dADD };

/* Costs of operations for interpreting MISHMASH bytecode on "a=-1" twisted
 * Edwards curves and Montgomery curves.
 * The cost are the number of modular multiplications and squarings.
 *
 * For "a=-1" twisted Edwards curves, we use 2 different coordinates system:
 *    - projective: (X:Y:Z)
 *    - extended: (X:Y:Z:T) and XY = ZT
 *      Note: if we forget T, we have projective coordinates
 * The costs corresponds to operations on different coordinates system:
 *    - DBL: doubling projective -> projective
 *    - DBLa: doubling projective -> extended
 *    - TPL: tripling projective -> projective
 *    - TPLa: tripling projective -> extended
 *    - ADD: addition extended,extended -> projective
 *    - ADDa: addition extended,extended -> extended
 *    - ADDd: addition extended,extended -> Montgomery curve
 *
 *  dADD and dDBL correspond to differential addition and doubling on a
 *  Montgomery curve.
 */
mishmash_cost_t ec_mixed_repr_opcost = {
    .DBL = EDWARDS_DBL, .DBLa = EDWARDS_DBLext,
    .TPL = EDWARDS_TPL, .TPLa = EDWARDS_TPLext,
    .ADD = EDWARDS_ADD, .ADDa = EDWARDS_ADDext, .ADDd = EDWARDS_ADDmontgomery,
    .dDBL = MONTGOMERY_DBL, .dADD = MONTGOMERY_dADD
  };

/***************************** P-1 ********************************************/

void
pm1_make_plan (pm1_plan_t *plan, const unsigned int B1, const unsigned int B2,
	       int verbose)
{
  mpz_t E;
  unsigned int p;
  size_t tmp_E_nrwords;

  verbose_output_print(0, 3, "# make plan for P-1 with B1=%u, B2=%u\n", B1, B2);

  /* Generate the exponent for stage 1 */
  plan->exp2 = 0;
  for (p = 1; p <= B1 / 2; p *= 2)
    plan->exp2++;

  plan->B1 = B1;
  mpz_init (E);
  mpz_set_ui (E, 1UL);
  prime_info pi;
  prime_info_init (pi);
  p = (unsigned int) getprime_mt (pi);
  ASSERT (p == 3);
  for ( ; p <= B1; p = (unsigned int) getprime_mt (pi))
    {
      unsigned long q;
      /* Uses p^k s.t. (p-1)p^(k-1) <= B1, except for p=2 because our
         base 2 is a QR for primes == 1 (mod 8) already */
      for (q = 1; q <= B1 / (p - 1); q *= p)
        mpz_mul_ui (E, E, p);
    }
  prime_info_clear (pi);

  if (verbose)
    gmp_printf ("pm1_make_plan: E = %Zd;\n", E);

  plan->E = mpz_export (NULL, &tmp_E_nrwords, -1, sizeof(unsigned long),
                        0, 0, E);
  plan->E_nrwords = (unsigned int) tmp_E_nrwords;
  mpz_clear (E);
  /* Find highest set bit in E. */
  ASSERT (plan->E[plan->E_nrwords - 1] != 0);
  plan->E_mask = ~0UL - (~0UL >> 1); /* Only MSB set */
  while ((plan->E[plan->E_nrwords - 1] & plan->E_mask) == 0UL)
    plan->E_mask >>= 1;

  /* stage2 is done with P+1 code */
  stage2_cost_t stage2_opcost = { .dadd = pp1_opcost.dADD,
                                  .dbl = pp1_opcost.DBL, .is_ecm = 0 };
  stage2_make_plan (&(plan->stage2), B1, B2, &stage2_opcost, verbose);
}


void
pm1_clear_plan (pm1_plan_t *plan)
{
  stage2_clear_plan (&(plan->stage2));

  free (plan->E);
  plan->E = NULL;
  plan->E_nrwords = 0;
  plan->B1 = 0;
}

/***************************** P+1 ********************************************/

/* Make byte code for addition chain for stage 1, and the parameters for
   stage 2 */

void
pp1_make_plan (pp1_plan_t *plan, const unsigned int B1, const unsigned int B2,
	       int verbose)
{
  verbose_output_print(0, 3, "# make plan for P+1 with B1=%u, B2=%u\n", B1, B2);

  /* Make bytecode for stage 1 */
  plan->exp2 = 0;
  for (unsigned int p = 1; p <= B1 / 2; p *= 2)
    plan->exp2++;

  plan->B1 = B1;

  /* Allocates plan->bc, fills it and returns bc_len */
  bytecode_prac_encode (&(plan->bc), B1, plan->exp2, 0, &pp1_opcost,
                        BYTECODE_COMPRESS, verbose);

  /* Make stage 2 plan */
  stage2_cost_t stage2_opcost = { .dadd = pp1_opcost.dADD,
                                  .dbl = pp1_opcost.DBL, .is_ecm = 0 };
  stage2_make_plan (&(plan->stage2), B1, B2, &stage2_opcost, verbose);
}

void
pp1_clear_plan (pp1_plan_t *plan)
{
  stage2_clear_plan (&(plan->stage2));
  free (plan->bc);
  plan->bc = NULL;
  plan->B1 = 0;

  /* Clear the cache for the chains computed by prac algorithm.
   * The first call to this function will free the memory, the other calls will
   * be nop.
   */
  bytecode_prac_cache_free ();
}

/***************************** ECM ********************************************/

/* Generate bytecode for stage 1, and the parameters for stage 2.
   The curve used in ECM will be computed using the given 'parameterization'
   with the given 'parameter' value.
   "extra_primes" controls whether some primes should be added (or left out!)
   on top of the primes and prime powers <= B1, for example to take into
   account the known factors in the group order. */

void
ecm_make_plan (ecm_plan_t *plan, const unsigned int B1, const unsigned int B2,
               const ec_parameterization_t parameterization,
               const unsigned long parameter, const int extra_primes,
               const int verbose)
{
  unsigned int pow3_extra;

  verbose_output_print(0, 3, "# make plan for ECM with B1=%u, B2=%u, "
                             "parameterization = %d, " "parameter = %lu, "
                             "extra primes = %d\n", B1, B2, parameterization,
                             parameter, extra_primes);

  /* If group order is divisible by 12 or 16, add two or four 2s to stage 1 */
  if (extra_primes)
  {
    if (parameterization & ECM_TORSION16)
      plan->exp2 = 4;
    else if (parameterization & ECM_TORSION12)
      plan->exp2 = 2;
  }
  else
    plan->exp2 = 0;
  for (unsigned int q = 1; q <= B1 / 2; q *= 2)
    plan->exp2++;
  if (verbose)
    printf ("Exponent of 2 in stage 1 primes: %u\n", plan->exp2);

  /* Make bytecode for stage 1 */
  plan->B1 = B1;
  plan->parameterization = parameterization;
  plan->parameter = parameter;

  /* If group order is divisible by 12, add another 3 to stage 1 primes */
  pow3_extra = (extra_primes && (parameterization & ECM_TORSION12)) ? 1 : 0 ;
  if (verbose && pow3_extra)
    printf ("Add another 3 to stage 1 primes\n");

  /* Allocates plan->bc, fills it and returns bc_len.
   * For Montgomery curves: we use Lucas chains computed with PRAC algorithm.
   */
  if (parameterization & FULLMONTY)
  {
    bytecode_prac_encode (&(plan->bc), B1, plan->exp2, pow3_extra,
        &ec_montgomery_opcost, BYTECODE_COMPRESS, verbose);
  }
  else if (parameterization & FULLMONTYTWED)
  {
    bytecode_mishmash_encode (&(plan->bc), B1, plan->exp2, pow3_extra,
                            &ec_mixed_repr_opcost, BYTECODE_COMPRESS, verbose);
  }
  else
    FATAL_ERROR_CHECK (1, "Unknown parameterization");

  /* Make stage 2 plan */
  stage2_cost_t stage2_opcost = { .dadd = ec_montgomery_opcost.dADD,
                                 .dbl = ec_montgomery_opcost.DBL, .is_ecm = 1 };
  stage2_make_plan (&(plan->stage2), B1, B2, &stage2_opcost, verbose);
}


void
ecm_clear_plan (ecm_plan_t *plan)
{
  stage2_clear_plan (&(plan->stage2));
  plan->B1 = 0;

  if (plan->bc != NULL)
  {
    free (plan->bc);
    plan->bc = NULL;
  }

  /* Clear the cache for the chains computed by prac algorithm.
   * The first call to this function will free the memory, the other calls will
   * be nop.
   */
  bytecode_prac_cache_free ();
}
