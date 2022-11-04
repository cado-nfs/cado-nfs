#ifndef ECM_H_
#define ECM_H_

#include "modredc_ul.h"
#include "modredc_15ul.h"
#include "modredc_2ul2.h"
#include "mod_mpz.h"
#include "bytecode.h"
#include "stage2.h"

typedef enum {
  BRENT12 = 1,
  MONTY12 = 2,
  MONTY16 = 4,
  MONTYTWED12 = 8,
  MONTYTWED16 = 16      // some code bits are missing
} ec_parameterization_t;

#define FULLMONTY (BRENT12 | MONTY12 | MONTY16)
#define FULLMONTYTWED (MONTYTWED12 | MONTYTWED16)
#define ECM_TORSION16 (MONTY16 | MONTYTWED16)
#define ECM_TORSION12 (BRENT12 | MONTY12 | MONTYTWED12)


typedef struct {
  bytecode bc;          /* Bytecode for stage 1 */
  unsigned int exp2;    /* Exponent of 2 in stage 1 primes */
  unsigned int B1;
  ec_parameterization_t parameterization;
  unsigned long parameter;  /* Used to compute curve coefficients with the
                             * above parameterization. */
  stage2_plan_t stage2;
} ecm_plan_t;


#ifdef __cplusplus
extern "C" {
#endif

int ecm_ul (modintredcul_t, const modulusredcul_t, const ecm_plan_t *);
int ecm_15ul (modintredc15ul_t, const modulusredc15ul_t, const ecm_plan_t *);
int ecm_2ul2 (modintredc2ul2_t, const modulusredc2ul2_t, const ecm_plan_t *);
int ecm_mpz (modintmpz_t, const modulusmpz_t, const ecm_plan_t *);

unsigned long ec_parameterization_point_order_ul (const ec_parameterization_t,
                                                  const unsigned long,
                                                  const unsigned long,
                                                  const unsigned long,
                                                  const modulusredcul_t,
                                                  const int);
unsigned long ec_parameterization_curve_order_ul (const ec_parameterization_t,
                                                  const unsigned long,
                                                  const modulusredcul_t);

void ecm_make_plan (ecm_plan_t *, const unsigned int, const unsigned int, 
		    const ec_parameterization_t, const unsigned long, const int, const int);

void ecm_clear_plan (ecm_plan_t *);

int ec_parameter_is_valid (ec_parameterization_t, const unsigned long);
unsigned long
ec_valid_parameter_from_sequence (ec_parameterization_t parameterization,
                       const unsigned long sequence_value);

#ifdef __cplusplus
}
#endif

#endif	/* ECM_H_ */
