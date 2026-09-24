#ifndef CADO_FACUL_ECM_H
#define CADO_FACUL_ECM_H

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

void ecm_make_plan (ecm_plan_t *, unsigned int, unsigned int, 
		    ec_parameterization_t, unsigned long, int, int);

void ecm_clear_plan (ecm_plan_t *);

int ec_parameter_is_valid (ec_parameterization_t, unsigned long);
unsigned long
ec_valid_parameter_from_sequence (ec_parameterization_t parameterization,
                       unsigned long sequence_value);

#ifdef __cplusplus
}
#endif

#endif	/* CADO_FACUL_ECM_H */
