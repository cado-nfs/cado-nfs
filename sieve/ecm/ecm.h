#include "modredc_ul.h"
#include "modredc_15ul.h"
#include "modredc_2ul2.h"
#include "modredc_3ul.h"
#include "stage2.h"

#define BRENT12 0
#define MONTY12 1
#define MONTY16 2

typedef struct {
  char *bc;             /* Bytecode for the Lucas chain for stage 1 */
  unsigned int bc_len;  /* Number of bytes in bytecode */
  unsigned int exp2;    /* Exponent of 2 in stage 1 primes */
  unsigned int B1;
  int parameterization; /* BRENT12 or MONTY12 */
  unsigned long sigma;  /* Sigma parameter for Brent curves, or
			   multiplier for Montgomery torsion-12 curves */
  stage2_plan_t stage2;
} ecm_plan_t;


int ecm_ul (modintredcul_t, const modulusredcul_t, const ecm_plan_t *);

unsigned long ell_pointorder_ul (const residueredcul_t, const int, 
                                 const unsigned long, const unsigned long,
                                 const modulusredcul_t, const int);
unsigned long ellM_curveorder_jacobi_ul (residueredcul_t, residueredcul_t, 
                                         modulusredcul_t);
int ecm_15ul (modintredc15ul_t, const modulusredc15ul_t, const ecm_plan_t *);
unsigned long ell_pointorder_15ul (const residueredc15ul_t, const int, 
                                   const unsigned long, const unsigned long,
                                   const modulusredc15ul_t, const int);
unsigned long ellM_curveorder_jacobi_15ul (residueredc15ul_t, residueredc15ul_t, 
                                           modulusredc15ul_t);
int ecm_2ul2 (modintredc2ul2_t, const modulusredc2ul2_t, const ecm_plan_t *);
unsigned long ell_pointorder_2ul2 (const residueredc2ul2_t, const int, 
                                   const unsigned long, const unsigned long,
                                   const modulusredc2ul2_t, const int);
unsigned long ellM_curveorder_jacobi_2ul2 (residueredc2ul2_t, residueredc2ul2_t, 
                                           modulusredc2ul2_t);
int ecm_3ul (modintredc3ul_t, const modulusredc3ul_t, const ecm_plan_t *);
unsigned long ell_pointorder_3ul (const residueredc3ul_t, const int, 
                                   const unsigned long, const unsigned long,
                                   const modulusredc3ul_t, const int);
unsigned long ellM_curveorder_jacobi_3ul (residueredc3ul_t, residueredc3ul_t, 
                                           modulusredc3ul_t);

void ecm_make_plan (ecm_plan_t *, const unsigned int, const unsigned int, 
		    const int, const unsigned long, const int, const int);
void ecm_clear_plan (ecm_plan_t *);
