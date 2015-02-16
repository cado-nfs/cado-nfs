#include "cado.h"
#include <stdint.h>     /* AIX wants it first (it's a bug) */
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "utils.h" /* for getprime() */
#include "pm1.h"
#include "pp1.h"
#include "ecm.h"
#include "prac_bc.h"
#include "portability.h"

#define PP1_DICT_NRENTRIES 6
static size_t pp1_dict_len[PP1_DICT_NRENTRIES] = {1, 1, 2, 2, 3, 4};
static literal_t *pp1_dict_entry[PP1_DICT_NRENTRIES] = 
  {"\xB", "\xA", "\xB\xA", "\x3\x0", "\x3\xB\xA", "\x3\x0\x3\x0"};
static code_t pp1_dict_code[PP1_DICT_NRENTRIES] = {0, 0, 10, 11, 13, 14};

static bc_dict_t pp1_dict = 
  {PP1_DICT_NRENTRIES, pp1_dict_len, pp1_dict_entry, pp1_dict_code};


#define ECM_DICT_NRENTRIES 4
static size_t ecm_dict_len[ECM_DICT_NRENTRIES] = {1, 1, 2, 2};
static literal_t *ecm_dict_entry[ECM_DICT_NRENTRIES] = 
{"\xB", "\xA", "\xB\xA", "\x3\x0"};
static code_t ecm_dict_code[ECM_DICT_NRENTRIES] = {0, 0, 10, 11};

static bc_dict_t ecm_dict = 
  {ECM_DICT_NRENTRIES, ecm_dict_len, ecm_dict_entry, ecm_dict_code};


void 
pm1_make_plan (pm1_plan_t *plan, const unsigned int B1, const unsigned int B2,
	       int verbose)
{
  mpz_t E;
  unsigned int p;
  size_t tmp_E_nrwords;
  
  if (verbose)
    printf("Making plan for P-1 with B1=%u, B2=%u\n", B1, B2);
  /* Generate the exponent for stage 1 */
  plan->exp2 = 0;
  for (p = 1; p <= B1 / 2; p *= 2)
    plan->exp2++;

  plan->B1 = B1;
  mpz_init (E);
  mpz_set_ui (E, 1UL);
  p = (unsigned int) getprime (p);
  ASSERT (p == 3);
  for ( ; p <= B1; p = (unsigned int) getprime (p))
    {
      unsigned long q;
      /* Uses p^k s.t. (p-1)p^(k-1) <= B1, except for p=2 because our 
         base 2 is a QR for primes == 1 (mod 8) already */
      for (q = 1; q <= B1 / (p - 1); q *= p)
        mpz_mul_ui (E, E, p);
    }
  
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
  getprime (0);
  
  stage2_make_plan (&(plan->stage2), B1, B2, verbose);
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


/* Make byte code for addition chain for stage 1, and the parameters for 
   stage 2 */
                    
void 
pp1_make_plan (pp1_plan_t *plan, const unsigned int B1, const unsigned int B2,
	       enum pp1_parameterization_t parameterization, int verbose)
{
  unsigned int p;
  const double addcost = 10, doublecost = 10, bytecost = 1, changecost = 1;
  const unsigned int compress = 1;
  bc_state_t *bc_state;
  const char *parameterization_names[2] = {"2/7", "6/5"};
  
  if (verbose)
    printf("Making plan for P+1 with B1=%u, B2=%u, starting point = %s\n",
           B1, B2, parameterization_names[parameterization]);

  plan->parameterization = parameterization;

  /* Make bytecode for stage 1 */
  plan->exp2 = 0;
  for (p = 1; p <= B1 / 2; p *= 2)
    plan->exp2++;

  plan->B1 = B1;
  bc_state = bytecoder_init (compress ? &pp1_dict : NULL);
  p = (unsigned int) getprime (p);
  ASSERT (p == 3);
  for ( ; p <= B1; p = (unsigned int) getprime (p))
    {
      unsigned long q;
      for (q = 1; q <= B1 / p; q *= p)
	prac_bytecode (p, addcost, doublecost, bytecost, changecost, bc_state);
    }
  bytecoder ((literal_t) 12, bc_state);
  bytecoder_flush (bc_state);
  plan->bc_len = bytecoder_size (bc_state);
  plan->bc = (char *) malloc (plan->bc_len);
  ASSERT (plan->bc != NULL);
  bytecoder_read (plan->bc, bc_state);
  bytecoder_clear (bc_state);
  getprime (0);

  if (!compress)
    {
      /* The very first chain init and very last chain end are hard-coded
	 in the stage 1 code and must be removed from the byte code. */
      size_t i;
      ASSERT (plan->bc[0] == 10); /* check that first code is chain init */
      ASSERT (plan->bc[plan->bc_len - 2] == 11); /* check that next-to-last 
						    code is chain end */
      /* check that last code is bytecode end */
      ASSERT (plan->bc[plan->bc_len - 1] == (literal_t) 12);
      /* Remove first code 10 and last code 11 */
      for (i = 1; i < plan->bc_len; i++)
	plan->bc[i - 1] = plan->bc[i];
      plan->bc[plan->bc_len - 3] = plan->bc[plan->bc_len - 2];
      plan->bc_len -= 2;
    }

  if (verbose)
    {
      int changes = 0;
      printf ("Byte code for stage 1: ");
      for (p = 0; p < plan->bc_len; p++)
        {
	  printf ("%s%d", (p == 0) ? "" : ", ", (int) (plan->bc[p]));
	  changes += (p > 0 && plan->bc[p-1] != plan->bc[p]);
        }
      printf ("\n");
      printf ("Length %d, %d code changes\n", plan->bc_len, changes);
    }
    
  /* Make stage 2 plan */
  stage2_make_plan (&(plan->stage2), B1, B2, verbose);
}

void 
pp1_clear_plan (pp1_plan_t *plan)
{
  stage2_clear_plan (&(plan->stage2));
  free (plan->bc);
  plan->bc = NULL;
  plan->bc_len = 0;
  plan->B1 = 0;
}


/* Make byte code for addition chain for stage 1, and the parameters for 
   stage 2. Parameterization chooses Brent-Suyama curves with order divisible
   by 12 (BRENT12), Montgomery with torsion 12 over Q (MONTY12) or Montgomery
   with torsion 16 over Q (MONTY16), sigma is the associated parameter.
   "extra_primes" controls whether some primes should be added (or left out!)
   on top of the primes and prime powers <= B1, for example to take into 
   account the known factors in the group order. */

void 
ecm_make_plan (ecm_plan_t *plan, const unsigned int B1, const unsigned int B2,
	       const int parameterization, const unsigned long sigma, 
	       const int extra_primes, const int verbose)
{
  unsigned int p, q;
  const double addcost = 6., doublecost = 5., /* TODO: find good ratio */
    bytecost = 1, changecost = 1;
  const unsigned int compress = 1;
  bc_state_t *bc_state;
  double totalcost = 0.;
  
  if (verbose)
    printf("Making plan for ECM with B1=%u, B2=%u, parameterization = %d, "
            "sigma=%lu, extra primes = %d\n", 
            B1, B2, parameterization, sigma, extra_primes);

  /* If group order is divisible by 12 or 16, add two or four 2s to stage 1 */
  if (extra_primes)
    plan->exp2 = (parameterization == MONTY16) ? 4 : 2;
  else
    plan->exp2 = 0;
  for (q = 1; q <= B1 / 2; q *= 2)
    plan->exp2++;
  totalcost += plan->exp2 * doublecost;
  
  /* Make bytecode for stage 1 */
  plan->B1 = B1;
  plan->parameterization = parameterization;
  if (parameterization & FULLMONTY)
    {
      plan->sigma = sigma;
      bc_state = bytecoder_init (compress ? &ecm_dict : NULL);
      p = (unsigned int) getprime (2UL);
      ASSERT (p == 3);
      /* If group order is divisible by 12, add another 3 to stage 1 primes */
      if (extra_primes && 
	  (parameterization == BRENT12 || parameterization == MONTY12))
	totalcost += prac_bytecode (3, addcost, doublecost, bytecost, 
				    changecost, bc_state);
      for ( ; p <= B1; p = (unsigned int) getprime (p))
	{
	  for (q = 1; q <= B1 / p; q *= p)
	    totalcost += prac_bytecode (p, addcost, doublecost, bytecost, 
					changecost, bc_state);
	}
      bytecoder ((literal_t) 12, bc_state);
      bytecoder_flush (bc_state);
      plan->bc_len = bytecoder_size (bc_state);
      plan->bc = (char *) malloc (plan->bc_len);
      ASSERT (plan->bc);
      bytecoder_read (plan->bc, bc_state);
      bytecoder_clear (bc_state);
      getprime (0);
      
      if (!compress)
	{
	  /* The very first chain init and very last chain end are hard-coded
	     in the stage 1 code and must be removed from the byte code. */
	  size_t i;
	  ASSERT (plan->bc[0] == 10); /* check that first code is chain init */
	  ASSERT (plan->bc[plan->bc_len - 2] == 11); /* check that next-to-last 
							code is chain end */
	  /* check that last code is bytecode end */
	  ASSERT (plan->bc[plan->bc_len - 1] == (literal_t) 12);
	  /* Remove first code 10 and last code 11 */
	  for (i = 1; i < plan->bc_len; i++)
	    plan->bc[i - 1] = plan->bc[i];
	  plan->bc[plan->bc_len - 3] = plan->bc[plan->bc_len - 2];
	  plan->bc_len -= 2;
	}

      if (verbose)
	{
	  int changes = 0;
	  printf ("Exponent of 2 in stage 1 primes: %u\n", plan->exp2);
	  printf ("Byte code for stage 1: ");
	  for (p = 0; p < plan->bc_len; p++)
	    {
	      printf ("%s%d", (p == 0) ? "" : ", ", (int) (plan->bc[p]));
	      changes += (p > 0 && plan->bc[p-1] != plan->bc[p]);
	    }
	  printf ("\n");
	  printf ("Length %d, %d code changes, total cost: %f\n", 
		  plan->bc_len, changes, totalcost);
	}
    }
  if (parameterization == TWED16)
    plan->E = &Ecurve14;
  
  /* Make stage 2 plan */
  stage2_make_plan (&(plan->stage2), B1, B2, verbose);
}
  

void 
ecm_clear_plan (ecm_plan_t *plan)
{
  stage2_clear_plan (&(plan->stage2));
  plan->B1 = 0;

  if (plan->parameterization & FULLMONTY)
    {
      free (plan->bc);
      plan->bc = NULL;
      plan->bc_len = 0;
    }
  
  if (plan->parameterization == TWED16)
    {
      plan->E = NULL;
    }
}


/* void */
/* ecmE_make_plan (ecmE_plan_t *plan, const unsigned int B1, const int parameterization) */
/* { */
/*   fprintf(stderr, "In ecmE_make_plan\n"); */

/*   plan->exp2 = 0; */
/*   plan->B1 = B1; */
/*   plan->parameterization = parameterization; */
/*   plan->E = &Ecurve14; */
/* } */

/* void  */
/* ecmE_clear_plan (ecmE_plan_t *plan MAYBE_UNUSED) */
/* { */
/* } */
