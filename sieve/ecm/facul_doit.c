#include "facul_doit.h"
#include "pm1.h"
#include "pp1.h"
#include "ecm.h"

extern unsigned long stats_called[];
extern unsigned long stats_found_n[];

int
primetest (const modulus_t m)
{
  residue_t one, r;
  int isprime;
  
  isprime = mod_sprp2 (m);
  if (isprime)
    {
      mod_init_noset0 (one, m);
      mod_init_noset0 (r, m);
      mod_set1 (one, m);
      mod_add (r, one, one, m);
      mod_add (r, r, one, m);   /* r = 3 */
      isprime = mod_sprp (r, m);
      mod_clear (one, m);
      mod_clear (r, m);
    }
  
  return isprime;
}

int
facul_doit (unsigned long *factors, const modulus_t m, 
	    const facul_strategy_t *strategy, const int method_start)
{
  residue_t r;
  modint_t n, f;
  modulusredc2ul2_t fm_2ul2, cfm_2ul2; /* Modulus for factor and cofactor */
  modulusredc15ul_t fm_15ul, cfm_15ul;
  modulusredcul_t fm_ul, cfm_ul;
  int i, found = 0, bt, fprime, cfprime, f_arith = 0, cf_arith = 0;
  
  mod_getmod_uls (n, m);
  mod_intset_ul (f, 1UL);
  mod_init (r, m);
  
  for (i = method_start; strategy->methods[i].method != 0; i++)
    {
      /* Simple-minded early abort for large input */
      if (i > 3 && mod_intbits (n) > LONG_BIT)
        break;
      
      if (i < STATS_LEN)
	  stats_called[i]++;
      
      if (strategy->methods[i].method == PM1_METHOD)
	bt = pm1 (f, m, (pm1_plan_t *) (strategy->methods[i].plan));
      else if (strategy->methods[i].method == PP1_METHOD)
	bt = pp1 (f, m, (pp1_plan_t *) (strategy->methods[i].plan));
      else if (strategy->methods[i].method == EC_METHOD)
	bt = ecm (f, m, (ecm_plan_t *) (strategy->methods[i].plan));
      else 
	{
	  /* A method value we don't know about. Something's wrong, bail out */
	  found = -1;
	  break;
	}
      
      /* The following possibilities exist:
	 bt:   Factor:    Cofactor:   Action:
	 0           1    composite   Try next method
	 1           1    composite   Could try again with careful bt
	 0    prime<lp     prime<lp   Store both, exit successfully
	 0    prime>lp     prime<lp   Not smooth, exit NOT_SMOOTH
	 0    prime<lp     prime>lp   Not smooth, exit NOT_SMOOTH
	 0    prime>lp     prime>lp   Not smooth, exit NOT_SMOOTH
	 1    prime<lp     prime<lp   Store both, exit successfully
	 1    prime>lp     prime<lp   Not smooth, exit NOT_SMOOTH
	 1    prime<lp     prime>lp   Not smooth, exit NOT_SMOOTH
	 1    prime>lp     prime>lp   Not smooth, exit NOT_SMOOTH
	 0    prime<lp    composite   Store prime, continue with cofactor
	 0    prime>lp    composite   Not smooth, exit NOT_SMOOTH
	 1    prime<lp    composite   Store prime, try same method with cofactor
	 1    prime>lp    composite   Not smooth, exit NOT_SMOOTH
	 0   composite     prime<lp   Store prime, continue next method
	 0   composite     prime>lp   Not smooth, exit NOT_SMOOTH
	 1   composite     prime<lp   Store prime, retry this method
	 1   composite     prime>lp   Not smooth, exit NOT_SMOOTH
	 0   composite            1   Could try with lower bounds
	 1   composite            1   Could try again with careful bt

	 Simplified:

	 bt:   Factor:    Cofactor:   Action:
	 0           1    composite   Try next method
	 1           1    composite   Could try again with careful bt
	 ?    prime<lp     prime<lp   Store both, exit successfully
	 ?    prime>lp            ?   Not smooth, exit NOT_SMOOTH
	 ?           ?     prime>lp   Not smooth, exit NOT_SMOOTH
	 0    prime<lp    composite   Store prime, continue with cofactor
	 1    prime<lp    composite   Store prime, same method with cofactor
	 0   composite     prime<lp   Store prime, continue next method
	 1   composite     prime<lp   Store prime, retry this method
	 0   composite            1   Could try with lower bounds
	 1   composite            1   Could try again with careful bt
	 
      */
      
      
      if (mod_intequal_ul (f, 1UL))
	{
	  if (bt == 0)
	    {
	      /* No factor found, no backtracking... this was a simple miss. */
	      continue;
	    }
	  else
	    {
	      /* Backtracking was used, so there was a point where all
		 factors had been found simultaneously, but backing up
		 to the previous checkpoint resulted in no factors being
		 found. We could try to do some more clever backtracking 
		 to discover the factors yet. TODO. For now, just continue
		 to the next method. */
	      continue;
	    }
	}
      
      if (mod_intequal (f, n))
	{
	  if (i < STATS_LEN)
	    stats_found_n[i]++;
	  if (bt == 0)
	    {
	      /* Input number was found without any backtracking happening?
		 Find out when this can occur and how to get a chance of
		 finding the factors yet. TODO. */
	      continue;
	    }
	  else
	    {
	      /* Backtracking was used, but could not separate the factors,
	         e.g. if both factors are found in stage 1 without 
		 multiplying/exponentiating by 2 at all. Better backtracking
		 might recover the factors yet. TODO. */
	      continue;
	    }
	}
      
      /* So we found a non-trivial factor. See if it is prime, if the 
	 cofactor is prime, and if one of them is, whether they are too
	 large for our smoothness bounds */
      
      /* A quick test if the factor is <= fbb^2 and >lpb */
      /* FIXME: must always use same width for comparison */
      fprime = (mod_intcmp (f, strategy->fbb2) <= 0); 
      if (fprime && mod_intcmp_ul (f, strategy->lpb) > 0)
	{
	  found = FACUL_NOT_SMOOTH; /* A prime > lpb, not smooth */
	  break;
	}
      
      /* Compute the cofactor */
      mod_intdivexact (n, n, f);
      
      /* See if cofactor is <= fbb^2 and > lpb */
      cfprime = (mod_intcmp (n, strategy->fbb2) <= 0);
      if (cfprime && mod_intcmp_ul (n, strategy->lpb) > 0)
	{
	  found = FACUL_NOT_SMOOTH; /* A prime > lpb, not smooth */
	  break;
	}
      
      /* Determine for certain if the factor is prime */
      if (!fprime)
	{
	  if (mod_intbits (f) <= MODREDCUL_MAXBITS)
	    {
	      f_arith = 0;
	      modredcul_initmod_uls (fm_ul, f);
	      fprime = primetest_ul (fm_ul);
            }
          else if (mod_intbits (f) <= MODREDC15UL_MAXBITS)
            {
              f_arith = 1;
	      modredc15ul_initmod_uls (fm_15ul, f);
	      fprime = primetest_15ul (fm_15ul);
            }
          else
            {
              f_arith = 2;
	      modredc2ul2_initmod_uls (fm_2ul2, f);
	      fprime = primetest_2ul2 (fm_2ul2);
            }

	  if (fprime && mod_intcmp_ul (f, strategy->lpb) > 0)
	    {
	      found = FACUL_NOT_SMOOTH; /* A prime > lpb, not smooth */
	      break;
	    }
	}
      
      /* Determine for certain if the cofactor is prime */
      if (!cfprime)
	{
	  if (mod_intbits (n) <= MODREDCUL_MAXBITS)
	    {
	      cf_arith = 0;
	      modredcul_initmod_uls (cfm_ul, n);
	      cfprime = primetest_ul (cfm_ul);
            }
	  else if (MOD_SIZE >= MODREDC15UL_SIZE && 
	           mod_intbits (n) <= MODREDC15UL_MAXBITS)
	    {
	      cf_arith = 1;
	      modredc15ul_initmod_uls (cfm_15ul, n);
	      cfprime = primetest_15ul (cfm_15ul);
            }
          else if (MOD_SIZE >= MODREDC2UL2_SIZE)
            {
              cf_arith = 2;
	      modredc2ul2_initmod_uls (cfm_2ul2, n);
	      cfprime = primetest_2ul2 (cfm_2ul2);
            }
          else
            abort ();

	  if (cfprime && mod_intcmp_ul (n, strategy->lpb) > 0)
	    {
	      found = FACUL_NOT_SMOOTH; /* A prime > lpb, not smooth */
	      break;
	    }
	}
      
      /* So each of factor and cofactor is either a prime < lpb, 
	 or is composite */

      if (fprime)
	factors[found++] = f[0]; /* f < lp, so it fits in 1 unsigned long */
      else
	{
	  int f2;
	  /* Factor the composite factor. Use the same method again so that
	     backtracking can separate the factors */
	  if (f_arith == 0)
	    f2 = facul_doit_ul (factors + found, fm_ul, strategy, i);
	  else if (f_arith == 1)
	    f2 = facul_doit_15ul (factors + found, fm_15ul, strategy, i);
          else
            f2 = facul_doit_2ul2 (factors + found, fm_2ul2, strategy, i);
          
	  if (f2 == FACUL_NOT_SMOOTH)
	    {
	      found = FACUL_NOT_SMOOTH;
	      break;
	    }
	  found += f2;
	}

      if (cfprime)
	factors[found++] = n[0]; /* n < lp, so it fits in 1 unsigned long */
      else
	{
	  int f2;
	  /* Factor the composite cofactor */
	  if (cf_arith == 0)
	    f2 = facul_doit_ul (factors + found, cfm_ul, strategy, i + 1);
	  else if (cf_arith == 1)
	    f2 = facul_doit_15ul (factors + found, cfm_15ul, strategy, i + 1);
          else
            f2 = facul_doit_2ul2 (factors + found, cfm_2ul2, strategy, i + 1);
          
	  if (f2 == FACUL_NOT_SMOOTH)
	    {
	      found = FACUL_NOT_SMOOTH;
	      break;
	    }
	  found += f2;
	}
      
      /* We found a non-trivial factorization and any composite 
	 factors/cofactors have been treated in recursive calls, 
	 so we can stop here */
      break;
    }
  
  mod_clear (r, m);
  return found;
}
