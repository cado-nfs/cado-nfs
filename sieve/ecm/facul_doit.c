#include "cado.h"
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
  modulusredcul_t fm_ul, cfm_ul;
#if     MOD_MAXBITS > MODREDCUL_MAXBITS
  modulusredc15ul_t fm_15ul, cfm_15ul;
#endif
#if     MOD_MAXBITS > MODREDC15UL_MAXBITS
  modulusredc2ul2_t fm_2ul2, cfm_2ul2; /* Modulus for factor and cofactor */
#endif
#if     MOD_MAXBITS > MODREDC2UL2_MAXBITS
  modulusredc3ul_t fm_3ul, cfm_3ul; /* Modulus for factor and cofactor */
#endif
  int i, found = 0, bt, fprime, cfprime;
  enum {
      CHOOSE_UL,
#if     MOD_MAXBITS > MODREDCUL_MAXBITS
      CHOOSE_15UL,
#endif
#if     MOD_MAXBITS > MODREDC15UL_MAXBITS
      CHOOSE_2UL2,
#endif
#if     MOD_MAXBITS > MODREDC2UL2_MAXBITS
      CHOOSE_3UL,
#endif
  } f_arith = CHOOSE_UL, cf_arith = CHOOSE_UL;
  
  mod_getmod_uls (n, m);
  mod_intset_ul (f, 1UL);
  mod_init (r, m);
  
  for (i = method_start; ; ++i)
    {
      if (strategy->methods[i].method == 0) {
          // We run out of methods in the strategy. 
          // We failed to factor this number !
          found = FACUL_NOT_FOUND;
//          printf("# Warning: Running out of methods\n");
          break;
      }

      /* Simple-minded early abort for large input.
         Note: before the test was "mod_intbits (n) > LONG_BIT" which was
         machine-dependent. However it would be better if the early abort
         test depends on the size of the number we are trying to factor,
         since for a large number we can invest more in cofactorization. */
//      if (i > 3 && mod_intbits (n) > 64)
//        break;

      // Basic early abort strategy:
      //   - let b_purged the bit-size of the largest primes that have
      //   been removed at this point with prob > 90%. 
      //   - If 2*b_purged > bit_size of remaining cofactor, then abort.
      //   - If a linear combination of cofac_size and b_purged > some_bound,
      //     then abort, where "some_bound" must be larger than
      //     mfb + fbb (start of cofactorization) and larger than
      //     lpb + ecmb (otherwise ecmb does not make sense). 
      //
      if (strategy->early_abort) {
          int nn = mod_intbits (n);
          int bb = strategy->purged_bits[i];

          if (2*bb > nn)
              break;

          int lpb = strategy->lpb_bits;
          int fbb = strategy->fbb_bits;
          int mfb = strategy->mfb;
          int ecmb = strategy->ecmb;
//          printf("lpb=%d fbb=%d mfb=%d ecmb=%d\n", lpb, fbb, mfb, ecmb);

          int delta_n = (ecmb - fbb) / 2;
          int delta_b = (mfb - lpb) / 2;
          // Ugly formula for the linear combination that fits the
          // conditions above with a margin of delta_n and delta_b
          // respetcively.
          // (ecmb-fbb-delta_n)*nn+(-lpb-delta_b+mfb)*bb
          //   +delta_n*(lpb+delta_b)+fbb*(lpb+delta_b) -mfb*ecmb
          int crit = (ecmb-fbb-delta_n)*nn+(-lpb-delta_b+mfb)*bb
              +delta_n*(lpb+delta_b)+fbb*(lpb+delta_b) -mfb*ecmb;
          if (crit >= 0) {
//              printf("Abort: nn=%d bb=%d\n", nn, bb);
              break;
          }
      }
      
      if (i < STATS_LEN)
	  stats_called[i]++;
      
      if (strategy->methods[i].method == PM1_METHOD)
	bt = pm1 (f, m, (pm1_plan_t *) (strategy->methods[i].plan));
      else if (strategy->methods[i].method == PP1_27_METHOD)
	bt = pp1_27 (f, m, (pp1_plan_t *) (strategy->methods[i].plan));
      else if (strategy->methods[i].method == PP1_65_METHOD)
	bt = pp1_65 (f, m, (pp1_plan_t *) (strategy->methods[i].plan));
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
      if (fprime && mod_intcmp (f, strategy->lpb) > 0)
	{
	  found = FACUL_NOT_SMOOTH; /* A prime > lpb, not smooth */
	  break;
	}
      
      /* Compute the cofactor */
      mod_intdivexact (n, n, f);
      
      /* See if cofactor is <= fbb^2 and > lpb */
      cfprime = (mod_intcmp (n, strategy->fbb2) <= 0);
      if (cfprime && mod_intcmp (n, strategy->lpb) > 0)
	{
	  found = FACUL_NOT_SMOOTH; /* A prime > lpb, not smooth */
	  break;
	}
      
      /* Determine for certain if the factor is prime */
      if (!fprime)
	{
	  if (mod_intbits (f) <= MODREDCUL_MAXBITS)
	    {
	      f_arith = CHOOSE_UL;
	      modredcul_initmod_uls (fm_ul, f);
	      fprime = primetest_ul (fm_ul);
            }
#if     MOD_MAXBITS > MODREDCUL_MAXBITS
          else if (mod_intbits (f) <= MODREDC15UL_MAXBITS)
            {
              f_arith = CHOOSE_15UL;
	      modredc15ul_initmod_uls (fm_15ul, f);
	      fprime = primetest_15ul (fm_15ul);
            }
#endif
#if     MOD_MAXBITS > MODREDC15UL_MAXBITS
	  else if (mod_intbits (f) <= MODREDC2UL2_MAXBITS)
            {
              f_arith = CHOOSE_2UL2;
	      modredc2ul2_initmod_uls (fm_2ul2, f);
	      fprime = primetest_2ul2 (fm_2ul2);
            }
#endif
#if     MOD_MAXBITS > MODREDC2UL2_MAXBITS
	  else if (mod_intbits (f) <= MODREDC3UL_MAXBITS)
            {
              f_arith = CHOOSE_3UL;
	      modredc3ul_initmod_uls (fm_3ul, f);
	      fprime = primetest_3ul (fm_3ul);
            }
#endif
          else
              abort();

	  if (fprime && mod_intcmp (f, strategy->lpb) > 0)
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
	      cf_arith = CHOOSE_UL;
	      modredcul_initmod_uls (cfm_ul, n);
	      cfprime = primetest_ul (cfm_ul);
            }
#if     MOD_MAXBITS > MODREDCUL_MAXBITS
	  else if (mod_intbits (n) <= MODREDC15UL_MAXBITS)
	    {
	      cf_arith = CHOOSE_15UL;
	      modredc15ul_initmod_uls (cfm_15ul, n);
	      cfprime = primetest_15ul (cfm_15ul);
            }
#endif
#if     MOD_MAXBITS > MODREDC15UL_MAXBITS
	  else if (mod_intbits (n) <= MODREDC2UL2_MAXBITS)
            {
              cf_arith = CHOOSE_2UL2;
	      modredc2ul2_initmod_uls (cfm_2ul2, n);
	      cfprime = primetest_2ul2 (cfm_2ul2);
            }
#endif
#if     MOD_MAXBITS > MODREDC2UL2_MAXBITS
	  else if (mod_intbits (n) <= MODREDC3UL_MAXBITS)
            {
              cf_arith = CHOOSE_3UL;
	      modredc3ul_initmod_uls (cfm_3ul, n);
	      cfprime = primetest_3ul (cfm_3ul);
            }
#endif
          else
            abort ();

	  if (cfprime && mod_intcmp (n, strategy->lpb) > 0)
	    {
	      found = FACUL_NOT_SMOOTH; /* A prime > lpb, not smooth */
	      break;
	    }
	}
      
      /* So each of factor and cofactor is either a prime < lpb, 
	 or is composite */

      if (fprime) {
        // We push only factors less than MAX_LONG.
        // We expect that large factors are just cofactors, and therefore
        // that there is only one which will be found by caller.
        // FIXME: return mpz_t's, here
#if MOD_MAXBITS > LONG_BIT
        if (f[1] == 0)
#endif
          factors[found++] = f[0]; 
      }
      else
	{
            int f2 = FACUL_NOT_SMOOTH;    /* placate gcc (!) */
	  /* Factor the composite factor. Use the same method again so that
	     backtracking can separate the factors */
          switch (f_arith) {
              case CHOOSE_UL:
                  f2 = facul_doit_ul (factors + found, fm_ul, strategy, i);
                  break;
#if     MOD_MAXBITS > MODREDCUL_MAXBITS
              case CHOOSE_15UL:
                  f2 = facul_doit_15ul (factors + found, fm_15ul, strategy, i);
                  break;
#endif
#if     MOD_MAXBITS > MODREDC15UL_MAXBITS
              case CHOOSE_2UL2:
                  f2 = facul_doit_2ul2 (factors + found, fm_2ul2, strategy, i);
                  break;
#endif
#if     MOD_MAXBITS > MODREDC2UL2_MAXBITS
              case CHOOSE_3UL:
                  f2 = facul_doit_3ul (factors + found, fm_3ul, strategy, i);
                  break;
#endif
          }
          
	  if (f2 == FACUL_NOT_SMOOTH)
	    {
	      found = FACUL_NOT_SMOOTH;
	      break;
	    }
	  found += f2;
	}

      if (cfprime) {
        // We push only factors less than MAX_LONG.
        // We expect that large factors are just cofactors, and therefore
        // that there is only one which will be found by caller.
        // FIXME: return mpz_t's, here
#if MOD_MAXBITS > LONG_BIT
        if (n[1] == 0)
#endif
          factors[found++] = n[0]; 
      }
      else
	{
	  int f2 = FACUL_NOT_SMOOTH;    /* placate gcc (!) */
	  /* Factor the composite cofactor */
          switch(cf_arith) {
              case CHOOSE_UL:
                  f2 = facul_doit_ul (factors + found, cfm_ul, strategy, i + 1);
                  break;
#if     MOD_MAXBITS > MODREDCUL_MAXBITS
              case CHOOSE_15UL:
                  f2 = facul_doit_15ul (factors + found, cfm_15ul, strategy, i + 1);
                  break;
#endif    
#if     MOD_MAXBITS > MODREDC15UL_MAXBITS
              case CHOOSE_2UL2:
                  f2 = facul_doit_2ul2 (factors + found, cfm_2ul2, strategy, i + 1);
                  break;
#endif
#if     MOD_MAXBITS > MODREDC2UL2_MAXBITS
              case CHOOSE_3UL:
                  f2 = facul_doit_3ul (factors + found, cfm_3ul, strategy, i + 1);
                  break;
#endif
          }
          
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
