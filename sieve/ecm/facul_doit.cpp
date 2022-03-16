/* This file is _NOT_ a standalone compilation unit. It is included by
 * facul_doit_blah.cpp ; therefore it's _normal_ if this file does not
 * seem to be self-contained.
 */

/* We don't even need to #include "cado.h", because the
 * facul_doit_blah.cpp has to do so already. However, our policy
 * detection script will insist on seeing this include, and it doesn't do
 * any harm
 */
#include "cado.h" // IWYU pragma: keep
#include <stdlib.h>

#ifndef FACUL_DOIT_READY_TO_INCLUDE_IMPL_CODE
#error "This file must not be used as a standalone compilation unit"
#endif

#include "facul_doit.hpp"
#include "modset.hpp"
#include "pm1.h"
#include "pp1.h"
#include "facul_ecm.h"
#include "mpqs.h"
#include "cxx_mpz.hpp"
#include "mod_ul.h"
#include "macros.h"

#ifdef ENABLE_UNSAFE_FACUL_STATS
extern unsigned long stats_called[];
extern unsigned long stats_found_n[];
extern int stats_current_index;
#endif

#define mod_intget_mpz MOD_RENAME(int_get_mpz)
void
mod_intget_mpz(mpz_t z, const modint_t x) {
#ifdef MOD_SIZE
    mpz_import(z, MOD_SIZE, -1, sizeof(unsigned long), 0, 0, x);
#else
    mpz_set(z, x);
#endif
}

#define mod_intget_cxx_mpz MOD_RENAME(int_get_cxx_mpz)
cxx_mpz mod_intget_cxx_mpz(const modint_t x) {
    cxx_mpz c;
    mod_intget_mpz(c, x);
    return c;
}

static inline const FaculModulusBase *
modset_init (modint_t m)
{
  return FaculModulusBase::MOD_APPEND_TYPE(init)(m);
}

int
facul_doit (std::vector<cxx_mpz> & factors, const modulus_t m, 
	    facul_strategy_oneside const & strategy, const int method_start)
{
  if (method_start >= (int) strategy.methods.size())
      return 0;

  modint_t n, f;
  const FaculModulusBase *fm = NULL, *cfm = NULL;
  int i, found = 0, bt, fprime, cfprime;
  
  mod_intinit (n);
  mod_intinit (f);
  mod_getmod_int (n, m);
  mod_intset_ul (f, 1UL);
  
  for (i = method_start; i < (int) strategy.methods.size() ; i++)
    {
      /* Simple-minded early abort for large input.
         Note: before the test was "mod_intbits (n) > LONG_BIT" which was
         machine-dependent. However it would be better if the early abort
         test depends on the size of the number we are trying to factor,
         since for a large number we can invest more in cofactorization. */
#if 0
      if (i > 3 && mod_intbits (n) > 64)
        break;
#endif

#ifdef ENABLE_UNSAFE_FACUL_STATS
      if (i < STATS_LEN)
	stats_called[i]++;
#endif
      
      if (strategy.methods[i].method == PM1_METHOD)
	bt = pm1 (f, m, (pm1_plan_t *) (strategy.methods[i].plan));
      else if (strategy.methods[i].method == PP1_27_METHOD)
	bt = pp1_27 (f, m, (pp1_plan_t *) (strategy.methods[i].plan));
      else if (strategy.methods[i].method == PP1_65_METHOD)
	bt = pp1_65 (f, m, (pp1_plan_t *) (strategy.methods[i].plan));
      else if (strategy.methods[i].method == EC_METHOD)
	bt = ecm (f, m, (ecm_plan_t *) (strategy.methods[i].plan));	
      else if (strategy.methods[i].method == MPQS_METHOD)
	bt = mpqs (f, m);
      else 
	{
	  /* A method value we don't know about. Something's wrong, bail out */
	  abort();
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
#ifdef ENABLE_UNSAFE_FACUL_STATS
	  if (i < STATS_LEN)
	    stats_found_n[i]++;
#endif  /* ENABLE_UNSAFE_FACUL_STATS */

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
      
      /* A quick test if the factor is <= fbb^2 and >2^lpb */
      double f_dbl = mod_intget_double (f);
      fprime = f_dbl < strategy.BB;
      if (fprime && mod_intbits (f) > strategy.lpb)
	{
	  found = FACUL_NOT_SMOOTH; /* A prime > 2^lpb, not smooth */
	  break;
	}

      /* if L^2 < f < B^3, it cannot be smooth */
      if (2 * strategy.lpb < mod_intbits (f) && f_dbl < strategy.BBB)
        {
          found = FACUL_NOT_SMOOTH;
          break;
        }
      
      /* Compute the cofactor */
      mod_intdivexact (n, n, f);
      
      /* See if cofactor is <= fbb^2 and > 2^lpb */
      double n_dbl = mod_intget_double (n);
      cfprime = n_dbl < strategy.BB;
      if (cfprime && mod_intbits (n) > strategy.lpb)
	{
	  found = FACUL_NOT_SMOOTH; /* A prime > 2^lpb, not smooth */
	  break;
	}

      if (2 * strategy.lpb < mod_intbits (n) &&
          n_dbl < strategy.BBB)
        {
          found = FACUL_NOT_SMOOTH;
          break;
        }

      /* Determine for certain if the factor is prime */
      if (!fprime)
	{
          fm = modset_init (f);
	  fprime = fm->isprime ();
          if (fprime) {
            delete fm;
            fm = NULL;
          }
	  if (fprime && mod_intbits (f) > strategy.lpb)
	    {
	      found = FACUL_NOT_SMOOTH; /* A prime > 2^lpb, not smooth */
	      break;
	    }
	}
      
      /* Determine for certain if the cofactor is prime */
      if (!cfprime)
	{
	  cfm = modset_init (n);
	  cfprime = cfm->isprime ();

          if (cfprime) {
            delete cfm;
            cfm = NULL;
          }
	  if (cfprime && mod_intbits (n) > strategy.lpb)
	    {
	      if (!fprime) {
	        delete fm;
                fm = NULL;
              }
	      found = FACUL_NOT_SMOOTH; /* A prime > 2^lpb, not smooth */
	      break;
	    }
	}
      
      /* So each of factor and cofactor is either a prime < 2^lpb,
	 or is composite */

      if (fprime) {
          factors.push_back(mod_intget_cxx_mpz(f));
          found++;
	}
      else
	{
          int found2 = FACUL_NOT_SMOOTH;    /* placate gcc (!) */
	  /* Factor the composite factor. Use the same method again so that
	     backtracking can separate the factors */
          found2 = fm->facul_doit (factors, strategy, i);
	  delete fm;
          fm = NULL;
	  if (found2 == FACUL_NOT_SMOOTH) {
            found = FACUL_NOT_SMOOTH;
            if (!cfprime) {
              delete cfm;
              cfm = NULL;
            }
            break;
          }
          found += found2;
	}
      
      if (cfprime) {
          factors.push_back(mod_intget_cxx_mpz(n));
          found++;
        }
      else
	{
	  int found2 = FACUL_NOT_SMOOTH;    /* placate gcc (!) */
	  /* Factor the composite cofactor */
	  found2 = cfm->facul_doit (factors, strategy, i + 1);
          delete cfm;
          cfm = NULL;
	  if (found2 == FACUL_NOT_SMOOTH)
	    {
	      found = FACUL_NOT_SMOOTH;
	      break;
	    }
	  found += found2;
	}
      /* We found a non-trivial factorization and any composite 
	 factors/cofactors have been treated in recursive calls, 
	 so we can stop here */
      ASSERT_ALWAYS(fm == NULL);
      ASSERT_ALWAYS(cfm == NULL);
      break;
    }
  
  ASSERT_ALWAYS(fm == NULL);
  ASSERT_ALWAYS(cfm == NULL);
  
  mod_intclear (n);
  mod_intclear (f);
  return found;
}



/*****************************************************************************/
/*                       STRATEGY BOOK                                       */
/*****************************************************************************/

/*
  This function applies one factoring method 'method' on an integer 'm'
  and returns: 
       FACUL_NOT_SMOOTH (-1) if m is not smooth.
       n if we found n prime factors and store them in 'factors'.

  Remark: if m has more than two factors, it's possible that 
  we need to try another factorization on f (or/and  m/f). So
  the values of our composite factor are stored in fm (or/and cfm).

  Note: BB was renamed BB from previously chosen name "assume_prime_thresh". I think it's correct.
*/
int
facul_doit_onefm (std::vector<cxx_mpz> & factors, const modulus_t m,
		  facul_method const & method,
		  const FaculModulusBase* &fm, const FaculModulusBase* &cfm,
                  unsigned long lpb, double BB, double BBB)
{
  residue_t r;
  modint_t n, f;
  int bt, fprime, cfprime;
  int found = 0;
  
  mod_intinit (n);
  mod_intinit (f);
  mod_getmod_int (n, m);
  mod_intset_ul (f, 1UL);
  mod_init (r, m);
  
  if (method.method == PM1_METHOD)
    bt = pm1 (f, m, (pm1_plan_t *) (method.plan));
  else if (method.method == PP1_27_METHOD)
    bt = pp1_27 (f, m, (pp1_plan_t *) (method.plan));
  else if (method.method == PP1_65_METHOD)
    bt = pp1_65 (f, m, (pp1_plan_t *) (method.plan));
  else if (method.method == EC_METHOD)
    bt = ecm (f, m, (ecm_plan_t *) (method.plan));
  else if (method.method == MPQS_METHOD)
    bt = mpqs (f, m);
  else 
    {
      /* A method value we don't know about. Something's wrong, bail out */
      ASSERT_ALWAYS(0);
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

  double f_dbl, n_dbl;

  if (mod_intequal_ul (f, 1UL))
    {
      if (bt == 0)
	{
	  /* No factor found, no backtracking... this was a simple miss. */
	  goto clean_up;
	}
      else
	{
	  /* Backtracking was used, so there was a point where all
	     factors had been found simultaneously, but backing up
	     to the previous checkpoint resulted in no factors being
	     found. We could try to do some more clever backtracking 
	     to discover the factors yet. TODO. For now, just continue
	     to the next method. */
	  goto clean_up;
	}
    }
      
  if (mod_intequal (f, n))
    {
#ifdef ENABLE_UNSAFE_FACUL_STATS
      if (stats_current_index < STATS_LEN)
	stats_found_n[stats_current_index]++;
#endif  /* ENABLE_UNSAFE_FACUL_STATS */
      if (bt == 0)
	{
	  /* Input number was found without any backtracking happening?
	     Find out when this can occur and how to get a chance of
	     finding the factors yet. TODO. */
	  goto clean_up;
	}
      else
	{
	  /* Backtracking was used, but could not separate the factors,
	     e.g. if both factors are found in stage 1 without 
	     multiplying/exponentiating by 2 at all. Better backtracking
	     might recover the factors yet. TODO. */
	  goto clean_up;
	}
    }

  /* So we found a non-trivial factor. See if it is prime, if the 
     cofactor is prime, and if one of them is, whether they are too
     large for our smoothness bounds */
  
  /* A quick test if the factor is <= fbb^2 and >2^lpb */
  f_dbl = mod_intget_double (f);
  fprime = f_dbl < BB;
  if (fprime && mod_intbits (f) > lpb)
    {
      found = FACUL_NOT_SMOOTH; /* A prime > 2^lpb, not smooth */
      goto clean_up;
    }

  /* if L^2 < f < B^3, it cannot be smooth */
  if (2 * lpb < mod_intbits (f) && f_dbl < BBB)
    {
      found = FACUL_NOT_SMOOTH;
      goto clean_up;
    }
    
  /* Compute the cofactor */
  mod_intdivexact (n, n, f);
      
  /* See if cofactor is <= fbb^2 and > 2^lpb */
  n_dbl = mod_intget_double (n);
  cfprime = n_dbl < BB;
  if (cfprime && mod_intbits (n) > lpb)
    {
      found = FACUL_NOT_SMOOTH; /* A prime > 2^lpb, not smooth */
      goto clean_up;
    }
  
  if (2 * lpb < mod_intbits (n) && n_dbl < BBB)
    {
      found = FACUL_NOT_SMOOTH;
      goto clean_up;
    }
  
  /* Determine for certain if the factor is prime */
  if (!fprime)
    {
      fm = modset_init (f);
      fprime = fm->isprime ();
      if (fprime) {
        delete fm;
        fm = NULL;
      }
      if (fprime && mod_intbits (f) > lpb)
	{
	  found = FACUL_NOT_SMOOTH; /* A prime > 2^lpb, not smooth */
	  goto clean_up;
	}
    }
      
  /* Determine for certain if the cofactor is prime */
  if (!cfprime)
    {
      cfm = modset_init (n);
      cfprime = cfm->isprime ();
      if (cfprime) {
        delete cfm;
        cfm = NULL;
      }
      if (cfprime && mod_intbits (n) > lpb)
        {
          if (!fprime) {
            delete fm;
            fm = NULL;
          }
          found = FACUL_NOT_SMOOTH; /* A prime > 2^lpb, not smooth */
          goto clean_up;
        }
    }
      
  /* So each of factor and cofactor is either a prime < 2^lpb, 
     or is composite */

  if (fprime) {
      factors.push_back(mod_intget_cxx_mpz(f));
      found++;
  }

  if (cfprime) {
      factors.push_back(mod_intget_cxx_mpz(n));
      found++;
  }

  /* if either f of cf is composite, it is returned in fm or cfm */

  //Free
 clean_up:
  mod_clear (r, m);
  mod_intclear (n);
  mod_intclear (f);
  return found;
}


