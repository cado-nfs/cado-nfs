/* Factors integers with P-1, P+1 and ECM. Input is in an mpz_t, 
   factors are unsigned long. Returns number of factors found, 
   or -1 in case of error. */

#include "cado.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "pm1.h"
#include "pp1.h"
#include "ecm.h"
#include "facul.h"
#include "facul_doit.h"
#include "mod_ul.h"

/* These global variables are only for statistics. In case of
 * multithreaded sieving, the stats might be wrong...
 */

unsigned long stats_called[STATS_LEN] = {
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
};
unsigned long stats_found_n[STATS_LEN] = {
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
};

// A chain constructed with basicstrat of length 60.
// Second and third integers are B1,B2.
// Fourth integers gives the bit size of the primes that have been purged
// at this stage of the chain, with prob > 90%.
#define CHAIN_LEN 60
static const int STRATEGY_CHAIN[CHAIN_LEN][4] = {
{ PM1_METHOD,   180,1000, 14 },
{ PP1_27_METHOD,        220,7665, 17 },
{ PM1_METHOD,   560,13965, 20 },
{ EC_METHOD,    125,6195, 22 },
{ EC_METHOD,    150,6195, 24 },
{ EC_METHOD,    220,9345, 25 },
{ EC_METHOD,    220,9345, 26 },
{ PP1_65_METHOD,        980,29925, 27 },
{ EC_METHOD,    385,16905, 28 },
{ EC_METHOD,    385,16905, 29 },
{ EC_METHOD,    385,16905, 30 },
{ EC_METHOD,    385,16905, 30 },
{ EC_METHOD,    385,16905, 31 },
{ PM1_METHOD,   2955,75705, 31 },
{ EC_METHOD,    465,20475, 32 },
{ EC_METHOD,    465,20475, 32 },
{ EC_METHOD,    465,20475, 33 },
{ EC_METHOD,    465,20475, 33 },
{ EC_METHOD,    675,29925, 34 },
{ EC_METHOD,    675,29925, 34 },
{ EC_METHOD,    675,29925, 35 },
{ EC_METHOD,    675,29925, 35 },
{ EC_METHOD,    675,29925, 35 },
{ EC_METHOD,    815,36015, 36 },
{ EC_METHOD,    815,36015, 36 },
{ EC_METHOD,    815,36015, 36 },
{ PP1_27_METHOD,        4260,131565, 37 },
{ EC_METHOD,    815,36015, 37 },
{ EC_METHOD,    815,36015, 37 },
{ EC_METHOD,    815,36015, 38 },
{ EC_METHOD,    815,36015, 38 },
{ EC_METHOD,    815,36015, 38 },
{ EC_METHOD,    815,36015, 39 },
{ EC_METHOD,    815,36015, 39 },
{ EC_METHOD,    815,36015, 39 },
{ EC_METHOD,    980,43365, 39 },
{ EC_METHOD,    980,43365, 39 },
{ EC_METHOD,    980,43365, 39 },
{ EC_METHOD,    980,43365, 40 },
{ EC_METHOD,    980,43365, 40 },
{ EC_METHOD,    980,43365, 40 },
{ EC_METHOD,    980,43365, 40 },
{ EC_METHOD,    980,43365, 40 },
{ EC_METHOD,    980,43365, 40 },
{ EC_METHOD,    980,43365, 41 },
{ EC_METHOD,    980,43365, 41 },
{ EC_METHOD,    980,43365, 41 },
{ EC_METHOD,    980,43365, 41 },
{ EC_METHOD,    980,43365, 41 },
{ EC_METHOD,    980,43365, 41 },
{ EC_METHOD,    980,43365, 42 },
{ EC_METHOD,    980,43365, 42 },
{ EC_METHOD,    980,43365, 42 },
{ EC_METHOD,    980,43365, 42 },
{ EC_METHOD,    980,43365, 42 },
{ EC_METHOD,    980,43365, 42 },
{ EC_METHOD,    980,43365, 42 },
{ EC_METHOD,    980,43365, 42 },
{ EC_METHOD,    980,43365, 42 },
{ EC_METHOD,    980,43365, 43 },
    };


/* Try a strategy for large lpb's. 
   In that case, n is ignored. This is not a bug. It is more the
   job of this function to decide the number of curves to run.
   */
facul_strategy_t *
facul_make_strategy(const unsigned long fbb, const unsigned int lpb, 
        const unsigned int mfb, const unsigned int ecmb)
{
  facul_strategy_t *strategy;
  facul_method_t *methods;
  int i;

  int nn;
  if ((int)ecmb < STRATEGY_CHAIN[CHAIN_LEN-1][3]) {
      // The ecm bound is smaller that what our strategy chain
      // can found. Select the appropriate subchain:
      nn = 0;
      while (STRATEGY_CHAIN[nn][3] < (int)ecmb) 
          nn++;
  } else {
      // Otherwise, take the full chain.
      nn = CHAIN_LEN;
  }
  printf("Chose nn = %d\n", nn);

  strategy = malloc (sizeof (facul_strategy_t));
  strategy->early_abort = 1;
  strategy->lpb_bits = lpb;
  strategy->mfb = mfb;
  strategy->ecmb = MIN(STRATEGY_CHAIN[nn-1][3], (int)ecmb);

  if (lpb < LONG_BIT) {
      strategy->lpb[0] = 1UL << lpb;
      strategy->lpb[1] = 0UL;
  } else {
      int ll = lpb - LONG_BIT;
      if (ll >= LONG_BIT) {
          fprintf(stderr, "Sorry, this large prime bound is too large\n");
          abort();
      }
      strategy->lpb[0] = 0UL;
      strategy->lpb[1] = 1UL << ll;
  }
  /* Store fbb^2 in fbb2 */
  ularith_mul_ul_ul_2ul (&(strategy->fbb2[0]), &(strategy->fbb2[1]), fbb, fbb);
  /* Store fbb as a number of bits */
  strategy->fbb_bits = modul_intbits(&fbb);

  methods = malloc ((1+nn) * sizeof (facul_method_t));
  strategy->methods = methods;
  unsigned int * purged = malloc ((1+nn) * sizeof (unsigned int));
  strategy->purged_bits = purged;

  for (i = 0; i < nn; ++i) {
      int type = STRATEGY_CHAIN[i][0];
      int B1   = STRATEGY_CHAIN[i][1];
      int B2   = STRATEGY_CHAIN[i][2];
      methods[i].method = type;
      if (type == PM1_METHOD) {
          methods[i].plan = malloc (sizeof (pm1_plan_t));
          pm1_make_plan (methods[i].plan, B1, B2, 0);
      } else if (type == PP1_27_METHOD || type == PP1_65_METHOD) {
          methods[i].plan = malloc (sizeof (pp1_plan_t));
          pp1_make_plan (methods[i].plan, B1, B2, 0);
      } else {
          methods[i].plan = malloc (sizeof (ecm_plan_t));
          ecm_make_plan (methods[i].plan, B1, B2, MONTY12, i+2, 1, 0);
      }
      purged[i] = MAX(STRATEGY_CHAIN[i][3], (int)strategy->fbb_bits);
  }

  // Use Brent sigma=11 for the first ECM
  i = 0;
  while (methods[i].method != EC_METHOD)
      i++;
  ecm_clear_plan(methods[i].plan);
  ecm_make_plan (methods[i].plan, STRATEGY_CHAIN[i][1],
          STRATEGY_CHAIN[i][2], BRENT12, 11, 1, 0);

  // Sentinel
  methods[nn].method = 0;
  methods[nn].plan = NULL;
  purged[nn] = 0;

  return strategy;
}


void 
facul_clear_strategy (facul_strategy_t *strategy)
{
  facul_method_t *methods = strategy->methods;
  int i = 0;

  for (i = 0; methods[i].method != 0; i++)
    {
      if (methods[i].method == PM1_METHOD)
        pm1_clear_plan (methods[i].plan);
      else if (methods[i].method == PP1_27_METHOD
              || methods[i].method == PP1_65_METHOD)
	pp1_clear_plan (methods[i].plan);
      else if (methods[i].method == EC_METHOD)
	ecm_clear_plan (methods[i].plan);
      methods[i].method = 0;
      free (methods[i].plan);
      methods[i].plan = NULL;
    }
  free (methods);
  methods = NULL;
  free (strategy->purged_bits);
  free (strategy);
}

static int
cmp_ul (const unsigned long *a, const unsigned long *b)
{
  if (*a < *b) return -1;
  if (*a == *b) return 0;
  return 1;
}


void facul_print_stats (FILE *stream)
{
  int i, notfirst;
  unsigned long sum;

  fprintf (stream, "# facul statistics.\n# histogram of methods called: ");
  notfirst = 0;
  sum = 0;
  for (i = 0; i < STATS_LEN; i++)
    {
      sum += stats_called[i];
      if (stats_called[i] > 0UL)
	fprintf (stream, "%s %d: %lu", 
		 (notfirst++) ? ", " : "", i, stats_called[i]);
    }
  fprintf (stream, ". Total: %lu\n", sum);

  fprintf (stream, "# histogram of input numbers found: ");
  notfirst = 0;
  sum = 0;
  for (i = 0; i < STATS_LEN; i++)
    {
      sum += stats_found_n[i];
      if (stats_found_n[i] > 0UL)
	fprintf (stream, "%s %d: %lu", 
		 (notfirst++) ? ", " : "", i, stats_found_n[i]);
    }
  fprintf (stream, ". Total: %lu\n", sum);
}


int
facul (unsigned long *factors, const mpz_t N, const facul_strategy_t *strategy)
{
  modintredc3ul_t n;
  int i, found = 0;
  
#ifdef PARI
  gmp_fprintf (stderr, "%Zd", N);
#endif

  if (mpz_sgn (N) <= 0)
    return -1;
  if (mpz_cmp_ui (N, 1UL) == 0)
    return 0;
  
  /* If the composite does not fit into our modular arithmetic, return
     no factor */
  if (mpz_sizeinbase (N, 2) > MODREDC3UL_MAXBITS)
    return 0;
  
  {
    size_t written;
    mpz_export (n, &written, -1, sizeof(unsigned long), 0, 0, N);
    for (i = written; i < MODREDC3UL_SIZE; i++)
      n[i] = 0UL;
  }
  
  /* Use the fastest modular arithmetic that's large enough for this input */
  i = modredc3ul_intbits (n);
  if (i <= MODREDCUL_MAXBITS)
    {
      modulusredcul_t m;
      modredcul_initmod_uls (m, n);
      found = facul_doit_ul (factors, m, strategy, 0);
      modredcul_clearmod (m);
    }
  else if (i <= MODREDC15UL_MAXBITS)
    {
      modulusredc15ul_t m;
      modredc15ul_initmod_uls (m, n);
      found = facul_doit_15ul (factors, m, strategy, 0);
      modredc15ul_clearmod (m);
    }
  else if (i <= MODREDC2UL2_MAXBITS)
    {
      modulusredc2ul2_t m;
      modredc2ul2_initmod_uls (m, n);
      found = facul_doit_2ul2 (factors, m, strategy, 0);
      modredc2ul2_clearmod (m);
    }
  else 
    {
      if (i < MODREDC3UL_MINBITS)  // 127 and 128 bits are not covered!
      {
//          printf("# Warning: can not factor numbers of 127 and 128 bits!\n");
          return 0; 
      }
      modulusredc3ul_t m;
      ASSERT (i <= MODREDC3UL_MAXBITS);
      modredc3ul_initmod_uls (m, n);
      found = facul_doit_3ul (factors, m, strategy, 0);
      modredc3ul_clearmod (m);
    }
  
  if (found > 1)
    {
      /* Sort the factors we found */
      qsort (factors, found, sizeof (unsigned long), 
	     (int (*)(const void *, const void *)) &cmp_ul);
    }

#ifdef PARI
  if (found > 1)
    {
      fprintf (stderr, " == ");
      for (i = 0; i < found; i++)
	fprintf (stderr, "%lu%s", factors[i], 
		 (i+1 < found) ? " * " : " /* PARI */\n");
    }
  else
    fprintf (stderr, "; /* PARI */\n");
#endif

  return found;
}
