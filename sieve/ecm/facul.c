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

/* A chain constructed with basicstrat.c, with fbb=18 */
#if 0
#define CHAIN_LEN 40
static const int STRATEGY_CHAIN[CHAIN_LEN][3] = {
  { PM1_METHOD    ,  125, 1785 }, 
  { PP1_27_METHOD ,  180, 3885 }, 
  { EC_METHOD    ,  125, 3045 }, 
  { PM1_METHOD    ,  385, 11445 }, 
  { EC_METHOD    ,  125, 6195 }, 
  { EC_METHOD    ,  150, 6195 }, 
  { PP1_65_METHOD ,  815, 20475 }, 
  { EC_METHOD    ,  220, 9345 }, 
  { EC_METHOD    ,  265, 11445 }, 
  { EC_METHOD    ,  385, 16905 }, 
  { EC_METHOD    ,  385, 16905 }, 
  { PM1_METHOD    ,  2460, 62895 }, 
  { EC_METHOD    ,  385, 16905 }, 
  { EC_METHOD    ,  385, 16905 }, 
  { EC_METHOD    ,  385, 16905 }, 
  { EC_METHOD    ,  465, 20475 }, 
  { EC_METHOD    ,  465, 20475 }, 
  { EC_METHOD    ,  465, 20475 }, 
  { EC_METHOD    ,  465, 20475 }, 
  { EC_METHOD    ,  675, 29925 }, 
  { EC_METHOD    ,  675, 29925 }, 
  { EC_METHOD    ,  675, 29925 }, 
  { PP1_27_METHOD ,  4260, 131565 }, 
  { EC_METHOD    ,  815, 36015 }, 
  { EC_METHOD    ,  815, 36015 }, 
  { EC_METHOD    ,  815, 36015 }, 
  { EC_METHOD    ,  815, 36015 }, 
  { EC_METHOD    ,  815, 36015 }, 
  { EC_METHOD    ,  815, 36015 }, 
  { EC_METHOD    ,  815, 36015 }, 
  { EC_METHOD    ,  815, 36015 }, 
  { EC_METHOD    ,  815, 36015 }, 
  { EC_METHOD    ,  815, 36015 }, 
  { EC_METHOD    ,  815, 36015 }, 
  { EC_METHOD    ,  815, 36015 }, 
  { EC_METHOD    ,  815, 36015 }, 
  { EC_METHOD    ,  980, 43365 }, 
  { EC_METHOD    ,  980, 43365 }, 
  { EC_METHOD    ,  980, 43365 }, 
  { EC_METHOD    ,  980, 43365 }, 
};
#else
#define CHAIN_LEN 60
static const int STRATEGY_CHAIN[CHAIN_LEN][3] = {
{ PM1_METHOD,   180,1000 },
{ PP1_27_METHOD,        220,7665 },
{ PM1_METHOD,   560,13965 },
{ EC_METHOD,    125,6195 },
{ EC_METHOD,    150,6195 },
{ EC_METHOD,    220,9345 },
{ EC_METHOD,    220,9345 },
{ PP1_65_METHOD,        980,29925 },
{ EC_METHOD,    385,16905 },
{ EC_METHOD,    385,16905 },
{ EC_METHOD,    385,16905 },
{ EC_METHOD,    385,16905 },
{ EC_METHOD,    385,16905 },
{ PM1_METHOD,   2955,75705 },
{ EC_METHOD,    465,20475 },
{ EC_METHOD,    465,20475 },
{ EC_METHOD,    465,20475 },
{ EC_METHOD,    465,20475 },
{ EC_METHOD,    675,29925 },
{ EC_METHOD,    675,29925 },
{ EC_METHOD,    675,29925 },
{ EC_METHOD,    675,29925 },
{ EC_METHOD,    675,29925 },
{ EC_METHOD,    815,36015 },
{ EC_METHOD,    815,36015 },
{ EC_METHOD,    815,36015 },
{ PP1_27_METHOD,        4260,131565 },
{ EC_METHOD,    815,36015 },
{ EC_METHOD,    815,36015 },
{ EC_METHOD,    815,36015 },
{ EC_METHOD,    815,36015 },
{ EC_METHOD,    815,36015 },
{ EC_METHOD,    815,36015 },
{ EC_METHOD,    815,36015 },
{ EC_METHOD,    815,36015 },
{ EC_METHOD,    980,43365 },
{ EC_METHOD,    980,43365 },
{ EC_METHOD,    980,43365 },
{ EC_METHOD,    980,43365 },
{ EC_METHOD,    980,43365 },
{ EC_METHOD,    980,43365 },
{ EC_METHOD,    980,43365 },
{ EC_METHOD,    980,43365 },
{ EC_METHOD,    980,43365 },
{ EC_METHOD,    980,43365 },
{ EC_METHOD,    980,43365 },
{ EC_METHOD,    980,43365 },
{ EC_METHOD,    980,43365 },
{ EC_METHOD,    980,43365 },
{ EC_METHOD,    980,43365 },
{ EC_METHOD,    980,43365 },
{ EC_METHOD,    980,43365 },
{ EC_METHOD,    980,43365 },
{ EC_METHOD,    980,43365 },
{ EC_METHOD,    980,43365 },
{ EC_METHOD,    980,43365 },
{ EC_METHOD,    980,43365 },
{ EC_METHOD,    980,43365 },
{ EC_METHOD,    980,43365 },
{ EC_METHOD,    980,43365 },
    };
#endif


/* Try a strategy for large lpb's. 
   In that case, n is ignored. This is not a bug. It is more the
   job of this function to decide the number of curves to run.
   */
static facul_strategy_t *
facul_make_strategy_large (const int MAYBE_UNUSED n, const unsigned long fbb, 
		     const unsigned int lpb)
{
  facul_strategy_t *strategy;
  facul_method_t *methods;
  int i;

  // the number of curves is nb_curves[lpb-35] 
  // TODO: Change that!!!!
  int nb_curves[60] = {
      12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17, 17, 18, 18,
      19, 19, 20, 20, 21, 21, 22, 22, 23, 23, 24, 24, 25, 25,
      26, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 32, 32,
  };

  int nn;
  if (lpb < 35) 
      nn = nb_curves[0];
  else if (lpb > 74)
      nn = nb_curves[39];
  else
      nn = nb_curves[lpb - 35];
  // FIXME: for the moment, activate always maximal number of curves.
  nn = CHAIN_LEN;

  strategy = malloc (sizeof (facul_strategy_t));
  strategy->lpb_bits = lpb;
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

  methods = malloc ((1+nn) * sizeof (facul_method_t));
  strategy->methods = methods;

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

  return strategy;
}


/* Make a simple minded strategy for factoring. We start with P-1 and
   P+1 (with x0=2/7), then an ECM curve with low bounds, then a bunch of
   ECM curves with larger bounds. How many methods to do in total is
   controlled by the n parameter: P-1, P+1 and the first ECM curve
   (with small bounds) are always done, then n ECM curves (with larger bounds)
*/

facul_strategy_t *
facul_make_strategy (const int n, const unsigned long fbb, 
		     const unsigned int lpb)
{
  if (lpb >= 35)
    return facul_make_strategy_large(n, fbb, lpb);
  facul_strategy_t *strategy;
  facul_method_t *methods;
  int i;
  
  strategy = malloc (sizeof (facul_strategy_t));
  strategy->lpb_bits = lpb;
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

  methods = malloc ((n + 4) * sizeof (facul_method_t));
  strategy->methods = methods;

  /* run one P-1 curve with B1=315 and B2=2205 */
  methods[0].method = PM1_METHOD;
  methods[0].plan = malloc (sizeof (pm1_plan_t));
  pm1_make_plan (methods[0].plan, 315, 2205, 0);

  /* run one P+1 curve with B1=525 and B2=3255 */
  methods[1].method = PP1_27_METHOD;
  methods[1].plan = malloc (sizeof (pp1_plan_t));
  pp1_make_plan (methods[1].plan, 525, 3255, 0);

  /* run one ECM curve with Montgomery parametrization, B1=105, B2=3255 */
  methods[2].method = EC_METHOD;
  methods[2].plan = malloc (sizeof (ecm_plan_t));
  ecm_make_plan (methods[2].plan, 105, 3255, MONTY12, 2, 1, 0);
  
  if (n > 0)
    {
      methods[3].method = EC_METHOD;
      methods[3].plan = malloc (sizeof (ecm_plan_t));
      ecm_make_plan (methods[3].plan, 315, 5355, BRENT12, 11, 1, 0);
    }

  for (i = 4; i < n + 3; i++)
    {
      methods[i].method = EC_METHOD;
      methods[i].plan = malloc (sizeof (ecm_plan_t));
      ecm_make_plan (methods[i].plan, 315, 5355, MONTY12, i - 1, 1, 0);
    }

  methods[n + 3].method = 0;
  methods[n + 3].plan = NULL;

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
