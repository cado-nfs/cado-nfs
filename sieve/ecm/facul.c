/* Factors integers with P-1, P+1 and ECM. Input is in an mpz_t, 
   factors are unsigned long. Returns number of factors found, 
   or -1 in case of error. */

#include "cado.h"
#include <stdint.h>	/* AIX wants it first (it's a bug) */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <regex.h>

#include "utils.h"      /* Happens to define ULONG_BITS, which we use */
#include "portability.h"
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
//for the auxiliary factorization.
unsigned long stats_called_aux[STATS_LEN] = {
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
};
int stats_current_index = 0; //only useful for stats_found_n!
unsigned long stats_found_n[STATS_LEN] = {
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
};

int
nb_curves (const unsigned int lpb)
{
  /* the following table, computed with the proba_cofactor() function in the
     facul.sage file, ensures a probability of at least about 90% to find a
     factor below 2^lpb with n = T[lpb] */
#define LPB_MAX 33
  int T[LPB_MAX+1] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 0-9 */
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 10-19 */
                      0, 0, 1 /*22:0.9074*/, 2 /*23:0.9059*/, 3 /*24:0.8990*/,
                      5 /*25:0.9194*/, 6 /*26:0.9065*/, 8 /*27:0.9053*/,
                      10 /*28:0.9010*/, 13 /*29:0.9091*/, 16 /*30:0.9134*/,
                      18 /*31:0.9039*/, 21 /*32:0.9076*/, 24/*33:0.8963*/};
  return (lpb <= LPB_MAX) ? T[lpb] : T[LPB_MAX];
#undef LPB_MAX
}

// ncurves = 100 corresponds to looking for 50-bit factors.
// ncurves = 200 corresponds to looking for 64-bit factors.


/* Make a simple minded strategy for factoring. We start with P-1 and
   P+1 (with x0=2/7), then an ECM curve with low bounds, then a bunch of
   ECM curves with larger bounds. How many methods to do in total is
   controlled by the n parameter: P-1, P+1 and the first ECM curve
   (with small bounds) are always done, then n ECM curves (with larger bounds)
*/

facul_strategy_t *
facul_make_strategy (const unsigned long fbb, const unsigned int lpb,
		     int n, const int verbose)
{
  facul_strategy_t *strategy;
  facul_method_t *methods;
  int i;

  if (n == -1)
    n = nb_curves (lpb);
  strategy = malloc (sizeof (facul_strategy_t));
  strategy->lpb = lpb;
  /* Store fbb^2 in assume_prime_thresh */
  strategy->assume_prime_thresh = (double) fbb * (double) fbb;
  strategy->BBB = (double) fbb * strategy->assume_prime_thresh;

  methods = malloc ((n + 4) * sizeof (facul_method_t));
  strategy->methods = methods;

  /* run one P-1 curve with B1=315 and B2=2205 */
  methods[0].method = PM1_METHOD;
  methods[0].plan = malloc (sizeof (pm1_plan_t));
  pm1_make_plan (methods[0].plan, 315, 2205, verbose);

  /* run one P+1 curve with B1=525 and B2=3255 */
  methods[1].method = PP1_METHOD;
  methods[1].plan = malloc (sizeof (pp1_plan_t));
  pp1_make_plan (methods[1].plan, 525, 3255, PP1_2_7, verbose);

  /* run one ECM curve with Montgomery parametrization, B1=105, B2=3255 */
  methods[2].method = EC_METHOD;
  methods[2].plan = malloc (sizeof (ecm_plan_t));
  
  ecm_make_plan (methods[2].plan, 105, 3255, MONTY12, 2, 1, verbose);

  if (n > 0)
    {
      methods[3].method = EC_METHOD;
      methods[3].plan = malloc (sizeof (ecm_plan_t));
      ecm_make_plan (methods[3].plan, 315, 5355, BRENT12, 11, 1, verbose);
    }
  /* heuristic strategy where B1 is increased by sqrt(B1) at each curve */
  double B1 = 105.0;
  for (i = 4; i < n + 3; i++)
    {
      double B2;
      unsigned int k;

      B1 += sqrt (B1);
      B2 = 50.0 * B1;
      /* we round B2 to (2k+1)*105, thus k is the integer nearest to
	 B2/210-0.5 */
      k = B2 / 210.0;
      methods[i].method = EC_METHOD;
      methods[i].plan = malloc (sizeof (ecm_plan_t));
      ecm_make_plan (methods[i].plan, (unsigned int) B1, (2 * k + 1) * 105,
		     MONTY12, i - 1, 1, 0);
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
      else if (methods[i].method == PP1_METHOD)
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

  fprintf (stream, "# histogram of auxiliary methods called: ");
  notfirst = 0;
  sum = 0;
  for (i = 0; i < STATS_LEN; i++)
    {
      sum += stats_called_aux[i];
      if (stats_called_aux[i] > 0UL)
	fprintf (stream, "%s %d: %lu", 
		 (notfirst++) ? ", " : "", i, stats_called_aux[i]);
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
  int found = 0;
  size_t bits;
  
#ifdef PARI
  gmp_fprintf (stderr, "%Zd", N);
#endif

  if (mpz_sgn (N) <= 0)
    return -1;
  if (mpz_cmp_ui (N, 1UL) == 0)
    return 0;
  
  /* If the composite does not fit into our modular arithmetic, return
     no factor */
  bits = mpz_sizeinbase (N, 2);
  if (bits > MODMPZ_MAXBITS)
    return 0;
  
  /* Use the fastest modular arithmetic that's large enough for this input */
  if (bits <= MODREDCUL_MAXBITS)
    {
      modulusredcul_t m;
      ASSERT(mpz_fits_ulong_p(N));
      modredcul_initmod_ul (m, mpz_get_ui(N));
      found = facul_doit_ul (factors, m, strategy, 0);
      modredcul_clearmod (m);
    }
  else if (bits <= MODREDC15UL_MAXBITS)
    {
      modulusredc15ul_t m;
      unsigned long t[2];
      modintredc15ul_t n;
      size_t written;
      mpz_export (t, &written, -1, sizeof(unsigned long), 0, 0, N);
      ASSERT_ALWAYS(written <= 2);
      modredc15ul_intset_uls (n, t, written);
      modredc15ul_initmod_int (m, n);
      found = facul_doit_15ul (factors, m, strategy, 0);
      modredc15ul_clearmod (m);
    }
  else if (bits <= MODREDC2UL2_MAXBITS)
    {
      modulusredc2ul2_t m;
      unsigned long t[2];
      modintredc2ul2_t n;
      size_t written;
      mpz_export (t, &written, -1, sizeof(unsigned long), 0, 0, N);
      ASSERT_ALWAYS(written <= 2);
      modredc2ul2_intset_uls (n, t, written);
      modredc2ul2_initmod_int (m, n);
      found = facul_doit_2ul2 (factors, m, strategy, 0);
      modredc2ul2_clearmod (m);
    } 
  else 
    {
      modulusmpz_t m;
      modmpz_initmod_int (m, N);
      found = facul_doit_mpz (factors, m, strategy, 0);
      modmpz_clearmod (m);
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


/*****************************************************************************/
/*                       STRATEGY BOOK                                       */
/*****************************************************************************/

/*
 * If the plan is already precomputed and stored at the index i, return
 * i. Otherwise return the last index of tab to compute and add this new
 * plan at this index!
 */
static int
get_index_method (facul_method_t* tab, unsigned int B1, unsigned int B2,
		  int method, int parameterization, unsigned long sigma)
{
  int i = 0;
  while (tab[i].method != 0){
    if (tab[i].plan == NULL){
      if(B1 == 0 && B2 == 0)
	//zero method!
	break;
      else
	{
	  i++;
	  continue;
	}
    }
    else if (tab[i].method == method) {
      if (method == PM1_METHOD)
	{
	  pm1_plan_t* plan = (pm1_plan_t*)tab[i].plan;
	  if (plan->B1 == B1 && plan->stage2.B2 == B2)
	    break;
	}
      else if (method == PP1_METHOD)
	{
	  pp1_plan_t* plan = (pp1_plan_t*)tab[i].plan;
	  if (plan->B1 == B1 && plan->stage2.B2 == B2 &&
	      (int) plan->parameterization == parameterization)
	    break;
	}
      else if (method == EC_METHOD)
	{
	  ecm_plan_t* plan = (ecm_plan_t*)tab[i].plan;
	  if (plan->B1 == B1 && plan->stage2.B2 == B2 &&
	      plan->parameterization == parameterization
	      && plan->sigma == sigma)
	    break;
	}
    }
    i++;
  }
  return i;
}


static void
return_data_ex (char** res, regmatch_t *pmatch, size_t nmatch,
		const char * str_process)
{
  //(re)init res
  for (size_t i = 0; i < nmatch;i++)
    res[i] = NULL;

  if ( pmatch[0].rm_so != pmatch[0].rm_eo)
    {
      for (size_t i = 1; i < nmatch; i++)
	{
	  int start = pmatch[i].rm_so;
	  int end = pmatch[i].rm_eo;
	  if (start == -1)
	    {
	      break;
	    }
	  else
	    {
	      int size = end-start;
	      char* el = malloc (sizeof(char)*(size+1));
	      assert (el != NULL);
	      strncpy (el, &str_process[start], size);
	      el[size] = '\0';
	      res[i-1] = el;
	    }
	}
    }
}

/*
 * process one line of our strategy file to collect our strategy book!
 */
static int
process_line (facul_strategies_t* strategies, unsigned int* index_st,
	      const char *str, const int verbose)
{
  int index_method = 0; /* this is the index of the current factoring
			   methods */
  unsigned int SIGMA[2] = {2,2};
  int is_first_brent12[2] = {true, true};

  regex_t preg_index, preg_fm;
  //regular expression for the sides
  const char *str_preg_index =
        "r0=([[:digit:]]+)"
        ",[[:space:]]*"
        "r1=([[:digit:]]+)";
  //regular expression for the strategy
  const char *str_preg_fm =
        "(S[[:alnum:]]+):[[:space:]]*"  /* side, like "S0: " */
        "([[:alnum:]-]+)"               /* method, like "PP1-65" or "ECM-M12" */
        ",[[:space:]]*([[:digit:]]+)"   /* B1, an integer */
        ",[[:space:]]*([[:digit:]]+)";  /* B2, an integer */
  regcomp (&preg_index, str_preg_index, REG_ICASE|REG_EXTENDED);
  regcomp (&preg_fm, str_preg_fm, REG_ICASE|REG_EXTENDED);

  //process the ligne
  const char * str_process = &str[0];
  int side = -1;
  while (str_process[0] != '\0' )
    {
      //init
      size_t nmatch = 5;
      regmatch_t *pmatch= calloc (sizeof(*pmatch), nmatch);
      char **res = malloc (sizeof(char*) * nmatch);
      /*TEST REGULAR EXPRESSION  'preg_index*/
      regexec (&preg_index, str_process, nmatch, pmatch, 0);
      return_data_ex (res, pmatch, nmatch, str_process);
      if (res[0] != NULL)
	{
	  /* changes the current strategy! */
	  index_st[0] = atoi(res[0]);
	  index_st[1] = atoi(res[1]);
	  /* re-init the value of sigma */
	  SIGMA[0] = 2;
	  SIGMA[1] = 2;
	  //todo: change it to add it in the parameters of our function.
	  //maybe unused if you use only one curve B12 by strategy.
	  is_first_brent12[0] = true;
	  is_first_brent12[1] = true;
	}

      /*else TEST REGULAR EXPRESSION  'preg_alg'*/
      else
	{
	  regexec (&preg_fm, str_process, nmatch, pmatch, 0);
	  return_data_ex (res, pmatch, nmatch, str_process);
	  if (res[0] != NULL)
	    {
	      /*add the new factoring method to the current strategy! */
	      if (index_st[0] > strategies->mfb[0] ||
		  index_st[1] > strategies->mfb[1])
		return 0;
	      
	      facul_method_side_t* methods =
		strategies->methods[index_st[0]][index_st[1]];

	      if (strcmp (res[0], "S1") == 0)
		side = 1;
	      else if (strcmp(res[0], "S0") == 0)
		side = 0;
	      else 
		side = atoi(res[0]);

	      unsigned int B1 = (unsigned int) atoi (res[2]);
	      unsigned int B2 = (unsigned int) atoi (res[3]);
	      if (B1 == 0 && B2 == 0)
		{ //zero method!
		  goto next_regex;
		}
	      unsigned long sigma = 0;
	      int curve = 0;
	      int method = 0;
	      //method
	      if (strcmp (res[1], "PM1") == 0)
		method = PM1_METHOD;
	      else if (strcmp (res[1], "PP1-27") == 0) {
		method = PP1_METHOD;
		curve = PP1_2_7;
              } else if (strcmp (res[1], "PP1-65") == 0) {
		method = PP1_METHOD;
		curve = PP1_6_5;
	      } else 
		{
		  method = EC_METHOD;
		  //curve
		  if (strcmp (res[1], "ECM-B12") == 0)
		    curve = BRENT12;
		  else if (strcmp (res[1], "ECM-M12") == 0)
		    curve = MONTY12;
		  else if (strcmp (res[1], "ECM-M16") == 0)
		    curve = MONTY16;
		  else
		    {
		      fprintf (stderr,
			       "error : the method '%s' is unknown!\n",
			       res[1]);
		      return -1;
		    }
		  //sigma!
		  if (curve == MONTY16)
		    sigma = 1;
		  else
		    {
		      if (curve == BRENT12 && is_first_brent12[side])
			{
			  sigma = 11;
			  is_first_brent12[side] = false;
			}
		      else
			sigma = SIGMA[side]++;
		    }
		}
	      //check if the method is already computed!
	      int index_prec_fm = 
		get_index_method(strategies->precomputed_methods, B1, B2,
				 method, curve, sigma);
	      if ( strategies->precomputed_methods[index_prec_fm].method == 0)
		{
		  /*
		   * The current method isn't already precomputed. So
		   * we will compute and store it.
		   */
		  void* plan = NULL;
		  if (method == PM1_METHOD)
		    {
		      plan = malloc (sizeof (pm1_plan_t));
		      pm1_make_plan (plan, B1, B2, verbose);
		    }
		  else if (method == PP1_METHOD)
		    {
		      plan = malloc (sizeof (pp1_plan_t));
		      pp1_make_plan (plan, B1, B2, curve, verbose);
		    }
		  else { //method == EC_METHOD
		    plan = malloc (sizeof (ecm_plan_t));
		    ecm_make_plan (plan,
				   B1, B2, curve, sigma, 1, verbose);
		  }
		  strategies->precomputed_methods[index_prec_fm].method =method;
		  strategies->precomputed_methods[index_prec_fm].plan = plan;
		  
		  ASSERT_ALWAYS (index_prec_fm+1 < NB_MAX_METHODS);
		  //to show the end of plan!
		  strategies->precomputed_methods[index_prec_fm+1].plan = NULL;
		  strategies->precomputed_methods[index_prec_fm+1].method = 0;
		}
	      /* 
	       * Add this method to the current strategy
	       * methods[index_st[0]][index_st[1]]. 
	       */
	      methods[index_method].method =
		&strategies->precomputed_methods[index_prec_fm];
	      methods[index_method].side = side;
	      index_method++;
	      ASSERT_ALWAYS (index_method < NB_MAX_METHODS);
	      
	      //to show the end of methods!
	      methods[index_method].method = NULL;
	    }
	  else//to end the while
	    {
	      pmatch[0].rm_eo = strlen(str_process);
	    }
	}
      next_regex:
      str_process = &str_process[pmatch[0].rm_eo];
      //free
      for (size_t i = 0; i < nmatch; i++)
	free(res[i]);
      free(res);
      free (pmatch);
    }
  
  regfree(&preg_index);
  regfree(&preg_fm);
  return 0;
}

/*
 * This function generates a chain of 'n' ecm used to the auxiliary
 * cofactorization.
 */
facul_method_t*
facul_make_aux_methods (int n, const int verbose)
{
  if (n <= 0)
    n = 30;//nb_curves

  facul_method_t *methods = malloc ((n+1) * sizeof (facul_method_t));

  /* run one ECM curve with Montgomery parametrization, B1=105, B2=3255 */
  methods[0].method = EC_METHOD;
  methods[0].plan = malloc (sizeof (ecm_plan_t));
  ecm_make_plan (methods[0].plan, 105, 3255, MONTY12, 2, 1, verbose);

  methods[1].method = EC_METHOD;
  methods[1].plan = malloc (sizeof (ecm_plan_t));
  ecm_make_plan (methods[1].plan, 315, 5355, BRENT12, 11, 1, verbose);

  /* heuristic strategy where B1 is increased by sqrt(B1) at each curve */
  double B1 = 105.0;
  for (int i = 2; i < n ; i++)
    {
      double B2;
      unsigned int k;

      B1 += sqrt (B1);
      B2 = 17.0 * B1;
      /* we round B2 to (2k+1)*105, thus k is the integer nearest to
	 B2/210-0.5 */
      k = B2 / 210.0;
      methods[i].method = EC_METHOD;
      methods[i].plan = malloc (sizeof (ecm_plan_t));
      ecm_make_plan (methods[i].plan, (unsigned int) B1, (2 * k + 1) * 105,
		     MONTY12, i + 1, 1, 0);
    }

  methods[n].method = 0;
  methods[n].plan = NULL;

  return methods;
}

/*
  Make a simple strategy for factoring. We start with
  P-1 and P+1 (with x0=2/7), then an ECM curve with low bounds, then
  a bunch of ECM curves with larger bounds. How many methods to do in
  total is controlled by the n parameter: P-1, P+1 and the first ECM
  curve (with small bounds) are always done, then n ECM curves (with
  larger bounds).
  This function is used when you don't give a strategy file.
*/

facul_method_t*
facul_make_default_strategy (int n, const int verbose)
{
  ASSERT_ALWAYS (n >= 0);  
  facul_method_t *methods = malloc ((n+4) * sizeof (facul_method_t));

  /* run one P-1 curve with B1=315 and B2=2205 */
  methods[0].method = PM1_METHOD;
  methods[0].plan = malloc (sizeof (pm1_plan_t));
  pm1_make_plan (methods[0].plan, 315, 2205, verbose);

  /* run one P+1 curve with B1=525 and B2=3255 */
  methods[1].method = PP1_METHOD;
  methods[1].plan = malloc (sizeof (pp1_plan_t));
  pp1_make_plan (methods[1].plan, 525, 3255, PP1_2_7, verbose);

  /* run one ECM curve with Montgomery parametrization, B1=105, B2=3255 */
  methods[2].method = EC_METHOD;
  methods[2].plan = malloc (sizeof (ecm_plan_t));
  ecm_make_plan (methods[2].plan, 105, 3255, MONTY12, 2, 1, verbose);

  if (n > 0)
    {
      methods[3].method = EC_METHOD;
      methods[3].plan = malloc (sizeof (ecm_plan_t));
      ecm_make_plan (methods[3].plan, 315, 5355, BRENT12, 11, 1, verbose);
    }

  /* heuristic strategy where B1 is increased by sqrt(B1) at each curve */
  double B1 = 105.0;
  for (int i = 4; i < n + 3; i++)
    {
      double B2;
      unsigned int k;

      B1 += sqrt (B1);
      B2 = 17.0 * B1;
      /* we round B2 to (2k+1)*105, thus k is the integer nearest to
	 B2/210-0.5 */
      k = B2 / 210.0;
      methods[i].method = EC_METHOD;
      methods[i].plan = malloc (sizeof (ecm_plan_t));
      ecm_make_plan (methods[i].plan, (unsigned int) B1, (2 * k + 1) * 105,
		     MONTY12, i - 1, 1, 0);
    }

#ifdef USE_MPQS
  /* replace last method by MPQS */
  if (n > 1)
    {
      ecm_clear_plan (methods[n+2].plan);
      /* the plan pointer will be freed in facul_clear_aux_methods() */
      methods[n+2].method = MPQS_METHOD;
    }
#endif

  methods[n+3].method = 0;
  methods[n+3].plan = NULL;
  return methods;
}


void 
facul_clear_aux_methods (facul_method_t *methods)
{
  if (methods == NULL)
    return;

  for (int i = 0; methods[i].method != 0; i++)
    {
      if (methods[i].method == PM1_METHOD)
	pm1_clear_plan (methods[i].plan);
      else if (methods[i].method == PP1_METHOD)
	pp1_clear_plan (methods[i].plan);
      else if (methods[i].method == EC_METHOD)
	ecm_clear_plan (methods[i].plan);
      methods[i].method = 0;
      free (methods[i].plan);
      methods[i].plan = NULL;
    }
  free (methods);
  methods = NULL;
}


/*
 * Create our strategy book from a file (if a file is given) and
 * otherwise from our default strategy.
 */
facul_strategies_t*
facul_make_strategies(const unsigned long rfbb, const unsigned int rlpb,
		      const unsigned int rmfb, const unsigned long afbb,
		      const unsigned int alpb, const unsigned int amfb,
		      int n0, int n1, FILE* file, const int verbose)
{
  //printf ("create strategies!\n");
  facul_strategies_t* strategies = malloc (sizeof(facul_strategies_t));
  ASSERT (strategies != NULL);
  strategies->mfb[0] = rmfb;
  strategies->mfb[1] = amfb;

  strategies->lpb[0] = rlpb;
  strategies->lpb[1] = alpb;
  /* Store fbb^2 in assume_prime_thresh */
  strategies->assume_prime_thresh[0] = (double) rfbb * (double) rfbb;
  strategies->assume_prime_thresh[1] = (double) afbb * (double) afbb;

  strategies->BBB[0] = (double) rfbb * strategies->assume_prime_thresh[0];
  strategies->BBB[1] = (double) afbb * strategies->assume_prime_thresh[1];

  //alloc methods!
  facul_method_side_t*** methods = malloc (sizeof (*methods) * (rmfb+1));
  unsigned int r, a;
  //init methods!
  for (r = 0; r <= rmfb; r++) {
    methods[r] = malloc (sizeof (*methods[r]) * (amfb+1));
    ASSERT (methods[r] != NULL);
    for (a = 0; a <= amfb; a++)
      {
	methods[r][a] = malloc (NB_MAX_METHODS * sizeof (facul_method_side_t));
	ASSERT (methods[r][a] != NULL);
	methods[r][a][0].method = NULL;
      }
  }
  strategies->methods = methods;
  
  /*Default strategy. */ 
  if (file == NULL)
    {//make_default_strategy
      verbose_output_print(0, 1, "# Using default strategy for the cofactorization\n");
      int ncurves[2];
      ncurves[0] = (n0 > -1) ? n0 : nb_curves (rlpb);
      ncurves[1] = (n1 > -1) ? n1 : nb_curves (alpb);
      int max_ncurves = ncurves[0] > ncurves[1]? ncurves[0]: ncurves[1];
      strategies->precomputed_methods =
	facul_make_default_strategy (max_ncurves,verbose);
      for (r = 0; r <= rmfb; r++)
	for (a = 0; a <= amfb; a++) {
	  int index = 0;
	  int first = (r < a);//begin by the largest cofactor.
	  for (int z  = 0; z < 2; z++)
	    {
	      int side = first ^ z;
	      for (int i = 0; i < ncurves[side] + 3; i++)
		{
		  methods[r][a][index].method =
		    &strategies->precomputed_methods[i];
		  methods[r][a][index].side = side;
		  index++;
		}
	    }
	  //add NULL to show the end of the strategy.
	  methods[r][a][index].method = NULL;
	}
    }
  else
    {/* to precompute our method from the file.*/
      verbose_output_print(0, 1, "# Read the cofactorization strategy file\n");
      facul_method_t* precomputed_methods =
	malloc (sizeof(facul_method_t) * NB_MAX_METHODS);
      precomputed_methods[0].method = 0;
      precomputed_methods[0].plan = NULL;
      
      strategies->precomputed_methods = precomputed_methods;
      unsigned int index_strategies[2] = {0,0};
      
      char line[10000];
      
      fseek (file, 0, SEEK_SET);
      while (fgets (line, sizeof(line), file) != NULL)
	{
	  //process each line of 'file'
	  int err = process_line (strategies, index_strategies,
				  line, verbose);
	  if (err == -1)
	    return NULL;
	}
    }
  /*
    For each strategy, one finds what is the last method used on each
    side.
  */
  for (r = 1; r <= rmfb; r++)
    for (a = 1; a <= amfb; a++){
      facul_method_side_t* methods =
	strategies->methods[r][a];
      if (methods == NULL)
	continue;
      int index_last_method[2] = {-1, -1};      //the default_values
      unsigned int i;
      for (i = 0; methods[i].method != 0; i++) {
	methods[i].is_the_last = 0;
	index_last_method[methods[i].side] = i;
      }
      //== -1 if it exists zero method for this side!
      if (index_last_method[0] != -1)
	  methods[index_last_method[0]].is_the_last = 1;
      if (index_last_method[1] != -1)
	  methods[index_last_method[1]].is_the_last = 1;
    }
  
  //Create the auxiliary methods!
  //add test to check if it's necessary to create our aux methods!
  //todo: choose a better choice of our number of curves!
  strategies->methods_aux = facul_make_aux_methods (30, verbose);

  return strategies;
}


void
facul_clear_strategies (facul_strategies_t *strategies)
{
  if (strategies == NULL)
    return ;

  //free precomputed_methods
  facul_method_t* fm = strategies->precomputed_methods;
  for (int j = 0; fm[j].method!=0; j++)
    {
      if (fm[j].plan == NULL)
	continue;
      if (fm[j].method == PM1_METHOD)
	pm1_clear_plan (fm[j].plan);
      else if (fm[j].method == PP1_METHOD)
	pp1_clear_plan (fm[j].plan);
      else if (fm[j].method == EC_METHOD)
	ecm_clear_plan (fm[j].plan);
      free (fm[j].plan);
      fm[j].method = 0;
      fm[j].plan = NULL;
    }
  free (fm);
  fm = NULL;

  //free methods
  unsigned int r;
  for (r = 0; r <= strategies->mfb[0]; r++) {
    unsigned int a;
    for (a = 0; a <= strategies->mfb[1]; a++)
      free (strategies->methods[r][a]);
    free (strategies->methods[r]);
  }
  free (strategies->methods);

  //free methods_aux
  facul_clear_aux_methods (strategies->methods_aux);

  //free strategies
  free (strategies);

}


int
facul_fprint_strategies (FILE* file, facul_strategies_t* strategies)
{
  if (file == NULL)
    return -1;
  unsigned int r;
  //print info lpb ...
  fprintf (file,
	  "(lpb = [%ld,%ld], as...=[%lf, %lf], BBB = [%lf, %lf])\n",
	  strategies->lpb[0], strategies->lpb[1],
	  strategies->assume_prime_thresh[0],
	  strategies->assume_prime_thresh[1],
	  strategies->BBB[0], strategies->BBB[1]);
  fprintf (file,
	  "mfb = [%d, %d]\n", strategies->mfb[0], strategies->mfb[1]);
  for (r = 0; r <= strategies->mfb[0]; r++) {
    unsigned int a;
    for (a = 0; a <= strategies->mfb[1]; a++) {
      fprintf (file, "[r = %d, a = %d]", r, a);
      facul_method_side_t* methods = strategies->methods[r][a];
      if (methods == NULL)
	continue;
      
      int i = 0;
      while (true)
	{
	  facul_method_t* fm = methods[i].method;
	  if (fm == NULL || fm->method == 0)
	    break;
	  if (fm->plan == NULL)//zero method!!!
	    fprintf (file,
		    "[side=%d, FM=%ld, B1=0, B2=0] ", methods[i].side,
		    fm->method);
	  else {
	    if (fm->method == PM1_METHOD)
	      {
		pm1_plan_t* plan = (pm1_plan_t*) fm->plan;
		fprintf (file,
			 "[side=%d, FM=%ld, B1=%d, B2=%d] ", methods[i].side,
			fm->method, plan->B1, plan->stage2.B2);
	      }
	    else if (fm->method == PP1_METHOD)
	      {
		pp1_plan_t* plan = (pp1_plan_t*) fm->plan;
		fprintf (file,
		       "[side=%d, FM=%ld, B1=%d, B2=%d] ", methods[i].side,
			fm->method, plan->B1, plan->stage2.B2);
	      }
	    else if (fm->method == EC_METHOD)
	      {
		ecm_plan_t* plan = (ecm_plan_t*) fm->plan;
		fprintf (file,
			"[side=%d, FM=%ld, B1=%d, B2=%d] ", methods[i].side,
			fm->method, plan->B1,plan->stage2.B2);
	      }
	    else
	      return -1;
	  }
	  i++;
	}
      fprintf (file, "\n");
    }
  }
  return 1;
}


/*
 * This is our auxiliary factorization.
 * It applies a bunch of ECM curves with larger bounds to find
 * a factor with high probability. It returns -1 if the factor
 * is not smooth, otherwise the number of
 * factors.
 */
static int
facul_aux (unsigned long *factors, const modset_t m,
	   const facul_strategies_t *strategies, int method_start, int side)
{
  int found = 0;
  facul_method_t* methods = strategies->methods_aux;
  if (methods == NULL)
    return found;

  int i = 0;
  for (i = method_start ;methods[i].method != 0; i++)
    {
      if (i < STATS_LEN)
	stats_called_aux[i]++;
      modset_t fm, cfm;

      int res_fac = 0;

      switch (m.arith) {
      case CHOOSE_UL:
	res_fac = facul_doit_onefm_ul(factors, m.m_ul,
				      methods[i], &fm, &cfm,
				      strategies->lpb[side],
				      strategies->assume_prime_thresh[side],
				      strategies->BBB[side]);
	break;
      case CHOOSE_15UL:
	res_fac = facul_doit_onefm_15ul(factors, m.m_15ul,
					methods[i], &fm, &cfm,
					strategies->lpb[side],
					strategies->assume_prime_thresh[side],
					strategies->BBB[side]);
	break;
      case CHOOSE_2UL2:
	res_fac = facul_doit_onefm_2ul2 (factors, m.m_2ul2,
					 methods[i], &fm, &cfm,
					 strategies->lpb[side],
					 strategies->assume_prime_thresh[side],
					 strategies->BBB[side]);
	break;
      case CHOOSE_MPZ:
	res_fac = facul_doit_onefm_mpz (factors, m.m_mpz,
					methods[i], &fm, &cfm,
					strategies->lpb[side],
					strategies->assume_prime_thresh[side],
					strategies->BBB[side]);
	break;
      default: abort();
      }
      //check our result!
      //res_fac contains the number of factors found!
      if (res_fac == -1)
	{
	  /*
	    The cofactor m is not smooth. So, one stops the
	    cofactorization.
	  */
	  found = FACUL_NOT_SMOOTH;
	  break;
	}
      if (res_fac == 0)
	{
	  /* Zero factor found. If it was the last method for this
	     side, then one stops the cofactorization. Otherwise, one
	     tries with an other method! */
	    continue;
	}

      found += res_fac;
      if (res_fac == 2)
	break;

      /*
	res_fac == 1!  Only one factor has been found. Hence, our
	factorization is not finished.
      */
      if (fm.arith != CHOOSE_NONE)
	{
	  int found2 = facul_aux (factors+res_fac, fm, strategies,
				  i+1, side);
	  if (found2 == FACUL_NOT_SMOOTH)
	    {
	      found = FACUL_NOT_SMOOTH;
	      modset_clear (&cfm);
	      modset_clear (&fm);
	      break;
	    }
	  else
	    found += found2;
	  modset_clear (&fm);
	}
      if (cfm.arith != CHOOSE_NONE)
	{
	  int found2 = facul_aux (factors+res_fac, cfm, strategies,
				  i+1, side);
	  if (found2 == FACUL_NOT_SMOOTH)
	    found = FACUL_NOT_SMOOTH;
	  else
	    found += found2;
	  modset_clear (&cfm);
	  break;
	}
      break;
    }

  return found;
}





/*
  This function tries to factor a pair of cofactors (m[0], m[1]) from
  strategies. It returns the number of factors found on each side, or
  -1 if the factor is not smooth.
  Remarks: - the values of factors found are stored in 'factors'.
           - the variable 'is_smooth' allows to know if a cofactor
             is already factored.
*/

static int*
facul_both_src (unsigned long **factors, const modset_t* m,
		const facul_strategies_t *strategies, int* cof,
		int* is_smooth)
{
  int* found = calloc(2, sizeof(int));

  facul_method_side_t* methods = strategies->methods[cof[0]][cof[1]];

  if (methods == NULL)
    return found;

  modset_t f[2][2];
  f[0][0].arith = CHOOSE_NONE;
  f[0][1].arith = CHOOSE_NONE;
  f[1][0].arith = CHOOSE_NONE;
  f[1][1].arith = CHOOSE_NONE;
  int stats_nb_side = 0, stats_index_transition = 0;
  for (int i = 0; methods[i].method != NULL; i++)
    {
      //{for the stats
      stats_current_index = i - stats_nb_side * stats_index_transition;
      if (methods[i].is_the_last)
	{
	  stats_nb_side = 1;
	  stats_index_transition = i+1;
	}
      //}
      int side = methods[i].side;
      if (is_smooth[side] != FACUL_MAYBE)
	continue;

      //{for the stats
      if (stats_current_index < STATS_LEN)
	stats_called[stats_current_index]++;
      //}
      int res_fac = 0;
      switch (m[side].arith) {
      case CHOOSE_UL:
	res_fac = facul_doit_onefm_ul(factors[side], m[side].m_ul,
				      *methods[i].method, &f[side][0], &f[side][1],
				      strategies->lpb[side],
				      strategies->assume_prime_thresh[side],
				      strategies->BBB[side]);
	break;
      case CHOOSE_15UL:
	res_fac = facul_doit_onefm_15ul(factors[side], m[side].m_15ul,
					*methods[i].method, &f[side][0], &f[side][1],
					strategies->lpb[side],
					strategies->assume_prime_thresh[side],
					strategies->BBB[side]);
	break;
      case CHOOSE_2UL2:
	res_fac = facul_doit_onefm_2ul2 (factors[side], m[side].m_2ul2,
					 *methods[i].method, &f[side][0], &f[side][1],
					 strategies->lpb[side],
					 strategies->assume_prime_thresh[side],
					 strategies->BBB[side]);
	break;
      case CHOOSE_MPZ:
	res_fac = facul_doit_onefm_mpz (factors[side], m[side].m_mpz,
					*methods[i].method, &f[side][0], &f[side][1],
					strategies->lpb[side],
					strategies->assume_prime_thresh[side],
					strategies->BBB[side]);
	break;
      default: abort();
      }
      //check our result!
      //res_fac contains the number of factors found!
      if (res_fac == -1)
	{
	  /*
	    The cofactor m[side] is not smooth. So, one stops the
	    cofactorization.
	  */
	  found[side] = -1;
	  break;
	}
      if (res_fac == 0)
	{
	  /* No factor found. If it was the last method for this
	     side, then one stops the cofactorization. Otherwise, one
	     tries with an other method.
	  */
	  if (methods[i].is_the_last)
            break;
	  else
	    continue;
	}
      found[side] = res_fac;
      
      if (res_fac == 2)
	{
	  /*
	    Indeed, using only one factoring method, you found two
	    primes factors of m (f, m/f) then m is factored and your
	    work is finished for this cofactor!!
	  */
	  is_smooth[side] = FACUL_SMOOTH;
	  continue;
	}
      /*
	res_fac == 1. Only one factor has been found. Hence, a
	auxiliary factorization will be necessary!
       */
      is_smooth[side] = FACUL_AUX;
    }
  //begin the auxiliary factorization!
  if (is_smooth[0] >= 1 && is_smooth[1] >= 1)
    for (int side = 0; side < 2; side++)
      if (is_smooth[side] == FACUL_AUX)
	for (int ind_cof = 0; ind_cof < 2; ind_cof++)
	  {
	    //factor f[side][0] or/and f[side][1]
	    if (f[side][ind_cof].arith != CHOOSE_NONE)
	      {
		int found2 = facul_aux (factors[side]+found[side],
					f[side][ind_cof], strategies, 0, side);
		if (found2 < 1)//FACUL_NOT_SMOOTH or FACUL_MAYBE
		  {
		    is_smooth[side] = FACUL_NOT_SMOOTH;
		    found[side] = found2;//FACUL_NOT_SMOOTH or FACUL_MAYBE
		    goto clean_up;
		  }
		else
		  {
		    is_smooth[side] = FACUL_SMOOTH;
		    found[side] += found2;
		  }
	      }
	  }

 clean_up:
  modset_clear (&f[0][0]);
  modset_clear (&f[0][1]);
  modset_clear (&f[1][0]);
  modset_clear (&f[1][1]);
  return found;
}


/*
  This function is like facul, but we will work with the both norms
  together.
  It returns the number of factors for each side!
*/
int*
facul_both (unsigned long **factors, mpz_t* N,
	    const facul_strategies_t *strategies, int* is_smooth)
{
  int cof[2];
  size_t bits;
  int* found = NULL;

  modset_t n[2];
  n[0].arith = CHOOSE_NONE;
  n[1].arith = CHOOSE_NONE;

#ifdef PARI
  gmp_fprintf (stderr, "(%Zd %Zd)", N[0], N[1]);
#endif

  if (mpz_sgn (N[0]) <= 0 || mpz_sgn (N[1]) <= 0)
      return found;

  if (mpz_cmp_ui (N[0], 1UL) == 0)
    is_smooth[0] = FACUL_SMOOTH;
  if (mpz_cmp_ui (N[1], 1UL) == 0)
    is_smooth[1] = FACUL_SMOOTH;

  for (int side = 0; side < 2; side++)
    {
      /* If the composite does not fit into our modular arithmetic, return
	 no factor */
      bits = mpz_sizeinbase (N[side], 2);
      cof[side] = bits;
      if (bits > MODMPZ_MAXBITS)
	return 0;

      /* Use the fastest modular arithmetic that's large enough for
	 this input */
      if (bits <= MODREDCUL_MAXBITS)
	{
	  ASSERT(mpz_fits_ulong_p(N[side]));
	  modredcul_initmod_ul (n[side].m_ul, mpz_get_ui(N[side]));
	  n[side].arith = CHOOSE_UL;
	}
      else if (bits <= MODREDC15UL_MAXBITS)
	{
	  unsigned long t[2];
	  modintredc15ul_t m;
	  size_t written;
	  mpz_export (t, &written, -1, sizeof(unsigned long), 0, 0, N[side]);
	  ASSERT_ALWAYS(written <= 2);
	  modredc15ul_intset_uls (m, t, written);
	  modredc15ul_initmod_int (n[side].m_15ul, m);
	  n[side].arith = CHOOSE_15UL;
	}
      else if (bits <= MODREDC2UL2_MAXBITS)
	{
	  unsigned long t[2];
	  modintredc2ul2_t m;
	  size_t written;
	  mpz_export (t, &written, -1, sizeof(unsigned long), 0, 0, N[side]);
	  ASSERT_ALWAYS(written <= 2);
	  modredc2ul2_intset_uls (m, t, written);
	  modredc2ul2_initmod_int (n[side].m_2ul2, m);
	  n[side].arith = CHOOSE_2UL2;
	}
      else
	{
	  modmpz_initmod_int (n[side].m_mpz, N[side]);
	  n[side].arith = CHOOSE_MPZ;
	}
      ASSERT (n[side].arith != CHOOSE_NONE);
    }

  found = facul_both_src (factors, n, strategies, cof, is_smooth);
  for (int side = 0; side < 2; side++)
    {
      if (found[side] > 1)
	{
	  /* Sort the factors we found */
	  qsort (factors[side], found[side], sizeof (unsigned long),
		 (int (*)(const void *, const void *)) &cmp_ul);
	}

#ifdef PARI
      if (found[side] > 1)
	{
	  fprintf (stderr, " == ");
	  for (int i = 0; i < found[side]; i++)
	    fprintf (stderr, "%lu%s", factors[side][i],
		     (i+1 < found[side]) ? " * " : " /* PARI */\n");
	}
      else
	fprintf (stderr, "; /* PARI */\n");
#endif
    }

  //Free
  modset_clear (&n[0]);
  modset_clear (&n[1]);

  return found;
}


/********************************************/
/*            modset_t                      */
/********************************************/


void
modset_clear (modset_t *modset)
{
  switch (modset->arith) {
  case CHOOSE_NONE: /* already clear */
    break;
  case CHOOSE_UL:
    modredcul_clearmod (modset->m_ul);
    break;
  case CHOOSE_15UL:
    modredc15ul_clearmod (modset->m_15ul);
    break;
  case CHOOSE_2UL2:
    modredc2ul2_clearmod (modset->m_2ul2);
    break;
  case CHOOSE_MPZ:
    modmpz_clearmod (modset->m_mpz);
    break;
  default:
    abort();
  }
  modset->arith = CHOOSE_NONE;
}
