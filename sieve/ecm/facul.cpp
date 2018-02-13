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
#include <vector>
#include <algorithm>

#include "utils.h"      /* verbose_ stuff */
#include "portability.h"
#include "pm1.h"
#include "pp1.h"
#include "facul_ecm.h"
#include "facul.hpp"
#include "facul_doit.hpp"

#ifdef ENABLE_UNSAFE_FACUL_STATS
/*
 * FIXME: the stats are not thread-safe!
 *
 * These global variables are only for statistics. In case of
 * multithreaded sieving, the stats might be wrong...
 */

unsigned long stats_called[STATS_LEN] = {
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
};
// for the auxiliary factorization.
unsigned long stats_called_aux[STATS_LEN] = {
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
};
int stats_current_index = 0; // only useful for stats_found_n
unsigned long stats_found_n[STATS_LEN] = {
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
};
#endif  /* ENABLE_UNSAFE_FACUL_STATS */

static int nb_curves90 (const unsigned int lpb);
#if 0
static int nb_curves95 (const unsigned int lpb);
static int nb_curves99 (const unsigned int lpb);
#endif

/* don't use nb_curves90, nb_curves95 and nb_curves99, only use nb_curves */
int nb_curves (const unsigned int lpb, const unsigned int mfb)
{
  /* for up to 2 large primes, we use a very large number of curves (say 100),
     since the more curves we need, the more likely the cofactor is smooth. */
  if (mfb <= 2 * lpb)
    return 100;
  return nb_curves90 (lpb);
}

static int
nb_curves90 (const unsigned int lpb)
{
  /* The following table, computed with do_table(10,33,ntries=10000) from the
     facul.sage file, ensures a probability of at least about 90% to find a
     factor below 2^lpb with n = T[lpb]. */
  int T[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 0-9 */
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 10-19 */
	     /* lpb=20 */ 1, /* 0:0.878100, 1:0.969600 */
	     /* lpb=21 */ 1, /* 1:0.940400 */
	     /* lpb=22 */ 1, /* 1:0.907400 */
	     /* lpb=23 */ 2, /* 1:0.859100, 2:0.903900 */
	     /* lpb=24 */ 4, /* 3:0.897000, 4:0.929200 */
	     /* lpb=25 */ 5, /* 4:0.884100, 5:0.916400 */
	     /* lpb=26 */ 6, /* 5:0.868600, 6:0.904200 */
	     /* lpb=27 */ 8, /* 7:0.873500, 8:0.901700 */
             /* lpb=28 */ 11, /* 10:0.896600, 11:0.918000 */
	     /* lpb=29 */ 13, /* 12:0.881700, 13:0.905600 */
	     /* lpb=30 */ 16, /* 15:0.889400, 16:0.909400 */
	     /* lpb=31 */ 18, /* 17:0.883500, 18:0.901300 */
	     /* lpb=32 */ 21, /* 20:0.884500, 21:0.903700 */
	     /* lpb=33 */ 25, /* 24:0.892800, 25:0.910000 */
        /* The extra ones below are computed with do_table(10,64,ntries=200) */
	     /* lpb=34 */ 24, /* 23:0.895000, 24:0.905000 */
	     /* lpb=35 */ 29, /* 28:0.895000, 29:0.905000 */
	     /* lpb=36 */ 30, /* 29:0.890000, 30:0.905000 */
	     /* lpb=37 */ 35, /* 34:0.890000, 35:0.900000 */
	     /* lpb=38 */ 37, /* 36:0.890000, 37:0.900000 */
	     /* lpb=39 */ 42, /* 41:0.890000, 42:0.900000 */
	     /* lpb=40 */ 45, /* 44:0.890000, 45:0.900000 */
	     /* lpb=41 */ 50, /* 49:0.885000, 50:0.905000 */
	     /* lpb=42 */ 56, /* 55:0.875000, 56:0.900000 */
	     /* lpb=43 */ 60, /* 59:0.895000, 60:0.900000 */
	     /* lpb=44 */ 68, /* 67:0.890000, 68:0.910000 */
	     /* lpb=45 */ 73, /* 72:0.895000, 73:0.905000 */
	     /* lpb=46 */ 82, /* 81:0.895000, 82:0.900000 */
	     /* lpb=47 */ 90, /* 89:0.895000, 90:0.900000 */
	     /* lpb=48 */ 93, /* 92:0.895000, 93:0.900000 */
	     /* lpb=49 */ 93, /* 93:0.900000 */
	     /* lpb=50 */ 111, /* 110:0.885000, 111:0.900000 */
	     /* lpb=51 */ 117, /* 116:0.890000, 117:0.905000 */
	     /* those below were computed with do_table(10,64,ntries=100) */
	     /* lpb=52 */ 130, /* 129:0.890000, 130:0.900000 */
	     /* lpb=53 */ 137, /* 136:0.890000, 137:0.900000 */
	     /* lpb=54 */ 142, /* 141:0.890000, 142:0.900000 */
	     /* lpb=55 */ 152, /* 151:0.890000, 152:0.920000 */
	     /* lpb=56 */ 167, /* 166:0.890000, 167:0.900000 */
	     /* lpb=57 */ 167, /* 167:0.910000 */
	     /* lpb=58 */ 186, /* 185:0.890000, 186:0.900000 */
	     /* lpb=59 */ 202, /* 201:0.880000, 202:0.900000 */
	     /* lpb=60 */ 227, /* 226:0.890000, 227:0.910000 */
	     /* lpb=61 */ 244, /* 243:0.890000, 244:0.910000 */
	     /* lpb=62 */ 244, /* 244:0.920000 */
	     /* lpb=63 */ 248, /* 247:0.880000, 248:0.900000 */
	     /* lpb=64 */ 294, /* 293:0.890000, 294:0.900000 */
  };
  const unsigned int nT = sizeof(T)/sizeof(int) - 1;
  return (lpb <= nT) ? T[lpb] : T[nT];
}

#if 0
static int
nb_curves95 (const unsigned int lpb)
{
    /* same, but with target probability 95% */
    /* do_table(10,64,ntries=500,target_prob=0.95)
     */
  int T[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 0-9 */
	     /* lpb=10 */ 0, /* 0:1.000000 */
	     /* lpb=11 */ 0, /* 0:1.000000 */
	     /* lpb=12 */ 0, /* 0:1.000000 */
	     /* lpb=13 */ 0, /* 0:1.000000 */
	     /* lpb=14 */ 0, /* 0:1.000000 */
	     /* lpb=15 */ 0, /* 0:0.998000 */
	     /* lpb=16 */ 0, /* 0:0.986000 */
	     /* lpb=17 */ 0, /* 0:0.972000 */
	     /* lpb=18 */ 0, /* 0:0.964000 */
	     /* lpb=19 */ 1, /* 0:0.926000, 1:0.986000 */
	     /* lpb=20 */ 1, /* 1:0.986000 */
	     /* lpb=21 */ 2, /* 1:0.946000, 2:0.976000 */
	     /* lpb=22 */ 3, /* 2:0.940000, 3:0.966000 */
	     /* lpb=23 */ 3, /* 3:0.952000 */
	     /* lpb=24 */ 5, /* 4:0.926000, 5:0.962000 */
	     /* lpb=25 */ 7, /* 6:0.934000, 7:0.956000 */
	     /* lpb=26 */ 8, /* 7:0.934000, 8:0.956000 */
	     /* lpb=27 */ 10, /* 9:0.936000, 10:0.956000 */
	     /* lpb=28 */ 13, /* 12:0.938000, 13:0.954000 */
	     /* lpb=29 */ 16, /* 15:0.946000, 16:0.956000 */
	     /* lpb=30 */ 18, /* 17:0.936000, 18:0.950000 */
	     /* lpb=31 */ 22, /* 21:0.934000, 22:0.958000 */
	     /* lpb=32 */ 26, /* 25:0.940000, 26:0.950000 */
	     /* lpb=33 */ 29, /* 28:0.942000, 29:0.952000 */
	     /* lpb=34 */ 33, /* 32:0.948000, 33:0.956000 */
	     /* lpb=35 */ 37, /* 36:0.938000, 37:0.952000 */
	     /* lpb=36 */ 42, /* 41:0.946000, 42:0.952000 */
	     /* lpb=37 */ 45, /* 44:0.946000, 45:0.958000 */
  };
  const unsigned int nT = sizeof(T)/sizeof(int) - 1;
  return (lpb <= nT) ? T[lpb] : T[nT];
}

static int
nb_curves99 (const unsigned int lpb)
{
    /* same, but with target probability 99% */
    /* do_table(10,64,ntries=100,target_prob=0.99)
     */
  int T[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 0-9 */
	     /* lpb=10 */ 0, /* 0:1.000000 */
	     /* lpb=11 */ 0, /* 0:1.000000 */
	     /* lpb=12 */ 0, /* 0:1.000000 */
	     /* lpb=13 */ 0, /* 0:1.000000 */
	     /* lpb=14 */ 0, /* 0:1.000000 */
	     /* lpb=15 */ 0, /* 0:1.000000 */
	     /* lpb=16 */ 0, /* 0:1.000000 */
	     /* lpb=17 */ 1, /* 0:0.980000, 1:1.000000 */
	     /* lpb=18 */ 1, /* 1:0.990000 */
	     /* lpb=19 */ 2, /* 1:0.980000, 2:1.000000 */
	     /* lpb=20 */ 2, /* 2:0.990000 */
	     /* lpb=21 */ 2, /* 2:0.990000 */
	     /* lpb=22 */ 4, /* 3:0.980000, 4:1.000000 */
	     /* lpb=23 */ 4, /* 4:0.990000 */
	     /* lpb=24 */ 5, /* 4:0.950000, 5:0.990000 */
	     /* lpb=25 */ 8, /* 7:0.980000, 8:0.990000 */
	     /* lpb=26 */ 12, /* 11:0.980000, 12:1.000000 */
	     /* lpb=27 */ 14, /* 13:0.980000, 14:1.000000 */
	     /* lpb=28 */ 14, /* 14:1.000000 */
	     /* lpb=29 */ 14, /* 14:0.990000 */
	     /* lpb=30 */ 17, /* 16:0.970000, 17:0.990000 */
	     /* lpb=31 */ 25, /* 24:0.980000, 25:1.000000 */
	     /* lpb=32 */ 25, /* 25:0.990000 */
	     /* lpb=33 */ 25, /* 25:0.990000 */
	     /* lpb=34 */ 32, /* 31:0.980000, 32:0.990000 */
	     /* lpb=35 */ 36, /* 35:0.980000, 36:0.990000 */
	     /* lpb=36 */ 39, /* 38:0.980000, 39:0.990000 */
	     /* lpb=37 */ 43, /* 42:0.980000, 43:0.990000 */
	     /* lpb=38 */ 49, /* 48:0.980000, 49:0.990000 */
	     /* lpb=39 */ 52, /* 51:0.970000, 52:0.990000 */
	     /* lpb=40 */ 56, /* 55:0.980000, 56:0.990000 */
	     /* lpb=41 */ 60, /* 59:0.980000, 60:0.990000 */
  };
  const unsigned int nT = sizeof(T)/sizeof(int) - 1;
  return (lpb <= nT) ? T[lpb] : T[nT];
}
#endif

/* Wrapper around facul_make_default_strategy */
facul_strategy_t *
facul_make_strategy (const unsigned long fbb, const unsigned int lpb,
		     const unsigned int mfb, int n, const int verbose)
{
  facul_strategy_t *strategy;

  if (n == -1)
    n = nb_curves (lpb, mfb);
  strategy = (facul_strategy_t*) malloc (sizeof (facul_strategy_t));
  strategy->lpb = lpb;
  /* Store fbb^2 in assume_prime_thresh */
  strategy->assume_prime_thresh = (double) fbb * (double) fbb;
  strategy->BBB = (double) fbb * strategy->assume_prime_thresh;

  strategy->methods = facul_make_default_strategy (n, verbose);
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
	pm1_clear_plan ((pm1_plan_t*) methods[i].plan);
      else if (methods[i].method == PP1_27_METHOD)
	pp1_clear_plan ((pp1_plan_t*) methods[i].plan);
      else if (methods[i].method == PP1_65_METHOD)
	pp1_clear_plan ((pp1_plan_t*) methods[i].plan);
      else if (methods[i].method == EC_METHOD)
	ecm_clear_plan ((ecm_plan_t*) methods[i].plan);
      methods[i].method = 0;
      free (methods[i].plan);
      methods[i].plan = NULL;
    }
  free (methods);
  methods = NULL;
  free (strategy);
}

struct cxx_mpz_cmp {
    inline bool operator()(cxx_mpz const& a, cxx_mpz const & b) {
        return mpz_cmp(a, b) < 0;
    }
};

void facul_print_stats (FILE *stream MAYBE_UNUSED)
{
#ifdef ENABLE_UNSAFE_FACUL_STATS
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
#endif  /* ENABLE_UNSAFE_FACUL_STATS */
}


int
facul (std::vector<cxx_mpz> & factors, cxx_mpz const & N, const facul_strategy_t *strategy)
{
  int found = 0;
  size_t bits;
  
    /* XXX ATTENTION: This function may be called recursively. In
     * particular it may happen that the factors[] vector is not empty. */

  size_t factors_previous_size = factors.size();

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
      std::sort(factors.begin() + factors_previous_size, factors.end(), cxx_mpz_cmp());
    }

  return found;
}


/*****************************************************************************/
/*                       STRATEGY BOOK                                       */
/*****************************************************************************/

/*
 * If the plan is already precomputed and stored at the index i, return
 * i. Otherwise return the last index of tab to compute and add this new
 * plan at this index.
 */
static int
get_index_method (facul_method_t* tab, unsigned int B1, unsigned int B2,
		  int method, int parameterization, unsigned long sigma)
{
  int i = 0;
  while (tab[i].method != 0){
    if (tab[i].plan == NULL){
      if(B1 == 0 && B2 == 0)
	// zero method
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
      else if (method == PP1_27_METHOD ||
	       method == PP1_65_METHOD)
	{
	  pp1_plan_t* plan = (pp1_plan_t*)tab[i].plan;
	  if (plan->B1 == B1 && plan->stage2.B2 == B2)
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
  // (re)init res
  for (size_t i = 0; i < nmatch;i++)
    res[i] = NULL;

  if ( pmatch[0].rm_so != pmatch[0].rm_eo)
    {
      for (size_t i = 1; i < nmatch; i++)
	{
	  int start = pmatch[i].rm_so;
	  int end = pmatch[i].rm_eo;
	  if (start == -1)
	      break;
	  else
	    {
	      int size = end-start;
	      char* el = (char*) malloc (size+1);
	      assert (el != NULL);
	      strncpy (el, &str_process[start], size);
	      el[size] = '\0';
	      res[i-1] = el;
	    }
	}
    }
}

/*
 * process one line of our strategy file to collect our strategy book.
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
  // regular expression for the sides
  const char *str_preg_index =
        "r0=([[:digit:]]+)"
        ",[[:space:]]*"
        "r1=([[:digit:]]+)";
  // regular expression for the strategy
  const char *str_preg_fm =
        "(S[[:alnum:]]+):[[:space:]]*"  /* side, like "S0: " */
        "([[:alnum:]-]+)"               /* method, like "PP1-65" or "ECM-M12" */
        ",[[:space:]]*([[:digit:]]+)"   /* B1, an integer */
        ",[[:space:]]*([[:digit:]]+)";  /* B2, an integer */
  regcomp (&preg_index, str_preg_index, REG_ICASE|REG_EXTENDED);
  regcomp (&preg_fm, str_preg_fm, REG_ICASE|REG_EXTENDED);

  // process the line
  const char * str_process = &str[0];
  int side = -1;
  while (str_process[0] != '\0' )
    {
      // init
      size_t nmatch = 5;
      regmatch_t *pmatch= (regmatch_t*) calloc (sizeof(*pmatch), nmatch);
      char **res = (char**) malloc (sizeof(char*) * nmatch);
      /*TEST REGULAR EXPRESSION  'preg_index*/
      regexec (&preg_index, str_process, nmatch, pmatch, 0);
      return_data_ex (res, pmatch, nmatch, str_process);
      if (res[0] != NULL)
	{
	  /* changes the current strategy. */
	  index_st[0] = atoi(res[0]);
	  index_st[1] = atoi(res[1]);
	  /* re-init the value of sigma */
	  SIGMA[0] = 2;
	  SIGMA[1] = 2;
	  // todo: change it to add it in the parameters of our function.
	  // maybe unused if you use only one curve B12 by strategy.
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
	      /*add the new factoring method to the current strategy. */
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
		{ // zero method
		  goto next_regex;
		}
	      unsigned long sigma = 0;
	      int curve = 0;
	      int method = 0;
	      // method
	      if (strcmp (res[1], "PM1") == 0)
		method = PM1_METHOD;
	      else if (strcmp (res[1], "PP1-27") == 0)
		method = PP1_27_METHOD;
	      else if (strcmp (res[1], "PP1-65") == 0)
		method = PP1_65_METHOD;
	      else 
		{
		  method = EC_METHOD;
		  // curve
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
	      // check if the method is already computed
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
		      plan = (pm1_plan_t*) malloc (sizeof (pm1_plan_t));
		      pm1_make_plan ((pm1_plan_t*) plan, B1, B2, verbose);
		    }
		  else if (method == PP1_27_METHOD ||
			   method == PP1_65_METHOD)
		    {
		      plan = (pp1_plan_t*) malloc (sizeof (pp1_plan_t));
		      pp1_make_plan ((pp1_plan_t*) plan, B1, B2, verbose);
		    }
		  else { // method == EC_METHOD
		    plan = (ecm_plan_t*) malloc (sizeof (ecm_plan_t));
		    ecm_make_plan ((ecm_plan_t*) plan,
				   B1, B2, curve, sigma, 1, verbose);
		  }
		  strategies->precomputed_methods[index_prec_fm].method =method;
		  strategies->precomputed_methods[index_prec_fm].plan = plan;
		  
		  ASSERT_ALWAYS (index_prec_fm+1 < NB_MAX_METHODS);
		  // to show the end of plan
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
	      
	      // to show the end of methods
	      methods[index_method].method = NULL;
	    }
	  else// to end the while
	    {
	      pmatch[0].rm_eo = strlen(str_process);
	    }
	}
      next_regex:
      str_process = &str_process[pmatch[0].rm_eo];
      // free
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
  facul_method_t *methods = (facul_method_t*) malloc ((n+4) * sizeof (facul_method_t));

  int i = 0;

#if 0
  /* This is relevant only for very weird experiments where we have small
   * factors that stick together.  */
  /* run one P-1 curve with B1=30 and B2=100 */
  methods[i].method = PM1_METHOD;
  methods[i].plan = (pm1_plan_t*) malloc (sizeof (pm1_plan_t));
  pm1_make_plan ((pm1_plan_t*) methods[i].plan, 30, 100, verbose);
  i++;
#endif

  /* run one P-1 curve with B1=315 and B2=2205 */
  methods[i].method = PM1_METHOD;
  methods[i].plan = (pm1_plan_t*) malloc (sizeof (pm1_plan_t));
  pm1_make_plan ((pm1_plan_t*) methods[i].plan, 315, 2205, verbose);
  i++;

  /* run one P+1 curve with B1=525 and B2=3255 */
  methods[i].method = PP1_27_METHOD;
  methods[i].plan = (pp1_plan_t*) malloc (sizeof (pp1_plan_t));
  pp1_make_plan ((pp1_plan_t*) methods[i].plan, 525, 3255, verbose);
  i++;

  /* run one ECM curve with Montgomery parametrization, B1=105, B2=3255 */
  methods[i].method = EC_METHOD;
  methods[i].plan = (ecm_plan_t*) malloc (sizeof (ecm_plan_t));
  ecm_make_plan ((ecm_plan_t*) methods[i].plan, 105, 3255, MONTY12, 2, 1, verbose);
  i++;

  if (n > 0)
    {
      methods[i].method = EC_METHOD;
      methods[i].plan = (ecm_plan_t*) malloc (sizeof (ecm_plan_t));
      ecm_make_plan ((ecm_plan_t*) methods[i].plan, 315, 5355, BRENT12, 11, 1, verbose);
      i++;
    }

  /* heuristic strategy where B1 is increased by c*sqrt(B1) at each curve */
  double B1 = 105.0;
  for (; i < n + 3; i++)
    {
      double B2;
      unsigned int k;

      B1 += sqrt (B1);
      /* The factor 50 was determined experimentally with testbench, to find
	 factors of 40 bits:
	 testbench -p -cof 1208925819614629174706189 -strat 549755813888 549755913888
	 This finds 1908 factors (out of 3671 input numbers) with n=24 curves
	 and 3.66s.
	 With B2=17*B1, and 29 curves, we find 1898 factors in 4.07s.
	 With B2=100*B1, and 21 curves we find 1856 factors in 3.76s.
	 Thus 50 seems close to optimal.
	 Warning: changing the value of B2 might break the sieving tests.
      */
      B2 = 50.0 * B1;
      /* we round B2 to (2k+1)*105, thus k is the integer nearest to
	 B2/210-0.5 */
      k = B2 / 210.0;
      methods[i].method = EC_METHOD;
      methods[i].plan = (ecm_plan_t*) malloc (sizeof (ecm_plan_t));
      ecm_make_plan ((ecm_plan_t*) methods[i].plan, (unsigned int) B1, (2 * k + 1) * 105,
		     MONTY12, i - 1, 1, 0);
    }

#ifdef USE_MPQS
  /* replace last method by MPQS */
  if (n > 1)
    {
      ecm_clear_plan ((ecm_plan_t*) methods[n+2].plan);
      /* the plan pointer will be freed in facul_clear_strategy() */
      methods[n+2].method = MPQS_METHOD;
    }
#endif

  methods[n+3].method = 0;
  methods[n+3].plan = NULL;
  return methods;
}


void 
facul_clear_methods (facul_method_t *methods)
{
  if (methods == NULL)
    return;

  for (int i = 0; methods[i].method != 0; i++)
    {
      if (methods[i].method == PM1_METHOD)
	pm1_clear_plan ((pm1_plan_t*) methods[i].plan);
      else if (methods[i].method == PP1_27_METHOD)
	pp1_clear_plan ((pp1_plan_t*) methods[i].plan);
      else if (methods[i].method == PP1_65_METHOD)
	pp1_clear_plan ((pp1_plan_t*) methods[i].plan);
      else if (methods[i].method == EC_METHOD)
	ecm_clear_plan ((ecm_plan_t*) methods[i].plan);
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
                      bool perfectly_sieved,
		      int n0, int n1, FILE* file, const int verbose)
{
  unsigned int max_curves_used_before_aux = 0;
  // printf ("create strategies\n");
  facul_strategies_t* strategies = (facul_strategies_t*) malloc (sizeof(facul_strategies_t));
  ASSERT_ALWAYS (strategies != NULL);
  strategies->mfb[0] = rmfb;
  strategies->mfb[1] = amfb;

  strategies->lpb[0] = rlpb;
  strategies->lpb[1] = alpb;
  /* Store fbb^2 in assume_prime_thresh */
  strategies->assume_prime_thresh[0] = (double) rfbb * (double) rfbb;
  strategies->assume_prime_thresh[1] = (double) afbb * (double) afbb;

  strategies->BBB[0] = (double) rfbb * strategies->assume_prime_thresh[0];
  strategies->BBB[1] = (double) afbb * strategies->assume_prime_thresh[1];
  if (!perfectly_sieved) {
      strategies->assume_prime_thresh[0] = 0.0;
      strategies->assume_prime_thresh[1] = 0.0;
  }

  // alloc methods
  facul_method_side_t*** methods = (facul_method_side_t***) malloc (sizeof (*methods) * (rmfb+1));
  ASSERT_ALWAYS (methods != NULL);

  if (file == NULL) {
      /* we have just one strategy, really. So we just allocate once. */
      strategies->uniform_strategy[0] = (facul_method_side_t*)  malloc (NB_MAX_METHODS * sizeof (facul_method_side_t));
      ASSERT_ALWAYS (strategies->uniform_strategy[0] != NULL);
      strategies->uniform_strategy[1] = (facul_method_side_t*)  malloc (NB_MAX_METHODS * sizeof (facul_method_side_t));
      ASSERT_ALWAYS (strategies->uniform_strategy[1] != NULL);
  } else {
      strategies->uniform_strategy[0] = NULL;
      strategies->uniform_strategy[1] = NULL;
  }

  // init methods
  for (unsigned int r = 0; r <= rmfb; r++) {
    methods[r] = (facul_method_side_t**) malloc (sizeof (*methods[r]) * (amfb+1));
    ASSERT_ALWAYS (methods[r] != NULL);
    if (file != NULL) {
        for (unsigned int a = 0; a <= amfb; a++)
        {
            methods[r][a] = (facul_method_side_t*)  malloc (NB_MAX_METHODS * sizeof (facul_method_side_t));
            ASSERT_ALWAYS (methods[r][a] != NULL);
            methods[r][a][0].method = NULL;
        }
    }
  }
  strategies->methods = methods;
  
  /*Default strategy. */ 
  if (file == NULL)
    {// make_default_strategy
      int ncurves[2];
      ncurves[0] = (n0 > -1) ? n0 : nb_curves (rlpb, rmfb);
      ncurves[1] = (n1 > -1) ? n1 : nb_curves (alpb, amfb);
      int max_ncurves = ncurves[0] > ncurves[1]? ncurves[0]: ncurves[1];
      max_curves_used_before_aux = max_ncurves + 4; // account for fixed methods.
      // There is an hardcoded bound on the number of methods.
      // If ncurves0 or ncurves1 passed by the user is too large,
      // we can not handle that.
      // TODO: is this really a limitation?
      ASSERT_ALWAYS (2*(max_ncurves + 4) <= NB_MAX_METHODS);

      verbose_output_print(0, 1, "# Using default strategy for the cofactorization: ncurves0=%d ncurves1=%d\n", ncurves[0], ncurves[1]);
      strategies->precomputed_methods =
	facul_make_default_strategy (max_ncurves, verbose);
      /* make two well-behaved facul_method_side_t* lists, based on which
       * side is first */
      for (int first = 0 ; first < 2 ; first++) {
          facul_method_side_t * next = strategies->uniform_strategy[first];
          for (int z = 0; z < 2; z++) {
              int side = first ^ z;
              for (int i = 0; i < ncurves[side] + 3; i++) {
                  next->method = strategies->precomputed_methods + i;
                  next->side = side;
                  next++;
              }
          }
          next->method = NULL;
      }
      /* now have all entries in our huge 2d array point to either of
       * these */
      for (unsigned int r = 0; r <= rmfb; r++) {
	for (unsigned int a = 0; a <= amfb; a++) {
	  int first = (r < a);// begin by the largest cofactor.
          methods[r][a] = strategies->uniform_strategy[first];
        }
      }
    }
  else
    {/* to precompute our method from the file.*/
      verbose_output_print(0, 1, "# Read the cofactorization strategy file\n");
      facul_method_t* precomputed_methods =
	(facul_method_t*) malloc (sizeof(facul_method_t) * NB_MAX_METHODS);
      precomputed_methods[0].method = 0;
      precomputed_methods[0].plan = NULL;
      
      strategies->precomputed_methods = precomputed_methods;
      unsigned int index_strategies[2] = {0,0};
      
      char line[10000];
      
      fseek (file, 0, SEEK_SET);
      while (fgets (line, sizeof(line), file) != NULL)
	{
	  // process each line of 'file'
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
  for (unsigned int r = 1; r <= rmfb; r++)
    for (unsigned int a = 1; a <= amfb; a++){
      facul_method_side_t* methods =
	strategies->methods[r][a];
      if (methods == NULL)
	continue;
      int index_last_method[2] = {-1, -1};      // the default_values
      unsigned int i;
      for (i = 0; methods[i].method != 0; i++) {
	methods[i].is_the_last = 0;
	index_last_method[methods[i].side] = i;
      }
      if (i > max_curves_used_before_aux)
          max_curves_used_before_aux = i;
      // == -1 if it exists zero method for this side
      if (index_last_method[0] != -1)
	  methods[index_last_method[0]].is_the_last = 1;
      if (index_last_method[1] != -1)
	  methods[index_last_method[1]].is_the_last = 1;
    }
  
  return strategies;
}


void
facul_clear_strategies (facul_strategies_t *strategies)
{
  if (strategies == NULL)
    return ;

  int uniform = (strategies->uniform_strategy[0] != NULL);

  // free precomputed_methods
  facul_method_t* fm = strategies->precomputed_methods;
  for (int j = 0; fm[j].method!=0; j++)
    {
      if (fm[j].plan == NULL)
	continue;
      if (fm[j].method == PM1_METHOD)
	pm1_clear_plan ((pm1_plan_t*) fm[j].plan);
      else if (fm[j].method == PP1_27_METHOD ||
	       fm[j].method == PP1_65_METHOD)
	pp1_clear_plan ((pp1_plan_t*) fm[j].plan);
      else if (fm[j].method == EC_METHOD)
	ecm_clear_plan ((ecm_plan_t*) fm[j].plan);
      free (fm[j].plan);
      fm[j].method = 0;
      fm[j].plan = NULL;
    }
  free (fm);
  fm = NULL;

  // free methods
  unsigned int r;
  for (r = 0; r <= strategies->mfb[0]; r++) {
    unsigned int a;
    if (!uniform) {
        for (a = 0; a <= strategies->mfb[1]; a++)
            free (strategies->methods[r][a]);
    }
    free (strategies->methods[r]);
  }
  free (strategies->methods);
  if (uniform) {
      free(strategies->uniform_strategy[0]);
      free(strategies->uniform_strategy[1]);
  }

  // free strategies
  free (strategies);

}


int
facul_fprint_strategies (FILE* file, facul_strategies_t* strategies)
{
  if (file == NULL)
    return -1;
  // print info lpb ...
  fprintf (file,
	  "(lpb = [%ld,%ld], as...=[%lf, %lf], BBB = [%lf, %lf])\n",
	  strategies->lpb[0], strategies->lpb[1],
	  strategies->assume_prime_thresh[0],
	  strategies->assume_prime_thresh[1],
	  strategies->BBB[0], strategies->BBB[1]);
  fprintf (file,
	  "mfb = [%d, %d]\n", strategies->mfb[0], strategies->mfb[1]);
  for (unsigned int r = 0; r <= strategies->mfb[0]; r++) {
    for (unsigned int a = 0; a <= strategies->mfb[1]; a++) {
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
	  if (fm->plan == NULL) // zero method
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
	    else if (fm->method == PP1_27_METHOD ||
		     fm->method == PP1_65_METHOD)
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

void
modset_get_z (mpz_t z, const struct modset_t *modset)
{
  ASSERT_ALWAYS(modset->arith != modset_t::CHOOSE_NONE);
  switch (modset->arith)
    {
    case modset_t::CHOOSE_UL:
      mpz_set_ui (z, modset->m_ul->m);
      break;
    case modset_t::CHOOSE_15UL:
      mpz_set_ui (z, modset->m_15ul->m[1]);
      mpz_mul_2exp (z, z, LONG_BIT);
      mpz_add_ui (z, z, modset->m_15ul->m[0]);
      break;
    case modset_t::CHOOSE_2UL2:
      mpz_set_ui (z, modset->m_2ul2->m[1]);
      mpz_mul_2exp (z, z, LONG_BIT);
      mpz_add_ui (z, z, modset->m_2ul2->m[0]);
      break;
    case modset_t::CHOOSE_MPZ:
      mpz_set (z, modset->m_mpz);
      break;
    default:
      ASSERT_ALWAYS(0);
    }
}

/*
 * This is our auxiliary factorization.
 * It applies a bunch of ECM curves with larger bounds to find
 * a factor with high probability. It returns -1 if the factor
 * is not smooth, otherwise the number of
 * factors.
 */
static int
facul_aux (std::vector<cxx_mpz> & factors, const modset_t m,
	   const facul_strategies_t *strategies,
	   facul_method_side_t *methods, int method_start, int side)
{
  int found = 0;
  if (methods == NULL)
    return found;

  for (int i = method_start; methods[i].method != NULL; i++)
    {
      if (methods[i].side != side)
	continue; /* this method is not for this side */

#ifdef ENABLE_UNSAFE_FACUL_STATS
      if (i < STATS_LEN)
	stats_called_aux[i]++;
#endif  /* ENABLE_UNSAFE_FACUL_STATS */
      modset_t fm, cfm;

      int res_fac = 0;

      switch (m.arith) {
          case modset_t::CHOOSE_UL:
	res_fac = facul_doit_onefm_ul(factors, m.m_ul,
				      *methods[i].method, &fm, &cfm,
				      strategies->lpb[side],
				      strategies->assume_prime_thresh[side],
				      strategies->BBB[side]);
	break;
          case modset_t::CHOOSE_15UL:
	res_fac = facul_doit_onefm_15ul(factors, m.m_15ul,
					*methods[i].method, &fm, &cfm,
					strategies->lpb[side],
					strategies->assume_prime_thresh[side],
					strategies->BBB[side]);
	break;
          case modset_t::CHOOSE_2UL2:
	res_fac = facul_doit_onefm_2ul2 (factors, m.m_2ul2,
					 *methods[i].method, &fm, &cfm,
					 strategies->lpb[side],
					 strategies->assume_prime_thresh[side],
					 strategies->BBB[side]);
	break;
          case modset_t::CHOOSE_MPZ:
	res_fac = facul_doit_onefm_mpz (factors, m.m_mpz,
					*methods[i].method, &fm, &cfm,
					strategies->lpb[side],
					strategies->assume_prime_thresh[side],
					strategies->BBB[side]);
	break;
      default:
	ASSERT_ALWAYS(0);
      }
      // check our result
      // res_fac contains the number of factors found
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
	     tries with an other method */
	    continue;
	}

      found += res_fac;
      if (res_fac == 2)
	break;

      /*
	res_fac == 1  Only one factor has been found. Hence, our
	factorization is not finished.
      */
      if (fm.arith != modset_t::CHOOSE_NONE)
	{
	  int found2 = facul_aux (factors, fm, strategies,
				  methods, i+1, side);
          if (found2 < 1)// FACUL_NOT_SMOOTH or FACUL_MAYBE
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
      if (cfm.arith != modset_t::CHOOSE_NONE)
	{
	  int found2 = facul_aux (factors, cfm, strategies,
				  methods, i+1, side);
          if (found2 < 1)// FACUL_NOT_SMOOTH or FACUL_MAYBE
          {
              found = FACUL_NOT_SMOOTH;
          }
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

static std::array<int, 2>
facul_both_src (std::array<std::vector<cxx_mpz>, 2> & factors, const modset_t* m,
		const facul_strategies_t *strategies, int* cof,
		int* is_smooth)
{
    std::array<int, 2> found;
    found.fill(0);

  facul_method_side_t* methods = strategies->methods[cof[0]][cof[1]];

  if (methods == NULL)
    return found;

  modset_t f[2][2];
  f[0][0].arith = modset_t::CHOOSE_NONE;
  f[0][1].arith = modset_t::CHOOSE_NONE;
  f[1][0].arith = modset_t::CHOOSE_NONE;
  f[1][1].arith = modset_t::CHOOSE_NONE;
#ifdef ENABLE_UNSAFE_FACUL_STATS
  int stats_nb_side = 0, stats_index_transition = 0;
#endif  /* ENABLE_UNSAFE_FACUL_STATS */
  int last_i[2]; /* last_i[s] is the index of last method tried on side s */
  for (int i = 0; methods[i].method != NULL; i++)
    {
#ifdef ENABLE_UNSAFE_FACUL_STATS
      stats_current_index = i - stats_nb_side * stats_index_transition;
      if (methods[i].is_the_last)
	{
	  stats_nb_side = 1;
	  stats_index_transition = i+1;
	}
#endif  /* ENABLE_UNSAFE_FACUL_STATS */
      int side = methods[i].side;
      if (is_smooth[side] != FACUL_MAYBE)
	{
	  /* If both sides are smooth, we can exit the loop,
	     otherwise we must continue with the next methods,
	     since methods might be interleaved between side 0 and 1,
	     thus we don't have an easy way to skip all methods for this side.
	     We could do this with another representation, say methods[0][i]
	     for side 0, 0 <= i < m, methods[1][j] for side 1, 0 <= j < n,
	     and which_method[k] = {0, 1} for 0 <= k < m+n. */
	  if (is_smooth[side] == FACUL_SMOOTH &&
              is_smooth[1-side] == FACUL_SMOOTH)
            break;
	  continue;
	}

#ifdef ENABLE_UNSAFE_FACUL_STATS
      if (stats_current_index < STATS_LEN)
	stats_called[stats_current_index]++;
#endif  /* ENABLE_UNSAFE_FACUL_STATS */
      int res_fac = 0;
      last_i[side] = i;
      switch (m[side].arith) {
          case modset_t::CHOOSE_UL:
	res_fac = facul_doit_onefm_ul(factors[side], m[side].m_ul,
				      *methods[i].method, &f[side][0], &f[side][1],
				      strategies->lpb[side],
				      strategies->assume_prime_thresh[side],
				      strategies->BBB[side]);
	break;
          case modset_t::CHOOSE_15UL:
	res_fac = facul_doit_onefm_15ul(factors[side], m[side].m_15ul,
					*methods[i].method, &f[side][0], &f[side][1],
					strategies->lpb[side],
					strategies->assume_prime_thresh[side],
					strategies->BBB[side]);
	break;
          case modset_t::CHOOSE_2UL2:
	res_fac = facul_doit_onefm_2ul2 (factors[side], m[side].m_2ul2,
					 *methods[i].method, &f[side][0], &f[side][1],
					 strategies->lpb[side],
					 strategies->assume_prime_thresh[side],
					 strategies->BBB[side]);
	break;
          case modset_t::CHOOSE_MPZ:
	res_fac = facul_doit_onefm_mpz (factors[side], m[side].m_mpz,
					*methods[i].method, &f[side][0], &f[side][1],
					strategies->lpb[side],
					strategies->assume_prime_thresh[side],
					strategies->BBB[side]);
	break;
      default:
	ASSERT_ALWAYS(0);
      }
      // check our result
      // res_fac contains the number of factors found, or -1 if not smooth
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
            Indeed, if using only one factoring method we found two
            prime factors of m (f, m/f) then m is factored and work is
            finished for this cofactor.
	  */
	  is_smooth[side] = FACUL_SMOOTH;
	  continue;
	}
      /*
	res_fac == 1. Only one factor has been found. Hence, an
	auxiliary factorization will be necessary.
       */
      is_smooth[side] = FACUL_AUX;
    }
  // begin the auxiliary factorization
  if (is_smooth[0] >= 1 && is_smooth[1] >= 1)
    for (int side = 0; side < 2; side++)
      if (is_smooth[side] == FACUL_AUX)
	for (int ind_cof = 0; ind_cof < 2; ind_cof++)
	  {
	    // factor f[side][0] or/and f[side][1]
	    if (f[side][ind_cof].arith != modset_t::CHOOSE_NONE)
	      {
		int found2 = facul_aux (factors[side],
					f[side][ind_cof], strategies,
					methods, last_i[side] + 1, side);
		if (found2 < 1)// FACUL_NOT_SMOOTH or FACUL_MAYBE
		  {
		    is_smooth[side] = FACUL_NOT_SMOOTH;
		    found[side] = found2;// FACUL_NOT_SMOOTH or FACUL_MAYBE
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
  This function is like facul, but we will work with both norms
  together.  It returns the number of factors for each side.
*/
std::array<int, 2>
facul_both (std::array<std::vector<cxx_mpz>, 2> & factors,
            std::array<cxx_mpz, 2> & N,
	    const facul_strategies_t *strategies, int* is_smooth)
{
  int cof[2];
  size_t bits;
  std::array<int, 2> found;
  found.fill(0);

  modset_t n[2];
  n[0].arith = modset_t::CHOOSE_NONE;
  n[1].arith = modset_t::CHOOSE_NONE;

#ifdef PARI
  gmp_fprintf (stderr, "(%Zd %Zd)", N[0], N[1]);
#endif

  /* cofactors should be positive */
  ASSERT (mpz_sgn (N[0]) > 0 && mpz_sgn (N[1]) > 0);

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
	return found;

      /* Use the fastest modular arithmetic that's large enough for
	 this input */
      if (bits <= MODREDCUL_MAXBITS)
	{
	  ASSERT(mpz_fits_ulong_p(N[side]));
	  modredcul_initmod_ul (n[side].m_ul, mpz_get_ui(N[side]));
	  n[side].arith = modset_t::CHOOSE_UL;
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
	  n[side].arith = modset_t::CHOOSE_15UL;
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
	  n[side].arith = modset_t::CHOOSE_2UL2;
	}
      else
	{
	  modmpz_initmod_int (n[side].m_mpz, N[side]);
	  n[side].arith = modset_t::CHOOSE_MPZ;
	}
      ASSERT (n[side].arith != modset_t::CHOOSE_NONE);
    }

  found = facul_both_src (factors, n, strategies, cof, is_smooth);
  for (int side = 0; side < 2; side++)
    {
      if (found[side] > 1)
	{
	  /* Sort the factors we found */
          ASSERT_ALWAYS(factors[side].size() == (size_t) found[side]);
          std::sort(factors[side].begin(), factors[side].end(), cxx_mpz_cmp());
	}
    }

  // Free
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
      case modset_t::CHOOSE_NONE: /* already clear */
    break;
      case modset_t::CHOOSE_UL:
    modredcul_clearmod (modset->m_ul);
    break;
      case modset_t::CHOOSE_15UL:
    modredc15ul_clearmod (modset->m_15ul);
    break;
      case modset_t::CHOOSE_2UL2:
    modredc2ul2_clearmod (modset->m_2ul2);
    break;
      case modset_t::CHOOSE_MPZ:
    modmpz_clearmod (modset->m_mpz);
    break;
  default:
    ASSERT_ALWAYS(0);
  }
  modset->arith = modset_t::CHOOSE_NONE;
}
