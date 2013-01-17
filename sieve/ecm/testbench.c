/* testbench - test program for different factoring strategies

   Examples:

   $ primes 3 1000000 | ./testbench -inp /dev/stdin -cof 4294967311 \
                                    -ecm 315 5355 6
   Strategy has 1 methods
   Of 78498 tries there were 64684 with a factor found
   Ratio: 0.824021
   Total time: 4.25 s, per call: 54.095060 usec, per factor: 65.647672 usec

   $ primes 3 1000000 | ./testbench -inp /dev/stdin -cof 4294967311 -pm1 \
                                    315 3000 -ecm 315 5355 6 -ecm 315 5355 -2
   Strategy has 3 methods
   Of 78498 tries there were 77224 with a factor found
   Ratio: 0.983770
   Total time: 3.10 s, per call: 39.472700 usec, per factor: 40.123899 usec
*/

#include "cado.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include "facul.h"
#include "pm1.h"
#include "pp1.h"
#include "ecm.h"
#include "modredc_ul.h"
#include "modredc_ul_default.h"
#include "getprime.h"
#include "timing.h"

#define MAX_METHODS 20

const char *method_name[] = {"P-1", "P+1", "ECM"};

void
print_pointorder (const unsigned long p, const unsigned long s, 
                  const unsigned long parameterization, const int verbose)
{
  residue_t sigma;
  modulus_t m;
  unsigned long o;

  modredcul_initmod_ul (m, p);
  modredcul_init (sigma, m);
  modredcul_set_ul (sigma, s, m);
      
  o = ell_pointorder_ul (sigma, parameterization, 12, 0, m, verbose);
  if (verbose)
    printf ("%lu %lu\n", p, o);

  modredcul_clear (sigma, m);  
  modredcul_clearmod (m);
}

static int
tryfactor (mpz_t N, const facul_strategy_t *strategy, 
           const int verbose, const int printfactors, const int printnonfactors)
{
  unsigned long f[16];
  int facul_code;

  if (verbose >= 2)
    gmp_printf ("Trying to factor %Zd\n", N);

  facul_code = facul (f, N, strategy);
  
  if (printfactors && facul_code > 0)
    {
      int j;
      for (j = 0; j < facul_code; j++)
        printf ("%lu ", f[j]);
      printf ("\n");
    }
  
  if (printnonfactors && facul_code == 0)
    {
      int j;
      for (j = 0; j < facul_code; j++)
        mpz_tdiv_q_ui (N, N, f[j]);
      gmp_printf ("%Zd\n", N);
    }
  
  return facul_code;
}

void print_help (char *programname)
{
  printf ("%s [options] [<start> <stop>]\n", programname);
  printf ("Run factoring method on numbers in interval [<start>, <stop>], or from file\n");
  printf ("<start> and <stop> must be given unless -inp or -inpraw is used\n");
  printf ("-pm1 <B1> <B2>      Run P-1 with B1, B2\n");
  printf ("-pp1 <B1> <B2>      Run P+1 with B1, B2\n");
  printf ("-ecm <B1> <B2> <s>  Run ECM with B1, B2 and parameter s, using BRENT12 curve\n");
  printf ("-ecmm12 <B1> <B2> <s>  Same, but using Montgomery torsion 12 curve\n");
  printf ("-ecmm16 <B1> <B2> <s>  Same, but using Montgomery torsion 16 curve\n");
  printf ("-strat   Use the facul default strategy. Don't use with -pm1, -pp1, -ecm\n");
  printf ("-fbb <n> Use <n> as factor base bound, e.g. for primality checks\n");
  printf ("-lpb <n> Use <n> as large prime bound, e.g. for early abort\n");
  printf ("-ep      Add certain extra primes in ECM stage 1 (e.g., 12 or 16)\n");
  printf ("-p       Try only primes in [<start>, <stop>] (default: all odd "
	  "numbers)\n");
  printf ("-q       Suppress normal output, output from -v, -vf and -vnf still appears\n");
  printf ("-v       More verbose output. Once: print statistics. Twice: print \n"
          "         numbers that are tried. Three: internal info from strategies\n");
  printf ("-vf      Print factors that are found\n");
  printf ("-vnf     Print input numbers where no factor was found\n");
  printf ("-cof <n> Multiply each number to test by <num> before calling facul.\n"
	  "         NOTE: facul does not report input numbers as factors,\n"
	  "         with -p you MUST use -cof or no factors will be found\n");
  printf ("-m <n>   Choose modulus, do separate counts of p %% modulus\n");
  printf ("-inp <f> Read decimal numbers to factor from file <f>\n");
  printf ("-inpraw <f>  Read numbers in GMP raw format from file <f>\n");
  printf ("-inpstop <n> Stop after reading <n> numbers\n");
}

int main (int argc, char **argv)
{
  unsigned long start, stop, i, mod = 0UL, inpstop = ULONG_MAX;
  unsigned long hits = 0, total = 0;
  unsigned long fbb = 0, lpb = ~(0UL);
  char *inp_fn = NULL;
  FILE *inp;
  mpz_t N, cof;
  facul_strategy_t *strategy;
  int nr_methods = 0;
  int only_primes = 0, verbose = 0, quiet = 0;
  int printfactors = 0;
  int printnonfactors = 0;
  int do_pointorder = 0;
  unsigned long po_sigma = 0, po_parameterization = 0;
  int inp_raw = 0;
  int strat = 0;
  int extra_primes = 0;
  unsigned long *primmod = NULL, *hitsmod = NULL;
  uint64_t starttime, endtime;

  strategy = malloc (sizeof(facul_strategy_t));
  strategy->methods = malloc ((MAX_METHODS + 1) * sizeof (facul_method_t));
  strategy->lpb = ~(0UL);
  strategy->fbb2[0] = 0UL;
  strategy->fbb2[1] = 0UL;

  /* Parse options */
  mpz_init (N);
  mpz_init (cof);
  mpz_set_ui (cof, 1UL);
  while (argc > 1 && argv[1][0] == '-')
    {
      if (argc > 1 && strcmp (argv[1], "-h") == 0)
	{
	  print_help (argv[0]);
	  return 0;
	}
      else if (argc > 3 && strcmp (argv[1], "-pm1") == 0 && 
	       nr_methods < MAX_METHODS)
	{
	  unsigned long B1, B2;
	  B1 = strtoul (argv[2], NULL, 10);
	  B2 = strtoul (argv[3], NULL, 10);
	  strategy->methods[nr_methods].method = PM1_METHOD;
	  strategy->methods[nr_methods].plan = malloc (sizeof (pm1_plan_t));
	  ASSERT (strategy->methods[nr_methods].plan != NULL);
	  pm1_make_plan (strategy->methods[nr_methods].plan, B1, B2, 
			 (verbose >= 3));
	  nr_methods++;
	  argc -= 3;
	  argv += 3;
	}
      else if (argc > 3 && strcmp (argv[1], "-pp1") == 0 && 
	       nr_methods < MAX_METHODS)
	{
	  unsigned long B1, B2;
	  B1 = strtoul (argv[2], NULL, 10);
	  B2 = strtoul (argv[3], NULL, 10);
	  strategy->methods[nr_methods].method = PP1_METHOD;
	  strategy->methods[nr_methods].plan = malloc (sizeof (pp1_plan_t));
	  ASSERT (strategy->methods[nr_methods].plan != NULL);
	  pp1_make_plan (strategy->methods[nr_methods].plan, B1, B2, 
			 (verbose >= 3));
	  nr_methods++;
	  argc -= 3;
	  argv += 3;
	}
      else if (argc > 4 && strncmp (argv[1], "-ecm", 4) == 0 && 
	       nr_methods < MAX_METHODS)
	{
	  unsigned long B1, B2;
	  long sigma;
	  int parameterization = BRENT12;
	  B1 = strtoul (argv[2], NULL, 10);
	  B2 = strtoul (argv[3], NULL, 10);
	  sigma = strtol (argv[4], NULL, 10);
	  if (strcmp (argv[1], "-ecmm12") == 0)
	    parameterization = MONTY12;
	  if (strcmp (argv[1], "-ecmm16") == 0)
	    parameterization = MONTY16;
          if (parameterization == MONTY12 && sigma == 1)
            {
              fprintf (stderr, "Parameter 1 does not lead to a valid curve. "
              "Use parameter > 1 for MONTY12 curves.\n");
              exit (EXIT_FAILURE);
            }
          if (parameterization == MONTY16 && (sigma < 0 || sigma > 1))
            {
              fprintf (stderr, "Only parameter 1 for MONTY16 curves so far.\n");
              exit (EXIT_FAILURE);
            }
	  strategy->methods[nr_methods].method = EC_METHOD;
	  strategy->methods[nr_methods].plan = malloc (sizeof (ecm_plan_t));
	  ASSERT (strategy->methods[nr_methods].plan != NULL);
	  ecm_make_plan (strategy->methods[nr_methods].plan, B1, B2, 
                         parameterization, labs (sigma), extra_primes, 
                         (verbose >= 3));
	  nr_methods++;
	  argc -= 4;
	  argv += 4;
	}
      else if (argc > 1 && strcmp (argv[1], "-strat") == 0 && 
	       nr_methods == 0)
        {
	  strat = 1;
	  argc -= 1;
	  argv += 1;
	}
     else if (argc > 2 && strncmp (argv[1], "-po", 3) == 0)
	{
	  do_pointorder = 1;
	  po_parameterization = BRENT12;
	  if (strcmp (argv[1], "-pom12") == 0)
	    po_parameterization = MONTY12;
	  if (strcmp (argv[1], "-pom16") == 0)
	    po_parameterization = MONTY16;
	  po_sigma = strtol (argv[2], NULL, 10);
	  argc -= 2;
	  argv += 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-fbb") == 0)
	{
	  fbb = strtoul (argv[2], NULL, 10);
	  ularith_mul_ul_ul_2ul (&(strategy->fbb2[0]), &(strategy->fbb2[1]), 
	                         fbb, fbb);
	  argc -= 2;
	  argv += 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-lpb") == 0)
	{
	  lpb = strtoul (argv[2], NULL, 10);
	  argc -= 2;
	  argv += 2;
	}
     else if (argc > 1 && strcmp (argv[1], "-ep") == 0)
	{
	  extra_primes = 1;
	  argc -= 1;
	  argv += 1;
	}
     else if (argc > 1 && strcmp (argv[1], "-p") == 0)
	{
	  only_primes = 1;
	  argc -= 1;
	  argv += 1;
	}
      else if (argc > 1 && strcmp (argv[1], "-v") == 0)
	{
	  verbose++;
	  argc -= 1;
	  argv += 1;
	}
      else if (argc > 1 && strcmp (argv[1], "-q") == 0)
	{
	  quiet++;
	  argc -= 1;
	  argv += 1;
	}
      else if (argc > 1 && strcmp (argv[1], "-vf") == 0)
	{
	  printfactors = 1;
	  argc -= 1;
	  argv += 1;
	}
      else if (argc > 1 && strcmp (argv[1], "-vnf") == 0)
	{
	  printnonfactors = 1;
	  argc -= 1;
	  argv += 1;
	}
      else if (argc > 2 && strcmp (argv[1], "-inpstop") == 0)
	{
	  inpstop = strtoul (argv[2], NULL, 10);
	  argc -= 2;
	  argv += 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-m") == 0)
	{
	  mod= strtoul (argv[2], NULL, 10);
	  hitsmod = (unsigned long *) malloc (mod * sizeof (unsigned long));
	  primmod = (unsigned long *) malloc (mod * sizeof (unsigned long));
	  for (i = 0; i < mod; i++)
	    hitsmod[i] = primmod[i]= 0;
	  argc -= 2;
	  argv += 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-cof") == 0)
	{
	  mpz_set_str (cof, argv[2], 10);
	  argc -= 2;
	  argv += 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-inp") == 0)
	{
	  inp_fn = argv[2];
	  argc -= 2;
	  argv += 2;
	}
      else if (argc > 2 && strcmp (argv[1], "-inpraw") == 0)
	{
	  inp_fn = argv[2];
	  inp_raw = 1;
	  argc -= 2;
	  argv += 2;
	}
      else
        {
	  printf ("Unrecoglized option: %s\n", argv[1]);
	  exit (EXIT_FAILURE);
        }
    }
  
  if (strat && nr_methods != 0)
    {
      printf ("Don't use -strat with -pm1, -pp1 or -ecm\n");
      exit (EXIT_FAILURE);
    }

  if (only_primes && inp_fn != NULL)
    fprintf (stderr, "-p has no effect with -inp or -inpraw\n");

  if (strat)
    {
      facul_clear_strategy (strategy);
      strategy = facul_make_strategy (15, fbb, (lpb == 0) ? 0 : 1UL << lpb);
    }
  else
    {
      if (!quiet) printf ("Strategy has %d methods\n", nr_methods);
      strategy->methods[nr_methods].method = 0;
    }

  if (inp_fn == NULL)
    {
      if (argc < 3)
	{
	  print_help (argv[0]);
	  return 1;
	}
      
      /* Get range to test */
      start = strtoul (argv[1], NULL, 10);
      if (start % 2UL == 0UL)
	start++;
      stop = strtoul (argv[2], NULL, 10);
      
      if (only_primes)
	{
	  do {
            i = getprime (1);
          } while (i < start);
	}
      else
	i = start;
      
      starttime = microseconds();

      /* The main loop */
      while (i <= stop)
	{
	  if (mod > 0)
	    primmod[i % mod]++;
	  
	  total++;
	  if (do_pointorder)
	    print_pointorder (i, po_sigma, po_parameterization, verbose);
          else
	    {
              mpz_mul_ui (N, cof, i);
              if (tryfactor (N, strategy, verbose, printfactors, printnonfactors))
                {
                  hits++;
                  if (mod > 0)
                    hitsmod[i % mod]++;
                }
            }

	  if (only_primes)
	    i = getprime (1);
	  else
	    i += 2;
	}
  } else {
    inp = fopen (inp_fn, "r");
    if (inp == NULL)
      {
	printf ("Could not open %s\n", inp_fn);
	exit (EXIT_FAILURE);
      }

    starttime = microseconds();

    /* Read lines from inp */
    while (!feof(inp) && inpstop-- > 0)
      {
        size_t read_ok;
	if (inp_raw)
	  read_ok = mpz_inp_raw (N, inp);
	else
	  read_ok = mpz_inp_str (N, inp, 0);
        if (read_ok == 0)
          break;
	if (mpz_sgn (N) <= 0)
	  continue;
	total++;

        if (do_pointorder && mpz_fits_ulong_p (N))
          print_pointorder (mpz_get_ui (N), po_sigma, po_parameterization, verbose);
        else
          {
            mpz_mul (N, N, cof);
            if (tryfactor (N, strategy, verbose, printfactors, printnonfactors))
              hits++;
          }
      }
    fclose (inp);
  }
  
  facul_clear_strategy (strategy);

  endtime = microseconds();

  if (!quiet)
    {
      printf ("Of %lu tries there were %lu with a factor found\n", 
              total, hits);
      printf ("Ratio: %f\n", (double)hits / (double)total);
    }
  
  if (!quiet && endtime > starttime)
    {
      double usrtime = endtime - starttime;
      printf ("Total time: %.2f s, per call: %f usec, per factor: %f usec\n",
              usrtime / 1000000., usrtime / (double) total, 
              usrtime / (double) hits);
    }
    
  if (!quiet && mod > 0)
    {
      printf ("Distribution of factors found over residue classes mod %lu:\n",
              mod);
      for (i = 0; i < mod; i++)
        if (hitsmod[i] > 0)
          printf ("%lu:%lu (%f)%s", 
                  i, hitsmod[i], (double)hitsmod[i] / (double) primmod[i],
                  (i < mod - 1) ? ", " : "\n");
    }

  mpz_clear (N);
  mpz_clear (cof);

  if (mod > 0)
    {
      free (hitsmod);
      free (primmod);
    }

  if (verbose >= 1)
    facul_print_stats (stdout);

  return 0;
}
