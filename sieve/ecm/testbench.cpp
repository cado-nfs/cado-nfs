/* testbench - test program for different factoring strategies

   Examples:

   $ primes 3 1000000 | ./testbench -inp /dev/stdin -cof 4294967311 \
                                    -ecm 315 5355 6
   Strategy has 1 method(s)
   Of 78498 tries there were 64684 with a factor found
   Ratio: 0.824021
   Total time: 4.25 s, per call: 54.095060 usec, per factor: 65.647672 usec

   $ primes 3 1000000 | ./testbench -inp /dev/stdin -cof 4294967311 -pm1 \
                                    315 3000 -ecm 315 5355 6 -ecm 315 5355 -2
   Strategy has 3 method(s)
   Of 78498 tries there were 77224 with a factor found
   Ratio: 0.983770
   Total time: 3.10 s, per call: 39.472700 usec, per factor: 40.123899 usec

To test the ul arithmetic on a 64-bit processor:

./testbench -cof 4294967311 -p -ecm 315 5355 11 3 10000000

To test the 15ul arithmetic on a 64-bit processor:

./testbench -cof 18446744073709551629 -p -ecm 315 5355 11 3 10000000

To test the 2ul2 arithmetic on a 64-bit processor:

./testbench -cof 79228162514264337593543950397 -p -ecm 315 5355 11 3 10000000

To test the mpz arithmetic on a 64-bit processor:

./testbench -cof 340282366920938463463374607431768211507 -p -ecm 315 5355 11 3 10000000

*/

#include "cado.h" // IWYU pragma: keep
// IWYU pragma: no_include <ext/alloc_traits.h>
#include <cstdint>     /* AIX wants it first (it's a bug) */
#include <climits>              // for ULONG_MAX
#include <cstdio>               // for printf, NULL, fprintf, stderr, fclose
#include <cstdlib>              // for strtoul, malloc, exit, free, strtol
#include <cstring>              // for strcmp, strncmp
#include <memory>                // for allocator_traits<>::value_type
#include <vector>                // for vector
#include <gmp.h>                 // for gmp_printf, mpz_srcptr, mpz_get_ui
#include "cxx_mpz.hpp"
#include "facul.hpp"             // for facul_strategy_t, facul, facul_clear...
#include "facul_strategies_stats.hpp"  // for facul_strategy_t, facul, facul_clear...
#include "facul_ecm.h"           // for ec_parameter_is_valid, BRENT12, ec_p...
#include "facul_fwd.hpp"         // for facul_method_t
#include "macros.h"              // for ASSERT
#include "modredc_ul.h"          // for modredcul_clearmod, modredcul_initmo...
#include "modredc_ul_default.h"  // for modulus_t
#include "pm1.h"                 // for pm1_make_plan, pm1_plan_t
#include "pp1.h"                 // for pp1_make_plan, pp1_plan_t
#include "timing.h"  // microseconds

#define MAX_METHODS 20

const char *method_name[] = {"P-1", "P+1", "ECM"};

static void
print_pointorder (const unsigned long p, const unsigned long parameter,
                  ec_parameterization_t parameterization, const int verbose)
{
  modulus_t m;
  unsigned long o, knownfac;

  modredcul_initmod_ul (m, p);
  
  if (parameterization & ECM_TORSION12)
    knownfac = 12;
  else if (parameterization & ECM_TORSION16)
    knownfac = 16;
  else
    knownfac = 1;
  
  o = ec_parameterization_point_order_ul (parameterization, parameter, knownfac,
                                          0, m, verbose);
  if (verbose)
    printf ("%lu %lu\n", p, o);

  modredcul_clearmod (m);
}


static int
tryfactor (cxx_mpz const & N, facul_strategy_oneside const & strategy, 
           const int verbose, const int printfactors, const int printnonfactors, 
           const int printcofactors)
{
  std::vector<cxx_mpz> f;
  int facul_code;

  if (verbose >= 2)
    gmp_printf ("Trying to factor %Zd\n", (mpz_srcptr) N);

  facul_code = facul (f, N, strategy);
  
  if (verbose >= 2)
    gmp_printf ("facul returned code %d\n", facul_code);

  if (printfactors && facul_code > 0)
    {
      int j;
      gmp_printf ("%Zd", (mpz_srcptr) f[0]);
      for (j = 1; j < facul_code; j++)
        gmp_printf (" %Zd", (mpz_srcptr) f[j]);
      if (printcofactors) {
        cxx_mpz c;
        mpz_set (c, N);
        for (j = 0; j < facul_code; j++)
          mpz_divexact (c, c, f[j]);
        if (mpz_cmp_ui (c, 1) != 0)
          gmp_printf (" %Zd", (mpz_srcptr) c);
      }
      printf ("\n");
    }
  
  if (printnonfactors && (facul_code == 0 || facul_code == FACUL_NOT_SMOOTH))
    {
      gmp_printf ("%Zd\n", (mpz_srcptr) N);
    }
  return facul_code;
}

static void print_help (char *programname)
{
  printf ("%s [options] [<start> <stop>]\n", programname);
  printf ("Run factoring method on numbers in interval [<start>, <stop>], or from file\n");
  printf ("<start> and <stop> must be given unless -inp or -inpraw is used\n");
  printf ("-pm1 <B1> <B2>      Run P-1 with B1, B2\n");
  printf ("-pp1_27 <B1> <B2>   Run P+1 with B1, B2, x0=2/7\n");
  printf ("-pp1_65 <B1> <B2>   Run P+1 with B1, B2, x0=6/5\n");
  printf ("-ecm <B1> <B2> <s>  Run ECM with B1, B2 and parameter s, using BRENT12 curve\n");
  printf ("-ecmm12 <B1> <B2> <s>  Same, but using Montgomery torsion 12 curve\n");
  printf ("-ecmm16 <B1> <B2> <s>  Same, but using Montgomery torsion 16 curve\n");
  printf ("-ecmem12 <B1> <B2> <s>  Same, but using Edwards torsion 12 curve\n");
  printf ("-strat   Use the facul default strategy. Don't use with -pm1, -pp1, -ecm\n");
  printf ("-fbb <n> Use <n> as factor base bound, e.g. for primality checks\n");
  printf ("-lpb <n> Use <n> as large prime bound, e.g. for early abort\n");
  printf ("-ep      Add certain extra primes in ECM stage 1 (e.g., 12 or 16)\n");
  printf ("         Affects only curves specified after the -ep parameter\n");
  printf ("-p       Try only primes in [<start>, <stop>] (default: all odd "
	  "numbers)\n");
  printf ("-q       Suppress normal output, output from -v, -vf and -vnf still appears\n");
  printf ("-v       More verbose output. Once: print statistics. Twice: print \n"
          "         numbers that are tried. Three: internal info from strategies\n");
  printf ("-vf      Print factors that are found\n");
  printf ("-vnf     Print input numbers where no factor was found\n");
  printf ("-vcf     Print cofactor if any factors were found\n");
  printf ("-cof <n> Multiply each number to test by <num> before calling facul.\n"
	  "         NOTE: facul does not report input numbers as factors,\n"
	  "         with -p you MUST use -cof or no factors will be found\n");
  printf ("-m <n>   Choose modulus, do separate counts of p %% modulus\n");
  printf ("-inp <f> Read decimal numbers to factor from file <f>\n");
  printf ("-inpraw <f>  Read numbers in GMP raw format from file <f>\n");
  printf ("-inpstop <n> Stop after reading <n> numbers\n");
  printf ("-po <s>, -pom12 <s>, -pom16 <s>, -poem12 <s> Compute order of starting point. Use -vf\n");
}

static unsigned long
next_prime (unsigned long start)
{
  mpz_t p;

  mpz_init_set_ui (p, start - 1);
  mpz_nextprime (p, p);
  start = mpz_get_ui (p);
  mpz_clear (p);
  return start;
}

// coverity[root_function]
int main (int argc, char **argv)
{
  char *argv0 = argv[0];
  unsigned long start, stop, i, mod = 0UL, inpstop = ULONG_MAX;
  unsigned long hits = 0, total = 0;
  unsigned long fbb = 0, lpb = ULONG_MAX;
  char *inp_fn = NULL;
  FILE *inp;
  cxx_mpz N, cof;
  int nr_methods = 0;
  int only_primes = 0, verbose = 0, quiet = 0;
  int printfactors = 0;
  int printnonfactors = 0;
  int printcofactors = 0;
  int do_pointorder = 0;
  unsigned long po_parameter = 0;
  ec_parameterization_t po_parameterization = BRENT12;
  int inp_raw = 0;
  int strat = 0;
  int extra_primes = 0;
  int ncurves = -1;
  unsigned long *primmod = NULL, *hitsmod = NULL;
  uint64_t starttime, endtime;
  std::vector<facul_method::parameters> method_params;

  /* Parse options */
  mpz_set_ui (cof, 1UL);
  for ( ; argc > 1 && argv[1][0] == '-' ; ) {
      if (argc > 1 && strcmp (argv[1], "-h") == 0) {
	  print_help (argv0);
	  return 0;
      } else if (argc > 3 && strcmp (argv[1], "-pm1") == 0) {
	  unsigned long B1, B2;
	  B1 = strtoul (argv[2], NULL, 10);
	  B2 = strtoul (argv[3], NULL, 10);
          method_params.emplace_back(PM1_METHOD, B1, B2);
	  argc -= 3;
	  argv += 3;
      } else if (argc > 3 && strcmp (argv[1], "-pp1_27") == 0) {
	  unsigned long B1, B2;
	  B1 = strtoul (argv[2], NULL, 10);
	  B2 = strtoul (argv[3], NULL, 10);
          method_params.emplace_back(PP1_27_METHOD, B1, B2);
	  argc -= 3;
	  argv += 3;
      } else if (argc > 3 && strcmp (argv[1], "-pp1_65") == 0) {
	  unsigned long B1, B2;
	  B1 = strtoul (argv[2], NULL, 10);
	  B2 = strtoul (argv[3], NULL, 10);
          method_params.emplace_back(PP1_65_METHOD, B1, B2);
	  argc -= 3;
	  argv += 3;
      } else if (argc > 4 && strncmp (argv[1], "-ecm", 4) == 0) {
	  unsigned long B1, B2;
	  unsigned long parameter;
	  ec_parameterization_t parameterization;
	  B1 = strtoul (argv[2], NULL, 10);
	  B2 = strtoul (argv[3], NULL, 10);
	  parameter = strtol (argv[4], NULL, 10);
	  if (strcmp (argv[1], "-ecm") == 0) {
	    parameterization = BRENT12;
          } else if (strcmp (argv[1], "-ecmm12") == 0) {
              parameterization = MONTY12;
          } else if (strcmp (argv[1], "-ecmm16") == 0) {
              parameterization = MONTY16;
          } else if (strcmp (argv[1], "-ecmem12") == 0) {
              parameterization = MONTYTWED12;
          } else {
              fprintf (stderr, "Unrecognized option: %s\n", argv[1]);
              exit (EXIT_FAILURE);
          }
          method_params.emplace_back(EC_METHOD, B1, B2,
                  parameterization, parameter, extra_primes);
	  argc -= 4;
	  argv += 4;
      } else if (argc > 2 && strcmp (argv[1], "-ncurves") == 0 && method_params.empty()) {
	  ncurves = strtoul (argv[2], NULL, 10);
	  argc -= 2;
	  argv += 2;
      } else if (argc > 1 && strcmp (argv[1], "-strat") == 0 && method_params.empty()) {
	  strat = 1;
	  argc -= 1;
	  argv += 1;
      } else if (argc > 2 && strncmp (argv[1], "-po", 3) == 0) {
	  do_pointorder = 1;
          if (strcmp (argv[1], "-po") == 0) {
              po_parameterization = BRENT12;
          } else if (strcmp (argv[1], "-pom12") == 0) {
              po_parameterization = MONTY12;
          } else if (strcmp (argv[1], "-pom16") == 0) {
              po_parameterization = MONTY16;
          } else if (strcmp (argv[1], "-poem12") == 0) {
              po_parameterization = MONTYTWED12;
          } else {
              fprintf (stderr, "Unrecognized option: %s\n", argv[1]);
              exit (EXIT_FAILURE);
          }
	  po_parameter = strtol (argv[2], NULL, 10);

          if (!ec_parameter_is_valid (po_parameterization, po_parameter))
          {
              fprintf (stderr, "Parameter %lu is not valid with parametrization '%s'\n",
                      po_parameter, argv[1]);
              exit (EXIT_FAILURE);
          }

          argc -= 2;
          argv += 2;
      } else if (argc > 2 && strcmp (argv[1], "-fbb") == 0) {
          fbb = strtoul (argv[2], NULL, 10);
          argc -= 2;
          argv += 2;
      } else if (argc > 2 && strcmp (argv[1], "-lpb") == 0) {
          lpb = strtoul (argv[2], NULL, 10);
          argc -= 2;
          argv += 2;
      } else if (argc > 1 && strcmp (argv[1], "-ep") == 0 && method_params.empty()) {
          extra_primes = 1;
          argc -= 1;
          argv += 1;
      } else if (argc > 1 && strcmp (argv[1], "-p") == 0) {
          only_primes = 1;
          argc -= 1;
          argv += 1;
      } else if (argc > 1 && strcmp (argv[1], "-v") == 0) {
          verbose++;
          argc -= 1;
          argv += 1;
      } else if (argc > 1 && strcmp (argv[1], "-q") == 0) {
          quiet++;
          argc -= 1;
          argv += 1;
      } else if (argc > 1 && strcmp (argv[1], "-vf") == 0) {
          printfactors = 1;
          argc -= 1;
          argv += 1;
      } else if (argc > 1 && strcmp (argv[1], "-vnf") == 0) {
          printnonfactors = 1;
          argc -= 1;
          argv += 1;
      } else if (argc > 1 && strcmp (argv[1], "-vcf") == 0) {
          printcofactors = 1;
          argc -= 1;
          argv += 1;
      } else if (argc > 2 && strcmp (argv[1], "-inpstop") == 0) {
          inpstop = strtoul (argv[2], NULL, 10);
          argc -= 2;
          argv += 2;
      } else if (argc > 2 && strcmp (argv[1], "-m") == 0) {
          mod= strtoul (argv[2], NULL, 10);
          hitsmod = (unsigned long *) malloc (mod * sizeof (unsigned long));
          primmod = (unsigned long *) malloc (mod * sizeof (unsigned long));
          for (i = 0; i < mod; i++)
              hitsmod[i] = primmod[i]= 0;
          argc -= 2;
          argv += 2;
      } else if (argc > 2 && strcmp (argv[1], "-cof") == 0) {
          mpz_set_str (cof, argv[2], 10);
          argc -= 2;
          argv += 2;
      } else if (argc > 2 && strcmp (argv[1], "-inp") == 0) {
          inp_fn = argv[2];
          argc -= 2;
          argv += 2;
      } else if (argc > 2 && strcmp (argv[1], "-inpraw") == 0) {
          inp_fn = argv[2];
          inp_raw = 1;
          argc -= 2;
          argv += 2;
      } else {
          printf ("Unrecognized option: %s\n", argv[1]);
          exit (EXIT_FAILURE);
      }
  }
  
  if (strat && !method_params.empty())
    {
      printf ("Don't use -strat with -pm1, -pp1 or -ecm\n");
      exit (EXIT_FAILURE);
    }

  facul_strategy_oneside strategy;

  if (only_primes && inp_fn != NULL)
    fprintf (stderr, "-p has no effect with -inp or -inpraw\n");

  if (strat) {
      /* we set mfb = 3*lpb to avoid a huge number of curves if ncurves is not
         given (case of 2 large primes) */
      strategy = facul_strategy_oneside (fbb, lpb, 3 * lpb, ncurves, (verbose / 3));
  } else {
      strategy = facul_strategy_oneside (fbb, lpb, 3 * lpb, method_params, (verbose / 3));
  }
  if (!quiet) printf ("Strategy has %d method(s)\n", nr_methods);

  if (inp_fn == NULL)
    {
      if (argc < 3)
	{
	  print_help (argv0);
	  return 1;
	}
      
      /* Get range to test */
      start = strtoul (argv[1], NULL, 10);
      if (start % 2UL == 0UL)
	start++;
      stop = strtoul (argv[2], NULL, 10);

      if (only_primes)
        i = next_prime (start);
      else
	i = start;
      
      starttime = microseconds();

      /* The main loop */
      while (i <= stop)
	{
	  if (mod > 0)
	    primmod[i % mod]++;
	  
	  total++;
	  if (do_pointorder) {
	      print_pointorder (i, po_parameter, po_parameterization, verbose+printfactors);
	      /* TODO: check point order */
          } else {
              mpz_mul_ui (N, cof, i);
              if (tryfactor (N, strategy, verbose, printfactors, printnonfactors, printcofactors)) {
                  hits++;
                  if (mod > 0)
                    hitsmod[i % mod]++;
              }
          }

	  if (only_primes)
            i = next_prime (i + 1);
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

        if (do_pointorder)
          {
            if (mpz_fits_ulong_p (N))
              print_pointorder (mpz_get_ui (N), po_parameter, po_parameterization, printfactors);
            else
              gmp_fprintf (stderr, "%Zd does not fit into an unsigned long, not computing group order\n", (mpz_srcptr) N);
          }
        else
          {
            mpz_mul (N, N, cof);
            if (tryfactor (N, strategy, verbose, printfactors, printnonfactors, printcofactors))
              hits++;
          }
      }
    fclose (inp);
  }
  
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
              usrtime / 1000000., usrtime / (double) (total + !total), 
              usrtime / (double) (hits + !hits));
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

  if (mod > 0)
    {
      free (hitsmod);
      free (primmod);
    }

  if (verbose >= 1)
    facul_print_stats (stdout);

  return 0;
}
