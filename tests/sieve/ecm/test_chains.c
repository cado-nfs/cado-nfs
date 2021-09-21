#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "sieve/ecm/bytecode.h"
#include "tests_common.h"
#include "getprime.h"   // prime_info
#include "macros.h"
#include "misc.h"

void
mpz_prod_primes_below_B1 (mpz_t E, unsigned int B1)
{
  mpz_set_ui (E, 1);
  prime_info pi;
  prime_info_init (pi);
  unsigned long p = getprime_mt (pi);
  ASSERT (p == 3);

  for ( ; p <= B1; p = (unsigned int) getprime_mt (pi))
  {
    for (unsigned long q = 1; q <= B1 / p; q *= p)
      mpz_mul_ui (E, E, p);
  }
  prime_info_clear (pi);
}

unsigned int
test_prac (unsigned int B1min, unsigned int B1max, int compress,
           const prac_cost_t *opcost, int verbose)
{
  unsigned int n = 0;
  /* Test all values of B1 in [B1min..B1max] */
  for (unsigned int B1 = B1min; B1 <= B1max; B1++)
  {
    if (verbose)
      printf ("##### %s with B1 = %u\n", __func__, B1);

    bytecode bc;
    mpz_t E;

    bytecode_prac_encode (&bc, B1, 0, 0, opcost, compress, verbose);

    mpz_init (E);
    mpz_prod_primes_below_B1 (E, B1);

    if (bytecode_prac_check (bc, E, verbose) != 0)
    {
      printf ("##### Test with B1 = %u failed with PRAC\n", B1);
      n++;
    }
    free (bc);
    mpz_clear (E);
  }
  return n;
}

unsigned int
test_mishmash (unsigned int B1min, unsigned int B1max, int compress,
               const mishmash_cost_t *opcost, int verbose)
{
  unsigned int n = 0;
  /* Test all values of B1 in [B1min..B1max] */
  for (unsigned int B1 = B1min; B1 <= B1max; B1++)
  {
    if (verbose)
      printf ("##### %s with B1 = %u\n", __func__, B1);

    bytecode bc;
    mpz_t E;

    bytecode_mishmash_encode (&bc, B1, 0, 0, opcost, compress, verbose);

    mpz_init (E);
    mpz_prod_primes_below_B1 (E, B1);

    if (bytecode_mishmash_check (bc, E, verbose) != 0)
    {
      printf ("##### Test with B1 = %u failed with MISHMASH\n", B1);
      n++;
    }
    free (bc);
    mpz_clear (E);
  }
  return n;
}

int main (int argc, const char **argv)
{
  unsigned int nerrors = 0;

  tests_common_cmdline (&argc, &argv, PARSE_SEED | PARSE_VERBOSE);

  int verbose = tests_common_get_verbose ();

  { /* PRAC test */
    /* Random opcost (between 0 and 16) [ 0x1p60 == double 2^60 ] */
    prac_cost_t cost;
    cost.DBL = (double) u64_random(state) / 0x1p60;
    cost.dADD = (double) u64_random(state) / 0x1p60;
    printf ("PRAC cost: DBL = %f; dADD = %f\n", cost.DBL, cost.dADD);

    /* compress prac chains */
    nerrors += test_prac (1, 1100, 1, &cost, verbose);
    /* uncompress prac chains */
    nerrors += test_prac (1, 1100, 0, &cost, verbose);
  }

  { /* mishmash test */
    /* Random opcost (between 0 and 16) [ 0x1p60 == double 2^60 ] */
    dbchain_cost_t dbchain_cost;
    memset (&dbchain_cost, 0, sizeof (dbchain_cost_t));
    precomp_cost_t precomp_cost;
    memset (&precomp_cost, 0, sizeof (precomp_cost_t));
    mishmash_cost_t cost;
    cost.DBL = (double) u64_random(state) / 0x1p60;
    cost.DBLa = (double) u64_random(state) / 0x1p60;
    cost.TPL = (double) u64_random(state) / 0x1p60;
    cost.TPLa = (double) u64_random(state) / 0x1p60;
    cost.ADD = (double) u64_random(state) / 0x1p60;
    cost.ADDa = (double) u64_random(state) / 0x1p60;
    cost.ADDd = (double) u64_random(state) / 0x1p60;
    cost.dDBL = (double) u64_random(state) / 0x1p60;
    cost.dADD = (double) u64_random(state) / 0x1p60;
    printf ("MISHMASH PRAC cost: DBL = %f ; DBLa = %f ; TPL = %f ; TPLa = %f ; "
            "ADD = %f ; ADDa = %f ; ADDd = %f ; dDBL = %f ; dADD = %f\n",
            cost.DBL, cost.DBLa, cost.TPL, cost.TPLa, cost.ADD, cost.ADDa,
            cost.ADDd, cost.dDBL, cost.dADD);

    /* compress mishmash chains */
    nerrors += test_mishmash (1, 1100, 1, &cost, verbose);
    /* uncompress mishmash chains */
    nerrors += test_mishmash (1, 1100, 0, &cost, verbose);
  }

  bytecode_prac_cache_free ();
  tests_common_clear ();
  return (nerrors == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
