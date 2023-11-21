#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include "tests_common.h"
#include "sieve/ecm/mpqs_doit.h"
#include "misc.h"

/* Set N to a crude approximation to 2^(b/2) */
void
set_sqrt_2exp(mpz_t N, const unsigned int e)
{
  mpz_set_ui (N, 1);
  if (e <= 1) {
    /* Return round(sqrt(2)) = 1 */
  } else if (e % 2 == 1) {
    /* N = 1.5 * 2^floor(e/2) */
    mpz_mul_2exp(N, N, (e - 1) / 2 - 1);
    mpz_mul_ui (N, N, 3);
  } else {
    /* N = 2^(e/2) */
    mpz_mul_2exp(N, N, e / 2);
  }
}


int
main (int argc, const char *argv[])
{
  mpz_t N, f, *primes;
  unsigned long bits, iter = 100, i, found = 0, n_primes;
  int verbose, quiet;

  tests_common_cmdline (&argc, &argv, PARSE_SEED | PARSE_ITER | PARSE_VERBOSE | PARSE_QUIET);
  tests_common_get_iter (&iter);
  quiet = tests_common_get_quiet();
  verbose = tests_common_get_verbose();

  if (argc > 1)
    bits = atoi (argv[1]);
  else
    /* generate a random size between one and two words */
#define BITS_PER_ULONG (8 * sizeof(unsigned long))
    bits = BITS_PER_ULONG + gmp_urandomm_ui(state, BITS_PER_ULONG + 1);
  if (!quiet)
    printf ("bits=%lu iter=%lu\n", bits, iter);

  mpz_init (N);
  mpz_init (f);
  /* With n primes, we can generate n*(n-1)/2 distinct products */
  /* We want n*(n-1)/2 >= iter  <==  n >= sqrt(iter*2) + 1 */
  n_primes = ceil(sqrt(2 * iter)) + 1;
  primes = malloc(n_primes * sizeof(mpz_t));

  set_sqrt_2exp(N, bits - 1);
  for(i = 0; i < n_primes; i++) {
    mpz_nextprime (N, N);
    mpz_init_set (primes[i], N);
  }

  size_t p1 = 1, p2 = 0;

  for (i = 0; i < iter; i++)
    {
      mpz_mul (N, primes[p1], primes[p2]);
      if (verbose)
        gmp_printf ("N=%Zd\n", N);
      mpqs_doit (f, N, verbose);
      if (mpz_cmp_ui (f, 1) > 0 && mpz_cmp (f, N) < 0)
        found ++;

      if (p1 == ++p2) {
        p1++;
        if (p1 == n_primes)
          abort();
        p2 = 0;
      }
    }
  mpz_clear (N);
  mpz_clear (f);

  if (!quiet)
    printf ("Found %lu factors out of %lu composites\n", found, iter);
  smooth_stat (-1);

  for(i = 0; i < n_primes; i++) {
      mpz_clear(primes[i]);
  }

  free(primes);
  tests_common_clear ();

  return 0;
}
