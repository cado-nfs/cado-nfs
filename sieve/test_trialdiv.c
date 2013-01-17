#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "ularith.h"
#include "modredc_ul.h"
#include "trialdiv.h"
#include "getprime.h"
#include <sys/time.h>

int main (int argc, char **argv)
{
  fbprime_t *primes;
  trialdiv_divisor_t *d;
  const size_t max_div = 16;
  unsigned long *factors;
  int i, len = 1, nr_primes = 1000, nr_N = 100000;
  unsigned long expect = 0, nr_div = 0;
  mpz_t M, N, pk, t1, t2;
  double usrtime;
  int verbose = 0;
  
  if (argc > 1 && argv[1][0] == '-' && argv[1][1] == 'v')
    {
      verbose = 1;
      argc--;
      argv++;
    }

  if (argc > 1)
    len = atoi (argv[1]);
  
  if (argc > 2)
    nr_N = atoi (argv[2]);

  if (argc > 3)
    nr_primes = atoi (argv[3]);
  
  if (len > TRIALDIV_MAXLEN)
    {
      printf ("Error, trialdiv not implemented for input sizes greater than "
	      "%d words\n", TRIALDIV_MAXLEN);
      exit (EXIT_FAILURE);
    }

  mpz_init (M);
  mpz_init (N);
  mpz_init (pk);
  mpz_init (t1);
  mpz_init (t2);
  mpz_set_ui (N, 1UL);
  mpz_mul_2exp (N, N, 8 * sizeof(unsigned long) * len);
  mpz_tdiv_q_ui (N, N, 3UL);
  mpz_nextprime (N, N);

  for (i = 0; i < 0; i++) /* To skip some primes */
    getprime(1UL);
  primes = malloc (nr_primes * sizeof (fbprime_t));

  mpz_add_ui (M, N, nr_N - 1); /* We'll divide integers in [N, M] */
  mpz_sub_ui (N, N, 1UL);
  for (i = 0; i < nr_primes; i++)
    {
      primes[i] = getprime(1UL);
      /* Estimate the number of divisors we'll find when trial dividing by
         this p and its powers. 
         The number of times d divides in [0, n] is floor(n/d)+1,
         so the number of times it divides in [N, M] is
         floor(M/d) - floor((N-1)/d).
      */
      mpz_set_ui (pk, 1UL);
      do
	{
	  mpz_mul_ui (pk, pk, primes[i]);
	  mpz_tdiv_q (t1, N, pk);
	  mpz_tdiv_q (t2, M, pk);
	  mpz_sub (t2, t2, t1);
	  if (!mpz_fits_ulong_p (t2))
	    abort ();
	  expect += mpz_get_ui (t2);
	} while (mpz_sgn (t2) != 0);
    }
  mpz_add_ui (N, N, 1UL);
  if (verbose)
    gmp_printf ("Trial dividing integers in [%Zd, %Zd] "
                "by the %d primes in [%lu, %lu]\n", 
                N, M, nr_primes, (unsigned long) primes[0], 
                (unsigned long) primes[nr_primes - 1]);
  getprime (0);

  d = trialdiv_init (primes, nr_primes);
  free (primes);
  factors = (unsigned long *) malloc (max_div * sizeof (unsigned long));
  
  for (i = 0; i < nr_N; i++)
    {
      size_t t;
      mpz_set (M, N);
      do {
	t = trialdiv (factors, M, d, max_div);
	nr_div += (t > max_div) ? max_div : t;
      } while (t > max_div); /* If more factors were found than fit in factor 
                                array, continue factoring on cofactor */
      mpz_add_ui (N, N, 1UL);
    }
  
  usrtime = microseconds();

  mpz_clear (t1);
  mpz_clear (t2);
  mpz_clear (pk);
  mpz_clear (M);
  mpz_clear (N);
  trialdiv_clear (d);
  free (factors);
  if (verbose)
    {
      printf ("Found %lu divisors, expected %lu\n", nr_div, expect);
      printf ("Time: %f s, per N: %f mus, per prime: %f ns\n", 
              usrtime / 1e6, usrtime / nr_N, 
              usrtime / nr_N / nr_primes * 1e3);
    }
  if (nr_div != expect)
    {
      printf ("Error: did not find the expected number of divisors!\n");
      exit (EXIT_FAILURE);
    }
  return 0;
}
