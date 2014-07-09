#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <sys/time.h>
#include "sieve/trialdiv.h"
#include "utils.h"
#include "test_iter.h"
#include "tests_common.h"

void
trialdiv_stdinput(const unsigned long pmax, const int verbose)
{
  fbprime_t *primes;
  trialdiv_divisor_t *d;
  unsigned long *factors;
  int nr_primes = 0;
  unsigned long p;
  const size_t max_div = 32;
  mpz_t N;

  for (p = getprime(1UL); p <= pmax; p = getprime(1UL))
    nr_primes++;
  primes = malloc (nr_primes * sizeof (fbprime_t));

  getprime(0UL);
  nr_primes = 0;
  for (p = getprime(1UL); p <= pmax; p = getprime(1UL))
    primes[nr_primes++] = p;

  d = trialdiv_init (primes, nr_primes);
  free (primes);
  factors = (unsigned long *) malloc (max_div * sizeof (unsigned long));
  
  mpz_init(N);
  while (!feof(stdin)) {
    size_t t, i;
    unsigned long bit; /* old versions of GMP do not have mp_bitcnt_t */
    
    if (mpz_inp_str (N, stdin, 10) == 0)
      break;
    
    if (mpz_sgn(N) == 0)
      break;
    
    if (mpz_sgn(N) < 0)
      mpz_neg (N, N);

    /* Now N is positive */
    bit = mpz_scan1(N, 0);
    if (bit > 0)
      mpz_tdiv_q_2exp (N, N, bit);

    /* Now N is positive and odd */
    if (mpz_cmp_ui(N, 1UL) == 0) {
      printf ("1\n");
      continue;
    }
    t = trialdiv (factors, N, d, max_div);
    ASSERT_ALWAYS (t <= max_div);

    if (verbose) {
      for (i = 0; i < bit; i++) {
        printf ("2 ");
      }
      for (i = 0; i < t; i++)
        printf ("%lu ", (unsigned long) factors[i]);
    }
    gmp_printf ("%Zd\n", N);
  }
  
  mpz_clear (N);
  trialdiv_clear (d);
}

/* performs iter random tests with a cofactor of n limbs */
void
test_trialdiv (int n, unsigned long iter)
{
  mpz_t N;
  trialdiv_divisor_t *d;
  fbprime_t f[1];
  unsigned long p, g[1], i, pmax;
  size_t s;
  int ret;

  mpz_init (N);
  pmax = trialdiv_get_max_p();
  for (i = 0; i < iter; i++)
    {
      if (i == 0)
        p = 3;
      else if (i == 1) {
        /* Find largest prime <= pmax */
        unsigned long r;
        for (r = pmax, p = 0; p == 0 || p > pmax; r -= 2)
          p = ulong_nextprime (r);
      } else {
          do p = ulong_nextprime (lrand48 () % pmax); while (p > pmax || p < 3);
      }
      f[0] = p;
      d = trialdiv_init (f, 1);

      mpz_urandomb (N, state, n * mp_bits_per_limb);
      ret = mpz_divisible_ui_p (N, p);
      s = trialdiv (g, N, d, 1);
      ASSERT_ALWAYS (s <= 2); /* s can be max_div+1, i.e., 2 */
      if (ret) {
        ASSERT_ALWAYS (s >= 1);
        ASSERT_ALWAYS (g[0] == p);
      } else
        ASSERT_ALWAYS (s == 0);

      /* now test a case where it should divide */
      mpz_sub_ui (N, N, mpz_fdiv_ui (N, p));
      if (mpz_sgn(N) == 0)
        mpz_set_ui(N, p);
      s = trialdiv (g, N, d, 1);
      ASSERT_ALWAYS (1 <= s && s <= 2); /* s can be max_div+1, i.e., 2 */
      ASSERT_ALWAYS (g[0] == p);
      trialdiv_clear(d);
    }
  mpz_clear (N);
}

int main (int argc, const char **argv)
{
  fbprime_t *primes;
  trialdiv_divisor_t *d;
  const size_t max_div = 16;
  unsigned long *factors;
  int i, len = 1, nr_primes = 1000, nr_N = 100000;
  unsigned long expect = 0, nr_div = 0;
  mpz_t M, N, pk, t1, t2;
  double usrtime;
  int verbose = 0, input = 0;
  unsigned long iter = 10000;

  tests_common_cmdline(&argc, &argv, PARSE_SEED | PARSE_ITER);
  tests_common_get_iter(&iter);

  while (1) {
    if (argc > 1 && argv[1][0] == '-' && argv[1][1] == 'v')
      {
        verbose = 1;
        argc--;
        argv++;
      }
    else if (argc > 1 && argv[1][0] == '-' && argv[1][1] == 'i')
      {
        input = 1;
        argc--;
        argv++;
      }
    else
      break;
  }

  if (argc > 1)
    len = atoi (argv[1]);
  
  if (argc > 2)
    nr_N = atoi (argv[2]);

  if (argc > 3)
    nr_primes = atoi (argv[3]);
  
  if (input) {
    /* First parameter is pmax, and is stored in len */
    trialdiv_stdinput (len, verbose);
    exit (EXIT_SUCCESS);
  }

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
	  ASSERT_ALWAYS (mpz_fits_ulong_p (t2));
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

  /* random tests */
  for (int i = 1; i <= TRIALDIV_MAXLEN; i++)
    test_trialdiv (i, iter);

  tests_common_clear();
  exit (EXIT_SUCCESS);
}
