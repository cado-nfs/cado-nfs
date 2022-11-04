#include "cado.h" // IWYU pragma: keep
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <gmp.h>               // for mpz_sgn, mpz_add_ui, mpz_set_ui, mpz_s...
#include <stdint.h>            // for uint64_t
#include <vector>              // for vector
#include "gmp_aux.h"           // for ulong_nextprime
#include "sieve/trialdiv.hpp"
#include "tests_common.h"
#include "portability.h" // IWYU pragma: keep
#include "cxx_mpz.hpp"  // cxx_mpz
#include "getprime.h"   // prime_info
#include "timing.h"     // microseconds
#include "macros.h"

void
trialdiv_stdinput(const unsigned long pmax, const int verbose)
{
    std::vector<unsigned long> primes;
  unsigned long p;
  cxx_mpz N;
  prime_info pi;

  prime_info_init (pi);
  for (p = getprime_mt (pi); p <= pmax; p = getprime_mt (pi))
    primes.push_back(p);
  prime_info_clear (pi);

  trialdiv_data d(primes);

  
  while (!feof(stdin)) {
    size_t i;
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

    std::vector<uint64_t> factors = d.trial_divide(N);

    if (verbose) {
      for (i = 0; i < bit; i++) {
        printf ("2 ");
      }
      for (auto p: factors)
        printf ("%" PRIu64 " ", p);
    }
    gmp_printf ("%Zd\n", (mpz_srcptr) N);
  }
}

/* performs iter random tests with a cofactor of n limbs */
void
test_trialdiv (int n, unsigned long iter)
{
  cxx_mpz N;
  unsigned long p, pmax;
  std::vector<uint64_t> g;
  int ret;

  pmax = trialdiv_data::max_p;
  for (unsigned long i = 0; i <iter; i++)
    {
      if (i == 0)
        p = 3;
      else if (i == 1) {
        /* Find largest prime <= pmax */
        unsigned long r;
        for (r = pmax, p = 0; p == 0 || p > pmax; r -= 2)
          p = ulong_nextprime (r);
      } else {
          do p = ulong_nextprime (gmp_urandomm_ui(state, pmax)); while (p > pmax || p < 3);
      }
      /* Trial divide candidate numbers with a trivial list of exactly 1
       * prime (p). This is really a dumb test */
      trialdiv_data d(std::vector<unsigned long>(1, p));

      mpz_urandomb (N, state, n * mp_bits_per_limb);
      ret = mpz_divisible_ui_p (N, p);
      g = d.trial_divide (N, 1);
      ASSERT_ALWAYS (g.size() <= 1); 
      if (ret) {
        ASSERT_ALWAYS (g.size() >= 1);
        ASSERT_ALWAYS (g[0] == p);
      } else {
        ASSERT_ALWAYS (g.size() == 0);
      }

      /* now test a case where it should divide */
      mpz_sub_ui (N, N, mpz_fdiv_ui (N, p));
      if (mpz_sgn(N) == 0)
        mpz_set_ui(N, p);
      g = d.trial_divide(N, 1);

      ASSERT_ALWAYS (g.size() == 1);
      ASSERT_ALWAYS (g[0] == p);
    }
}

// coverity[root_function]
int main (int argc, const char **argv)
{
  int i, len = 1, nr_primes = 1000, nr_N = 100000;
  unsigned long expect = 0, nr_div = 0;
  cxx_mpz M, N, pk, t1, t2;
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

  if (input) {
    if (argc <= 1) {
        fprintf(stderr, "-i requires a value for pmax on the command line\n");
        exit (EXIT_FAILURE);
    }
    trialdiv_stdinput (atol(argv[1]), verbose);
    exit (EXIT_SUCCESS);
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

  mpz_set_ui (N, 1UL);
  mpz_mul_2exp (N, N, 8 * sizeof(unsigned long) * len);
  mpz_tdiv_q_ui (N, N, 3UL);
  mpz_nextprime (N, N);

  std::vector<unsigned long> primes;
  prime_info pi;
  prime_info_init (pi);
  for (i = 0; i < 0; i++) /* To skip some primes */
    getprime_mt (pi);
  primes.reserve(nr_primes);

  mpz_add_ui (M, N, nr_N - 1); /* We'll divide integers in [N, M] */
  mpz_sub_ui (N, N, 1UL);
  for (i = 0; i < nr_primes; i++)
    {
      primes.push_back(getprime_mt (pi));
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
                (mpz_srcptr) N,
                (mpz_srcptr) M,
                nr_primes, (unsigned long) primes[0], 
                (unsigned long) primes[nr_primes - 1]);
  prime_info_clear (pi);

  trialdiv_data d(primes);

  for (i = 0; i < nr_N; i++)
    {
      mpz_set (M, N);
      nr_div += d.trial_divide(M).size();
      mpz_add_ui (N, N, 1UL);
    }
  
  usrtime = microseconds();

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
