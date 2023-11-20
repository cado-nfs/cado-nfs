#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "rootfinder.h"
#include "tests_common.h"
#include "mpz_poly.h"
#include "macros.h"
#include "gmp_aux.h"
#include "portability.h" // IWYU pragma: keep

int cmp(mpz_t * a, mpz_t * b)
{
    return mpz_cmp(*a, *b);
}

#define MAX_BITS 40

/* Check roots of polynomial ff[d]*x^d + ... + ff[1]*x + ff[0] mod pp.
   If pp is the empty string, generate a random integer (and then generate
   random coefficients for ff[]).
   If nroots <> -1, it is the expected number of roots. */
void
test (int d, const char *pp, const char *ff[], int nroots)
{
  mpz_t p, *f, *r, v;
  int i, n, n1;
  mpz_poly F;

  mpz_init (p);
  if (strlen (pp) > 0)
    mpz_set_str (p, pp, 0);
  else
    {
      do mpz_urandomb (p, state, MAX_BITS); while (mpz_even_p (p));
    }
  mpz_init (v);
  f = (mpz_t *) malloc ((d + 1) * sizeof(mpz_t));
  F->coeff = f;
  F->deg = d;
  for (i = 0; i <= d; i++)
    {
      mpz_init (f[i]);
    }
 retry:
  mpz_set (v, p);
  for (i = 0; i <= d; i++)
    {
      if (strlen (pp) > 0)
        mpz_set_str (f[i], ff[d - i], 0);
      else
        mpz_urandomb (f[i], state, MAX_BITS);
      mpz_gcd (v, v, f[i]);
    }
  if (mpz_cmp_ui (v, 1) > 0) /* f vanishes modulo some prime factor of p */
    {
      if (strlen (pp) > 0)
        {
          gmp_fprintf (stderr, "Error, given polynomial vanishes mod %Zd\n",
                       v);
          exit (1);
        }
      else
        goto retry;
    }
  n = mpz_poly_roots_gen (&r, F, p, state);

  if (tests_common_get_verbose()) {
      for(int i = 0 ; i < n ; i++) {
          gmp_printf("%Zd\n", r[i]);
      }
  }

  if (mpz_probab_prime_p (p, 5) && mpz_sizeinbase (p, 2) <= 64)
    {
      n1 = mpz_poly_roots_uint64 (NULL, F, mpz_get_uint64 (p), state);
      ASSERT_ALWAYS(n1 == n);
    }
  if (nroots != -1)
    ASSERT_ALWAYS(n == nroots);
  for (i = 0; i < n; i++)
    {
      mpz_poly_eval (v, F, r[i]);
      ASSERT_ALWAYS(mpz_divisible_p (v, p));
    }
  for (i = 0; i <= d; i++)
    mpz_clear (f[i]);
  for (i = 0; i < n; i++)
    mpz_clear (r[i]);
  free (f);
  free (r);
  mpz_clear (p);
  mpz_clear (v);
}

int
main (int argc, const char *argv[])
{
    int d;
    const char* test0[] = {"4294967291", "1", "0", "-3"};
    const char* test1[] = {"18446744073709551557", "1", "2", "3", "5"};
    const char* test2[] = {"18446744073709551629", "1", "-1", "7", "-1"};
    const char* test3[] = {"246089", "1", "2", "3", "4"};
    const char* test4[] = {"9", "1", "0", "0", "1"};
    unsigned long iter = 100;

    tests_common_cmdline (&argc, &argv, PARSE_SEED | PARSE_ITER | PARSE_VERBOSE);
    tests_common_get_iter (&iter);

    d = argc - 3;

    if (d < 1 && d != -2) {
	fprintf(stderr, "Usage: test_rootfinder [-seed nnn] [-iter nnn] p a_d [...] a_1 a_0\n");
	exit(1);
    }

    if (d >= 1) {
      test (d, argv[1], argv + 2, -1);
    } else {
        test (2, test0[0], test0 + 1, 2);
        test (3, test1[0], test1 + 1, 1);
        test (3, test2[0], test2 + 1, 0);
        test (3, test3[0], test3 + 1, 27);
        test (3, test4[0], test4 + 1, 3);

        while (iter--)
          {
            const char * foo[] = {
                "10315", "28842", "29254", "36865",
                "48751", "56362", "63592", "71203" };
            d = 1 + gmp_urandomm_ui(state, 7);
            test (d, "", foo, -1);
          }
    }

    tests_common_clear ();
    exit (EXIT_SUCCESS);
}

#if 0
// magma code for producing test cases.
s:=1.2;            
p:=10;
while p lt 2^200 do
    for i in [1..100] do
        p:=NextPrime(p);
        d:=Random([2..7]);
        coeffs:=[Random(GF(p)):i in [0..d]];
        F:=PolynomialRing(GF(p))!coeffs;
        printf "in %o", p;
        for c in Reverse(coeffs) do printf " %o", c; end for;
        printf "\n";
        r:=Sort([x[1]: x in Roots(F) ]);
        printf "out";
        for c in r do printf " %o", c; end for;
        printf "\n";
    end for;
    p := Ceiling(p*s);
end while;


#endif

