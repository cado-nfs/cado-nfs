#include "cado.h" // IWYU pragma: keep
#include "macros.h"
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "mod_ul.h"        // for modul_initmod_ul, modul_clearmod, modul_se...
#include "mpz_poly.h"      // for mpz_poly
#include "modul_poly.h"
#include "gmp_aux.h"
#include "tests_common.h"
#include "cado_poly.h"
#include "portability.h" // IWYU pragma: keep

void
test_modul_poly_is_irreducible (unsigned long iter)
{
  modul_poly_t f;
  modulusul_t p;
  int d, i, irred, n;
  unsigned long q;
  residueul_t r[MAX_DEGREE];

  /* first try some hard-coded polynomials */
  modul_poly_init (f, MAX_DEGREE);

  modul_initmod_ul (p, 3);
  modul_poly_set_immediate (f, 0, p, 1);
  ASSERT_ALWAYS (modul_poly_is_irreducible (f, p) != 0);
  ASSERT_ALWAYS (modul_poly_is_squarefree (f, p) == 0);

  modul_initmod_ul (p, 3);
  modul_poly_set_immediate (f, -1, p);
  ASSERT_ALWAYS (modul_poly_is_irreducible (f, p) != 0);
  ASSERT_ALWAYS (modul_poly_is_squarefree (f, p) == 0);

  modul_initmod_ul (p, 3);
  modul_poly_set_immediate (f, 2, p, 2, 1, 2);
  ASSERT_ALWAYS (modul_poly_is_irreducible (f, p) == 0);
  ASSERT_ALWAYS (modul_poly_is_squarefree (f, p) == 0);

  modul_initmod_ul (p, 5);
  modul_poly_set_immediate (f, 2, p, 2, 1, 2);
  ASSERT_ALWAYS (modul_poly_is_irreducible (f, p) == 0);

  modul_initmod_ul (p, 7);
  modul_poly_set_immediate (f, 2, p, 2, 1, 2);
  ASSERT_ALWAYS (modul_poly_is_irreducible (f, p) != 0);
  ASSERT_ALWAYS (modul_poly_is_squarefree (f, p) != 0);

  modul_initmod_ul (p, 11);
  modul_poly_set_immediate (f, 2, p, 2, 1, 2);
  ASSERT_ALWAYS (modul_poly_is_irreducible (f, p) != 0);
  ASSERT_ALWAYS (modul_poly_is_squarefree (f, p) != 0);

  /* Has distinct factors: (x + 3) * (x + 6) */
  modul_initmod_ul (p, 17);
  modul_poly_set_immediate (f, 2, p, 2, 1, 2);
  ASSERT_ALWAYS (modul_poly_is_irreducible (f, p) == 0);
  ASSERT_ALWAYS (modul_poly_is_squarefree (f, p) != 0);

  /* 20*x^2 + 18*x + 19 == (x + 20)^2 (mod 23) */
  modul_initmod_ul (p, 23);
  modul_poly_set_immediate (f, 2, p, 19, 18, 20);
  ASSERT_ALWAYS (modul_poly_is_irreducible (f, p) == 0);
  ASSERT_ALWAYS (modul_poly_is_squarefree (f, p) == 0);

  while (iter--)
    {
      d = 1 + gmp_urandomm_ui(state, MAX_DEGREE - 1);
      q = gmp_urandomb_ui(state, 31);
      q = ulong_nextprime (q);
      /* modul_poly_cantor_zassenhaus only works for odd primes */
      q += (q == 2);
      modul_initmod_ul (p, q);
      for (i = 0; i <= d; i++)
        modul_set_ul (f->coeff[i], gmp_urandomb_ui(state, 31), p);
      while (modul_is0 (f->coeff[d], p))
        modul_set_ul (f->coeff[d], gmp_urandomb_ui(state, 31), p);
      f->degree = d;
      irred = modul_poly_is_irreducible (f, p);
      int squarefree = modul_poly_is_squarefree (f, p);
      if (d == 1) {
        ASSERT_ALWAYS(irred != 0);
        ASSERT_ALWAYS(squarefree);
      }
      if (irred) {
        ASSERT_ALWAYS(squarefree);
      }
      if (d == 1 || (d == 2 && irred == 0 && squarefree))
        {
          n = modul_poly_cantor_zassenhaus (r, f, p, state);
          ASSERT_ALWAYS(n == d);
        }
      modul_clearmod (p);
    }
  modul_poly_clear (f);
}

void
test_modul_poly_roots_ulong (unsigned long iter)
{
  unsigned long r[MAX_DEGREE];
  mpz_t f[MAX_DEGREE + 1];
  int d, i, n;
  modulusul_t p;
  residueul_t y, x;
  modul_poly_t fp;

  for (i = 0; i <= MAX_DEGREE; i++)
    mpz_init (f[i]);

  /* hard-coded examples */
  mpz_set_ui (f[0], 0);
  mpz_set_ui (f[1], 0);
  mpz_set_ui (f[2], 1);
  d = 2;
  mpz_poly F;
  F->coeff =  f;
  F->deg = d;
  modul_initmod_ul (p, 113);
  n = modul_poly_roots_ulong (r, F, p, state);
  ASSERT_ALWAYS(n == 1 && r[0] == 0);
  modul_clearmod (p);

  modul_initmod_ul (p, 13);
  mpz_set_ui (f[0], 4);
  mpz_set_ui (f[1], 10);
  mpz_set_ui (f[2], 12);
  mpz_set_ui (f[3], 0);
  mpz_set_ui (f[4], 9);
  mpz_set_ui (f[5], 3);
  mpz_set_ui (f[6], 1);
  F->deg = 6;
  n = modul_poly_roots (NULL, F, p, state);
  ASSERT_ALWAYS (n == 5);

  while (iter--)
    {
      d = 1 + gmp_urandomm_ui(state, MAX_DEGREE - 1);
      F->deg = d;
      for (i = 0; i <= d; i++)
        mpz_urandomb (f[i], state, 64);
      modul_initmod_ul (p, ulong_nextprime (gmp_urandomb_ui(state, 31)));
      while (mpz_divisible_ui_p (f[d], modul_getmod_ul (p)))
        mpz_urandomb (f[d], state, 64);
      n = modul_poly_roots_ulong (r, F, p, state);
      ASSERT_ALWAYS(0 <= n && n <= d);
      modul_poly_init (fp, d);
      modul_poly_set_mod (fp, F, p);
      if (n > 0 && fp->degree > 1)
        ASSERT_ALWAYS(modul_poly_is_irreducible (fp, p) == 0);
      /* if n=0, f might be irreducible or not mod p (product of two
         degree-2 factors for example),
         if d=1, f is irreducible */
      if (d == 1 || (d <= 3 && n == 0))
        ASSERT_ALWAYS(modul_poly_is_irreducible(fp, p) != 0);
      for (i = 0; i < n; i++)
        {
          modul_set_ul (x, r[i], p);
          modul_poly_eval (y, fp, x, p);
          ASSERT_ALWAYS(modul_is0 (y, p));
        }
      modul_poly_clear (fp);
      modul_clearmod (p);
    }
  for (i = 0; i <= MAX_DEGREE; i++)
    mpz_clear (f[i]);
}

int
main (int argc, const char *argv[])
{
  unsigned long iter = 1000;

  tests_common_cmdline (&argc, &argv, PARSE_SEED | PARSE_ITER);
  tests_common_get_iter (&iter);
  test_modul_poly_is_irreducible (iter);
  test_modul_poly_roots_ulong (iter);
  tests_common_clear ();
  exit (EXIT_SUCCESS);
}

