/* test fb_root_in_qlattice_31bits */

#include "cado.h"
#include "macros.h"
#include "gmp.h"
#include "las-qlattice.hpp"

/* return (R*b1-a1)/(a0-R*b0) % p */
static fbprime_t
ref_fb_root_in_qlattice_31bits (fbprime_t p, fbprime_t R, qlattice_basis basis)
{
  mpz_t t, P;
  fbprime_t r;
  int64_t num, den;

  num = (int64_t) R * basis.b1 - basis.a1;
  den = basis.a0 - (int64_t) R * basis.b0;
  mpz_init_set_int64 (t, den);
  mpz_init_set_ui (P, p);
  mpz_invert (t, t, P);
  mpz_mul_int64 (t, t, num);
  mpz_mod (t, t, P);
  r = mpz_get_ui (t);
  mpz_clear (t);
  mpz_clear (P);
  return r;
}

static void
test_fb_root_in_qlattice_31bits (unsigned long N)
{
  fbprime_t *p, *R, r, rref;
  uint32_t *invp;
  gmp_randstate_t rstate;
  unsigned long i, j;
  mpz_t t, u;
  qlattice_basis *basis;
  double st;
  unsigned long Nsmall = N / 10; /* use a smaller set for correctness test */

  gmp_randinit_default (rstate);
  mpz_init (t);
  mpz_init (u);

  p = (fbprime_t*) malloc (N * sizeof (fbprime_t));
  ASSERT_ALWAYS(p != NULL);
  R = (fbprime_t*) malloc (N * sizeof (fbprime_t));
  ASSERT_ALWAYS(R != NULL);
  invp = (uint32_t*) malloc (N * sizeof (uint32_t));
  ASSERT_ALWAYS(invp != NULL);
  basis = (qlattice_basis*) malloc (N * sizeof (qlattice_basis));

  /* generate p[i], R[i], invp[i] for 0 <= i < N */
  for (i = 0; i < N; i++)
    {
      do
	{
	  mpz_urandomb (t, rstate, 31);
	  mpz_nextprime (t, t);
	}
      while (mpz_sizeinbase (t, 2) != 31);
      p[i] = mpz_get_ui (t);
      mpz_urandomm (u, rstate, t);
      R[i] = mpz_get_ui (u);
      ASSERT_ALWAYS (R[i] < p[i]);
      /* invp[i] is -1/p[i] mod 2^32 */
      mpz_ui_pow_ui (u, 2, 32);
      mpz_neg (t, t);
      mpz_invert (t, t, u);
      invp[i] = mpz_get_ui (t);
    }

  /* generate basis[j] for 0 <= j < N */
  for (j = 0; j < N; j++)
    {
      /* R*b1-a1 and a0-R*b0 should fit in an int64_t, with R < 2^31,
	 and additionally |a0-R*b0|,|R*b1-a1| < 2^32*p:
	 we take 0 <= b0,b1 < 2^31, and |a0|, |a1| < 2^61 */
      mpz_urandomb (t, rstate, 62);
      basis[j].a0 = mpz_get_ui (t) - 0x2000000000000000UL;
      mpz_urandomb (t, rstate, 62);
      basis[j].a1 = mpz_get_ui (t) - 0x2000000000000000UL;
      mpz_urandomb (t, rstate, 31);
      basis[j].b0 = mpz_get_ui (t);
      mpz_urandomb (t, rstate, 31);
      basis[j].b1 = mpz_get_ui (t);
    }

  /* correctness test */
  for (i = 0; i < Nsmall; i++)
    for (j = 0; j < Nsmall; j++)
      {
	r = fb_root_in_qlattice_31bits (p[i], R[i], invp[i], basis[j]);
	rref = ref_fb_root_in_qlattice_31bits (p[i], R[i], basis[j]);
	if (r != rref)
	  {
	    fprintf (stderr, "Error for p=%u R=%u a0=%ld b0=%ld a1=%ld b1=%ld\n",
		     p[i], R[i], basis[j].a0, basis[j].b0, basis[j].a1, basis[j].b1);
	    fprintf (stderr, "fb_root_in_qlattice_31bits gives %u\n", r);
	    fprintf (stderr, "ref_fb_root_in_qlattice_31bits gives %u\n", rref);
	    exit (1);
	  }
      }

  /* efficiency test */
  st = seconds ();
  r = 0;
  for (j = 0; j < N; j++)
    for (i = 0; i < N; i++)
      r += fb_root_in_qlattice_31bits (p[i], R[i], invp[i], basis[j]);
  st = seconds () - st;
  printf ("%lu tests took %.2fs (r=%u)\n", N * N, st, r);

  free (p);
  free (R);
  free (invp);
  free (basis);

  gmp_randclear (rstate);
  mpz_clear (t);
  mpz_clear (u);
}

int
main (int argc, char *argv[])
{
  unsigned long N;

  ASSERT_ALWAYS (argc == 1 || argc == 2);

  N = (argc == 1) ? 1000 : strtoul (argv[1], NULL, 10);

  test_fb_root_in_qlattice_31bits (N);

  return 0;
}



