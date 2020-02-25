/* test fb_root_in_qlattice_31bits */

#include "cado.h"
#include "macros.h"
#include "gmp.h"
#include "las-qlattice.hpp"

/* return (R*b1-a1)/(a0-R*b0) % p */
static fbprime_t
ref_fb_root_in_qlattice (fbprime_t p, fbprime_t R, qlattice_basis basis)
{
  mpz_t t, P;
  fbprime_t r;
  unsigned long tt;

  mpz_init_set_int64 (t, -basis.b0);
  mpz_mul_ui (t, t, R);
  mpz_add_int64 (t, t, basis.a0);
  mpz_init_set_ui (P, p);
  mpz_invert (t, t, P);
  tt = mpz_get_ui (t);
  mpz_set_int64 (t, basis.b1);
  mpz_mul_ui (t, t, R);
  mpz_sub_int64 (t, t, basis.a1);
  mpz_mul_ui (t, t, tt);
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
	  mpz_urandomb (t, rstate, 32);
	  mpz_nextprime (t, t);
	}
      while (mpz_sizeinbase (t, 2) != 32);
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
      /* R*b1-a1 and a0-R*b0 should fit in an int64_t, with R < p < 2^32,
	 and additionally |a0-R*b0|,|R*b1-a1| < 2^32*p (which holds here
	 since p has 32 bits, thus 2^32*p > 2^63);
	 we take 0 <= b0,b1 < 2^30, and |a0|, |a1| < 2^61 */
      mpz_urandomb (t, rstate, 63);
      basis[j].a0 = mpz_get_ui (t) - 0x4000000000000000UL;
      mpz_urandomb (t, rstate, 63);
      basis[j].a1 = mpz_get_ui (t) - 0x4000000000000000UL;
      mpz_urandomb (t, rstate, 30);
      basis[j].b0 = mpz_get_ui (t);
      mpz_urandomb (t, rstate, 30);
      basis[j].b1 = mpz_get_ui (t);
    }

  /* efficiency test */
  st = seconds ();
  r = 0;
  for (j = 0; j < N; j++)
    for (i = 0; i < N; i++)
      r += fb_root_in_qlattice_31bits (p[i], R[i], invp[i], basis[j]);
  st = seconds () - st;
  printf ("fb_root_in_qlattice_31bits: %lu tests took %.2fs (r=%u)\n",
	  N * N, st, r);

  /* correctness test */
  for (i = 0; i < Nsmall; i++)
    for (j = 0; j < Nsmall; j++)
      {
	r = fb_root_in_qlattice_31bits (p[i], R[i], invp[i], basis[j]);
	rref = ref_fb_root_in_qlattice (p[i], R[i], basis[j]);
	if (r != rref)
	  {
	    fprintf (stderr, "Error for p=%u R=%u a0=%ld b0=%ld a1=%ld b1=%ld\n",
		     p[i], R[i], basis[j].a0, basis[j].b0, basis[j].a1, basis[j].b1);
	    fprintf (stderr, "fb_root_in_qlattice_31bits gives %u\n", r);
	    fprintf (stderr, "ref_fb_root_in_qlattice gives %u\n", rref);
	    exit (1);
	  }
      }

  free (p);
  free (R);
  free (invp);
  free (basis);

  gmp_randclear (rstate);
  mpz_clear (t);
  mpz_clear (u);
}

static void
test_fb_root_in_qlattice_127bits (unsigned long N)
{
  fbprime_t *p, *R, r, r31, rref;
  uint32_t *invp32;
  uint64_t *invp64;
  gmp_randstate_t rstate;
  unsigned long i, j;
  mpz_t t, u;
  qlattice_basis *basis;
  double st;
  unsigned long Nsmall = N; /* use a smaller set for correctness test */

  gmp_randinit_default (rstate);
  mpz_init (t);
  mpz_init (u);

  p = (fbprime_t*) malloc (N * sizeof (fbprime_t));
  ASSERT_ALWAYS(p != NULL);
  R = (fbprime_t*) malloc (N * sizeof (fbprime_t));
  ASSERT_ALWAYS(R != NULL);
  invp32 = (uint32_t*) malloc (N * sizeof (uint32_t));
  ASSERT_ALWAYS(invp32 != NULL);
  invp64 = (uint64_t*) malloc (N * sizeof (uint64_t));
  ASSERT_ALWAYS(invp64 != NULL);
  basis = (qlattice_basis*) malloc (N * sizeof (qlattice_basis));

  /* generate p[i], R[i], invp32[i], invp64[i] for 0 <= i < N */
  for (i = 0; i < N; i++)
    {
      do
	{
	  mpz_urandomb (t, rstate, 32);
	  mpz_nextprime (t, t);
	}
      while (mpz_sizeinbase (t, 2) != 32);
      p[i] = mpz_get_ui (t);
      mpz_urandomm (u, rstate, t);
      R[i] = mpz_get_ui (u);
      ASSERT_ALWAYS (R[i] < p[i]);
      /* invp32[i] is -1/p[i] mod 2^32 */
      mpz_ui_pow_ui (u, 2, 32);
      mpz_neg (t, t);
      mpz_invert (t, t, u);
      invp32[i] = mpz_get_ui (t);
      /* invp64[i] is -1/p[i] mod 2^64 */
      mpz_ui_pow_ui (u, 2, 64);
      mpz_set_ui (t, p[i]);
      mpz_neg (t, t);
      mpz_invert (t, t, u);
      invp64[i] = mpz_get_ui (t);
    }

  /* generate basis[j] for 0 <= j < N */
  for (j = 0; j < N; j++)
    {
      mpz_urandomb (t, rstate, 63);
      basis[j].a0 = mpz_get_ui (t) - 0x4000000000000000UL;
      mpz_urandomb (t, rstate, 63);
      basis[j].a1 = mpz_get_ui (t) - 0x4000000000000000UL;
      mpz_urandomb (t, rstate, 30);
      basis[j].b0 = mpz_get_ui (t);
      mpz_urandomb (t, rstate, 30);
      basis[j].b1 = mpz_get_ui (t);
    }

  /* efficiency test */
  st = seconds ();
  r = 0;
  for (j = 0; j < N; j++)
    for (i = 0; i < N; i++)
      r += fb_root_in_qlattice_127bits (p[i], R[i], invp64[i], basis[j]);
  st = seconds () - st;
  printf ("fb_root_in_qlattice_127bits: %lu tests took %.2fs (r=%u)\n",
          N * N, st, r);

  /* correctness test */
  for (i = 0; i < Nsmall; i++)
    for (j = 0; j < Nsmall; j++)
      {
	r = fb_root_in_qlattice_127bits (p[i], R[i], invp64[i], basis[j]);
	r31 = fb_root_in_qlattice_31bits (p[i], R[i], invp32[i], basis[j]);
	if (r != r31)
	  {
	    fprintf (stderr, "Error for p=%u R=%u a0=%ld b0=%ld a1=%ld b1=%ld\n",
		     p[i], R[i], basis[j].a0, basis[j].b0, basis[j].a1, basis[j].b1);
	    fprintf (stderr, "fb_root_in_qlattice_127bits gives %u\n", r);
	    fprintf (stderr, "fb_root_in_qlattice_31bits  gives %u\n", r31);
	    exit (1);
	  }
	rref = ref_fb_root_in_qlattice (p[i], R[i], basis[j]);
	if (r != rref)
	  {
	    fprintf (stderr, "Error for p=%u R=%u a0=%ld b0=%ld a1=%ld b1=%ld\n",
		     p[i], R[i], basis[j].a0, basis[j].b0, basis[j].a1, basis[j].b1);
	    fprintf (stderr, "fb_root_in_qlattice_127bits gives %u\n", r);
	    fprintf (stderr, "ref_fb_root_in_qlattice gives %u\n", rref);
	    exit (1);
	  }
      }

  free (p);
  free (R);
  free (invp32);
  free (invp64);
  free (basis);

  gmp_randclear (rstate);
  mpz_clear (t);
  mpz_clear (u);
}

/* exercise bugs in invmod_redc_64 and fb_root_in_qlattice_127bits */
static void
bug20200225 (void)
{
  uint64_t a = 811915144;
  uint64_t b = 2190486847;
  uint64_t expected = 2164754518;
  uint64_t got = invmod_redc_64 (a, b);
  if (got != expected)
    {
      fprintf (stderr, "Error in invmod_redc_64 for a=%lu b=%lu\n", a, b);
      fprintf (stderr, "Expected %lu\n", expected);
      fprintf (stderr, "Got      %lu\n", got);
      exit (1);
    }

  fbprime_t p = 3628762957;
  fbprime_t R = 1702941053;
  uint64_t invp = 5839589727713490555UL;
  qlattice_basis basis[1];
  basis[0].a0 = -2503835703516628395L;
  basis[0].b0 = 238650852;
  basis[0].a1 = -3992552824749287692L;
  basis[0].b1 = 766395543;
  expected = 987485779;
  got = fb_root_in_qlattice_127bits (p, R, invp, basis[0]);
  if (got != expected)
    {
      fprintf (stderr, "Error in fb_root_in_qlattice_127bits for p=%u R=%u "
	       "a0=%ld b0=%ld a1=%ld b1=%ld\n",
	       p, R, basis[0].a0, basis[0].b0, basis[0].a1, basis[0].b1);
      fprintf (stderr, "Expected %lu\n", expected);
      fprintf (stderr, "Got      %lu\n", got);
      exit (1);
    }

  p = 3725310689;
  R = 2661839516;
  invp = 1066179678986106591UL;
  basis[0].a0 = 3008222006914909739L;
  basis[0].b0 = 877054135;
  basis[0].a1 = 3170231873717741170L;
  basis[0].b1 = 932375769;
  expected = 2956728450;
  got = fb_root_in_qlattice_127bits (p, R, invp, basis[0]);
  if (got != expected)
    {
      fprintf (stderr, "Error in fb_root_in_qlattice_127bits for p=%u R=%u "
	       "a0=%ld b0=%ld a1=%ld b1=%ld\n",
	       p, R, basis[0].a0, basis[0].b0, basis[0].a1, basis[0].b1);
      fprintf (stderr, "Expected %lu\n", expected);
      fprintf (stderr, "Got      %lu\n", got);
      exit (1);
    }
}

int
main (int argc, char *argv[])
{
  unsigned long N;

  ASSERT_ALWAYS (argc == 1 || argc == 2);

  N = (argc == 1) ? 1000 : strtoul (argv[1], NULL, 10);

  bug20200225 ();

  test_fb_root_in_qlattice_31bits (N);
  test_fb_root_in_qlattice_127bits (N);

  return 0;
}



