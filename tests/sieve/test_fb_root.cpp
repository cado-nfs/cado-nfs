/* test fb_root_in_qlattice_31bits */

#include "cado.h" // IWYU pragma: keep
#include <cinttypes>               // for PRId64, PRIu64
#include <cstdint>                 // for uint32_t, uint64_t
#include <cstring>                 // for strcmp
#include <cstdio>
#include <cstdlib>
#include <memory>                   // for allocator_traits<>::value_type
#include <vector>
#include <gmp.h>
#include "fb-types.h"               // for fbprime_t, FBPRIME_FORMAT
#include "gmp_aux.h"                // for mpz_add_int64, mpz_init_set_int64
#include "las-arith.hpp"            // for invmod_redc_32
#include "macros.h"
#include "las-qlattice.hpp"
#include "las-fbroot-qlattice.hpp"
#include "timing.h"  // seconds

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
test_fb_root_in_qlattice_31bits (gmp_randstate_t rstate, int test_speed, unsigned long N)
{
  fbprime_t *p, *R, r, rref;
  uint32_t *invp;
  unsigned long i, j;
  mpz_t t, u;
  std::vector<qlattice_basis> basis;
  double st;

  mpz_init (t);
  mpz_init (u);

  p = (fbprime_t*) malloc (N * sizeof (fbprime_t));
  ASSERT_ALWAYS(p != NULL);
  R = (fbprime_t*) malloc (N * sizeof (fbprime_t));
  ASSERT_ALWAYS(R != NULL);
  invp = (uint32_t*) malloc (N * sizeof (uint32_t));
  ASSERT_ALWAYS(invp != NULL);
  basis.assign(N, qlattice_basis());

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

  if (test_speed) {
      /* efficiency test */
      st = seconds ();
      r = 0;
      for (j = 0; j < N; j++)
          for (i = 0; i < N; i++)
              r += fb_root_in_qlattice_31bits (p[i], R[i], invp[i], basis[j]);
      st = seconds () - st;
      printf ("fb_root_in_qlattice_31bits: %lu tests took %.2fs (r=%" FBPRIME_FORMAT ")\n",
              N * N, st, r);
  } else {
      /* correctness test */
      for (i = 0; i < N; i++)
          for (j = 0; j < N; j++)
          {
              r = fb_root_in_qlattice_31bits (p[i], R[i], invp[i], basis[j]);
              rref = ref_fb_root_in_qlattice (p[i], R[i], basis[j]);
              if (r != rref)
              {
                  fprintf (stderr, "Error for p:=%" FBPRIME_FORMAT "; R:=%" FBPRIME_FORMAT "; a0:=%" PRId64 "; b0:=%" PRId64 "; a1:=%" PRId64 "; b1:=%" PRId64 ";\n",
                          p[i], R[i], basis[j].a0, basis[j].b0, basis[j].a1, basis[j].b1);
                  fprintf (stderr, "fb_root_in_qlattice_31bits gives %" FBPRIME_FORMAT "\n", r);
                  fprintf (stderr, "ref_fb_root_in_qlattice gives %" FBPRIME_FORMAT "\n", rref);
                  exit (1);
              }
          }
  }

  free (p);
  free (R);
  free (invp);

  mpz_clear (t);
  mpz_clear (u);
}

static void
test_fb_root_in_qlattice_127bits (gmp_randstate_t rstate, int test_speed, unsigned long N)
{
  fbprime_t *p, *R, r, r31, rref;
  uint32_t *invp32;
  uint64_t *invp64;
  unsigned long i, j;
  mpz_t t, u;
  std::vector<qlattice_basis> basis;
  double st;

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
  basis.assign(N, qlattice_basis());

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

  if (test_speed) {
      /* efficiency test */
      st = seconds ();
      r = 0;
      for (j = 0; j < N; j++)
          for (i = 0; i < N; i++)
              r += fb_root_in_qlattice_127bits (p[i], R[i], invp64[i], basis[j]);
      st = seconds () - st;
      printf ("fb_root_in_qlattice_127bits: %lu tests took %.2fs (r=%" FBPRIME_FORMAT ")\n",
              N * N, st, r);
  } else {
      /* correctness test */
      for (i = 0; i < N; i++)
          for (j = 0; j < N; j++)
          {
              r = fb_root_in_qlattice_127bits (p[i], R[i], invp64[i], basis[j]);
              r31 = fb_root_in_qlattice_31bits (p[i], R[i], invp32[i], basis[j]);
              if (r != r31)
              {
                  fprintf (stderr, "Error for p:=%" FBPRIME_FORMAT "; R:=%" FBPRIME_FORMAT "; a0:=%" PRId64 "; b0:=%" PRId64 "; a1:=%" PRId64 "; b1:=%" PRId64 ";\n",
                          p[i], R[i], basis[j].a0, basis[j].b0, basis[j].a1, basis[j].b1);
                  fprintf (stderr, "fb_root_in_qlattice_127bits gives %" FBPRIME_FORMAT "\n", r);
                  fprintf (stderr, "fb_root_in_qlattice_31bits  gives %" FBPRIME_FORMAT "\n", r31);
                  exit (1);
              }
              rref = ref_fb_root_in_qlattice (p[i], R[i], basis[j]);
              if (r != rref)
              {
                  fprintf (stderr, "Error for p:=%" FBPRIME_FORMAT "; R:=%" FBPRIME_FORMAT "; a0:=%" PRId64 "; b0:=%" PRId64 "; a1:=%" PRId64 "; b1:=%" PRId64 ";\n",
                          p[i], R[i], basis[j].a0, basis[j].b0, basis[j].a1, basis[j].b1);
                  fprintf (stderr, "fb_root_in_qlattice_127bits gives %" FBPRIME_FORMAT "\n", r);
                  fprintf (stderr, "ref_fb_root_in_qlattice gives %" FBPRIME_FORMAT "\n", rref);
                  exit (1);
              }
          }
  }

  free (p);
  free (R);
  free (invp32);
  free (invp64);

  mpz_clear (t);
  mpz_clear (u);
}

/* exercise bugs in fb_root_in_qlattice_127bits */
static void
bug20200225 (void)
{
  uint64_t got, expected;
  fbprime_t p, R;
  uint64_t invp;
  unsigned long a, b;

  /* exercises bug in assembly part of invmod_redc_32 (starting around line
     326 with 2nd #ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM) */
  a = 8088625;
  b = 2163105767;
  expected = 2062858318;
  got = invmod_redc_32 (a, b);
  if (got != expected)
    {
      fprintf (stderr, "Error in invmod_redc_32 for a=%lu b=%lu\n", a, b);
      fprintf (stderr, "Expected %" PRIu64 "\n", expected);
      fprintf (stderr, "Got      %" PRIu64 "\n", got);
      exit (1);
    }

  a = 76285;
  b = 2353808591;
  expected = 2102979166;
  got = invmod_redc_32(a, b);
  if (got != expected)
    {
      fprintf (stderr, "Error in invmod_redc_32 for a:=%lu; b:=%lu;\n",
	       a, b);
      fprintf (stderr, "Expected %" PRIu64 "\n", expected);
      fprintf (stderr, "Got      %" PRIu64 "\n", got);
      exit (1);
    }

  p = 3628762957;
  R = 1702941053;
  invp = 5839589727713490555UL;
  qlattice_basis basis[1];
  basis[0].a0 = -2503835703516628395L;
  basis[0].b0 = 238650852;
  basis[0].a1 = -3992552824749287692L;
  basis[0].b1 = 766395543;
  expected = 987485779;
  got = fb_root_in_qlattice_127bits (p, R, invp, basis[0]);
  if (got != expected)
    {
      fprintf (stderr, "Error in fb_root_in_qlattice_127bits for p:=%" FBPRIME_FORMAT "; R:=%" FBPRIME_FORMAT "; "
	       "a0:=%" PRId64 "; b0:=%" PRId64 "; a1:=%" PRId64 "; b1:=%" PRId64 ";\n",
	       p, R, basis[0].a0, basis[0].b0, basis[0].a1, basis[0].b1);
      fprintf (stderr, "Expected %" PRIu64 "\n", expected);
      fprintf (stderr, "Got      %" PRIu64 "\n", got);
      exit (1);
    }

  /* exercises bug in invmod_redc_32, already tested directly above */
  p = 2163105767;
  R = 1743312141;
  invp = 3235101737;
  basis[0].a0 = -30118114923155082L;
  basis[0].b0 = 749622022;
  basis[0].a1 = 2851499432479966615L;
  basis[0].b1 = 443074848;
  expected = 1879080852;
  got = fb_root_in_qlattice_31bits (p, R, invp, basis[0]);
  if (got != expected)
    {
      fprintf (stderr, "Error in fb_root_in_qlattice_31bits for p=%" FBPRIME_FORMAT " R=%" FBPRIME_FORMAT " "
	       "a0=%" PRId64 " b0=%" PRId64 " a1=%" PRId64 " b1=%" PRId64 "\n",
	       p, R, basis[0].a0, basis[0].b0, basis[0].a1, basis[0].b1);
      fprintf (stderr, "Expected %" PRIu64 "\n", expected);
      fprintf (stderr, "Got      %" PRIu64 "\n", got);
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
      fprintf (stderr, "Error in fb_root_in_qlattice_127bits for p:=%" FBPRIME_FORMAT "; R:=%" FBPRIME_FORMAT "; "
	       "a0:=%" PRId64 "; b0:=%" PRId64 "; a1:=%" PRId64 "; b1:=%" PRId64 ";\n",
	       p, R, basis[0].a0, basis[0].b0, basis[0].a1, basis[0].b1);
      fprintf (stderr, "Expected %" PRIu64 "\n", expected);
      fprintf (stderr, "Got      %" PRIu64 "\n", got);
      exit (1);
    }
}

int
main (int argc, char *argv[])
{
  unsigned long N = 1000;
  unsigned long seed = 1;
  int correctness_only = 0;

  setbuf(stdout, NULL);
  setbuf(stderr, NULL);

  int wild = 0;
  for( ; argc > 1 ; argc--,argv++) {
      if (strcmp(argv[1], "-s") == 0) {
          seed = strtoul (argv[2], NULL, 10);
          argc--,argv++;
      } else if (strcmp(argv[1], "-c") == 0) {
          correctness_only = 1;
      } else if (!wild++) {
          N = strtoul (argv[1], NULL, 10);
      } else {
          fprintf(stderr, "unparsed arg: %s\n", argv[1]);
          exit (EXIT_FAILURE);
      }
  }

  gmp_randstate_t rstate;
  gmp_randinit_default(rstate);

  /* do correctness tests first */
  bug20200225 ();
  gmp_randseed_ui(rstate, seed);
  test_fb_root_in_qlattice_31bits (rstate, 0, N / 10);
  gmp_randseed_ui(rstate, seed);
  test_fb_root_in_qlattice_127bits (rstate, 0, N / 10);

  if (!correctness_only) {
      gmp_randseed_ui(rstate, seed);
      test_fb_root_in_qlattice_31bits (rstate, 1, N);
      gmp_randseed_ui(rstate, seed);
      test_fb_root_in_qlattice_127bits (rstate, 1, N);
  }

  gmp_randclear(rstate);

  return 0;
}



