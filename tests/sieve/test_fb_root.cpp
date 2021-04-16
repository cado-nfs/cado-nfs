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

/*
 * R encodes (p:1) if R<p, or (1:R-p) if R >= p
 *
 * when R < p (affine root),
 *      return num/den%p with num=R*b1-a1 den=a0-R*b0
 * when R >= p (projective root), 
 *      return num/den%p with num=b1-(R-p)*a1 den=(R-p)*a0-b0
 *
 * It is assumed that R*b1-a1 and a0-R*b0 cannot be both multiples of p,
 * since the matrix (a0,b0,a1,b1) has determinant coprime to p, and
 * gcd(1,R) == 1
 *
 * if den in either expression is divisile by p, then the
 * resulting root is projective. Therefore, we return p + den/num
 * */
static cxx_mpz
ref_fb_root_in_qlattice (fbprime_t p, fbprime_t R, qlattice_basis basis)
{
  cxx_mpz num, den;

  if (R < p) {
      den = -basis.b0;
      mpz_mul_ui (den, den, R);
      mpz_add_int64 (den, den, basis.a0);

      num = basis.b1;
      mpz_mul_ui (num, num, R);
      mpz_sub_int64 (num, num, basis.a1);
  } else {
      den = basis.a0;
      mpz_mul_ui (den, den, R-p);
      mpz_sub_int64 (den, den, basis.b0);

      num = -basis.a1;
      mpz_mul_ui (num, num, R-p);
      mpz_add_int64 (num, num, basis.b1);
  }

  if (mpz_gcd_ui(NULL, den, p) != 1) {
      /* root is projective */
      mpz_invert (num, num, cxx_mpz(p));
      mpz_mul(num, num, den);
      mpz_mod_ui(num, num, p);
      mpz_add_ui(num, num, p);
  } else {
      mpz_invert (den, den, cxx_mpz(p));
      mpz_mul(num, num, den);
      mpz_mod_ui(num, num, p);
  }
  return num;
}

static void
test_fb_root_in_qlattice_31bits (gmp_randstate_t rstate, int test_speed, unsigned long N)
{
  fbprime_t *p, *R, r;
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
              cxx_mpz rref = ref_fb_root_in_qlattice (p[i], R[i], basis[j]);
              if (mpz_cmp_ui(rref, r) != 0)
              {
                  fprintf (stderr, "Error for p:=%" FBPRIME_FORMAT "; R:=%" FBPRIME_FORMAT "; a0:=%" PRId64 "; b0:=%" PRId64 "; a1:=%" PRId64 "; b1:=%" PRId64 ";\n",
                          p[i], R[i], basis[j].a0, basis[j].b0, basis[j].a1, basis[j].b1);
                  fprintf (stderr, "fb_root_in_qlattice_31bits gives %" FBPRIME_FORMAT "\n", r);
                  gmp_fprintf (stderr, "ref_fb_root_in_qlattice gives %Zd\n", mpz_srcptr(rref));
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
  fbprime_t *p, *R, r, r31;
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
              cxx_mpz rref = ref_fb_root_in_qlattice (p[i], R[i], basis[j]);
              if (mpz_cmp_ui(rref, r) != 0)
              {
                  fprintf (stderr, "Error for p:=%" FBPRIME_FORMAT "; R:=%" FBPRIME_FORMAT "; a0:=%" PRId64 "; b0:=%" PRId64 "; a1:=%" PRId64 "; b1:=%" PRId64 ";\n",
                          p[i], R[i], basis[j].a0, basis[j].b0, basis[j].a1, basis[j].b1);
                  fprintf (stderr, "fb_root_in_qlattice_127bits gives %" FBPRIME_FORMAT "\n", r);
                  gmp_fprintf (stderr, "ref_fb_root_in_qlattice gives %Zd\n", mpz_srcptr(rref));
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
  uint64_t got;
  cxx_mpz expected;
  fbprime_t p, R;
  uint64_t invp;
  unsigned long a, b;

  /* exercises bug in assembly part of invmod_redc_32 (starting around line
     326 with 2nd #ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM) */
  a = 8088625;
  b = 2163105767;
  expected = 2062858318;
  got = invmod_redc_32 (a, b);
  if (mpz_cmp_ui(expected, got) != 0)
    {
      fprintf (stderr, "Error in invmod_redc_32 for a=%lu b=%lu\n", a, b);
      gmp_fprintf (stderr, "Expected %Zd\n", mpz_srcptr(expected));
      fprintf (stderr, "Got      %" PRIu64 "\n", got);
      exit (1);
    }

  a = 76285;
  b = 2353808591;
  expected = 2102979166;
  got = invmod_redc_32(a, b);
  if (mpz_cmp_ui(expected, got) != 0)
    {
      fprintf (stderr, "Error in invmod_redc_32 for a:=%lu; b:=%lu;\n",
	       a, b);
      gmp_fprintf (stderr, "Expected %Zd\n", mpz_srcptr(expected));
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
  if (mpz_cmp_ui(expected, got) != 0)
    {
      fprintf (stderr, "Error in fb_root_in_qlattice_127bits for p:=%" FBPRIME_FORMAT "; R:=%" FBPRIME_FORMAT "; "
	       "a0:=%" PRId64 "; b0:=%" PRId64 "; a1:=%" PRId64 "; b1:=%" PRId64 ";\n",
	       p, R, basis[0].a0, basis[0].b0, basis[0].a1, basis[0].b1);
      gmp_fprintf (stderr, "Expected %Zd\n", mpz_srcptr(expected));
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
  if (mpz_cmp_ui(expected, got) != 0)
    {
      fprintf (stderr, "Error in fb_root_in_qlattice_31bits for p=%" FBPRIME_FORMAT " R=%" FBPRIME_FORMAT " "
	       "a0=%" PRId64 " b0=%" PRId64 " a1=%" PRId64 " b1=%" PRId64 "\n",
	       p, R, basis[0].a0, basis[0].b0, basis[0].a1, basis[0].b1);
      gmp_fprintf (stderr, "Expected %Zd\n", mpz_srcptr(expected));
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
  if (mpz_cmp_ui(expected, got) != 0)
    {
      fprintf (stderr, "Error in fb_root_in_qlattice_127bits for p:=%" FBPRIME_FORMAT "; R:=%" FBPRIME_FORMAT "; "
	       "a0:=%" PRId64 "; b0:=%" PRId64 "; a1:=%" PRId64 "; b1:=%" PRId64 ";\n",
	       p, R, basis[0].a0, basis[0].b0, basis[0].a1, basis[0].b1);
      gmp_fprintf (stderr, "Expected %Zd\n", mpz_srcptr(expected));
      fprintf (stderr, "Got      %" PRIu64 "\n", got);
      exit (1);
    }

  /* This is playing with the 32-bit limit */
  p = 3486784401;       // 3^20
  R = 2009510725;
  invp = 898235023;
  basis[0].a0=-1353180941;
  basis[0].b0=-5;
  basis[0].a1=-223660881;
  basis[0].b1=8;
  expected = ref_fb_root_in_qlattice (p, R, basis[0]);
  got = fb_root_in_qlattice_127bits (p, R, invp, basis[0]);
  if (mpz_cmp_ui(expected, got) != 0)
    {
      fprintf (stderr, "Error in fb_root_in_qlattice_127bits for p:=%" FBPRIME_FORMAT "; R:=%" FBPRIME_FORMAT "; "
	       "a0:=%" PRId64 "; b0:=%" PRId64 "; a1:=%" PRId64 "; b1:=%" PRId64 ";\n",
	       p, R, basis[0].a0, basis[0].b0, basis[0].a1, basis[0].b1);
      gmp_fprintf (stderr, "Expected %Zd\n", mpz_srcptr(expected));
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



