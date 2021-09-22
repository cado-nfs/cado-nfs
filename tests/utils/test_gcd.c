#include "cado.h" // IWYU pragma: keep
#include <stdint.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "gcd.h"
#include "macros.h"
#include "gmp_aux.h"
#include "test_iter.h"
#include "tests_common.h"
#include "misc.h"

uint64_t B = 0;

static void
cmp_mpz_gcd_i64(const int64_t a, const int64_t b, const uint64_t g)
{
  mpz_t ma, mb, mg;
  mpz_init (ma);
  mpz_init (mb);
  mpz_init (mg);

  mpz_set_int64 (ma, a);
  mpz_set_int64 (mb, b);
  mpz_gcd (mg, ma, mb);
  /* mpz_gcd() is specified to be non-negative for all argument pairs */
  uint64_t r = mpz_get_uint64 (mg);
  
  if (r != g) {
    fprintf (stderr, "GCD(%" PRId64 ", %" PRId64 ") = %" PRIu64
             " is incorrect, GMP has %" PRIu64 "\n",
             a, b, g, r);
    abort();
  }

  mpz_clear (ma);
  mpz_clear (mb);
  mpz_clear (mg);
}

static void
cmp_mpz_gcd_ui64(const uint64_t a, const uint64_t b, const uint64_t g)
{
  mpz_t ma, mb, mg;
  mpz_init (ma);
  mpz_init (mb);
  mpz_init (mg);

  mpz_set_uint64 (ma, a);
  mpz_set_uint64 (mb, b);
  mpz_gcd (mg, ma, mb);
  
  if (mpz_get_uint64 (mg) != g) {
    fprintf (stderr, "GCD(%" PRIu64 ", %" PRIu64 ") = %" PRIu64
             " is incorrect, GMP has %" PRIu64 "\n",
             a, b, g, mpz_get_uint64 (mg));
    abort();
  }

  mpz_clear (ma);
  mpz_clear (mb);
  mpz_clear (mg);
}

static void
cmp_mpz_xgcd_ui64(const uint64_t u, const uint64_t a, const uint64_t b, const uint64_t g)
{
    cmp_mpz_gcd_ui64(a, b, g);

    if (g == 0) return;

    mpz_t ma, mb, mu;
    mpz_init (ma);
    mpz_init (mb);
    mpz_init (mu);

    mpz_set_uint64 (ma, a);
    mpz_set_uint64 (mb, b);
    mpz_divexact_uint64(mb, mb, g);
    mpz_divexact_uint64(ma, ma, g);
    mpz_invert(mu, ma, mb);

    if (mpz_cmp_uint64(mu, u) != 0) {
        fprintf (stderr, "GCD(a=%" PRIu64 ", b=%" PRIu64 ") = g=%" PRIu64
                " but u=%" PRIu64 " is not an inverse of a/g mod b/g"
                ", GMP has %" PRIu64 "\n",
                a, b, g, u, mpz_get_uint64 (mu));
        abort();
    }

    mpz_clear (ma);
    mpz_clear (mb);
    mpz_clear (mu);
}

void
test_gcd_int64 (const unsigned long iter)
{
  int64_t a, b;
  uint64_t g;
  unsigned long i;
  
  for (i = 0; i < iter; i++)
    {
      a = (i == 0 || i == 1) ? 0 : u64_random (state);
      b = (i == 0 || i == 2) ? 0 : u64_random (state);
      g = gcd_int64 (a, b);
      cmp_mpz_gcd_i64(a, b, g);
    }
}

void
test_gcd_uint64 (const unsigned long iter)
{
  uint64_t a, b, g;
  unsigned long i;
  
  for (i = 0; i < iter; i++)
    {
      a = (i == 0 || i == 1) ? 0 : u64_random (state);
      b = (i == 0 || i == 2) ? 0 : u64_random (state);
      g = gcd_uint64 (a, b);
      cmp_mpz_gcd_ui64(a, b, g);
    }
}

void
test_gcd_ul (const unsigned long iter)
{
  unsigned long a, b, g;
  unsigned long i;

  ASSERT_ALWAYS (sizeof(unsigned long) <= sizeof(uint64_t));
  for (i = 0; i < iter; i++)
    {
      a = (unsigned long) (i == 0 || i == 1) ? 0 : u64_random (state);
      b = (unsigned long) (i == 0 || i == 2) ? 0 : u64_random (state);
      g = gcd_ul (a, b);
      ASSERT_ALWAYS(sizeof(unsigned long) <= sizeof(uint64_t));
      cmp_mpz_gcd_ui64(a, b, g);
    }
}

void
test_xgcd_ul (const unsigned long iter)
{
  unsigned long a, b, g;
  unsigned long i;

  ASSERT_ALWAYS (sizeof(unsigned long) <= sizeof(uint64_t));
  for (i = 0; i < iter; i++)
    {
      a = (unsigned long) (i == 0 || i == 1) ? 0 : u64_random (state);
      b = (unsigned long) (i == 0 || i == 2) ? 0 : u64_random (state);
      unsigned long u;
      g = xgcd_ul (&u, a, b);
      ASSERT_ALWAYS(sizeof(unsigned long) <= sizeof(uint64_t));
      cmp_mpz_xgcd_ui64(u, a, b, g);
    }
}

void
test_bin_gcd_int64_safe_ab (const int64_t a, const int64_t b)
{
  const uint64_t g = bin_gcd_int64_safe (a, b);
  cmp_mpz_gcd_i64 (a, b, g);
}

void
test_bin_gcd_int64_safe (const unsigned long iter)
{
  test_iter_t iter_a, iter_b;
  unsigned long i;

  test_iter_init(iter_a, 100, test_iter_int64_next);
  while (!test_iter_isdone(iter_a)) {
    int64_t a;
    test_iter_next(&a, iter_a);
    test_iter_init(iter_b, 100, test_iter_int64_next);
    while (!test_iter_isdone(iter_b)) {
      int64_t b;
      test_iter_next(&b, iter_b);
      test_bin_gcd_int64_safe_ab(a, b);
    }
  }

  for (i = 0; i < iter; i++)
    {
      int64_t a = (i == 0 || i == 1) ? 0 : u64_random (state);
      int64_t b = (i == 0 || i == 2) ? 0 : u64_random (state);
      test_bin_gcd_int64_safe_ab(a,b);
    }
}

int
main (int argc, const char *argv[])
{
  unsigned long iter = 200000;
  tests_common_cmdline(&argc, &argv, PARSE_SEED | PARSE_ITER);
  tests_common_get_iter(&iter);
  test_gcd_int64 (iter);
  test_gcd_uint64 (iter);
  test_gcd_ul (iter);
  test_xgcd_ul (iter);
  test_bin_gcd_int64_safe (iter);
  tests_common_clear();
  exit (EXIT_SUCCESS);
}
