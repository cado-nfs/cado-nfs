#include "cado.h"
#include <stdint.h>
#include <inttypes.h>
#include <time.h>
#include <gmp.h>
#include "gcd.h"
#include "macros.h"
#include "gmp_aux.h"
#include "test_iter.h"
#include "tests_common.h"

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

void
test_gcd_int64 (const unsigned long iter)
{
  int64_t a, b;
  uint64_t g;
  unsigned long i;
  
  for (i = 0; i < iter; i++)
    {
      a = (i == 0 || i == 1) ? 0 : random_int64 ();
      b = (i == 0 || i == 2) ? 0 : random_int64 ();
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
      a = (i == 0 || i == 1) ? 0 : random_uint64 ();
      b = (i == 0 || i == 2) ? 0 : random_uint64 ();
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
      a = (unsigned long) (i == 0 || i == 1) ? 0 : random_uint64 ();
      b = (unsigned long) (i == 0 || i == 2) ? 0 : random_uint64 ();
      g = gcd_ul (a, b);
      ASSERT_ALWAYS(sizeof(unsigned long) <= sizeof(uint64_t));
      cmp_mpz_gcd_ui64(a, b, g);
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
      int64_t a = (i == 0 || i == 1) ? 0 : random_int64 ();
      int64_t b = (i == 0 || i == 2) ? 0 : random_int64 ();
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
  test_bin_gcd_int64_safe (iter);
  tests_common_clear();
  exit (EXIT_SUCCESS);
}
