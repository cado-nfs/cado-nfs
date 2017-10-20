#include "cado.h"
#include <stdint.h>
#include <time.h>
#include <limits.h>
#include "gmp_aux.h"
#include "macros.h"
#include "tests_common.h"

static void
test_mpz_set_uint64 ()
{
  uint64_t q;
  mpz_t z;

  mpz_init (z);

  q = 0;
  mpz_set_uint64 (z, q);
  if (mpz_cmp_ui (z, 0) != 0)
    abort();

  q = 4294967295UL;
  mpz_set_uint64 (z, q);
  if (mpz_cmp_d (z, 4294967295.0) != 0)
    abort ();

  q ++;
  mpz_set_uint64 (z, q);
  if (mpz_cmp_d (z, 4294967296.0) != 0)
    abort ();

  q = (q-1) * (q-1) + 2 * (q-1); /* 2^64-1 */
  mpz_set_uint64 (z, q);
  mpz_add_ui (z, z, 1);
  if (mpz_cmp_d (z, 18446744073709551616.0) != 0)
    abort ();
  mpz_clear (z);
}

static void
test_mpz_set_int64 ()
{
  int64_t q;
  mpz_t z;

  mpz_init (z);

  q = 0;
  mpz_set_int64 (z, q);
  if (mpz_cmp_ui (z, 0) != 0)
    abort();

  q = -1;
  mpz_set_int64 (z, q);
  if (mpz_cmp_si (z, -1) != 0)
    abort();

  q = 2147483647L;
  mpz_set_int64 (z, q);
  if (mpz_cmp_d (z, 2147483647.0) != 0)
    abort ();

  q = -2147483648L;
  mpz_set_int64 (z, q);
  if (mpz_cmp_d (z, -2147483648.0) != 0)
    abort ();

  q = 2147483648L;
  mpz_set_int64 (z, q);
  if (mpz_cmp_d (z, 2147483648.0) != 0)
    abort ();

  q = ((q-1) << 32) + 2 * (q-1) + 1; /* 2^63-1 */
  mpz_set_int64 (z, q);
  mpz_add_ui (z, z, 1);
  if (mpz_cmp_d (z, 9223372036854775808.0) != 0)
    abort ();

  q = -q-1; /* -2^63 */
  mpz_set_int64 (z, q);
  if (mpz_cmp_d (z, -9223372036854775808.0) != 0)
    abort ();

  mpz_clear (z);
}

void
test_mpz_get_uint64 (const unsigned long iter)
{
  uint64_t q, r;
  mpz_t z;
  unsigned long i;

  mpz_init (z);
  for (i = 0; i < iter; i++)
    {
      q = i & 3; /* two bits */
      q = (q << 31) + lrand48 (); /* 33 bits */
      q = (q << 31) + lrand48 (); /* 64 bits */
      mpz_set_uint64 (z, q);
      r = mpz_get_uint64 (z);
      if (r != q)
        abort();
    }
  mpz_clear (z);
}

void
test_mpz_get_int64 (const unsigned long iter)
{
  int64_t q, r;
  mpz_t z;
  unsigned long i;

  mpz_init (z);
  for (i = 0; i < iter; i++)
    {
      q = mrand48 (); /* [-2^31, 2^31] */
      q = (q << 32) + mrand48 (); /* 63 bits */
      mpz_set_int64 (z, q);
      r = mpz_get_int64 (z);
      if (r != q)
        abort();
    }
  mpz_clear (z);
}

void
test_mpz_fits_int64_p ()
{
  mpz_t z;


  mpz_init (z);
  /* check 2^63-1 fits but not 2^63 */
  mpz_set_ui (z, 1);
  mpz_mul_2exp (z, z, 63);
  if (mpz_fits_int64_p (z))
    abort();
  mpz_sub_ui (z, z, 1);
  if (mpz_fits_int64_p (z) == 0)
    abort();
  /* check -2^63 fits but not -2^63-1 */
  mpz_add_ui (z, z, 1);
  mpz_neg (z, z);
  if (mpz_fits_int64_p (z) == 0)
    abort();
  mpz_sub_ui (z, z, 1);
  if (mpz_fits_int64_p (z))
    abort();
  mpz_clear (z);
}

void
test_mpz_mul_uint64 (const unsigned long iter)
{
  mpz_t a, b, aa, cc;
  uint64_t c;
  unsigned long i;

  mpz_init (a);
  mpz_init (b);
  mpz_init (aa);
  mpz_init (cc);
  for (i = 0; i < iter; i++)
    {
      mpz_urandomb (b, state, 128);
      c = i & 3;
      c = (c << 31) + lrand48 ();
      c = (c << 31) + lrand48 ();
      mpz_mul_uint64 (a, b, c);
      mpz_set_uint64 (cc, c);
      mpz_mul (aa, b, cc);
      if (mpz_cmp (aa, a) != 0)
        abort ();
    }
  mpz_clear (a);
  mpz_clear (b);
  mpz_clear (aa);
  mpz_clear (cc);
}

void
test_mpz_mul_int64 (const unsigned long iter)
{
  mpz_t a, b, aa, cc;
  int64_t c;
  unsigned long i;

  mpz_init (a);
  mpz_init (b);
  mpz_init (aa);
  mpz_init (cc);
  for (i = 0; i < iter; i++)
    {
      mpz_urandomb (b, state, 128);
      c = mrand48 ();
      c = (c << 32) + mrand48 ();
      mpz_mul_int64 (a, b, c);
      mpz_set_int64 (cc, c);
      mpz_mul (aa, b, cc);
      if (mpz_cmp (aa, a) != 0)
        abort ();
    }
  mpz_clear (a);
  mpz_clear (b);
  mpz_clear (aa);
  mpz_clear (cc);
}

void
test_mpz_addmul_int64 (const unsigned long iter)
{
  mpz_t a, b, aa, cc;
  int64_t c;
  unsigned long i;

  mpz_init (a);
  mpz_init (b);
  mpz_init (aa);
  mpz_init (cc);
  for (i = 0; i < iter; i++)
    {
      mpz_urandomb (b, state, 128);
      c = mrand48 ();
      c = (c << 32) + mrand48 ();
      mpz_addmul_int64 (a, b, c);
      mpz_set_int64 (cc, c);
      mpz_addmul (aa, b, cc);
      if (mpz_cmp (aa, a) != 0)
        abort ();
    }
  mpz_clear (a);
  mpz_clear (b);
  mpz_clear (aa);
  mpz_clear (cc);
}

/* this function tests both ulong_nextprime and ulong_isprime */
void
test_ulong_nextprime (const unsigned long iter)
{
  unsigned long q, r, s;
  unsigned long i;

  q = lrand48 () % 300000000;
  for (i = 0; i < iter && q < 300000000; i++)
    {
      for (s = q + 1; s != 0 && ulong_isprime (s) == 0; s++);
      r = ulong_nextprime (q);
      if (r != s)
        abort ();
      q = r;
    }
  ASSERT_ALWAYS(ulong_nextprime (ULONG_MAX) == 0);
}

void
test_nbits ()
{
  uintmax_t p;
  int i;

  p = i = 1;
  while (p < 2*p)
    {
      if (nbits (p) != i)
        abort ();
      if (nbits (p-1) != i-1)
        abort ();
      p *= 2;
      i += 1;
    }
}

void
test_mpz_ndiv_q (const unsigned long iter)
{
  mpz_t q, n, d, t;
  unsigned long i;

  mpz_init (q);
  mpz_init (n);
  mpz_init (d);
  mpz_init (t);

  for (i = 0; i < iter; i++)
    {
      mpz_urandomb (n, state, 128);
      do mpz_urandomb (d, state, 64); while (mpz_cmp_ui (d, 0) == 0);
      mpz_ndiv_q (q, n, d);
      /* check |n-q*d| <= |d|/2 */
      mpz_mul (t, q, d);
      mpz_sub (t, n, t);
      mpz_mul_2exp (t, t, 1);
      if (mpz_cmpabs (t, d) > 0)
        abort ();
    }

  mpz_clear (q);
  mpz_clear (n);
  mpz_clear (d);
  mpz_clear (t);
}

static void
test_ulong_nextcomposite (void)
{
  unsigned long q;

  q = ulong_nextcomposite (2, 2);
  ASSERT_ALWAYS(q == 4);

  q = ulong_nextcomposite (2, 3);
  ASSERT_ALWAYS(q == 9);

  q = ulong_nextcomposite (2, 4);
  ASSERT_ALWAYS(q == 25);

  q = ulong_nextcomposite (1000000, 1000);
  ASSERT_ALWAYS(q == 1018081);

  q = 1000000-1;
  for (int i = 0; i < 1000; i++)
    q = ulong_nextcomposite (q+1, 1000);
  ASSERT_ALWAYS(q == 1475069);
}

static void
test_next_mpz_with_factor_constraints (void)
{
  mpz_t r, s;
  mpz_init(r);
  mpz_init(s);
  unsigned long fact[3];
  int ret;

  mpz_set_ui(s, 25631719);
  ret = next_mpz_with_factor_constraints(r, fact, s, 0, 256, 33554432);
  // should get 25631719 (it is a prime less than 2^25)
  printf("%d %lu\n", ret, fact[0]);
  ASSERT_ALWAYS(ret == 1 && fact[0] == 25631719);
  ASSERT_ALWAYS(mpz_cmp_ui(r, 25631719) == 0);
  
  ret = next_mpz_with_factor_constraints(r, fact, s, 1, 256, 33554432);
  // should get 25631731 (it is a prime, and all intermediate values have
  // small factors)
  ASSERT_ALWAYS(ret == 1 && fact[0] == 25631731);
  ASSERT_ALWAYS(mpz_cmp_ui(r, 25631731) == 0);

  ret = next_mpz_with_factor_constraints(r, fact, s, 1, 20, 33554432);
  // should get 25631729, whose smallest factor is 23.
  ASSERT_ALWAYS(ret == 2 && fact[0] == 23 && fact[1] == 1114423);
  ASSERT_ALWAYS(mpz_cmp_ui(r, 25631729) == 0);

  ret = next_mpz_with_factor_constraints(r, fact, s, 0, 20, 1114425);
  // should get 25631729 again, whose largest prime factor fits.
  ASSERT_ALWAYS(ret == 2 && fact[0] == 23 && fact[1] == 1114423);
  ASSERT_ALWAYS(mpz_cmp_ui(r, 25631729) == 0);

  mpz_set_ui(s, 25631719);
  ret = next_mpz_with_factor_constraints(r, fact, s, 0, 25, 100000);
  // should get 25631737 = 29*307*2879
  ASSERT_ALWAYS(ret == 3 && fact[0] == 29 && fact[1] == 307 && fact[2] == 2879);
  ASSERT_ALWAYS(mpz_cmp_ui(r, 25631737) == 0);

  mpz_set_ui(s, 5426767); // this is 31^2*5647; must be skipped
  ret = next_mpz_with_factor_constraints(r, fact, s, 0, 25, 100000);
  ASSERT_ALWAYS(mpz_cmp_ui(r, 5426777) == 0);

  mpz_set_ui(s, 25030009); // this is 5003^2 must be skipped
  ret = next_mpz_with_factor_constraints(r, fact, s, 0, 100, 100000);
  ASSERT_ALWAYS(mpz_cmp_ui(r, 25030039) == 0);

  mpz_clear(r);
  mpz_clear(s);
}




int
main (int argc, const char *argv[])
{
  unsigned long iter = 1000000;
  tests_common_cmdline(&argc, &argv, PARSE_SEED | PARSE_ITER);
  tests_common_get_iter(&iter);
  test_mpz_set_uint64 ();
  test_mpz_set_int64 ();
  test_mpz_get_uint64 (iter);
  test_mpz_get_int64 (iter);
  test_mpz_fits_int64_p ();
  test_mpz_mul_uint64 (iter);
  test_mpz_mul_int64 (iter);
  test_mpz_addmul_int64 (iter);
  test_ulong_nextprime (iter / 20);
  test_nbits ();
  test_mpz_ndiv_q (iter);
  test_ulong_nextcomposite ();
  test_next_mpz_with_factor_constraints ();
  tests_common_clear ();
  exit (EXIT_SUCCESS);
}
