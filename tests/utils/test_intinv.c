#include "cado.h"
#include "modredc_ul.h"
#include "test_iter.h"
#include "tests_common.h"
#include "gcd.h"
#include "gmp_aux.h"

static void
test_modredcul_intinv (unsigned long iter)
{
  unsigned long x, y, g;
  residueredcul_t r, A;
  modulusredcul_t m;
  int ret;
  mpz_t xx, yy;

  mpz_init (xx);
  mpz_init (yy);
  while (iter--)
    {
      x = random_uint64 ();
      y = random_uint64 ();
      if (x < y)
        {
          g = x;
          x = y;
          y = x;
        }
      if ((x & 1) == 0)
        x ++; /* x should be odd */
      if (x == y)
        continue;
#if ULONG_BITS==64
      if (iter == 0)
        {
          x = 14170321801878169673UL;
          y = 8493495769022353815UL;
        }
#endif
      modredcul_initmod_ul (m, x); /* x should be odd */
      modredcul_init (A, m);
      modredcul_init (r, m);
      modredcul_intset_ul (A, y);
      ret = modredcul_intinv (r, A, m);
      g = gcd_uint64 (x, y);
      if (g != 1)
        ASSERT_ALWAYS (ret == 0);
      else
        {
          ASSERT_ALWAYS (ret != 0);
          mpz_set_uint64 (xx, x);
          mpz_set_uint64 (yy, y);
          mpz_invert (yy, yy, xx);
          g = mpz_get_uint64 (yy);
          ASSERT_ALWAYS (g == modredcul_intget_ul (r));
        }
      modredcul_clear (A, m);
      modredcul_clear (r, m);
      modredcul_clearmod (m);
    }
  mpz_clear (xx);
  mpz_clear (yy);
}

int
main (int argc, const char *argv[])
{
  unsigned long iter = 200000;

  tests_common_cmdline (&argc, &argv, PARSE_SEED | PARSE_ITER);
  tests_common_get_iter (&iter);
  test_modredcul_intinv (iter);
  tests_common_clear ();
  exit (EXIT_SUCCESS);
}


