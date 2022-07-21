#include "cado.h" // IWYU pragma: keep
#include <cinttypes>        // for PRId64, PRIu64
#include <cstdint>          // for uint64_t, uint8_t
#include <cstring>          // for strcmp
#include <cstdio>
#include <cstdlib>
#include <memory>            // for allocator_traits<>::value_type
#include <vector>            // for vector
#include <gmp.h>
#include "cxx_mpz.hpp"       // for cxx_mpz
#include "gmp_aux.h"         // for mpz_set_int64, mpz_set_uint64
#include "portability.h" // IWYU pragma: keep
#include "relation-tools.h"
#include "relation.hpp"
#include "tests_common.h"
#include "typedefs.h"        // for PRpr
#include "misc.h"        // for PRpr

unsigned long
mpz_compute_r (mpz_t a, mpz_t b, mpz_t p)
{
  if (mpz_invert (b, b, p) == 0)
    return mpz_get_ui(p);
  else
  {
    mpz_mul (b, a, b);
    mpz_mod (b, b, p);

    return mpz_get_ui (b);
  }
 }

int
test_compute_r (unsigned int nb)
{
  int err = 0;
  mpz_t tp, ta, tb;
  mpz_init (tp);
  mpz_init (ta);
  mpz_init (tb);

  for (unsigned int i = 0; i < nb; i++)
  {
    uint64_t a;
    uint64_t b;
    unsigned long p;

    /* 5% of tests are for the case where p = 2^k */
    if (i % 20 == 0)
    {
      unsigned int exp = gmp_urandomb_ui(state, 5);
      exp = (exp == 0) ? 1 : exp;
      p = 1UL << exp;
      mpz_set_ui (tp, p);
    }
    else
    {
      mpz_set_ui (tp, gmp_urandomb_ui(state, 31));
      mpz_nextprime(tp, tp);
      p = mpz_get_ui (tp);
    }

    a = i64_random (state);
    /* 5% of tests are for the case where b = 0 mod p (with b > 0)
     * We do not need to test for free relations as they never go throught
     * relation_compute_r
     */
    if (i < (nb / 20))
      b = gmp_urandomb_ui(state, 31) * p;
    else {
      b = u64_random (state);
      b += !b;
    }
    mpz_set_int64 (ta, a);
    mpz_set_uint64 (tb, b);

    unsigned long r = relation_compute_r (a, b, p);

    unsigned long r2 = mpz_compute_r (ta, tb, tp);
    if (r != r2)
    {
      gmp_fprintf (stderr, "ERROR: a=%" PRId64 " b=%" PRIu64" p=%" PRpr "\n"
                   "Got r=%" PRpr " instead of %" PRpr "\n", a, b, p, r, r2);
      err++;
    }
  }
  mpz_clear (tp);
  mpz_clear (ta);
  mpz_clear (tb);
  return err;
}

int test_compute_all_r (unsigned int nb)
{
  int err = 0;
  mpz_t tp, ta, tb;
  mpz_init(tp);
  mpz_init(ta);
  mpz_init(tb);

  for (unsigned int i = 0; i < nb; i++)
  {
      int64_t a = i64_random(state);
      uint64_t b = u64_random(state); b += !b;
      relation t1(a, b);

    for (uint8_t k = 0; k <= gmp_urandomm_ui(state, 5); k++)
    {
      mpz_set_ui (tp, gmp_urandomb_ui(state, 31));
      mpz_nextprime(tp, tp);
      t1.add(0, tp, NULL);
    }
    for (uint8_t k = 0; k <= gmp_urandomm_ui(state, 5); k++)
    {
      mpz_set_ui (tp, gmp_urandomb_ui(state, 31));
      mpz_nextprime(tp, tp);
      t1.add(1, tp, NULL);
    }

    relation t2 = t1;
    t1.fixup_r();

    for (uint8_t k = 0; k < t2.sides[1].size() ; k++)
    {
      mpz_set_int64 (ta, t2.a);
      mpz_set_uint64 (tb, t2.b);
      mpz_set (tp, t2.sides[1][k].p);
      unsigned long r = mpz_compute_r (ta, tb, tp);
      if (r != mpz_get_ui(t1.sides[1][k].r))
      {
        gmp_fprintf (stderr, "ERROR: a=%" PRId64 " b=%" PRIu64" p=%" PRpr "\n"
                     "Got r=%" PRpr " instead of %" PRpr "\n",
                     t2.a, t2.b,
                     (mpz_srcptr) t2.sides[1][k].p,
                     (mpz_srcptr) t1.sides[1][k].r, r);
        err++;
      }
    }
  }

  mpz_clear (ta);
  mpz_clear (tb);
  mpz_clear (tp);

  return err;
}

int
check_str_err (const char * s1, const char * s2, mpz_t t)
{
  if (strcmp(s1, s2) != 0)
  {
    gmp_fprintf(stderr, "ERROR with integer %Zd: got \"%s\" instead of "
                        "\"%s\"\n", t, s2, s1);
    return 1;
  }
  else
    return 0;

}

int
test_conversion (unsigned int nb)
{
  int err = 0;
  char *s1, *s2, *tmp;
  s1 = (char *) malloc (25 * sizeof(char));
  s2 = (char *) malloc (25 * sizeof(char));
  mpz_t t;
  mpz_init(t);


  for (unsigned int i = 0; i < nb; i++)
  {
    uint64_t a = u64_random(state);
    int64_t b = i64_random(state);

    mpz_set_uint64 (t, a);

    mpz_get_str(s1, 10, t);
    tmp = u64toa10 (s2, a);
    *tmp = '\0';
    err += check_str_err (s1, s2, t);

    mpz_get_str(s1, 16, t);
    tmp = u64toa16 (s2, a);
    *tmp = '\0';
    err += check_str_err (s1, s2, t);


    mpz_set_int64 (t, b);

    mpz_get_str(s1, 10, t);
    tmp = d64toa10 (s2, b);
    *tmp = '\0';
    err += check_str_err (s1, s2, t);

    mpz_get_str(s1, 16, t);
    tmp = d64toa16 (s2, b);
    *tmp = '\0';
    err += check_str_err (s1, s2, t);
  }
  mpz_clear(t);
  free(s1);
  free(s2);

  return err;
}

int
main (int argc, const char *argv[])
{
  int err = 0;
  unsigned long iter = 10000;

  tests_common_cmdline(&argc, &argv, PARSE_SEED | PARSE_ITER);
  tests_common_get_iter(&iter);

  err += test_compute_r (iter);
  err += test_compute_all_r (iter / 10);
  err += test_conversion (iter);

  if (err)
    fprintf (stderr, "# %d erro%s found\n", err, (err == 1) ? "r" : "rs");
  tests_common_clear();
  return (err) ? EXIT_FAILURE : EXIT_SUCCESS;
}
