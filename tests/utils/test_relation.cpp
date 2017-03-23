#include "cado.h"
#include "utils.h"
#include "relation.hpp"
#include "relation-tools.h"
#include "tests_common.h"

unsigned long
mpz_compute_r (mpz_t a, mpz_t b, mpz_t p)
{
  if (mpz_invert (b, b, p) == 0)
    return mpz_get_ui(p);
  else
  {
    mpz_mul (b, a, b);
    mpz_neg (b, b);
    mpz_mod (b, b, p);

    return mpz_get_ui (b);
  }
 }

int
test_compute_r (unsigned int nb)
{
  int err = 0;
  cxx_mpz tp, ta, tb;

  for (unsigned int i = 0; i < nb; i++)
  {
    uint64_t a;
    uint64_t b;
    unsigned long p;

    /* 5% of tests are for the case where p = 2^k */
    if (i % 20 == 0)
    {
      unsigned int exp = (lrand48() & 0x0000001FUL);
      exp = (exp == 0) ? 1 : exp;
      p = 1UL << exp;
      mpz_set_ui (tp, p);
    }
    else
    {
      mpz_set_ui (tp, lrand48());
      mpz_nextprime(tp, tp);
      p = mpz_get_ui (tp);
    }

    a = random_uint64 ();
    /* 5% of tests are for the case where b = 0 mod p (with b > 0)
     * We do not need to test for free relations as they never go through
     * relation_compute_r
     */
    if (i < (nb / 20))
      b = lrand48() * p;
    else
      b = random_int64 () + 1; /* b > 0 */
    mpz_set_uint64 (ta, a);
    mpz_set_int64 (tb, b);

    unsigned long r = relation_compute_r (a, b, p);

    unsigned long r2 = mpz_compute_r (ta, tb, tp);
    if (r != r2)
    {
      gmp_fprintf (stderr, "ERROR: a=%" PRId64 " b=%" PRIu64" p=%" PRpr "\n"
                   "Got r=%" PRpr " instead of %" PRpr "\n", a, b, p, r, r2);
      err++;
    }
  }
  return err;
}

int test_compute_all_r (unsigned int nb)
{
  int err = 0;
  cxx_mpz tp, ta, tb;

  for (unsigned int i = 0; i < nb; i++)
  {
    relation t1(random_int64(), random_uint64 () + 1); // b > 0

    for (uint8_t k = 0; k <= lrand48() % 5; k++)
    {
      mpz_set_ui (tp, lrand48());
      mpz_nextprime(tp, tp);
      t1.add(0, tp, NULL);
    }
    for (uint8_t k = 0; k <= lrand48() % 5; k++)
    {
      mpz_set_ui (tp, lrand48());
      mpz_nextprime(tp, tp);
      t1.add(1, tp, NULL);
    }

    relation t2 = t1;
    t1.fixup_r();

    for (uint8_t k = 0; k < t2.sides[1].size() ; k++)
    {
      mpz_set(ta, t2[0]);
      mpz_set(tb, t2[1]);
      mpz_set(tp, t2.sides[1][k].p);
      unsigned long r = mpz_compute_r (ta, tb, tp);
      if (r != mpz_get_ui(t1.sides[1][k].r))
      {
        gmp_fprintf (stderr, "ERROR: a=%Zx b=%Zx p=%" PRpr "\n"
                     "Got r=%" PRpr " instead of %" PRpr "\n",
                     t2[0], t2[1], t2.sides[1][k].p, t1.sides[1][k].r, r);
        err++;
      }
    }
  }

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
    uint64_t a = random_uint64 ();
    uint64_t b = random_int64 ();

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
