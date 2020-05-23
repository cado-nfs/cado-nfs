#include "cado.h" // IWYU pragma: keep
#include <time.h>
#include "tests_common.h"
#include "cado_poly.h"

void
test_cado_poly_set ()
{
  cado_poly p, q;

  cado_poly_init (p);
  cado_poly_init (q);

  const double s = 3.1415;

  mpz_set_ui (q->n, 1000000007);
  q->skew = s;
  q->pols[0]->deg = 1;
  mpz_set_si (q->pols[0]->coeff[0], -123128869);
  mpz_set_ui (q->pols[0]->coeff[1], 1000000008);
  q->pols[1]->deg = 2;
  mpz_set_ui (q->pols[1]->coeff[0], 228868283);
  mpz_set_ui (q->pols[1]->coeff[1], 887036294);
  mpz_set_ui (q->pols[1]->coeff[2], 429156742);

  cado_poly_set (p, q);

  ASSERT_ALWAYS (mpz_cmp_ui (p->n, 1000000007) == 0);
  ASSERT_ALWAYS (p->skew == s);
  ASSERT_ALWAYS (mpz_cmp_si (p->pols[0]->coeff[0], -123128869) == 0);
  ASSERT_ALWAYS (mpz_cmp_ui (p->pols[0]->coeff[1], 1000000008) == 0);
  ASSERT_ALWAYS (mpz_cmp_ui (p->pols[1]->coeff[0], 228868283) == 0);
  ASSERT_ALWAYS (mpz_cmp_ui (p->pols[1]->coeff[1], 887036294) == 0);
  ASSERT_ALWAYS (mpz_cmp_ui (p->pols[1]->coeff[2], 429156742) == 0);
  mpz_t m;
  mpz_init(m);
  cado_poly_getm(m, p, p->n);
  ASSERT_ALWAYS (mpz_cmp_ui (m, 123128869) == 0);
  mpz_clear(m);
  cado_poly_clear (p);
  cado_poly_clear (q);
}

int
main (int argc, const char *argv[])
{
  tests_common_cmdline(&argc, &argv, 0);
  test_cado_poly_set ();
  tests_common_clear();
  exit (EXIT_SUCCESS);
}
