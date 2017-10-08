#include <stdint.h>
#include <stdio.h>
#include <inttypes.h>
#include "cado.h"
#include "utils.h"
#include "macros.h"
#include "ideal.h"

int main()
{
  mpz_poly h;
  mpz_poly_init(h, 0);
  mpz_poly_setcoeff_si(h, 1, 1);
  mpz_poly_setcoeff_si(h, 0, 3);
  uint64_t r = 5;
  unsigned int t = 4;

  ideal_t ideal;
  ideal_init(ideal);
  ideal_set_part(ideal, r, h);
  ASSERT_ALWAYS(ideal->r == r);
  ASSERT_ALWAYS(mpz_poly_cmp(ideal->h, h) == 0);
  ideal_clear(ideal);

  ideal_1_t ideal_1;
  ideal_1_init(ideal_1);
  ideal_1_set_part(ideal_1, r, h, t);
  ASSERT_ALWAYS(ideal_1->ideal->r == r);
  ASSERT_ALWAYS(mpz_poly_cmp(ideal_1->ideal->h, h) == 0);
  // We skip this test on i386: due to the x87 extended precision,
  // comparing doubles is a nightmare.
#ifndef __i386
  ASSERT_ALWAYS(ideal_1->log == log2((double)r));
#endif
  //Test Tr.
  mpz_t Tr;
  mpz_init(Tr);
  mpz_t coeff;
  mpz_init(coeff);
  mpz_poly_getcoeff(coeff, 0, h);
  mpz_set(Tr, coeff);
  mpz_mul_si(coeff, coeff, -1);
  ASSERT_ALWAYS(mpz_cmp(ideal_1->Tr[0], Tr) == 0);
  for (unsigned int i = 1; i < t - 1; i++) {
    mpz_mul(Tr, Tr, coeff);
    ASSERT_ALWAYS(mpz_cmp(ideal_1->Tr[i], Tr) == 0);
  }
  mpz_clear(Tr);
  mpz_clear(coeff);
  ideal_1_clear(ideal_1, t);

  ideal_u_t ideal_u;
  ideal_u_init(ideal_u);
  mpz_poly_setcoeff_si(h, 4, 1);
  ideal_u_set_part(ideal_u, r, h, t);
  ASSERT_ALWAYS(ideal_u->ideal->r == r);
  ASSERT_ALWAYS(mpz_poly_cmp(ideal_u->ideal->h, h) == 0);
  // We skip this test on i386: due to the x87 extended precision,
  // comparing doubles is a nightmare.
#ifndef __i386
  ASSERT_ALWAYS(ideal_u->log == log2(pow((double)r,(double)h->deg)));
#endif
  //TODO: Test Tr.

  ideal_u_clear(ideal_u, t);

  ideal_pr_t ideal_pr;
  ideal_pr_init(ideal_pr);
  ideal_pr_set_part(ideal_pr, r, t);
  ASSERT_ALWAYS(ideal_pr->ideal->r == r);
  ASSERT_ALWAYS(ideal_pr->ideal->h->deg == 1);
  ideal_pr_clear(ideal_pr, t);

  ideal_spq_t ideal_spq;
  ideal_spq_init(ideal_spq);
  h->deg = 1;
  ideal_spq_set_part(ideal_spq, r, h, t, 0);
  ASSERT_ALWAYS(ideal_spq_get_q(ideal_spq) == r);
  ASSERT_ALWAYS(ideal_spq_get_deg_g(ideal_spq) == 1);
  ideal_spq_clear(ideal_spq, t);

  mpz_poly_clear(h);
  return 0;
}
