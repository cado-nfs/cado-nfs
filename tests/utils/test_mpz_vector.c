#include "cado.h"
#include "macros.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <gmp.h>
#include "mpz_poly.h"
#include "mpz_vector.h"

void tests_common()
{
  unsigned int dim = 3;

  mpz_vector_t c0;
  mpz_vector_init(c0, dim);

  for (unsigned int i = 0; i < c0->dim; i++) {
    ASSERT_ALWAYS(mpz_cmp_ui(c0->c[i], 0) == 0);
  }

  mpz_t z;
  mpz_init(z);
  mpz_set_str(z, "17", 10);
  mpz_vector_setcoordinate(c0, 0, z);
  mpz_clear(z);

  uint64_t u = 42;
  mpz_vector_setcoordinate_ui(c0, 1, UINT_MAX);
  mpz_vector_setcoordinate_uint64(c0, 2, u);

  ASSERT_ALWAYS(mpz_cmp_ui(c0->c[0], 17) == 0);
  ASSERT_ALWAYS(!mpz_vector_is_coordinate_zero(c0, 0));
  ASSERT_ALWAYS(mpz_cmp_ui(c0->c[1], UINT_MAX) == 0);
  ASSERT_ALWAYS(mpz_cmp_ui(c0->c[2], 42) == 0);

  mpz_vector_t c1;
  mpz_vector_init(c1, dim);
  mpz_vector_swap(c0, c1);

  for (unsigned int i = 0; i < c0->dim; i++) {
    ASSERT_ALWAYS(mpz_cmp_ui(c0->c[i], 0) == 0);
  }
  ASSERT_ALWAYS(mpz_cmp_ui(c1->c[0], 17) == 0);
  ASSERT_ALWAYS(mpz_cmp_ui(c1->c[1], UINT_MAX) == 0);
  ASSERT_ALWAYS(mpz_cmp_ui(c1->c[2], 42) == 0);

  mpz_vector_set(c1, c0);

  for (unsigned int i = 0; i < c0->dim; i++) {
    ASSERT_ALWAYS(mpz_cmp_ui(c0->c[i], 0) == 0);
    ASSERT_ALWAYS(mpz_vector_is_coordinate_zero(c0, i));
  }

  mpz_vector_clear(c0);
  mpz_vector_clear(c1);
}

void test_dot_product()
{
  unsigned int dim = 3;

  mpz_vector_t c0;
  mpz_vector_init(c0, dim);
  mpz_vector_t c1;
  mpz_vector_init(c1, dim);

  mpz_vector_setcoordinate_ui(c0, 1, 2);
  mpz_vector_setcoordinate_ui(c0, 2, 5);

  mpz_vector_setcoordinate_ui(c1, 0, 5);

  mpz_t tmp;
  mpz_init(tmp);
  
  mpz_vector_dot_product(tmp, c0, c1);

  ASSERT_ALWAYS(mpz_cmp_ui(tmp, 0) == 0);
 
  mpz_clear(tmp);

  mpz_vector_clear(c0);
  mpz_vector_clear(c1);
}

void test_norm()
{
  mpz_vector_t c;
  mpz_vector_init(c, 3);

  mpz_vector_setcoordinate_ui(c, 1, 2);
  mpz_vector_setcoordinate_ui(c, 2, 5);

  mpz_t tmp;
  mpz_init(tmp);
  
  mpz_vector_norm(tmp, c);

  ASSERT_ALWAYS(mpz_cmp_ui(tmp, 29) == 0);
 
  mpz_clear(tmp);

  mpz_vector_clear(c);
}

void test_mpz_poly()
{
  mpz_vector_t c;
  mpz_vector_init(c, 3);
  mpz_poly_t p;
  mpz_poly_init(p, -1);

  mpz_vector_setcoordinate_ui(c, 1, 2);
  mpz_vector_setcoordinate_ui(c, 2, 5);

  mpz_vector_get_mpz_poly(p, c);
  for (unsigned int i = 0; i < c->dim; i++) {
    ASSERT_ALWAYS(mpz_cmp(p->coeff[i], c->c[i]) == 0);
  }
  ASSERT_ALWAYS(p->deg == (int)(c->dim - 1));

  mpz_vector_setcoordinate_ui(c, 2, 0);
  mpz_vector_get_mpz_poly(p, c);
  for (unsigned int i = 0; i < c->dim; i++) {
    ASSERT_ALWAYS(mpz_cmp(p->coeff[i], c->c[i]) == 0);
  }
  ASSERT_ALWAYS(p->deg == 1);

  mpz_poly_clear(p);
  mpz_vector_clear(c);
}

void test_submul()
{
  unsigned int dim = 3;

  mpz_vector_t c0;
  mpz_vector_init(c0, dim);
  mpz_vector_t c1;
  mpz_vector_init(c1, dim);

  mpz_vector_setcoordinate_ui(c0, 1, 2);
  mpz_vector_setcoordinate_ui(c0, 2, 5);

  mpz_vector_setcoordinate_ui(c1, 0, 5);

  mpz_t z;
  mpz_init(z);
  mpz_set_str(z, "-17", 10);

  mpz_vector_submul(c0, z, c1);

  ASSERT_ALWAYS(mpz_cmp_ui(c0->c[0], 85) == 0);
  ASSERT_ALWAYS(mpz_cmp_ui(c0->c[1], 2) == 0);
  ASSERT_ALWAYS(mpz_cmp_ui(c0->c[2], 5) == 0);

  mpz_clear(z);
  mpz_vector_clear(c0);
  mpz_vector_clear(c1);
}

int main()
{
  tests_common();
  test_dot_product();
  test_norm();
  test_mpz_poly();
  test_submul();
  exit (EXIT_SUCCESS);
}
