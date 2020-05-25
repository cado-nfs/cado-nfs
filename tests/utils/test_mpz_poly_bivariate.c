#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "macros.h"
#include "mpz_poly.h"
#include "mpz_poly_bivariate.h"

void test_mpz_poly_bivariate_trivialities()
{
  mpz_poly tmp;
  mpz_poly_init(tmp, 2);

  mpz_poly_bivariate_t f;
  mpz_poly_bivariate_init(f, 2);

  mpz_poly_setcoeff_int64(tmp, 0, 11);
  mpz_poly_setcoeff_int64(tmp, 1, 17);
  mpz_poly_setcoeff_int64(tmp, 2, 42);
  mpz_poly_bivariate_setcoeff(f, 0, tmp);

  mpz_poly_setcoeff_int64(tmp, 0, 12);
  mpz_poly_setcoeff_int64(tmp, 1, 18);
  mpz_poly_setcoeff_int64(tmp, 2, 43);
  mpz_poly_bivariate_setcoeff(f, 2, tmp);

  //(11+17*x+42*x^2)+(12+18*x+43*x^2)*y^2
  mpz_poly_bivariate_fprintf(stdout, f);

  mpz_poly_bivariate_clear(f);

  mpz_poly_bivariate_init_y_x(f, 2, 3);
  mpz_poly_setcoeff_int64(tmp, 0, 11);
  mpz_poly_setcoeff_int64(tmp, 1, 17);
  mpz_poly_setcoeff_int64(tmp, 2, 42);
  mpz_poly_bivariate_setcoeff(f, 0, tmp);

  mpz_poly_setcoeff_int64(tmp, 0, 12);
  mpz_poly_setcoeff_int64(tmp, 1, 18);
  mpz_poly_setcoeff_int64(tmp, 2, 43);
  mpz_poly_bivariate_setcoeff(f, 2, tmp);

  mpz_poly_setcoeff_int64(tmp, 0, 12);
  mpz_poly_setcoeff_int64(tmp, 1, 18);
  mpz_poly_setcoeff_int64(tmp, 5, 1);
  mpz_poly_bivariate_setcoeff(f, 5, tmp);

  //(11+17*x+42*x^2)+(12+18*x+43*x^2)*y^2+(12+18*x+x^5)*y^5
  ASSERT_ALWAYS(f->deg_y == 5);
  ASSERT_ALWAYS(f->deg_x == 5);
  ASSERT_ALWAYS(f->coeff[1]->deg == -1);
  ASSERT_ALWAYS(f->coeff[5]->deg == 5);
  ASSERT_ALWAYS(mpz_cmp_ui(f->coeff[2]->coeff[1], 18) == 0);

  tmp->deg = -1;
  mpz_poly_setcoeff_int64(tmp, 0, 0);
  mpz_poly_bivariate_setcoeff(f, 5, tmp);

  //(11+17*x+42*x^2)+(12+18*x+43*x^2)*y^2
  ASSERT_ALWAYS(f->deg_y == 2);
  ASSERT_ALWAYS(f->deg_x == 2);
  ASSERT_ALWAYS(f->coeff[1]->deg == -1);
  ASSERT_ALWAYS(f->coeff[0]->deg == 2);
  ASSERT_ALWAYS(mpz_cmp_ui(f->coeff[0]->coeff[2], 42) == 0);

  mpz_t y;
  mpz_init(y);
  mpz_set_si(y, 42);
  //tmp = 75894*x^2 + 31769*x + 21179
  mpz_poly_bivariate_eval_y(tmp, f, y);
  ASSERT_ALWAYS(tmp->deg == 2);
  ASSERT_ALWAYS(mpz_cmp_ui(tmp->coeff[0], 21179) == 0);
  ASSERT_ALWAYS(mpz_cmp_ui(tmp->coeff[1], 31769) == 0);
  ASSERT_ALWAYS(mpz_cmp_ui(tmp->coeff[2], 75894) == 0);

  mpz_set_si(y, 0);
  //tmp = 42*x^2 + 17*x + 11
  mpz_poly_bivariate_eval_y(tmp, f, y);
  ASSERT_ALWAYS(tmp->deg == 2);
  ASSERT_ALWAYS(mpz_cmp_ui(tmp->coeff[0], 11) == 0);
  ASSERT_ALWAYS(mpz_cmp_ui(tmp->coeff[1], 17) == 0);
  ASSERT_ALWAYS(mpz_cmp_ui(tmp->coeff[2], 42) == 0);

  mpz_t x;
  mpz_init(x);
  mpz_set_si(x, 42);
  //tmp = 76620*y^2 + 74813
  mpz_poly_bivariate_eval_x(tmp, f, x);
  ASSERT_ALWAYS(tmp->deg == 2);
  ASSERT_ALWAYS(mpz_cmp_ui(tmp->coeff[0], 74813) == 0);
  ASSERT_ALWAYS(mpz_cmp_ui(tmp->coeff[1], 0) == 0);
  ASSERT_ALWAYS(mpz_cmp_ui(tmp->coeff[2], 76620) == 0);

  mpz_set_si(x, 0);
  //tmp = 12*y^2 + 11
  mpz_poly_bivariate_eval_x(tmp, f, x);
  ASSERT_ALWAYS(tmp->deg == 2);
  ASSERT_ALWAYS(mpz_cmp_ui(tmp->coeff[0], 11) == 0);
  ASSERT_ALWAYS(mpz_cmp_ui(tmp->coeff[1], 0) == 0);
  ASSERT_ALWAYS(mpz_cmp_ui(tmp->coeff[2], 12) == 0);

  f->deg_y = -1;
  mpz_poly_bivariate_eval_y(tmp, f, y);
  ASSERT_ALWAYS(tmp->deg == -1);
  mpz_clear(y);

  mpz_poly_bivariate_eval_x(tmp, f, x);
  ASSERT_ALWAYS(tmp->deg == -1);
  mpz_clear(x);


  mpz_poly_bivariate_clear(f);
  mpz_poly_clear(tmp);
}

void test_mpz_poly_bivariate_resultant()
{
  mpz_poly tmp;
  mpz_poly_init(tmp, 2);

  mpz_poly_bivariate_t f;
  mpz_poly_bivariate_init(f, 2);

  mpz_poly_setcoeff_int64(tmp, 0, 11);
  mpz_poly_setcoeff_int64(tmp, 1, 17);
  mpz_poly_setcoeff_int64(tmp, 2, 42);
  mpz_poly_bivariate_setcoeff(f, 0, tmp);

  mpz_poly_setcoeff_int64(tmp, 0, 12);
  mpz_poly_setcoeff_int64(tmp, 1, 18);
  mpz_poly_setcoeff_int64(tmp, 2, 43);
  mpz_poly_bivariate_setcoeff(f, 2, tmp);

  mpz_poly_bivariate_t g;
  mpz_poly_bivariate_init(g, 2);

  mpz_poly_setcoeff_int64(tmp, 0, 10);
  mpz_poly_setcoeff_int64(tmp, 1, 15);
  mpz_poly_setcoeff_int64(tmp, 2, 100);
  mpz_poly_bivariate_setcoeff(g, 0, tmp);

  mpz_poly_setcoeff_int64(tmp, 0, 20);
  mpz_poly_setcoeff_int64(tmp, 1, 41);
  mpz_poly_setcoeff_int64(tmp, 2, 17);
  mpz_poly_bivariate_setcoeff(g, 2, tmp);

  mpz_poly_bivariate_resultant_y(tmp, f, g);
  ASSERT_ALWAYS(tmp->deg == 8);
  ASSERT_ALWAYS(mpz_cmp_ui(tmp->coeff[0], 10000) == 0);
  ASSERT_ALWAYS(mpz_cmp_ui(tmp->coeff[1], 86200) == 0);
  ASSERT_ALWAYS(mpz_cmp_ui(tmp->coeff[2], 150561) == 0);
  ASSERT_ALWAYS(mpz_cmp_si(tmp->coeff[3], -238512) == 0);
  ASSERT_ALWAYS(mpz_cmp_si(tmp->coeff[4], -1060332) == 0);
  ASSERT_ALWAYS(mpz_cmp_si(tmp->coeff[5], -2938364) == 0);
  ASSERT_ALWAYS(mpz_cmp_ui(tmp->coeff[6], 1450628) == 0);
  ASSERT_ALWAYS(mpz_cmp_ui(tmp->coeff[7], 3112648) == 0);
  ASSERT_ALWAYS(mpz_cmp_ui(tmp->coeff[8], 12859396) == 0);

  mpz_poly_bivariate_resultant_x(tmp, f, g);
  ASSERT_ALWAYS(tmp->deg == 8);
  ASSERT_ALWAYS(mpz_cmp_ui(tmp->coeff[0], 467750) == 0);
  ASSERT_ALWAYS(mpz_cmp_ui(tmp->coeff[2], 38960) == 0);
  ASSERT_ALWAYS(mpz_cmp_si(tmp->coeff[4], -996138) == 0);
  ASSERT_ALWAYS(mpz_cmp_ui(tmp->coeff[6], 44919) == 0);
  ASSERT_ALWAYS(mpz_cmp_ui(tmp->coeff[8], 622660) == 0);

  mpz_poly_setcoeff_int64(tmp, 1, 11);
  mpz_poly_setcoeff_int64(tmp, 3, 17);
  mpz_poly_setcoeff_int64(tmp, 4, 42);
  mpz_poly_bivariate_setcoeff(f, 0, tmp);

  mpz_poly_setcoeff_int64(tmp, 0, 100);
  mpz_poly_setcoeff_int64(tmp, 2, 46);
  mpz_poly_setcoeff_int64(tmp, 6, 43);
  mpz_poly_bivariate_setcoeff(f, 4, tmp);

  mpz_poly_setcoeff_int64(tmp, 0, 12);
  mpz_poly_setcoeff_int64(tmp, 2, 18);
  mpz_poly_setcoeff_int64(tmp, 6, 43);
  mpz_poly_bivariate_setcoeff(f, 5, tmp);

  mpz_poly_setcoeff_int64(tmp, 0, 10);
  mpz_poly_setcoeff_int64(tmp, 1, 15);
  mpz_poly_setcoeff_int64(tmp, 2, 100);
  mpz_poly_bivariate_setcoeff(g, 1, tmp);

  mpz_t Z_tmp;
  mpz_init(Z_tmp);
  mpz_poly_setcoeff_int64(tmp, 0, 20);
  mpz_poly_setcoeff_int64(tmp, 1, 41);
  mpz_poly_setcoeff_int64(tmp, 2, 17);
  mpz_poly_bivariate_setcoeff(g, 4, tmp);

  /*
     f=(467750+11*x+38960*x^2+17*x^3+42*x^4+44919*x^6+622660*x^8)
     +(12+18*x+43*x^2)*y^2
     +(100+11*x+46*x^2+17*x^3+42*x^4+43*x^6+622660*x^8)*y^4
     +(12+11*x+18*x^2+17*x^3+42*x^4+43*x^6+622660*x^8)*y^5

     g=(10+15*x+100*x^2)
     +(10+15*x+100*x^2+17*x^3+42*x^4+43*x^6+622660*x^8)*y
     +(20+41*x+17*x^2)*y^2
     +(20+41*x+17*x^2+17*x^3+42*x^4+43*x^6+622660*x^8)*y^4
   */
  mpz_poly_bivariate_resultant_y(tmp, f, g);
  ASSERT_ALWAYS(tmp->deg == 72);
  mpz_set_str(Z_tmp, "1494774748820383618990632550811343575455342810669737",
      10);
  ASSERT_ALWAYS(mpz_cmp(tmp->coeff[42], Z_tmp) == 0);

  mpz_poly_bivariate_resultant_x(tmp, f, g);
  ASSERT_ALWAYS(tmp->deg == 72);
  mpz_set_str(Z_tmp,
      "392307271360601863563100893262849685478724955452018428897254854297277064613707763200",
      10);
  ASSERT_ALWAYS(mpz_cmp(tmp->coeff[42], Z_tmp) == 0);
  mpz_clear(Z_tmp);

  mpz_poly_bivariate_clear(f);
  mpz_poly_bivariate_clear(g);
  mpz_poly_clear(tmp);
}

int main()
{
  test_mpz_poly_bivariate_trivialities();
  test_mpz_poly_bivariate_resultant();
  exit (EXIT_SUCCESS);
}
