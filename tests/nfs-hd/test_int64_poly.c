#include <gmp.h>
#include "cado.h" // IWYU pragma: keep
#include "int64_poly.h"
#include "mpz_poly.h"
#include "macros.h"

int main()
{
  int64_poly_t pol0;
  int64_poly_init(pol0, 2);

  for (int i = 0; i < pol0->alloc; i++) {
    ASSERT_ALWAYS(int64_poly_getcoeff(i, pol0) == 0);
  }

  int64_poly_setcoeff(pol0, 2, 1);
  int64_poly_setcoeff(pol0, 6, 2);

  for (int i = 3; i < pol0->deg; i++) {
    ASSERT_ALWAYS(int64_poly_getcoeff(i, pol0) == 0);
  }

  int64_poly_t pol1;
  int64_poly_init(pol1, 0);
  int64_poly_copy(pol1, pol0);

  ASSERT_ALWAYS(pol0->deg == pol1->deg);
  for (int i = 0; i < pol0->deg; i++) {
    ASSERT_ALWAYS(int64_poly_getcoeff(i, pol0) == int64_poly_getcoeff(i, pol1));
  }
  ASSERT_ALWAYS(int64_poly_equal(pol0, pol1));

  int64_poly_set_xi(pol0, 4);
  for (int i = 0; i < pol0->deg; i++) {
    ASSERT_ALWAYS(int64_poly_getcoeff(i, pol0) == 0);
  }
  ASSERT_ALWAYS(int64_poly_getcoeff(pol0->deg, pol0) == 1);

  int64_poly_set_bxi(pol0, 5, 4);
  for (int i = 0; i < pol0->deg; i++) {
    ASSERT_ALWAYS(int64_poly_getcoeff(i, pol0) == 0);
  }
  ASSERT_ALWAYS(int64_poly_getcoeff(pol0->deg, pol0) == 4);
  ASSERT_ALWAYS(!int64_poly_equal(pol0, pol1));

  ASSERT_ALWAYS((uint64_t) int64_poly_max(pol1) ==
      int64_poly_infinity_norm(pol1));

  mpz_poly pol2;
  mpz_poly_init(pol2, 0);
  int64_poly_to_mpz_poly(pol2, pol1);

  ASSERT_ALWAYS(pol2->deg == pol1->deg);
  mpz_t tmp;
  mpz_init(tmp);
  for (int i = 0; i <= pol2->deg; i++) {
    mpz_poly_getcoeff(tmp, i, pol2);
    ASSERT_ALWAYS(0 == mpz_cmp_ui(tmp, int64_poly_getcoeff(i, pol1)));
  }
  mpz_clear(tmp);

  int64_poly_clear(pol0);
  int64_poly_clear(pol1);
  mpz_poly_clear(pol2);
  
  return 0;
}
