#include "cado.h"
#include "utils.h"
#include "macros.h"
#include "alpha3d.h"
#include "mpz_poly.h"

//From test_mpz_poly
static void mpz_poly_setcoeffs_si_var(mpz_poly f, int d, ...)
{
    va_list ap;
    va_start(ap, d);
    mpz_poly_realloc(f, d + 1);
    for(int i = 0 ; i <= d ; i++) {
        mpz_set_si(f->coeff[i], va_arg(ap, int));
    }
    mpz_poly_cleandeg(f, d);
    va_end(ap);
}

int
main (int argc, char *argv[])
{
  mpz_poly f;
  gmp_randstate_t rstate;
  unsigned int N = 10000;

  gmp_randinit_default(rstate);

  if (argc >= 2)
    N = atoi (argv[1]);

  unsigned long seed = 0;
  if (argc >= 3)
    seed = atoi (argv[2]);
  gmp_randseed_ui(rstate, seed);

  mpz_poly_init(f, 6);

  mpz_poly_setcoeffs_si_var(f, 6, 1, -91348, -228385, -20, 228370, 91354, 1);
  double alpha = alpha3d(f, 2000, rstate, N);
  printf ("alpha = %f\n", alpha);
  ASSERT_ALWAYS(-2.0 < alpha && alpha < -1.6);

  mpz_poly_setcoeffs_si_var(f, 6, 23667000, 135452818, -16372955, -473340000,
      -338632045, 6549182, 23667000);
  alpha = alpha3d(f, 2000, rstate, N);
  printf ("alpha = %f\n", alpha);
  ASSERT_ALWAYS(-11.7 < alpha && alpha < -11.3);

  gmp_randclear(rstate);
  mpz_poly_clear(f);

  return 0;
}
