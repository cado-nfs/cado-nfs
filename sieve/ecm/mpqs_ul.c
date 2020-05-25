#include "cado.h" // IWYU pragma: keep
#include "mod_ul_default.h"
#include <gmp.h>
#include "mpqs.h"
#include "mpqs_doit.h"

int
mpqs_ul (modint_t f, const modulus_t m)
{
  mpz_t fz, mz;

  mpz_init (fz);
  mpz_init_set_ui (mz, m[0]);
  mpqs_doit (fz, mz, 0);
  f[0] = mpz_get_ui (fz);
  mpz_clear (fz);
  mpz_clear (mz);

  return 0; /* no back-tracking */
}
