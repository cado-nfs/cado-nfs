#include "cado.h" // IWYU pragma: keep
#include <gmp.h>
#include "mpqs_doit.h"

int
mpqs_mpz (mpz_t f, const mpz_t m)
{
  mpqs_doit (f, m, 0);

  return 0; /* no back-tracking */
}
