#include "cado.h" // IWYU pragma: keep
#include "utils_mpz.h"
#include <stdint.h>
#include <stdlib.h>
#include <gmp.h>
#include "macros.h"

int mpz_invert_ui(mpz_ptr rop, mpz_srcptr op1, const uint64_t op2)
{
  mpz_t op3;
  mpz_init(op3);
  mpz_set_ui(op3, op2);
  int res = mpz_invert(rop, op1, op3);
  mpz_clear(op3);

  return res;
}
