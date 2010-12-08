#include "cado.h"
#include <stdint.h>
#include "gcd_uint64.h"
#include "macros.h"

uint64_t
gcd_uint64 (uint64_t a, uint64_t b)
{
  uint64_t t;

  ASSERT (b != 0);

  if (a >= b)
    a %= b;

  while (a > 0)
    {
      /* Here 0 < a < b */
      t = b % a; /* 0 <= t < a */
      b = a;
      a = t; /* 0 <= a < b */
    }

  return b;
}
