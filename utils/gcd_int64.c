#include "cado.h"
#include <stdint.h>
#include "gcd_int64.h"
#include "macros.h"

int64_t
gcd_int64 (int64_t a, int64_t b)
{
  int64_t t;

  ASSERT (b != 0);

  if (a < 0)
    a = -a;
  if (b < 0)
    b = -b;
  
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
