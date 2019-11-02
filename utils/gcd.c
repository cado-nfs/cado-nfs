#include "cado.h"
#include <stdint.h>
#include "gcd.h"
#include "mod_ul.h"
#include "macros.h"
#include "misc.h" /* for cado_ctzl */

uint64_t
gcd_int64 (const int64_t a, const int64_t b)
{
  return gcd_uint64(safe_abs64(a), safe_abs64(b));
}

uint64_t
gcd_uint64 (uint64_t a, uint64_t b)
{
  uint64_t t;

  if (b == 0)
    return a;

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

unsigned long
gcd_ul (unsigned long a, unsigned long b)
{
  unsigned long t;

  if (b == 0)
    return a;

  if (a >= b)
    a %= b;

  while (a > 0)
    {
      /* Here 0 < a < b */
      t = b % a;
      b = a;
      a = t;
    }

  return b;
}

/* return 1/a mod p, assume 0 <= a < p */
unsigned long
invert_ul (unsigned long a, unsigned long p)
{
  modulusul_t q;
  residueul_t b;

  modul_initmod_ul (q, p);
  modul_init (b, q);
  assert (a < p);
  modul_set_ul_reduced (b, a, q);
  modul_inv (b, b, q);
  a = modul_get_ul (b, q);
  modul_clear (b, q);
  modul_clearmod (q);
  return a;
}

/* Binary gcd; any input allowed. */
uint64_t
bin_gcd_int64_safe (const int64_t a, const int64_t b)
{
  uint64_t ua = safe_abs64(a), ub = safe_abs64(b);
  int s, t;

  if (ua == 0)
    return ub;

  if (ub == 0)
    return ua;

  s = cado_ctz64 (ua);
  t = cado_ctz64 (ub);
  ua >>= s;
  ub >>= t;
  if (t < s)
    s = t;
  /* Here ua, ub > 0 and both odd */

  while (1)
    {
      while (ua >= ub) {
        /* Here, ua >= ub > 0, and ua, ub are both odd */
        ua -= ub;
        if (ua == 0)
          return ub << s;
        ua >>= cado_ctz64 (ua);
      }
      while (ub >= ua) {
       /* Here, ub >= ua > 0, and ua, ub are both odd */
       ub -= ua;
       if (ub == 0)
         return ua << s;
       ub >>= cado_ctz64 (ub);
      }
    }
}

