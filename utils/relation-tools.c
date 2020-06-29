#include "cado.h" // IWYU pragma: keep
#include "typedefs.h"   // p_r_values_t
#include "mod_ul.h"     // modulusul_t
#include "relation-tools.h"

/*  a few conversion relations. We happen to export these */
static inline void
sswap(char *pd, char *pe)
{
  char c;
  while (pe > pd)
    {
      c = *pe; *pe-- = *pd; *pd++ = c;
    }
}

char *
u64toa16 (char *p, uint64_t m)
{
  char *op;
  static char v[] = "0123456789abcdef";

  op = p;
  do
    {
      *p++ = v[m & 0xf];
      m >>= 4;
    }
  while (m);
  sswap(op, &(p[-1]));
  return (p);
}

char *
u64toa10 (char *p, uint64_t m)
{
  char *op;
  uint64_t n;

  op = p;
  do
    {
      n = m / 10;
      *p++ = (unsigned char) ((m - n * 10) + '0');
      m = n;
    }
  while (m);
  sswap(op, &(p[-1]));
  return (p);
}

char *
d64toa10 (char *p, int64_t m)
{
  if (m < 0)
    {
      *p++ = '-';
      m = -m;
    }
  return (u64toa10 (p, (uint64_t) m));
}

char *
d64toa16 (char *p, int64_t m)
{
  if (m < 0)
    {
      *p++ = '-';
      m = -m;
    }
  return (u64toa16 (p, (uint64_t) m));
}

/* return a/b mod p, and p when gcd(b,p) <> 1: this corresponds to a
   projective root */
/* We use the fact that sizeof(p_r_values_t) <= sizeof(unsigned long);
   this is a condition that is checked in typedefs.h.
   So we can use mod_initmod_ul function and mod_get_ul without risk. */
p_r_values_t
relation_compute_r (int64_t a, uint64_t b, p_r_values_t p)
{
  int inv;
  unsigned long root;
  modulusul_t m;
  residueul_t r, t, pa, br;

  modul_initmod_ul (m, p);
  modul_init (t, m);
  modul_init (br, m);

  modul_set_uint64 (br, b, m); /* Does reduction mod p */
  if (p & 1UL)
    inv = modul_inv_odd(t, br, m);
  else
    inv = modul_inv_powerof2(t, br, m);
  if (inv) /* if inv = 1 then t = 1/b mod p */
  {
    modul_init (pa, m);
    modul_init (r, m);

    modul_set_int64 (pa, a, m); /* Does reduction mod p */

    modul_mul(r, t, pa, m);
    root = modul_get_ul (r, m);

    modul_clear (pa, m); /* No-ops. Here for the sake of pedantry */
    modul_clear (r, m);
  }
  else /* if inv = 0 then p divides b */
    root = p;
  
  modul_clear (t, m); /* No-ops. Here for the sake of pedantry */
  modul_clear (br, m);
  modul_clearmod (m);
  return root;
}


