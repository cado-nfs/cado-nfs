#include "cado.h" // IWYU pragma: keep
#include "mod_ul.h"
#include "mod_ul_default.h"

#include "mod_ul_common.c"

#include "modredc_ul.h"
#include "macros.h"


/* Put 1/s (mod t) in r and return 1 if s is invertible, 
   or set r to 0 and return 0 if s is not invertible */

int
mod_inv (residue_t r, const residue_t sp, const modulus_t m)
{
  long u1, v1;
  unsigned long u2, v2, s, t;
#ifndef NDEBUG
  residue_t tmp;
#endif

#ifndef NDEBUG
  /* Remember input in case r overwrites it */
  mod_init_noset0 (tmp, m);
  mod_set (tmp, sp, m);
#endif

  s = mod_get_ul (sp, m);
  t = mod_getmod_ul (m);

  ASSERT (t > 0UL);
  ASSERT (s < t);

  if (s == 0UL)
    {
      r[0] = 0UL; /* Not invertible */
#ifndef NDEBUG
      mod_clear (tmp, m);
#endif
      return 0;
    }

  if (s == 1UL)
    {
      r[0] = 1UL;
#ifndef NDEBUG
      mod_clear (tmp, m);
#endif
      return 1;
    }

  u1 = 1L;
  u2 = s;
  v1 = - (long) (t / s); /* No overflow, since s >= 2 */
  v2 = t % s;

  if (v2 == 1UL)
    {
       u1 = v1 + t;
    }
  else 
    {
      while (v2 != 0UL)
	{
	  unsigned long q;
	  /* unroll twice and swap u/v */
	  q = u2 / v2;
	  ASSERT_EXPENSIVE (q <= (unsigned long) LONG_MAX);
	  u1 = u1 - (long) q * v1;
	  u2 = u2 - q * v2;
	  
	  if (u2 == 0UL)
	    {
	      u1 = v1;
	      u2 = v2;
	      break;
	    }
	  
	  q = v2 / u2;
	  ASSERT_EXPENSIVE (q <= (unsigned long) LONG_MAX);
	  v1 = v1 - (long) q * u1;
	  v2 = v2 - q * u2;
	}
  
      if (u2 != 1UL)
	{
	  /* printf ("s=%lu t=%lu found %lu\n", s[0], t[0], u2); */
	  r[0] = 0UL; /* non-trivial gcd */
#ifndef NDEBUG
          mod_clear (tmp, m);
#endif
	  return 0;
	}

      if (u1 < 0L)
        u1 = u1 + t;
    }

  ASSERT ((unsigned long) u1 < t);
    
  mod_set_ul (r, (unsigned long) u1, m);

#ifndef NDEBUG
  mod_mul (tmp, tmp, r, m);
  ASSERT(mod_is1 (tmp, m));
  mod_clear (tmp, m);
#endif

  return 1;
}

/* even_inv_lookup_table[i] is 1/(2*i+1) mod 128 */
static unsigned long even_inv_lookup_table[64] = {
  1, 43, 77, 55, 57, 35, 69, 111, 113, 27, 61, 39, 41, 19, 53, 95, 97, 11, 45,
  23, 25, 3, 37, 79, 81, 123, 29, 7, 9, 115, 21, 63, 65, 107, 13, 119, 121, 99,
  5, 47, 49, 91, 125, 103, 105, 83, 117, 31, 33, 75, 109, 87, 89, 67, 101, 15,
  17, 59, 93, 71, 73, 51, 85, 127 } ;


/* Faster modul_inv for the case where m = 2^k */
int
modul_inv_powerof2 (residue_t r, const residue_t A, const modulus_t m)
{
  unsigned long x = m[0], y = A[0];

  ASSERT (!(x & (x-1))); /* assert that x is a power of 2 */
  ASSERT (y < x);
  if (!(y & 1UL))
    return 0;
  else
  {
    if (!(x >> 4)) /* x = 2, 4 or 8 */
      r[0] = y;
    else if (!(x >> 8)) /* x = 16, 32, 64, or 128 */
      r[0] = even_inv_lookup_table[(y-1) >> 1] & (x-1);
    else
    {
      modulusul_t m2;
      residueul_t B;
      unsigned long h = x >> (ularith_ctz(x) >> 1);
      modul_initmod_ul (m2, h);
      modul_init_noset0 (B, m2);
      modul_set_ul_reduced (B, (y & (h-1)), m2);

      modul_inv_powerof2 (r, B, m2);
      unsigned long t = (r[0] * r[0]) & (x-1);
      t = (t * y) & (x-1);
      r[0] = (r[0] << 1) & (x-1);
      r[0] = (r[0] - t) & (x-1);

      modul_clear (B, m2);
      modul_clearmod (m2);
    }
    return 1;
  }
}

/* Faster modul_inv for the case where m is odd */
int
modul_inv_odd (residue_t r, const residue_t A, const modulus_t m)
{
  modulusredcul_t mm;
  modredcul_initmod_ul_raw (mm, m[0]);
  int ret = modredcul_intinv (r, A, mm);
  modredcul_clearmod(mm);
  return ret;
}
