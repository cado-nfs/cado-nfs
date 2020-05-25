#include "cado.h" // IWYU pragma: keep
#include "mod64.hpp"
#include "modredc64.hpp" // IWYU pragma: keep

typedef Modulus64 Modulus;
#include "mod64_common.cpp"
#include "macros.h"

/* Put 1/s (mod t) in r and return 1 if s is invertible, 
   or set r to 0 and return 0 if s is not invertible */

bool
Modulus64::inv (Residue &r, const Residue &sp) const
{
  int64_t u1, v1;
  uint64_t u2, v2, s, t;
#ifndef NDEBUG
  Residue tmp(*this);
#endif

#ifndef NDEBUG
  /* Remember input in case r overwrites it */
  set (tmp, sp);
#endif

  s = get_u64 (sp);
  t = getmod_u64 ();

  ASSERT (t > 0);
  ASSERT (s < t);

  if (s == 0)
    {
      set0(r); /* Not invertible */
      return 0;
    }

  if (s == 1)
    {
      set1(r);
      return 1;
    }

  u1 = 1;
  u2 = s;
  v1 = - (int64_t) (t / s); /* No overflow, since s >= 2 */
  v2 = t % s;

  if (v2 == 1)
    {
       u1 = v1 + t;
    }
  else 
    {
      while (v2 != 0)
	{
	  uint64_t q;
	  /* unroll twice and swap u/v */
	  q = u2 / v2;
	  ASSERT_EXPENSIVE (q <= (uint64_t) INT64_MAX);
	  u1 = u1 - (int64_t) q * v1;
	  u2 = u2 - q * v2;
	  
	  if (u2 == 0)
	    {
	      u1 = v1;
	      u2 = v2;
	      break;
	    }
	  
	  q = v2 / u2;
	  ASSERT_EXPENSIVE (q <= (uint64_t) INT64_MAX);
	  v1 = v1 - (int64_t) q * u1;
	  v2 = v2 - q * u2;
	}
  
      if (u2 != 1)
	{
	  /* printf ("s=%lu t=%lu found %lu\n", s[0], t[0], u2); */
	  set0(r); /* non-trivial gcd */
	  return 0;
	}

      if (u1 < 0)
        u1 = u1 + t;
    }

  ASSERT ((uint64_t) u1 < t);
    
  set (r, (uint64_t) u1);

#ifndef NDEBUG
  mul (tmp, tmp, r);
  ASSERT(is1 (tmp));
#endif

  return 1;
}

/* even_inv_lookup_table[i] is 1/(2*i+1) mod 128 */
static uint64_t even_inv_lookup_table[64] = {
  1, 43, 77, 55, 57, 35, 69, 111, 113, 27, 61, 39, 41, 19, 53, 95, 97, 11, 45,
  23, 25, 3, 37, 79, 81, 123, 29, 7, 9, 115, 21, 63, 65, 107, 13, 119, 121, 99,
  5, 47, 49, 91, 125, 103, 105, 83, 117, 31, 33, 75, 109, 87, 89, 67, 101, 15,
  17, 59, 93, 71, 73, 51, 85, 127 } ;


/* Faster modul_inv for the case where m = 2^k */
bool
Modulus64::inv_powerof2 (Residue &r, const Residue &A) const
{
  uint64_t x = getmod_u64(), y = get_u64(A);

  ASSERT (!(x & (x-1))); /* assert that x is a power of 2 */
  ASSERT (y < x);
  if (!(y & 1))
    return 0;
  else
  {
    if (!(x >> 4)) /* x = 2, 4 or 8 */
      set(r, y);
    else if (!(x >> 8)) /* x = 16, 32, 64, or 128 */
      set(r, even_inv_lookup_table[(y-1) >> 1] & (x-1));
    else
    {
      uint64_t h = x >> (u64arith_ctz(x) >> 1);
            Modulus64 m2(h);
      Residue B(*this);
      m2.set_reduced (B, (y & (h-1)));

      m2.inv_powerof2 (r, B);
      uint64_t t = get_u64(r);
      set(r,  (2 * t - t*t*y) & (x-1));
    }
    return 1;
  }
}

/* Faster modul_inv for the case where m is odd */
bool
Modulus64::inv_odd (Residue &r, const Residue &A) const
{
    /* TODO: do REDC trick like in ModulusREDC64 */
    return inv(r, A);
}
