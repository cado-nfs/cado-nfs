/*****************************************************************
 *       Some basic number theory functions for inlining         *
 *****************************************************************/

#include <stdint.h>
#include "utils/misc.h" /* for ctzl */

static inline unsigned long
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

static inline unsigned long 
eulerphi_ul (unsigned long n)
{
  unsigned long p, r = 1UL;
  
  if (n == 0UL) /* Undefined, we return 0 */
    return 0UL;

  if (n % 2UL == 0UL)
    {
      n /= 2UL;
      while (n % 2UL == 0UL)
        {
          n /= 2UL;
          r *= 2UL;
        }
    }

  for (p = 3UL; p*p <= n; p += 2UL)
    {
      if (n % p == 0UL)
        {
          n /= p;
          r *= p - 1UL;
          while (n % p == 0UL)
            {
              n /= p;
              r *= p;
            }
        }
    }
  /* Now n is either 1 or a prime */
  if (n > 1UL)
    r *= n - 1UL;
  
  return r;
}

/* Returns 0 if n is prime, otherwise the smallest prime factor of n */
static inline unsigned long
iscomposite (const unsigned long n)
{
  unsigned long i, i2;

  if (n % 2 == 0)
    return (n == 2) ? 0 : 2;

  /* (i + 2)^2 = i^2 + 4*i + 4 */
  for (i = 3, i2 = 9; i2 <= n; i2 += (i+1) * 4, i += 2)
    if (n % i == 0)
	return i;

  return 0;
}

static inline uint32_t
signed_mod_longto32 (long a, uint32_t p)
{
  uint32_t amodp;
  if (a < 0)
    {
      amodp = ((unsigned long)(-a)) % p;
      if (amodp > 0)
	amodp = p - amodp;
    }
  else
    amodp = ((unsigned long) a) % p;

  return amodp;
}


#if 0
/* Square root. Returns the largest integer x so that x^2 <= n.
   If e != NULL, stores n - x^2 in *e. Should be compiled with 
   -funroll-loops for best performance. */

static inline unsigned long  
sqrtint_ul (const unsigned long n, unsigned long *e)
{
  int i;   
  unsigned long xs, c, d, s2;

  d = n; /* d = n - x^2 */
  xs = 0UL;
  s2 = 1UL << (sizeof (unsigned long) * 8 - 2);

  for (i = sizeof (unsigned long) * 4 - 1; i != 0; i--)
    {
      /* Here, s2 = 1 << (2*i) */
      /* xs = x << (i + 1), the value of x shifted left i+1 bits */

      c = xs + s2; /* c = (x + 2^i) ^ 2 - x^2 = 2^(i+1) * x + 2^(2*i) */
      xs >>= 1; /* Now xs is shifted only i positions */
      if (d >= c)
        {
          d -= c;
          xs |= s2; /* x |= 1UL << i <=> xs |= 1UL << (2*i) */
        }
      s2 >>= 2;
    }

  c = xs + s2; 
  xs >>= 1;
  if (d >= c)
    {
      d -= c;   
      xs |= s2;
    }
 
  if (e != NULL)
    *e = d;
  return xs;
}
#endif

/* Binary gcd; a can be odd or even, b can be zero. */
static inline int64_t
bin_gcd (int64_t a, int64_t b)
{
  ASSERT (b == 0 || b % 2 == 1);
  while (b != 0)
    {
      /* if a is odd, reduce a wrt b, i.e., cancel the two low bits of a,
         so that a + q*b = 0 (mod 4) */
      b >>= ctzl (b);
      a = ((a ^ b) & 2) ? a + b : a - b;
      /* if a was even, then since b is now odd, the new a is odd */
      if (a == 0)
        return (b > 0) ? b : -b;
      a >>= ctzl (a);
      /* from here on, a and b are odd (or zero) */
      ASSERT(a & 1);
      /* reduce b wrt a */
      b = ((b ^ a) & 2) ? b + a : b - a;
    }
  return (a > 0) ? a : -a;
}

/* Binary gcd; any input allowed. */
static inline int64_t
bin_gcd_safe (int64_t a, int64_t b)
{
  int s, t;

  if (a == 0)
    return b;

  if (b == 0)
    return a;

  s = ctzl (a);
  t = ctzl (b);
  a >>= s;
  b >>= t;
  if (t < s)
    s = t;
  
  while (b != 0)
    {
      /* if a is odd, reduce a wrt b, i.e., cancel the two low bits of a,
         so that a + q*b = 0 (mod 4) */
      b >>= ctzl (b);
      a = ((a ^ b) & 2) ? a + b : a - b;
      /* if a was even, then since b is now odd, the new a is odd */
      if (a == 0)
        return (b > 0) ? b : -b;
      a >>= ctzl (a);
      /* from here on, a and b are odd (or zero) */
      ASSERT(a & 1);
      /* reduce b wrt a */
      b = ((b ^ a) & 2) ? b + a : b - a;
    }
  return (a > 0) ? (a << s) : ((-a) << s);
}
