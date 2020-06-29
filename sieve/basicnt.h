#ifndef BASICNT_H_
#define BASICNT_H_


/*****************************************************************
 *       Some basic number theory functions for inlining         *
 *****************************************************************/

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

static inline unsigned long iscomposite (const unsigned long n);
static inline uint32_t signed_mod_longto32 (long a, uint32_t p);
#ifdef __cplusplus
}
#endif

/* Returns 0 if n is prime, otherwise the smallest prime factor of n */
static inline unsigned long iscomposite (const unsigned long n)
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

static inline uint32_t signed_mod_longto32 (long a, uint32_t p)
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


#endif	/* BASICNT_H_ */
