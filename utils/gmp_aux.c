/* auxiliary routines on GMP data-types */

#include "cado.h"
#include <stdint.h>
#include "gmp_aux.h"
#include "macros.h"

/* Set z to q. Warning: on 32-bit machines, we cannot use mpz_set_ui! */
void
mpz_set_uint64 (mpz_t z, uint64_t q)
{
  if (sizeof (unsigned long int) == 8)
    mpz_set_ui (z, (unsigned long int) q);
  else
    {
      ASSERT_ALWAYS (sizeof (unsigned long int) == 4);
      mpz_set_ui (z, (unsigned long int) (q >> 32));
      mpz_mul_2exp (z, z, 32);
      mpz_add_ui (z, z, (unsigned long int) (q & 4294967295UL));
    }
}

/* this functions discards z */
uint64_t
mpz_get_uint64 (mpz_t z)
{
  uint64_t q;

  if (sizeof (unsigned long int) == 8)
    q = mpz_get_ui (z);
  else
    {
      ASSERT_ALWAYS (sizeof (unsigned long int) == 4);
      q = mpz_get_ui (z); /* get the low word of z */
      mpz_div_2exp (z, z, 32);
      q += (uint64_t) mpz_get_ui (z) << 32;
    }
  return q;
}

/* returns the smallest prime > q */
uint64_t
uint64_nextprime (uint64_t q)
{
  mpz_t p;

  mpz_init (p);
  mpz_set_uint64 (p, q);
  mpz_nextprime (p, p);
  q = mpz_get_uint64 (p);
  mpz_clear (p);
  return q;
}

/* same as above, for an unsigned long */
unsigned long
ulong_nextprime (unsigned long q)
{
  mpz_t p;

  mpz_init (p);
  mpz_set_ui (p, q);
  mpz_nextprime (p, p);
  ASSERT_ALWAYS (mpz_fits_ulong_p (p));
  q = mpz_get_ui (p);
  mpz_clear (p);
  return q;
}

#define REPS 1 /* number of Miller-Rabin tests in isprime */

int
isprime (unsigned long p)
{
  mpz_t P;
  int res;
  
  mpz_init_set_ui (P, p);
  res = mpz_probab_prime_p (P, REPS);
  mpz_clear (P);
  return res;
}

/* return the number of bits of p, counting from the least significant end */
int nbits (uintmax_t p)
{
  int k;

  for (k = 0; p != 0; p >>= 1, k ++);
  return k;
}

