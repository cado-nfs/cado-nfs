#include "cado.h"
#include "modredc_ul.h"
#include "modredc_ul_default.h"
#include "mod_ul_common.c"

int
modredcul_inv (residueredcul_t r, const residueredcul_t A,
               const modulusredcul_t m)
{
  unsigned long x = m[0].m, y, u, v;
  int t, lsh;

  ASSERT (A[0] < x);
  ASSERT (x & 1UL);

  if (A[0] == 0UL)
    return 0;

  /* Let A = a*2^w, so we want the Montgomery representation of 1/a,
     which is 2^w/a. We start by getting y = a */
  y = modredcul_get_ul (A, m);

  /* We simply set y = a/2^w and t=0. The result before
     correction will be 2^(w+t)/a so we have to divide by t, which
     may be >64, so we may have to do a full and a variable width REDC. */
  y = modredcul_get_ul (&y, m);
  /* Now y = a/2^w */
  t = 0;

  u = 1UL; v = 0UL;

  // make y odd
  lsh = ularith_ctz(y);
  y >>= lsh;
  t += lsh;
  /* v <<= lsh; ??? v is 0 here */

  // Here y and x are odd, and y < x
  do {
    /* Here, y and x are odd, 0 < y < x, u is odd and v is even */
    do {
      x -= y; v += u;
      if (x == 0)
        break;
      lsh = ularith_ctz(x);
      ASSERT_EXPENSIVE (lsh > 0);
      x >>= lsh;
      t += lsh;
      u <<= lsh;
    } while (x > y); /* ~50% branch taken :( */
    /* Here, y and x are odd, 0 < x =< y, u is even and v is odd */

    /* x is the one that got reduced, test if we're done */

    if (x <= 1)
      break;

    /* Here, y and x are odd, 0 < x < y, u is even and v is odd */
    do {
      y -= x; u += v;
      if (y == 0)
        break;
      lsh = ularith_ctz(y);
      ASSERT_EXPENSIVE (lsh > 0);
      y >>= lsh;
      t += lsh;
      v <<= lsh;
    } while (x < y); /* about 50% branch taken :( */
    /* Here, y and x are odd, 0 < y =< x, u is odd and v is even */
    /* y is the one that got reduced, test if we're done */
  } while (y > 1);

  if ((x & y) == 0UL) /* Non-trivial GCD */
    return 0;

  if (y != 1)
    {
      /* So x is the one that reached 1.
	 We maintained ya == u2^t (mod m) and xa = -v2^t (mod m).
	 So 1/a = -v2^t.
       */
      u = m[0].m - v;
      /* Now 1/a = u2^t */
    }

  /* Here, u = 2^w * 2^t / a. We want 2^w / a. */

  /* Here, the inverse of y is u/2^t mod x. To do the division by t,
     we use a variable-width REDC. We want to add a multiple of m to u
     so that the low t bits of the sum are 0 and we can right-shift by t. */
  if (t >= LONG_BIT)
    {
      unsigned long tlow, thigh;
      tlow = u * m[0].invm; /* tlow <= 2^w-1 */
      ularith_mul_ul_ul_2ul (&tlow, &thigh, tlow, m[0].m); /* thigh:tlow <= (2^w-1)*m */
      u = thigh + ((u != 0UL) ? 1UL : 0UL);
      /* thigh:tlow + u < (2^w-1)*m + m < 2^w*m. No correction necesary */
      t -= LONG_BIT;
    }

  ASSERT (t < LONG_BIT);
  if (t > 0)
    {
      unsigned long tlow, thigh;
      /* Necessarily t < LONG_BIT, so the shift is ok */
      /* Doing a left shift first and then a full REDC needs a modular addition
	 at the end due to larger summands and thus is probably slower */
      tlow = ((u * m[0].invm) & ((1UL << t) - 1UL)); /* tlow <= 2^t-1 */
      ularith_mul_ul_ul_2ul (&tlow, &thigh, tlow, m[0].m); /* thigh:tlow <= m*(2^t-1) */
      ularith_add_ul_2ul (&tlow, &thigh, u); /* thigh:tlow <= m*2^t-1 (since u<m) */
      /* Now the low t bits of tlow are 0 */
      ASSERT_EXPENSIVE ((tlow & ((1UL << t) - 1UL)) == 0UL);
      ularith_shrd (&tlow, thigh, t);
      u = tlow;
      ASSERT_EXPENSIVE ((thigh >> t) == 0UL && u < m[0].m);
    }

  r[0] = u;
  return 1;
}

/* same as modredcul_inv, but for classical representation (not Montgomery) */
int
modredcul_intinv (residueredcul_t r, const residueredcul_t A,
               const modulusredcul_t m)
{
  unsigned long x = m[0].m, y, u, v;
  int t, lsh;

  ASSERT (A[0] < x);
  ASSERT (x & 1UL);

  if (A[0] == 0UL)
    return 0;

  y = A[0];
  t = 0;

  u = 1UL; v = 0UL;

  // make y odd
  lsh = ularith_ctz(y);
  y >>= lsh;
  t += lsh;
  /* v <<= lsh; ??? v is 0 here */

  // Here y and x are odd, and y < x
  do {
    /* Here, y and x are odd, 0 < y < x, u is odd and v is even */
    do {
      x -= y; v += u;
      lsh = ularith_ctz(x);
      ASSERT_EXPENSIVE (lsh > 0);
      x >>= lsh;
      t += lsh;
      u <<= lsh;
    } while (x > y); /* ~50% branch taken :( */
    /* Here, y and x are odd, 0 < x =< y, u is even and v is odd */

    /* x is the one that got reduced, test if we're done */

    if (x <= 1)
      break;

    /* Here, y and x are odd, 0 < x < y, u is even and v is odd */
    do {
      y -= x; u += v;
      lsh = ularith_ctz(y);
      ASSERT_EXPENSIVE (lsh > 0);
      y >>= lsh;
      t += lsh;
      v <<= lsh;
    } while (x < y); /* about 50% branch taken :( */
    /* Here, y and x are odd, 0 < y =< x, u is odd and v is even */
    /* y is the one that got reduced, test if we're done */
  } while (y > 1);

  if ((x & y) == 0UL) /* Non-trivial GCD */
    return 0;

  if (y != 1)
    {
      /* So x is the one that reached 1.
	 We maintained ya == u2^t (mod m) and xa = -v2^t (mod m).
	 So 1/a = -v2^t.
       */
      u = m[0].m - v;
      /* Now 1/a = u2^t */
    }

  /* Here, u = 2^w * 2^t / a. We want 2^w / a. */

  /* Here, the inverse of y is u/2^t mod x. To do the division by t,
     we use a variable-width REDC. We want to add a multiple of m to u
     so that the low t bits of the sum are 0 and we can right-shift by t. */
  if (t >= LONG_BIT)
    {
      unsigned long tlow, thigh;
      tlow = u * m[0].invm; /* tlow <= 2^w-1 */
      ularith_mul_ul_ul_2ul (&tlow, &thigh, tlow, m[0].m); /* thigh:tlow <= (2^w-1)*m */
      u = thigh + ((u != 0UL) ? 1UL : 0UL);
      /* thigh:tlow + u < (2^w-1)*m + m < 2^w*m. No correction necesary */
      t -= LONG_BIT;
    }

  ASSERT (t < LONG_BIT);
  if (t > 0)
    {
      unsigned long tlow, thigh;
      /* Necessarily t < LONG_BIT, so the shift is ok */
      /* Doing a left shift first and then a full REDC needs a modular addition
	 at the end due to larger summands and thus is probably slower */
      tlow = ((u * m[0].invm) & ((1UL << t) - 1UL)); /* tlow <= 2^t-1 */
      ularith_mul_ul_ul_2ul (&tlow, &thigh, tlow, m[0].m); /* thigh:tlow <= m*(2^t-1) */
      ularith_add_ul_2ul (&tlow, &thigh, u); /* thigh:tlow <= m*2^t-1 (since u<m) */
      /* Now the low t bits of tlow are 0 */
      ASSERT_EXPENSIVE ((tlow & ((1UL << t) - 1UL)) == 0UL);
      ularith_shrd (&tlow, thigh, t);
      u = tlow;
      ASSERT_EXPENSIVE ((thigh >> t) == 0UL && u < m[0].m);
    }

  r[0] = u;
  return 1;
}
