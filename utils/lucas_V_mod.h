#include "macros.h"
#ifndef LUCAS_V_MOD_H_
#define LUCAS_V_MOD_H_

#ifndef mod_init
  #error "One of the mod*_default.h headers must be included before this file"
#endif

/* Function to compute the Lucas sequence V_n (x, 1)
 * The sequence is written V_n (x) in the following.
 * It it defined by:
 *    V_0 (x) = 2
 *    V_1 (x) = x
 *    V_{n+2} (x) = x * V_{n+1} (x) - V_n (x)  for n >= 0
 *    V_{-n} (x) = V_n (x)
 * Properties:
 *    V_{mn} (x) = V_m (V_n (x))
 *    V_{m+n} (x) = V_m (x) * V_n (x) - V_{m-n} (x)
 *    V_n (x+1/x) = x^n + 1/x^n
 *
 * The second properties is used to evaluate V_n (x) with a Montgomery Ladder or
 * the PRAC algorithm.
 */

#ifdef __cplusplus
extern "C" {
#endif

/* Given a = V_n (x), b = V_m (x) and d = V_{n-m} (x), compute V_{m+n} (x).
 * r can be the same variable as a or b but must not be the same variable as d.
 */
static inline void
mod_V_dadd (residue_t r, const residue_t a, const residue_t b,
            const residue_t d, const modulus_t m)
{
  ASSERT (r != d);
  mod_mul (r, a, b, m);
  mod_sub (r, r, d, m);
}

/* Given a = V_n (x) and two = 2, compute V_{2n} (x).
 * r can be the same variable as a but must not be the same variable as two.
 */
static inline void
mod_V_dbl (residue_t r, const residue_t a, const residue_t two,
           const modulus_t m)
{
  ASSERT (r != two);
  mod_sqr (r, a, m);
  mod_sub (r, r, two, m);
}


/* Compute r = V_k (b) and rp1 = V_{k+1} (b) if rp1 != NULL
 * r or rp1 can be the same variable as b.
 */
static void
mod_V_eval_ul (residue_t r, residue_t rp1, const residue_t b,
               const unsigned long k, const modulus_t m)
{
  residue_t t0, t1, two;

  mod_init (t0, m);
  mod_init (t1, m);
  mod_init (two, m);

  mod_set1 (two, m);
  mod_add (two, two, two, m);

  if (k == 0UL)
  {
    mod_set (r, two, m);
    if (rp1)
      mod_set (rp1, b, m);
  }
  else if (k == 1UL)
  {
    mod_set (r, b, m);
    if (rp1)
      mod_V_dbl (rp1, b, two, m);
  }
  else if (k == 2UL)
  {
    if (rp1)
    {
      mod_V_dbl (t1, b, two, m);
      mod_V_dadd (rp1, t1, b, b, m);
      mod_set (r, t1, m);
    }
    else
      mod_V_dbl (r, b, two, m);
  }
  else /* k >= 3 */
  {
    /* Montgomery Ladder */
    unsigned long mask;

    mask = ~(0UL);
    mask -= mask/2;   /* Now the most significant bit of i is set */
    while ((mask & k) == 0)
      mask >>= 1;

    /* Most significant bit of k is 1, do it outside the loop */
    mod_set (t0, b, m);         /* starting value t0 = V_1 (b) = b */
    mod_V_dbl (t1, b, two, m);  /* starting value t1 = V_2 (b) */
    mask >>= 1;

    /* If the second most significant bit of k is 0, then we do the iteration
     * manually (to avoid to compute again V_2 (b))
     * As k >= 3, we know that in this case k has at least 3 bits.
     */
    if (!(k & mask)) /* (t0,t1) <- (V_2 (b), V_3 (b)) */
    {
      mod_set (t0, t1, m);
      mod_V_dadd (t1, t1, b, b, m);
      mask >>= 1;
    }

    for ( ; mask > 1; mask >>= 1) /* t0 = V_j (b) and t1 = V_{j+1} (b) */
    {
      if (k & mask) /* j -> 2*j+1. Compute V_{2j+1} and V_{2j+2} */
      {
        mod_V_dadd (t0, t1, t0, b, m);
        mod_V_dbl (t1, t1, two, m);
      }
      else /* j -> 2*j. Compute V_{2j} and V_{2j+1} */
      {
        mod_V_dadd (t1, t1, t0, b, m);
        mod_V_dbl (t0, t0, two, m);
      }
    }

    /* Deal with least significant bit outside the loop */
    if (k & mask)
    {
      mod_V_dadd (t0, t1, t0, b, m); /* cannot have r instead of t0, if r is the
                                      * same variable as b, the assert in
                                      * mod_V_dadd would fail */
      mod_set (r, t0, m);
      if (rp1)
        mod_V_dbl (rp1, t1, two, m);
    }
    else
    {
      mod_V_dbl (r, t0, two, m);
      if (rp1)
      {
        mod_V_dadd (t1, t1, t0, b, m); /* same as above */
        mod_set (rp1, t1, m);
      }
    }
  }

  mod_clear (t0, m);
  mod_clear (t1, m);
  mod_clear (two, m);
}

#ifdef __cplusplus
}
#endif


#endif /* LUCAS_V_MOD_H_ */
