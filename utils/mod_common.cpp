#include "cado.h" // IWYU pragma: keep
#include "u64arith.h"

/* Some methods that are implemented entirely via existing class methods,
 * and which thus can share source code for all arithmetic implentations */

#ifndef MOD_NO_SHARED_MOD_POW_UL
/* Compute r = b^e. Here, e is a uint64_t */
void
Modulus::pow (Residue &r, const Residue &b, const uint64_t e) const
{
  uint64_t mask;
  Residue t(*this);
  
  if (e == 0)
    {
      set1 (r);
      return;
    }

  /* Find highest set bit in e. */
  mask = (UINT64_C(1) << 63) >> u64arith_clz (e);

  set (t, b);

  while ((mask >>= 1) > 0) {
      sqr (t, t);
      if (e & mask) {
        mul (t, t, b);
      }
  }
  /* Now e = 0, mask = 1, and r^mask * b^0 = r^mask is the result we want */
  set (r, t);
}
#endif /* MOD_NO_SHARED_MOD_POW_UL */


#ifndef MOD_NO_SHARED_MOD_POW_MP
/* Compute r = b^e. Here e is a multiple precision integer 
   sum_{i=0}^{e_nrwords-1} e[i] * (machine word base)^i.
   Assume e[e_nrwords-1] is not zero when e_nrwords > 0.
*/
void
Modulus::pow (Residue &r, const Residue &b, const uint64_t *e, 
	    const size_t e_nrwords) const
{
  uint64_t mask;
  Residue t(*this);
  int i = e_nrwords;

  while (i > 0 && e[i - 1] == 0)
      i--;
  
  if (i == 0) {
      set1(r);
      return;
  }
  i--;

  /* Find highest set bit in e[i]. */
  mask = (UINT64_C(1) << 63) >> u64arith_clz (e[i]);
  /* t = 1, so t^(mask/2) * b^e = t^mask * b^e  */

  /* Exponentiate */

  set (t, b);       /* (r*b)^mask * b^(e-mask) = r^mask * b^e */
  mask >>= 1;

    for ( ; i >= 0; i--) {
        while (mask > 0) {
            sqr (t, t);
            if (e[i] & mask)
                mul (t, t, b);
            mask >>= 1;            /* (r^2)^(mask/2) * b^e = r^mask * b^e */
        }
        mask = UINT64_C(1) << 63;
    }
    set (r, t);
}
#endif /* MOD_NO_SHARED_MOD_POW_MP */

#ifndef MOD_NO_SHARED_MOD_POW_INT
void
Modulus::pow (Residue &r, const Residue &b, const Integer &e) const
{
  size_t i = e.getWordCount();
  const int bits = e.getWordSize();
  const Integer::WordType msb = (Integer::WordType) 1 << (bits-1);
  Integer::WordType word, mask;
  Residue t(*this);

  while (i > 0 && e.getWord(i - 1) == 0)
      i--;
  
  if (i == 0) {
      set1(r);
      return;
  }

  word = e.getWord(i - 1);
  mask = msb >> u64arith_clz (word);

  /* Exponentiate */

  set (t, b);
  mask >>= 1;

    for ( ; i > 0; i--) {
        word = e.getWord(i - 1);
        while (mask > 0) {
            sqr (t, t);
            if (word & mask)
                mul (t, t, b);
            mask >>= 1;
        }
        mask = msb;
    }
    set (r, t);
}
#endif

/* Returns 1 if r1 == 1 (mod m) or if r1 == -1 (mod m) or if
   one of r1^(2^1), r1^(2^2), ..., r1^(2^(po2-1)) == -1 (mod m),
   zero otherwise. Requires -1 (mod m) in minusone. */

bool
Modulus::find_minus1 (Residue &r1, const Residue &minusone, const int po2) const
{
    int i;

    if (is1 (r1) || equal (r1, minusone))
        return true;

    for (i = 1 ; i < po2; i++) {
        sqr (r1, r1);
        if (equal (r1, minusone))
            break;
    }

    return i < po2;
}

/* Compute r = V_k (b) and rp1 = V_{k+1} (b) if rp1 != NULL
 * r or rp1 can be the same variable as b.
 */
void
Modulus::V (Residue &r, Residue *rp1, const Residue &b,
            const uint64_t k) const
{
  Residue t0(*this), t1(*this), two(*this);

  set1 (two);
  add1 (two, two);

  if (k == 0UL)
  {
    set (r, two);
    if (rp1)
      set (*rp1, b);
  }
  else if (k == 1UL)
  {
    set (r, b);
    if (rp1)
      V_dbl (*rp1, b, two);
  }
  else if (k == 2UL)
  {
    if (rp1)
    {
      V_dbl (t1, b, two);
      V_dadd (*rp1, t1, b, b);
      set (r, t1);
    }
    else
      V_dbl (r, b, two);
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
    set (t0, b);         /* starting value t0 = V_1 (b) = b */
    V_dbl (t1, b, two);  /* starting value t1 = V_2 (b) */
    mask >>= 1;

    /* If the second most significant bit of k is 0, then we do the iteration
     * manually (to avoid to compute again V_2 (b))
     * As k >= 3, we know that in this case k has at least 3 bits.
     */
    if (!(k & mask)) /* (t0,t1) <- (V_2 (b), V_3 (b)) */
    {
      set (t0, t1);
      V_dadd (t1, t1, b, b);
      mask >>= 1;
    }

    for ( ; mask > 1; mask >>= 1) /* t0 = V_j (b) and t1 = V_{j+1} (b) */
    {
      if (k & mask) /* j -> 2*j+1. Compute V_{2j+1} and V_{2j+2} */
      {
        V_dadd (t0, t1, t0, b);
        V_dbl (t1, t1, two);
      }
      else /* j -> 2*j. Compute V_{2j} and V_{2j+1} */
      {
        V_dadd (t1, t1, t0, b);
        V_dbl (t0, t0, two);
      }
    }

    /* Deal with least significant bit outside the loop */
    if (k & mask)
    {
      V_dadd (t0, t1, t0, b); /* cannot have r instead of t0, if r is the
                               * same variable as b, the assert in
                               * mod_V_dadd would fail */
      set (r, t0);
      if (rp1)
        V_dbl (*rp1, t1, two);
    }
    else
    {
      V_dbl (r, t0, two);
      if (rp1)
      {
        V_dadd (t1, t1, t0, b); /* same as above */
        set (*rp1, t1);
      }
    }
  }
}

/* Compute modular inverses for n input residues. If c is not NULL,
   computes r[i] = c*a[i]^-1.
   If any of the residues is not invertible, returns 0 and contents of r are
   undefined.
   a and r must be non-overlapping. */
bool
Modulus::batchinv (Residue *r, const Residue *a, const size_t n,
              const Residue *c) const
{
    Residue R(*this);

    if (n == 0)
        return true;

    set(r[0], a[0]);
    for (size_t i = 1; i < n; i++) {
        mul(r[i], r[i-1], a[i]);
    }

    if (!inv(R, r[n-1]))
        return false;

    if (c != NULL) {
        mul(R, R, *c);
    }

    for (size_t i = n-1; i > 0; i--) {
        mul(r[i], R, r[i-1]);
        mul(R, R, a[i]);
    }
    set(r[0], R);
    return 1;
}
