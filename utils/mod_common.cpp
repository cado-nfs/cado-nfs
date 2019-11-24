#include "cado.h"
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
