#ifndef LAS_QLATTICE_H_
#define LAS_QLATTICE_H_

#include <stdint.h>
#include "las-types.h"
#include "fb.h"         /* fbprime_t */

#ifdef __cplusplus
extern "C" {
#endif

int SkewGauss (sieve_info_ptr si, double skewness);

static inline fbprime_t
fb_root_in_qlattice (const fbprime_t p, const fbprime_t R,
        const uint32_t invp, sieve_info_srcptr si);


#ifdef __cplusplus
}
#endif

/* implementations for inlines */
#include "las-arith.h"

static inline fbprime_t
fb_root_in_qlattice (const fbprime_t p, const fbprime_t R,
        const uint32_t invp, sieve_info_srcptr si)
{
  int64_t aux1, aux2;
  uint64_t u, v;
  fbprime_t add;

    /* Handle powers of 2 separately, REDC doesn't like them */
    if (UNLIKELY(p % 2 == 0))
      {
	fbprime_t u, v;
	ASSERT(p == (p & -p)); /* Test that p is power of 2 */
	if (R < p) /* Root in a,b-plane is non-projective */
	  {
	    u = R * si->b1 - si->a1;
	    v = si->a0 - R * si->b0;
	  }
	else /* Root in a,b-plane is projective */
	  {
	    u = si->b1 - (R - p) * si->a1;
	    v = (R - p) * si->a0 - si->b0;
	  }

	if (v % 2 != 0)
	  {
	    /* root in i,j-plane is non-projective */
	    v = invmod_po2 (v);
	    return (u * v) & (p - 1);
	  }
	else
	  {
	    /* root in i,j-plane is projective */
	    u = invmod_po2 (u);
	    return p + ((u * v) & (p - 1));
	  }
      }

    // Use Signed Redc for the computation:
    // Numerator and denominator will get divided by 2^32, but this does
    // not matter, since we take their quotient.

    if (LIKELY(R < p)) /* Root in a,b-plane is affine */
      {
	aux1 = ((int64_t)R)*si->b1 - si->a1;
	aux2 = si->a0 - ((int64_t)R)*si->b0;
      }
    else /* Root in a,b-plane is projective */
      {
	aux1 = si->b1 - ((int64_t)(R - p))*si->a1;
	aux2 = ((int64_t)(R - p))*si->a0 - si->b0;
      }
    u = redc_32(aux1, p, invp); /* 0 <= u < p */
    v = redc_32(aux2, p, invp); /* 0 <= den < p */

    add = 0;
    if (UNLIKELY(!invmod_REDC(&v, p)))
      {
	/* root in i,j-plane is projective */
	if (UNLIKELY(!invmod_REDC(&u, p)))
          {
            fprintf (stderr, "Error, root in (i,j)-plane is projective\n");
            exit (EXIT_FAILURE); /* Should never happen! */
          }
	add = p;
      }
    u *= v;
    return (fbprime_t) redc_u32 (u, p, invp) + add;
}

#endif	/* LAS_QLATTICE_H_ */
