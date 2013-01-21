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

/* The version fb_root_in_qlattice_31bits mandates that the coordinates
 * of the q-lattice are at most 31 bits, so that combinations such as
 * Rb1-a1 always fit within the interval ]-2^32p, +2^32p[
 */

/* This helper function is used for powers of 2. See below */
static inline fbprime_t
fb_root_in_qlattice_po2 (const fbprime_t p, const fbprime_t R,
        sieve_info_srcptr si);

static inline fbprime_t
fb_root_in_qlattice_31bits (const fbprime_t p, const fbprime_t R,
        const uint32_t invp, sieve_info_srcptr si)
{
  int64_t aux1, aux2;
  uint64_t u, v;
  fbprime_t add;

    /* Handle powers of 2 separately, REDC doesn't like them */
    if (UNLIKELY(p % 2 == 0))
        return fb_root_in_qlattice_po2(p, R, si);

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
    if (UNLIKELY(!invmod_redc_32(&v, p)))
      {
	/* root in i,j-plane is projective */
	if (UNLIKELY(!invmod_redc_32(&u, p)))
          {
            fprintf (stderr, "Error, root in (i,j)-plane is projective\n");
            exit (EXIT_FAILURE); /* Should never happen! */
          }
	add = p;
      }
    u *= v;
    return (fbprime_t) redc_u32 (u, p, invp) + add;
}

/* This one is dog slow, but should be correct in all cases */
static inline fbprime_t
fb_root_in_qlattice_generic (const fbprime_t p, const fbprime_t R,
        const uint64_t invp, sieve_info_srcptr si)
{
    /* Handle powers of 2 separately, REDC doesn't like them */
    if (UNLIKELY(p % 2 == 0))
        return fb_root_in_qlattice_po2(p, R, si);

    int64_t aux1, aux2;
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

    /* The root in the (i,j) plane is (aux1:aux2). Now let's put it
     * in of the two forms:
     * (x:1) -> we return x=aux1/aux2.
     * (1:x), with p|x -> we return p+x = p+aux2/aux1
     *
     * (note that p is not necessarily a prime, it may be a prime power
     */
    /* Do a full 64-bit redc */
    uint64_t u = redc_64(aux1, p, invp); /* 0 <= u < p */
    uint64_t v = redc_64(aux2, p, invp); /* 0 <= den < p */

    fbprime_t add = 0;
    if (UNLIKELY(!invmod_redc_64(&v, p)))
      {
	/* root in i,j-plane is projective */
	if (UNLIKELY(!invmod_redc_64(&u, p)))
          {
            fprintf (stderr, "Error, root in (i,j)-plane is projective\n");
            exit (EXIT_FAILURE); /* Should never happen! */
          }
	add = p;
      }
    u *= v;
    u = redc_64(u, p, invp);
    return (fbprime_t) u + add;
}

/* This is just for powers of 2, and is used by both versions above */

static inline fbprime_t fb_root_in_qlattice_po2 (const fbprime_t p, const fbprime_t R,
        sieve_info_srcptr si)
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

static inline fbprime_t
fb_root_in_qlattice(const fbprime_t p, const fbprime_t R,
        const uint32_t invp, sieve_info_srcptr si)
{
    return fb_root_in_qlattice_31bits(p, R, invp, si);
}

#endif	/* LAS_QLATTICE_H_ */
