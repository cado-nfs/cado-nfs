#ifndef LAS_QLATTICE_HPP_
#define LAS_QLATTICE_HPP_

#include <cstdint>
#include <gmp.h>
#include "fb-types.h"         /* fbprime_t */
#include "las-base.hpp"
#include "portability.h"
#include "las-todo-entry.hpp"

/* implementations for inlines */
#include "las-arith.hpp"

struct qlattice_basis {
    las_todo_entry doing;

    int64_t a0=0, b0=0, a1=0, b1=0;
    unsigned long q_ulong=0;
    // q (== doing.p) itself or 0 if q is too large to fit

    // handy to have here.
    sublat_t sublat;

    qlattice_basis() = default;
    inline double skewed_norm0(double s) const { return a0*a0/s+b0*b0*s; }
    inline double skewed_norm1(double s) const { return a1*a1/s+b1*b1*s; }

    // Assumes ell is prime.
    bool is_coprime_to(unsigned long ell) const {
        if (doing.is_prime()) {
            return (ell != q_ulong);
        } else {
            for (auto const & p : doing.prime_factors)
                if (p == ell)
                    return false;
            return true;
        }
    }

    inline bool fits_31bits() const {
        return !(
                 a0 <= INT64_C(-2147483648) ||
                 a0 >= INT64_C( 2147483648) ||
                 a1 <= INT64_C(-2147483648) ||
                 a1 >= INT64_C( 2147483648) ||
                 b0 <= INT64_C(-2147483648) ||
                 b0 >= INT64_C( 2147483648) ||
                 b1 <= INT64_C(-2147483648) ||
                 b1 >= INT64_C( 2147483648)
                 );
    }

    struct too_skewed : public std::exception { };

    qlattice_basis(las_todo_entry const & doing, double skew);
};

static inline fbprime_t
fb_root_in_qlattice_31bits (const fbprime_t p, const fbprime_t R,
        const uint32_t invp, const qlattice_basis &basis);
static inline fbprime_t
fb_root_in_qlattice_127bits (const fbprime_t p, const fbprime_t R,
        const uint64_t invp, const qlattice_basis &basis);


/* fb_root_in_qlattice returns (R*b1-a1)/(a0-R*b0) mod p */
#if defined(SUPPORT_LARGE_Q)
/* The reason why the special-q is constrained to some limit is quite
 * clearly linked to the fb_root_in_qlattice variant being used. However,
 * it does not seem to be exactly 31 or 127 bits. This should be
 * investigated */
#define MAX_SPECIALQ_BITSIZE    126
static inline fbprime_t
fb_root_in_qlattice(const fbprime_t p, const fbprime_t R,
        const redc_invp_t invp, const qlattice_basis &basis);
static inline fbprime_t
fb_root_in_qlattice(const fbprime_t p, const fbprime_t R,
        const redc_invp_t invp, const qlattice_basis &basis)
{
    return fb_root_in_qlattice_127bits(p, R, invp, basis);
}
#else

#define MAX_SPECIALQ_BITSIZE    30
static inline fbprime_t
fb_root_in_qlattice(const fbprime_t p, const fbprime_t R,
        const redc_invp_t invp, const qlattice_basis &basis);
static inline fbprime_t
fb_root_in_qlattice(const fbprime_t p, const fbprime_t R,
        const redc_invp_t invp, const qlattice_basis &basis)
{
    return fb_root_in_qlattice_31bits(p, R, invp, basis);
}
#endif

/* This helper function is used for powers of 2. See below */
static inline fbprime_t
fb_root_in_qlattice_po2 (const fbprime_t p, const fbprime_t R,
        const qlattice_basis &basis);

/* The version fb_root_in_qlattice_31bits mandates that the coordinates
 * of the q-lattice are at most 31 bits, so that combinations such as
 * Rb1-a1 always fit within the interval ]-2^32p, +2^32p[.
 * It makes 3 calls to redc_32 and 1 to invmod_redc_32.
 */
static inline fbprime_t
fb_root_in_qlattice_31bits (const fbprime_t p, const fbprime_t R,
        const uint32_t invp, const qlattice_basis &basis)
{
  int64_t aux1, aux2;
  uint32_t u, v;

    /* Handle powers of 2 separately, REDC doesn't like them */
  if (UNLIKELY(!(p & 1)))
    return fb_root_in_qlattice_po2(p, R, basis);

    // Use Signed Redc for the computation:
    // Numerator and denominator will get divided by 2^32, but this does
    // not matter, since we take their quotient.

  if (LIKELY(R < p)) /* Root in a,b-plane is affine */
    {
      aux1 = (int64_t)R * basis.b1 - basis.a1;
      aux2 = basis.a0 - (int64_t)R *basis.b0;
    }
  else /* Root in a,b-plane is projective */
    {
      aux1 = basis.b1 - (int64_t)(R - p) * basis.a1;
      aux2 = (int64_t)(R - p) * basis.a0 - basis.b0;
    }
  /* USE_NATIVE_MOD is slightly slower on Intel i5-4590 with gcc 9.2.1:
   * test_fb_root 10000 reports 14.49s instead of 13.26s
   * (same pattern on i7-8550U)
   */
//#define USE_NATIVE_MOD
#ifdef USE_NATIVE_MOD
  u = (aux1 >= 0) ? aux1 % p : p - ((-aux1) % p);
  v = (aux2 >= 0) ? aux2 % p : p - ((-aux2) % p);
#else
  u = redc_32(aux1, p, invp); /* 0 <= u < p */
  v = redc_32(aux2, p, invp); /* 0 <= v < p */
#endif

  aux2 = invmod_redc_32(v, p);
  if (LIKELY(aux2)) {
    aux1 = 0;
    aux2 *= u;
  }
  else 
    {
      /* root in i,j-plane is projective */
      aux2 = invmod_redc_32(u, p);
      if (UNLIKELY(!aux2))
	{
	  fprintf (stderr, "Error, root in (i,j)-plane is projective\n");
          ASSERT_ALWAYS(0);
	}
      aux1 = p;
      aux2 *= v;
    }
  return (fbprime_t) (redc_u32 (aux2, p, invp) + aux1);
}

/* This one is slower, but should be correct under the relaxed condition
 * that q be at most 127 bits or so, so that the coordinates of the
 * Q-lattice can be as large as 63 bits. We call redc 7 times here, instead
 * of 3 for fb_root_in_qlattice_31bits.
 */
static inline fbprime_t
fb_root_in_qlattice_127bits (const fbprime_t p, const fbprime_t R,
        const uint64_t invp, const qlattice_basis &basis)
{
  int64_t aux1, aux2;
  uint64_t u, v;
  
    /* Handle powers of 2 separately, REDC doesn't like them */
  if (UNLIKELY(!(p & 1 )))
    return fb_root_in_qlattice_po2(p, R, basis);
  
  if (LIKELY(R < p)) /* Root in a,b-plane is affine */
    {
        /* The products R*redc(b1) must not overflow. Therefore we must
         * be extra careful when p exceeds NextPrime(Floor(2^31.5))==3037000507
      aux1 = ((int64_t)R)*(int64_t) redc_32(basis.b1, p, invp) - (int64_t) redc_32(basis.a1, p, invp);
      aux2 = (int64_t) redc_32(basis.a0, p, invp) - ((int64_t)R)*(int64_t) redc_32(basis.b0, p, invp);
         */
        uint64_t Rl = R;
        uint64_t b1l = redc_32(basis.b1, p, invp);
        uint64_t b0l = redc_32(basis.b0, p, invp);
        aux1 = Rl*b1l;
        aux2 = Rl*b0l;
        /* If we have an overflow in the products above, replace by
         * something which is in the same congruence class mod p, and is
         * also within the correct range.
         *
         * aux1 < 0 can happen if 2^31.5 <= p <= 2^32-1 < 2^32, which
         * implies in particular
         *      2^63 <= (p-1)^2 <= 2^64-2*2^32+1
         *                      <= 2^64-2*p
         * we also have
         *      2^63.5 <= p*2^32
         *             <= (2^32-1)*2^32
         *             <= 2^64-p
         * If aux1 < 0, then it is a representative of something in the range
         * [2^63..2^64-2*p[
         * when we subtract a correction equal to p*2^32, we have;
         *      2^63-2^64+p <= aux1 - correction <= 2^64-2*p-2^63.5
         *          -2^63+p <= aux1 - correction <= 0.58*2^63-2*p
         *
         * (we could have used p<=2^32-5 as well, since p is prime -- but
         * assuming p odd is sufficient)
         */
        if (aux1 < 0) aux1 -= ((uint64_t)p)<<32;
        if (aux2 < 0) aux2 -= ((uint64_t)p)<<32;
        /* We don't want any of the two subtractions below to overflow
         * either.  The bounds above show that it can't happen if we
         * added a correction, because we are bounded away from the type
         * limits.
         *
         * If we did _not_ have to add a correction, then we're happy as
         * well, since aux1 is positive in that case, and the range
         * [0..2^63-1[ is safe for both subtractions below.
         */
        aux1 = aux1 - redc_32(basis.a1, p, invp);
        aux2 = redc_32(basis.a0, p, invp) - aux2;
    }
  else /* Root in a,b-plane is projective */
    {
        /*
      aux1 = (int64_t) redc_32(basis.b1, p, invp) - ((int64_t)(R - p))*(int64_t) redc_32(basis.a1, p, invp);
      aux2 = ((int64_t)(R - p))*(int64_t) redc_32(basis.a0, p, invp) - (int64_t) redc_32(basis.b0, p, invp);
      */
        uint64_t Rpl = R - p;
        uint64_t a1l = redc_32(basis.a1, p, invp);
        uint64_t a0l = redc_32(basis.a0, p, invp);
        aux1 = Rpl*a1l;
        aux2 = Rpl*a0l;
        /* same analysis as above */
        if (aux1 < 0) aux1 -= ((uint64_t)p)<<32;
        if (aux2 < 0) aux2 -= ((uint64_t)p)<<32;
        aux1 = aux1 - redc_32(basis.b1, p, invp);
        aux2 = redc_32(basis.b0, p, invp) - aux2;
    }
  
  /* The root in the (i,j) plane is (aux1:aux2). Now let's put it
   * in of the two forms:
   * (x:1) -> we return x=aux1/aux2.
   * (1:x), with p|x -> we return p+x = p+aux2/aux1
   *
   * (note that p is not necessarily a prime, it may be a prime power
   */
#ifdef USE_NATIVE_MOD
  u = (aux1 >= 0) ? aux1 % p : p - ((-aux1) % p);
  v = (aux2 >= 0) ? aux2 % p : p - ((-aux2) % p);
#else
  u = redc_32(aux1, p, invp); /* 0 <= u < p */
  v = redc_32(aux2, p, invp); /* 0 <= v < p */
#endif
  
  aux2 = invmod_redc_32(v, p);
  if (LIKELY(aux2)) {
    aux1 = 0;
    /* Warning: since 0 <= u < p and 0 <= aux2 < p, we have
       0 <= aux2 < p^2, which might overflow the int64_t range
       if p >= 3037000507. To avoid this, we subtract p if aux2 >= 2^31:
       * if aux2 < 2^31, then aux2 * u < 2^31*p < 2^63
       * if aux2 >= 2^31, this implies that p >= 2^31 since aux2 < p,
       then (aux2 - p) * u > (2^31-p) * u > -2^31*p > -2^63 */
    aux2 = (aux2 < 2147483648L) ? aux2 * u : (aux2 - p) * u;
  }
  else 
    {
      /* root in i,j-plane is projective */
      aux2 = invmod_redc_32(u, p);
      if (UNLIKELY(!aux2))
	{
	  fprintf (stderr, "Error, root in (i,j)-plane is projective\n");
	  ASSERT_ALWAYS(0);
	}
      aux1 = p;
      /* Warning: we have the same overflow problem as above. */
      aux2 *= v;
    }
  return (fbprime_t) (redc_32 (aux2, p, invp) + aux1);
}

/* This is just for powers of 2, and is used by both versions above */

static inline fbprime_t fb_root_in_qlattice_po2 (const fbprime_t p, const fbprime_t R, const qlattice_basis &basis)
{
    fbprime_t u, v;
    ASSERT(p == (p & -p)); /* Test that p is power of 2 */
    if (R < p) /* Root in a,b-plane is non-projective */
      {
	u = (int64_t)R * (basis.b1 % p) - basis.a1;
	v = basis.a0 - (int64_t)R * (basis.b0 % p);
      }
    else /* Root in a,b-plane is projective */
      {
        u = basis.b1 - (int64_t)(R - p) * (basis.a1 % p);
        v = (int64_t)(R - p) * (basis.a0 % p) - basis.b0;
      }
    
    if (v & 1)
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

std::ostream& operator<<(std::ostream& os, qlattice_basis const & Q);

#endif	/* LAS_QLATTICE_HPP_ */
