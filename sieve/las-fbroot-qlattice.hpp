#ifndef CADO_LAS_FBROOT_QLATTICE_HPP
#define CADO_LAS_FBROOT_QLATTICE_HPP

#include <cstdio>
#include <cstdint>

#include <gmp.h>

#include "fb-types.hpp"
#include "las-arith.hpp"
#include "las-qlattice.hpp"
#include "macros.h"

/* CARRYCHECK is a ternary value here: 0 means no carry check, 1 means carry
   check, and 2 means that the function should choose the value of CARRYCHECK
   depending on the value of p. */
template <int CARRYCHECK = 2>
static inline fb_root_p1
fb_root_in_qlattice_31bits (fbprime_t p, fb_root_p1 R,
        redc_invp_t invp, const qlattice_basis &basis);
template <int CARRYCHECK = 2>
static inline fb_root_p1
fb_root_in_qlattice_127bits (fbprime_t p, fb_root_p1 R,
        redc_invp_t invp, const qlattice_basis &basis);

/* batch calls expect affine inputs and affine outputs */
template <int CARRYCHECK = 2>
static inline bool
fb_root_in_qlattice_31bits_batch (fbroot_t *r_ij, fbprime_t p,
        const fbroot_t *r_ab, redc_invp_t invp,
        const qlattice_basis &basis, size_t n_roots);
#if 0
template <int CARRYCHECK = 2>
static inline bool
fb_root_in_qlattice_127bits_batch (fbroot_t *r_ij, fbprime_t p,
        const fbroot_t *r_ab, redc_invp_t invp,
        const qlattice_basis &basis, size_t n_roots);
#endif
/* Specialize the case CARRYCHECK=2 and call the template instance
   with CARRYCHECK=0 or 1, depending on the value of p. */
template<> inline fb_root_p1
fb_root_in_qlattice_31bits<2>(const fbprime_t p, const fb_root_p1 R,
        const redc_invp_t invp, const qlattice_basis &basis)
{
    if (redc_no_carry(p))
        return fb_root_in_qlattice_31bits<0>(p, R, invp, basis);
    else
        return fb_root_in_qlattice_31bits<1>(p, R, invp, basis);
}

template<> inline  bool
fb_root_in_qlattice_31bits_batch<2> (fbroot_t *r_ij, const fbprime_t p,
        const fbroot_t *r_ab, const redc_invp_t invp, const qlattice_basis &basis,
        const size_t n_roots)
{
    if (redc_no_carry(p))
        return fb_root_in_qlattice_31bits_batch<0> (r_ij, p, r_ab, invp, basis,
                                                 n_roots);
    else
        return fb_root_in_qlattice_31bits_batch<1> (r_ij, p, r_ab, invp, basis,
                                                 n_roots);
}

template<> inline fb_root_p1
fb_root_in_qlattice_127bits<2>(const fbprime_t p, const fb_root_p1 R,
        const redc_invp_t invp, const qlattice_basis &basis)
{
    if (redc_no_carry(p))
        return fb_root_in_qlattice_127bits<0>(p, R, invp, basis);
    else
        return fb_root_in_qlattice_127bits<1>(p, R, invp, basis);
}

#if 0
template<> inline bool
fb_root_in_qlattice_127bits_batch<2> (fbroot_t *r_ij, const fbprime_t p,
        const fbroot_t *r_ab, const redc_invp_t invp, const qlattice_basis &basis,
        const size_t n_roots)
{
    if (redc_no_carry(p))
        return fb_root_in_qlattice_127bits_batch<0> (r_ij, p, r_ab, invp, basis,
                                              n_roots);
    else
        return fb_root_in_qlattice_127bits_batch<1> (r_ij, p, r_ab, invp, basis,
                                              n_roots);
}
#endif

/* fb_root_in_qlattice returns (R*b1-a1)/(a0-R*b0) mod p */
static inline fb_root_p1
fb_root_in_qlattice(const fbprime_t p, const fb_root_p1 R,
        const redc_invp_t invp, const qlattice_basis &basis)
{
#if defined(SUPPORT_LARGE_Q)
/* Here, we only assume that the q-lattice basis entries are int64_t's.
 * This means that redc_32 is a bit sub-optimal.
 */
    return fb_root_in_qlattice_127bits(p, R, invp, basis);
#else
/* We want the q-lattice entries to fit within 31 bits */
    return fb_root_in_qlattice_31bits(p, R, invp, basis);
#endif
}


static inline bool
fb_root_in_qlattice_batch (fbroot_t *r_ij MAYBE_UNUSED,
        const fbprime_t p MAYBE_UNUSED,
        const fbroot_t *r_ab MAYBE_UNUSED,
        const redc_invp_t invp MAYBE_UNUSED,
        const qlattice_basis &basis MAYBE_UNUSED,
        size_t n_roots MAYBE_UNUSED)
{
#if defined(SUPPORT_LARGE_Q)
    /* See #30012 and merge request !198 -- !196 was not complete! */
    /* By returning false, we force all roots to be processed one by one.
     */
    return false;
    /*
    return fb_root_in_qlattice_127bits_batch (r_ij, p, r_ab, invp, basis,
                                              n_roots);
                                              */
#else
    return fb_root_in_qlattice_31bits_batch (r_ij, p, r_ab, invp, basis,
                                             n_roots);
#endif
}

/* This helper function is used for powers of 2. See below */
static inline fb_root_p1
fb_root_in_qlattice_po2 (fbprime_t p, fb_root_p1 R,
        const qlattice_basis &basis);

/* This version fb_root_in_qlattice_31bits mandates that the coordinates
 * of the q-lattice are at most 31 bits, so that combinations such as
 * R*b1-a1 always fit within the interval ]-2^32p, +2^32p[.
 * It makes 3 calls to redc_32 and 1 to invmod_redc_32.
 * If CARRYCHECK is 0, we require that p satisfies redc_no_carry(p).
 */
template <int CARRYCHECK>
static inline fb_root_p1
fb_root_in_qlattice_31bits (const fbprime_t p, const fb_root_p1 R,
        const redc_invp_t invp, const qlattice_basis &basis)
{
  int64_t aux1, aux2;
  uint32_t u, v;

  ASSERT_EXPENSIVE(basis.fits_31bits());
  static_assert (CARRYCHECK == 0 || CARRYCHECK == 1 , "Invalid CARRYCHECK value" );
  
  /* Handle powers of 2 separately, REDC doesn't like them */
  if (UNLIKELY(!(p & 1)))
    return fb_root_in_qlattice_po2(p, R, basis);

    // Use Signed Redc for the computation:
    // Numerator and denominator will get divided by 2^32, but this does
    // not matter, since we take their quotient.

  if (R.is_affine()) {
      /* With 0 <= R <= 2^32-1, -2^31 <= a1, b1 <= 2^31-1, we have
       * -2^63 + 1 <= aux1 <= 2^63 - 2^32 + 1
       * With 0 <= R <= p-1, -2^31 <= a1, b1 <= 2^31-1, we have
       * -2^31*p + 1 <= aux1 <= (p-1)*2^31 + 1, thus -2^32*p < aux1 < 2^32*p
       */
      aux1 = (int64_t)R.r * basis.b1 - basis.a1;
      /* With 0 <= R <= 2^32-1, -2^31 <= a0, b0 <= 2^31-1, we have
       * -2^63 + 2^32 - 1 <= aux2 <= 2^63 - 1
       * With 0 <= R <= p-1, -2^31 <= a0, b0 <= 2^31-1, we have
       * -p(2^31 - 1) - 1 <= aux2 <= p*2^31 - 1, thus -2^32*p < aux2 < 2^32*p
       */
      aux2 = basis.a0 - (int64_t)R.r *basis.b0;
    }
  else /* Root in a,b-plane is projective */
    {
      aux1 = basis.b1 - (int64_t)(R.r) * basis.a1;
      aux2 = (int64_t)(R.r) * basis.a0 - basis.b0;
    }
  /* USE_NATIVE_MOD is slightly slower on Intel i5-4590 with gcc 9.2.1:
   * test_fb_root 10000 reports 14.49s instead of 13.26s
   * (same pattern on i7-8550U)
   */
// #define USE_NATIVE_MOD 1
#ifdef USE_NATIVE_MOD
  u = (aux1 >= 0) ? aux1 % p : p - ((-aux1) % p);
  v = (aux2 >= 0) ? aux2 % p : p - ((-aux2) % p);
#else
  u = redc_32<CARRYCHECK>(aux1, p, invp); /* 0 <= u < p */
  v = redc_32<CARRYCHECK>(aux2, p, invp); /* 0 <= v < p */
#endif

  aux2 = invmod_redc_32(v, p, invp);
  if (LIKELY(aux2)) {
      return fb_root_p1::affine_root(redc_u32<CARRYCHECK> (aux2 * u, p, invp));
  } else {
      /* root in i,j-plane is projective */
      aux2 = invmod_redc_32(u, p, invp);
      if (UNLIKELY(!aux2))
	{
	  fprintf (stderr, "fb_root_in_qlattice_31bits(%" FBPRIME_FORMAT ", %"
                   FBROOT_FORMAT ", %" PRIu32 "): Error, root %" PRIu32 "/%"
                   PRIu32 " in (i,j)-plane is projective\n", p, R.to_old_format(p), invp, u, v);
          ASSERT_ALWAYS(0);
	}
      return fb_root_p1::projective_root(redc_u32<CARRYCHECK> (aux2 * v, p, invp));
  }
}

/** Transforms roots r_ab (mod p) according to a lattice basis.
 *
 * Each transformed root r_ij[i] is
 * (r_ab[i] * b1 - a1) / (r_ab[i] * a0 - b0) mod p, where a0, a1, b0, b1
 * are the lattice basis coordinates.
 * This version mandates that each lattice coordinate is in [-2^31, 2^31-1].
 * It makes 3 calls to redc_32 per root and 1 to batchinvredc_u32 for the
 * whole batch.
 *
 * \param [out] r_ij    If all transformed roots are affine, contains the
 *                      transformed roots, i.e., in i,j-coordinates in the
 *                      q-lattice. Otherwise, contents undefined
 * \param [in]  p       The modulus for the root transform
 * \param [in]  r_ab    The roots in a,b-coordinates (not in q-lattice)
 * \param [in]  invp    Must satisfy p*invp == -1 (mod p)
 * \param [in]  basis   The basis of the q-lattice
 * \param [in]  n_roots The number of roots in r_ab and r_ij
 * \return true if all inverses exist, and false otherwise.
 */

template <int CARRYCHECK>
static inline bool
fb_root_in_qlattice_31bits_batch (fbroot_t *r_ij, const fbprime_t p, 
        const fbroot_t *r_ab, const redc_invp_t invp, const qlattice_basis &basis,
        const size_t n_roots)
{
  /* p must be odd for REDC to work */
  ASSERT(p % 2 == 1);
  ASSERT_EXPENSIVE(basis.fits_31bits());
  static_assert (CARRYCHECK == 0 || CARRYCHECK == 1 , "Invalid CARRYCHECK value" );

  for (size_t i_root = 0; i_root < n_roots; i_root++) {
      const int64_t den = basis.a0 - (int64_t)r_ab[i_root] *basis.b0;
      /* USE_NATIVE_MOD is slightly slower on Intel i5-4590 with gcc 9.2.1:
       * test_fb_root 10000 reports 14.49s instead of 13.26s
       * (same pattern on i7-8550U)
       */
//#define USE_NATIVE_MOD
#ifdef USE_NATIVE_MOD
      r_ij[i_root] = (den >= 0) ? den % p : p - ((-den) % p);
#else
      // Use Signed Redc for the computation:
      // Numerator and denominator will get divided by 2^32, but this does
      // not matter, since we take their quotient.
      // We use r_ij[] as temp storage
      r_ij[i_root] = redc_32<CARRYCHECK>(den, p, invp); /* 0 <= v < p */
#endif
  }

  uint32_t inverses[n_roots];
  // If any transformed root is projective, return false.
  if (batchinvredc_u32<CARRYCHECK>(inverses, r_ij, n_roots, p, invp) == 0) {
      return false;
  }

  for (size_t i_root = 0; i_root < n_roots; i_root++) {
      int64_t num = (int64_t)r_ab[i_root] * basis.b1 - basis.a1;
#ifdef USE_NATIVE_MOD
      uint32_t u = (num >= 0) ? num % p : p - ((-num) % p);
#else
      uint32_t u = redc_32<CARRYCHECK>(num, p, invp); /* 0 <= u < p */
#endif
      num = (int64_t) u * (int64_t) inverses[i_root];
      r_ij[i_root] = (fbroot_t) (redc_u32<CARRYCHECK>(num, p, invp));
  }

  return true;
}

static fb_root_p1_t<cxx_mpz>
reference_fb_root_in_qlattice (fbprime_t p, fb_root_p1 R, qlattice_basis const & basis)
{
  cxx_mpz num, den;

  if (R.is_affine()) {
      den = -basis.b0;
      mpz_mul_ui (den, den, R.r);
      mpz_add_int64 (den, den, basis.a0);

      num = basis.b1;
      mpz_mul_ui (num, num, R.r);
      mpz_sub_int64 (num, num, basis.a1);
  } else {
      den = basis.a0;
      mpz_mul_ui (den, den, R.r);
      mpz_sub_int64 (den, den, basis.b0);

      num = -basis.a1;
      mpz_mul_ui (num, num, R.r);
      mpz_add_int64 (num, num, basis.b1);
  }

  if (mpz_gcd_ui(nullptr, den, p) != 1) {
      /* root is projective */
      mpz_invert (num, num, cxx_mpz(p));
      mpz_mul(num, num, den);
      mpz_mod_ui(num, num, p);
      return { num, true };
  } else {
      mpz_invert (den, den, cxx_mpz(p));
      mpz_mul(num, num, den);
      mpz_mod_ui(num, num, p);
      return { num, false };
  }
}


/* This one is slower, but should be correct under the relaxed condition
 * that q be at most 127 bits or so, so that the coordinates of the
 * Q-lattice can be as large as 63 bits. We call redc many more times here.
 */
template <int CARRYCHECK>
static inline fb_root_p1
fb_root_in_qlattice_127bits (const fbprime_t p, const fb_root_p1 R,
        redc_invp_t MAYBE_UNUSED, const qlattice_basis &basis)
{
    auto ref = reference_fb_root_in_qlattice(p, R, basis);

    return { (fbprime_t) mpz_get_ui(ref.r), ref.proj };

#if 0
  int64_t aux1, aux2;
  uint64_t u, v;
  
  static_assert (CARRYCHECK == 0 || CARRYCHECK == 1 , "Invalid CARRYCHECK value" );
    /* Handle powers of 2 separately, REDC doesn't like them */
  if (UNLIKELY(!(p & 1 )))
    return fb_root_in_qlattice_po2(p, R, basis);
  
  if (R.is_affine()) {
        /* The products R*redc(b1) must not overflow. Therefore we must
         * be extra careful when p exceeds NextPrime(Floor(2^31.5))==3037000507
      aux1 = ((int64_t)R)*(int64_t) redc_32(basis.b1, p, invp) - (int64_t) redc_32(basis.a1, p, invp);
      aux2 = (int64_t) redc_32(basis.a0, p, invp) - ((int64_t)R)*(int64_t) redc_32(basis.b0, p, invp);
         */
#if 0
      /* we decompose the basis entries as hi*2^32+lo with lo and hi
       * signed. We form these as 32-bit integers, but cast them to
       * 64-bit ones because it's needed for the multiplications that we
       * do.
       */
      const int64_t b0_lo = (int32_t) basis.b0;
      const int64_t b0_hi = (int32_t) (basis.b0 >> 32) + (b0_lo < 0);
      const int64_t b1_lo = (int32_t) basis.b1;
      const int64_t b1_hi = (int32_t) (basis.b1 >> 32) + (b1_lo < 0);
      const int64_t a0_lo = (int32_t) basis.a0;
      const int64_t a0_hi = (int32_t) (basis.a0 >> 32) + (a0_lo < 0);
      const int64_t a1_lo = (int32_t) basis.a1;
      const int64_t a1_hi = (int32_t) (basis.a1 >> 32) + (a1_lo < 0);
      const int64_t r = R.r;
      /* compute b1*r-a1 and a0*r-b0. We do this separately for both words.
       * We get resuls of reach redc_32<> in [0,p).
       */
      const uint32_t aux1_lo = redc_32<CARRYCHECK>(r * b1_lo - a1_lo);
      const uint32_t aux1_hi = redc_32<CARRYCHECK>(r * b1_hi - a1_hi);
      const uint32_t aux2_lo = redc_32<CARRYCHECK>(r * a0_lo - b0_lo);
      const uint32_t aux2_hi = redc_32<CARRYCHECK>(r * a0_hi - b0_hi);
#endif


        const uint64_t Rl = R.r;
        const uint64_t b1l = redc_32<CARRYCHECK>(basis.b1, p, invp);
        const uint64_t b0l = redc_32<CARRYCHECK>(basis.b0, p, invp);
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
        aux1 = aux1 - redc_32<CARRYCHECK>(basis.a1, p, invp);
        aux2 = redc_32<CARRYCHECK>(basis.a0, p, invp) - aux2;
    }
  else /* Root in a,b-plane is projective */
    {
        /*
      aux1 = (int64_t) redc_32(basis.b1, p, invp) - ((int64_t)(R - p))*(int64_t) redc_32(basis.a1, p, invp);
      aux2 = ((int64_t)(R - p))*(int64_t) redc_32(basis.a0, p, invp) - (int64_t) redc_32(basis.b0, p, invp);
      */
        const uint64_t Rpl = R.r;
        const uint64_t a1l = redc_32<CARRYCHECK>(basis.a1, p, invp);
        const uint64_t a0l = redc_32<CARRYCHECK>(basis.a0, p, invp);
        aux1 = Rpl*a1l;
        aux2 = Rpl*a0l;
        /* same analysis as above */
        if (aux1 < 0) aux1 -= ((uint64_t)p)<<32;
        if (aux2 < 0) aux2 -= ((uint64_t)p)<<32;
        aux1 = aux1 - redc_32<CARRYCHECK>(basis.b1, p, invp);
        aux2 = redc_32<CARRYCHECK>(basis.b0, p, invp) - aux2;
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
  u = redc_32<CARRYCHECK>(aux1, p, invp); /* 0 <= u < p */
  v = redc_32<CARRYCHECK>(aux2, p, invp); /* 0 <= v < p */
#endif
  
  aux2 = invmod_redc_32(v, p, invp);
  if (LIKELY(aux2)) {
    /* Warning: since 0 <= u < p and 0 <= aux2 < p, we have
       0 <= aux2 < p^2, which might overflow the int64_t range
       if p >= 3037000507. To avoid this, we subtract p if aux2 >= 2^31:
       * if aux2 < 2^31, then aux2 * u < 2^31*p < 2^63
       * if aux2 >= 2^31, this implies that p >= 2^31 since aux2 < p,
       then (aux2 - p) * u > (2^31-p) * u > -2^31*p > -2^63 */

      /* gcc-10.2.1 turns this into a cmov:
          db4a:       48 3d ff ff ff 7f       cmp    $0x7fffffff,%rax
          db50:       48 0f 4f c2             cmovg  %rdx,%rax
       */
      if (aux2 >= 2147483648L)
          aux2 -= p;
      return fb_root_p1::affine_root(redc_32<CARRYCHECK> (aux2 * u, p, invp));
  } else {
      /* root in i,j-plane is projective */
      aux2 = invmod_redc_32(u, p, invp);
      if (UNLIKELY(!aux2))
	{
	  fprintf (stderr, "Error, root in (i,j)-plane is projective\n");
	  ASSERT_ALWAYS(0);
	}
      /* Warning: we have the same overflow problem as above. */
      return fb_root_p1::projective_root(redc_32<CARRYCHECK> (aux2 * v, p, invp));
    }
#endif

}

#if 0
/* see #30112 / !198 */
/** Transforms roots r_ab (mod p) according to a lattice basis.
 *
 * Each transformed root r_ij[i] is
 * (r_ab[i] * b1 - a1) / (r_ab[i] * a0 - b0) mod p, where a0, a1, b0, b1
 * are the lattice basis coordinates.
 * This version allows the lattice coordinates in [-2^63, 2^63-1].
 * It makes 7 calls to redc_32 per root and 1 to batchinvredc_u32 for the
 * whole batch.
 *
 * \param [out] r_ij    If all transformed roots are affine, contains the
 *                      transformed roots, i.e., in i,j-coordinates in the
 *                      q-lattice. Otherwise, contents undefined
 * \param [in]  p       The modulus for the root transform
 * \param [in]  r_ab    The roots in a,b-coordinates (not in q-lattice)
 * \param [in]  invp    Must satisfy p*invp == -1 (mod p)
 * \param [in]  basis   The basis of the q-lattice
 * \param [in]  n_roots The number of roots in r_ab and r_ij
 * \return true if all inverses exist, and false otherwise.
 */
template <int CARRYCHECK>
static inline bool
fb_root_in_qlattice_127bits_batch (fbroot_t *r_ij, const fbprime_t p,
        const fbroot_t *r_ab, const redc_invp_t invp, const qlattice_basis &basis,
        const size_t n_roots)
{
    ASSERT(p % 2 == 1);
    static_assert (CARRYCHECK == 0 || CARRYCHECK == 1 , "Invalid CARRYCHECK value" );

    for (size_t i_root = 0; i_root < n_roots; i_root++) {
        int64_t den;
        const uint64_t Rl = r_ab[i_root];
        const uint64_t b0l = redc_32<CARRYCHECK>(basis.b0, p, invp);
        den = Rl*b0l;
        if (den < 0) den -= ((uint64_t)p)<<32;
        den = redc_32<CARRYCHECK>(basis.a0, p, invp) - den;

      // We use r_ij[] as temp storage
#ifdef USE_NATIVE_MOD
        r_ij[i_root] = (den >= 0) ? den % p : p - ((-den) % p);
#else
        r_ij[i_root] = redc_32<CARRYCHECK>(den, p, invp); /* 0 <= r_ij[i_root]v < p */
#endif
    }

    uint32_t inverses[n_roots];
    // If any transformed root is projective, return false.
    if (batchinvredc_u32<CARRYCHECK>(inverses, r_ij, n_roots, p, invp) == 0) {
        return false;
    }

    for (size_t i_root = 0; i_root < n_roots; i_root++) {
        int64_t aux1;
        uint32_t u;
        const uint64_t Rl = r_ab[i_root];
        const uint64_t b1l = redc_32<CARRYCHECK>(basis.b1, p, invp);
        aux1 = Rl*b1l;
        if (aux1 < 0) aux1 -= ((uint64_t)p)<<32;
        aux1 = aux1 - redc_32<CARRYCHECK>(basis.a1, p, invp);

#ifdef USE_NATIVE_MOD
        u = (aux1 >= 0) ? aux1 % p : p - ((-aux1) % p);
#else
        u = redc_32<CARRYCHECK>(aux1, p, invp); /* 0 <= u < p */
#endif
        r_ij[i_root] = (fbprime_t) mulmodredc_u32<CARRYCHECK>(u, inverses[i_root], p, invp);
    }

    return true;
}
#endif

static inline bool
fb_root_in_qlattice_127bits_batch (fbroot_t *r_ij, fbprime_t p,
        const fbroot_t *r_ab, redc_invp_t invp,
        const qlattice_basis &basis, size_t n_roots)
{
    // tests definitely want to have an fb_root_in_qlattice_127bits_batch
    // function, so we'll make them one. However this function is not
    // used in production.
    // Batch transform failed: do roots one at a time.
    for (unsigned char i = 0; i != n_roots; i++) {
        auto R = fb_root_in_qlattice_127bits(p, r_ab[i], invp, basis);
        if (R.proj) return false;
        r_ij[i] = R.r;
    }
    return true;
}


/* This is just for powers of 2, and is used by both versions above */

static inline fb_root_p1 fb_root_in_qlattice_po2 (const fbprime_t p, const fb_root_p1 R, const qlattice_basis &basis)
{
    fbprime_t u, v;
    ASSERT(p == (p & -p)); /* Test that p is power of 2 */
    if (R.is_affine()) {
	u = (int64_t)R.r * (basis.b1 % p) - basis.a1;
	v = basis.a0 - (int64_t)R.r * (basis.b0 % p);
    } else {
        u = basis.b1 - (int64_t)(R.r) * (basis.a1 % p);
        v = (int64_t)(R.r) * (basis.a0 % p) - basis.b0;
    }
    
    if (v & 1) {
        /* root in i,j-plane is non-projective */
        v = invmod_po2 (v);
        return fb_root_p1::affine_root((u * v) & (p - 1));
    } else {
        /* root in i,j-plane is projective */
        u = invmod_po2 (u);
        return fb_root_p1::projective_root((u * v) & (p - 1));
    }
}

#endif	/* CADO_LAS_FBROOT_QLATTICE_HPP */
