#ifndef UTILS_ARITHXX_API_IMPL_HPP_
#define UTILS_ARITHXX_API_IMPL_HPP_

#include <cstdint>
#include <cstddef>

#include "arithxx_common.hpp"
#include "u64arith.h"

/* This provides the _default_ instantiations for the non-inline functions
 * in arithxx_details::api.
 *
 * This include file is used only by the .cpp files. The minimal version
 * is when these file just emit the specialization, and are done with it.
 *
 * Optionally, a .cpp file may also override one of the specializations.
 * This is done in mod_mpz_new.cpp
 */

/* Compute r = b^e. Here e is a multiple precision integer
   sum_{i=0}^{e_nrwords-1} e[i] * (machine word base)^i.
   Assume e[e_nrwords-1] is not zero when e_nrwords > 0.
*/
template <typename layer>
void arithxx_details::api<layer>::pow(Residue & r,
                                                   Residue const & b,
                                                   uint64_t const * e,
                                                   size_t const e_nrwords) const
{
    auto const & me = downcast();
    size_t i = e_nrwords;
    auto const msb = uint64_t(1) << 63;
    uint64_t mask;
    Residue t(me);

    while (i > 0 && e[i - 1] == 0)
        i--;

    if (i == 0) {
        me.set1(r);
        return;
    }

    /* Find highest set bit in e[i]. */
    mask = msb >> u64arith_clz(e[i - 1]);
    /* t = 1, so t^(mask/2) * b^e = t^mask * b^e  */

    me.set(t, b); /* (r*b)^mask * b^(e-mask) = r^mask * b^e */
    mask >>= 1;

    for (; i > 0; i--) {
        auto const word = e[i - 1];
        for ( ; mask > 0 ; mask >>=1) {
            me.sqr(t, t);
            if (word & mask)
                me.mul(t, t, b);
        }
        mask = msb;
    }
    me.set(r, t);
}

template <typename layer>
inline void arithxx_details::api<layer>::pow2(Residue & r, uint64_t const * e,
                        size_t const e_nrwords) const
{
    auto const & me = downcast();
    size_t i = e_nrwords;
    auto const msb = uint64_t(1) << 63;
    uint64_t mask;
    Residue t(me);

    while (i > 0 && e[i - 1] == 0)
        i--;

    if (i == 0) {
        me.set1(r);
        return;
    }

    mask = msb >> u64arith_clz(e[i - 1]);
    mask >>= 1;

    me.set1(t);
    me.add(t, t, t);

    for (; i > 0; i--) {
        auto const word = e[i - 1];
        for ( ; mask > 0 ; mask >>=1) {
            me.sqr (t, t);
            if (word & mask) {
                me.add (t, t, t);
            }
        }
        mask = msb;
    }
    me.set(r, t);
}

template <typename layer>
inline void arithxx_details::api<layer>::pow(Residue & r, Residue const & b, uint64_t e) const
{
    pow(r, b, &e, 1);
}

template <typename layer>
inline void arithxx_details::api<layer>::pow2(Residue &r, const uint64_t e) const {
    pow2(r, &e, 1);
}

/* this does not work with mpz! */
template <typename layer>
inline void arithxx_details::api<layer>::pow(Residue & r, Residue const & b, Integer const & e) const {
    pow(r, b, e.data(), e.size_in_words());
}
template <typename layer>
inline void arithxx_details::api<layer>::pow2(Residue &r, const Integer &e) const {
    pow2(r, e.data(), e.size_in_words());
}

/* Compute r = V_k (b) and rp1 = V_{k+1} (b) if rp1 != NULL
 * r or rp1 can be the same variable as b.
 */
template <typename layer>
void arithxx_details::api<layer>::V(Residue & r, Residue * rp1,
                                                 Residue const & b,
                                                 Integer const & k) const
{
    auto const & me = downcast();
    Residue t0(me), t1(me), two(me);

    me.set1(two);
    me.add1(two, two);

    if (k == 0UL) {
        me.set(r, two);
        if (rp1)
            me.set(*rp1, b);
    } else if (k == 1UL) {
        me.set(r, b);
        if (rp1)
            me.V_dbl(*rp1, b, two);
    } else if (k == 2UL) {
        if (rp1) {
            me.V_dbl(t1, b, two);
            me.V_dadd(*rp1, t1, b, b);
            me.set(r, t1);
        } else
            me.V_dbl(r, b, two);
    } else /* k >= 3 */
    {
        /* Montgomery Ladder */
        Integer mask = Integer(1) << (k.bits()-1);

        /* Most significant bit of k is 1, do it outside the loop */
        me.set(t0, b);        /* starting value t0 = V_1 (b) = b */
        me.V_dbl(t1, b, two); /* starting value t1 = V_2 (b) */


        // t_0 = V_{k // mask}(b)
        // t_1 = V_{k // mask}(b)


        /* XXX explicit conversions to bool are mandatory for cxx_mpz's!
         */
        /* If the second most significant bit of k is 0, then we do the
         * iteration manually (to avoid to compute again V_2 (b)) As k >= 3, we
         * know that in this case k has at least 3 bits.
         */
        mask >>= 1;
        if (!bool(k & mask)) /* (t0,t1) <- (V_2 (b), V_3 (b)) */
        {
            me.set(t0, t1);
            me.V_dadd(t1, t1, b, b);
            // t_0 = V_{k // mask}(b)
            // t_1 = V_{k // mask}(b)
            mask >>= 1;
        }

        for (; mask > 1; mask >>= 1) /* t0 = V_j (b) and t1 = V_{j+1} (b) */
        {
            if (bool(k & mask)) /* j -> 2*j+1. Compute V_{2j+1} and V_{2j+2} */
            {
                me.V_dadd(t0, t1, t0, b);
                me.V_dbl(t1, t1, two);
            } else /* j -> 2*j. Compute V_{2j} and V_{2j+1} */
            {
                me.V_dadd(t1, t1, t0, b);
                me.V_dbl(t0, t0, two);
            }
            // t_0 = V_{k // mask}(b)
            // t_1 = V_{k // mask}(b)
        }

        /* Deal with least significant bit outside the loop */
        if (bool(k & mask)) {
            me.V_dadd(t0, t1, t0, b); /* cannot have r instead of t0, if r is the
                                    * same variable as b, the assert in
                                    * mod_V_dadd would fail */
            me.set(r, t0);
            if (rp1)
                me.V_dbl(*rp1, t1, two);
        } else {
            me.V_dbl(r, t0, two);
            if (rp1) {
                me.V_dadd(t1, t1, t0, b); /* same as above */
                me.set(*rp1, t1);
            }
        }
    }
}

/* Compute modular inverses for n input residues. If c is not NULL,
   computes r[i] = c*a[i]^-1.
   If any of the residues is not invertible, returns 0 and contents of r are
   undefined.
   a and r must be non-overlapping. */
template <typename layer>
bool arithxx_details::api<layer>::batchinv(Residue * r,
                                                        Residue const * a,
                                                        size_t const n,
                                                        Residue const * c) const
{
    auto const & me = downcast();
    Residue R(me);

    if (n == 0)
        return true;

    me.set(r[0], a[0]);
    for (size_t i = 1; i < n; i++) {
        me.mul(r[i], r[i - 1], a[i]);
    }

    if (!me.inv(R, r[n - 1]))
        return false;

    if (c)
        me.mul(R, R, *c);

    for (size_t i = n - 1; i > 0; i--) {
        me.mul(r[i], R, r[i - 1]);
        me.mul(R, R, a[i]);
    }
    me.set(r[0], R);
    return true;
}

template <typename layer>
bool arithxx_details::api<layer>::intinv(Integer & r, Integer const & a) const
{
    auto const & me = downcast();
    Residue R(me);
    bool const b = downcast().inv(R, downcast()(a));
    r = me.get(R);
    return b;
}
template <typename layer>
bool arithxx_details::api<layer>::inv_odd(Residue & r, Residue const & a) const
{
    return downcast().inv(r, a);
}
template <typename layer>
bool arithxx_details::api<layer>::inv_powerof2(Residue & r, Residue const & a) const
{
    return downcast().inv(r, a);
}

/* Returns 1 if m is a strong probable prime wrt base 2, 0 otherwise.
   Assumes m > 1 and is odd.
 */
template <typename layer>
bool arithxx_details::api<layer>::is_strong_pseudoprime_base2() const
{
    auto const & me = downcast();
    Residue r(me), minusone(me);
    int po2 = 1;

    /* If m == 1,7 (mod 8), then 2 is a quadratic residue, and we must find
       -1 with one less squaring. This does not reduce the number of
       pseudo-primes because strong pseudo-primes are also Euler pseudo-primes,
       but makes identifying composites a little faster on average. */
    auto const mod8 = uint64_t(Integer(me.m) & 7);
    if (mod8 == 1 || mod8 == 7)
        po2--;

    /* Set mm1 to the odd part of m-1 */
    auto mm1 = (Integer(me.m) - 1) >> 1;
    int k = mm1.ctz();
    po2 += k;
    mm1 >>= k;
    /* Hence, m-1 = mm1 * 2^po2 */

    me.set1(minusone);
    me.neg(minusone, minusone);

    /* Exponentiate */
    me.pow2(r, mm1);

    /* Now r == 2^mm1 (mod m) */

    /* Returns true if r1 == 1 (mod m) or if r1 == -1 (mod m) or if
       one of r1^(2^1), r1^(2^2), ..., r1^(2^(po2-1)) == -1 (mod m),
       zero otherwise. Requires -1 (mod m) in minusone. */
    if (me.is1(r) || me.equal(r, minusone))
        return true;

    for (int i = 1; i < po2; i++) {
        me.sqr(r, r);
        if (me.equal(r, minusone))
            return true;
    }

    return false;
}

/* This implements the "strong quadratic "V" Lucas pseudoprimality test",
 * with Q=1. I checked manually that none of the Miller-Rabin strong
 * pseudoprimes also pass this test up to 2^64
 */
template <typename layer>
bool arithxx_details::api<layer>::is_strong_lucas_pseudoprime() const
{
    auto const & me = downcast();

    /* Find the sequence parameter P */
    if (Integer(me.m) == 2) return true;
    if (!(Integer(me.m) & 1)) return false;

    Residue P = me(3);
    Residue D = me(5);
    for(int i = 0 ; me.jacobi(D) != -1 ; i++) {
        me.add1(P, P);
        me.add(D, D, P);
        me.add(D, D, P);
        me.add(D, D, P);
        me.add(D, D, P);
        me.add1(P, P);
        if (i == 20) {
            /* N might be a square, in which case we're going to loop
             * forever here. */
            if (mpz_perfect_square_p(cxx_mpz(Integer(me.m))))
                return false;
            // return mpz_probab_prime_p(N, 2);
        }
    }

    /* Compute the sequence and decide */
    Integer ell = me.getmod() + 1;
    int k = 0;
    for( ; (ell & 1) == 0 ; k++, ell>>=1) ;
    ASSERT_ALWAYS(k);
    Residue v(me);
    me.V(v, nullptr, P, ell);
    /* if v is already -2 or +2, things are not going to change */
    Residue const two = me(2);
    Residue mtwo(me);
    me.neg(mtwo, two);
    if (me.equal(v, two) || me.equal(v, mtwo))
        return true;
    for(int i = 0 ; i < k ; i++) {
        if (me.is0(v))
            return true;
        /* invariant: v = v_n = alpha^n + beta^n is neither 2 nor -2,
        */
        me.V_dbl(v, v, two); // v_{2n} = v_n^2-2
        /* Here, if N is prime then v_{2n} can't be +2 because that would
         * mean that v_n^2=4 and thus v_n == +2 or -2, which should both
         * have been checked earlier. So v_{2n} == 2 is a sign of
         * compositeness. Likewise, v_{2n}=-2 can only happen if
         * v_n^2=0. Since we checked v_n==0 before, that would mean
         * that n is not squarefree, and thus not prime.
         */
        if (me.equal(v, two) || me.equal(v, mtwo))
            return false;
    }
    return me.equal(v, two);
}

template <typename layer>
bool arithxx_details::api<layer>::is_prime() const
{
    auto const & me = downcast();

    if (!me.is_strong_pseudoprime_base2())
        return false;

    if (me.sprp2_is_enough())
        return true;

    return me.is_strong_lucas_pseudoprime();
}


#endif	/* UTILS_ARITHXX_API_IMPL_HPP_ */
