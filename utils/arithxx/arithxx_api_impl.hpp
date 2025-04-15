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

/* Returns 1 if r1 == 1 (mod m) or if r1 == -1 (mod m) or if
   one of r1^(2^1), r1^(2^2), ..., r1^(2^(po2-1)) == -1 (mod m),
   zero otherwise. Requires -1 (mod m) in minusone. */

template <typename layer>
bool arithxx_details::api<layer>::find_minus1(
    Residue & r1, Residue const & minusone, int const po2) const
{
    auto const & me = downcast();
    int i;

    if (me.is1(r1) || me.equal(r1, minusone))
        return true;

    for (i = 1; i < po2; i++) {
        me.sqr(r1, r1);
        if (me.equal(r1, minusone))
            break;
    }

    return i < po2;
}

/* Compute r = b^e. Here, e is a uint64_t */

template <typename layer>
void arithxx_details::api<layer>::pow(Residue & r,
                                                   Residue const & b,
                                                   uint64_t const e) const
{
    auto const & me = downcast();
    uint64_t mask;
    Residue t(me);

    if (e == 0) {
        me.set1(r);
        return;
    }

    /* Find highest set bit in e. */
    mask = (uint64_t(1) << 63) >> u64arith_clz(e);

    me.set(t, b);

    while ((mask >>= 1) > 0) {
        me.sqr(t, t);
        if (e & mask) {
            me.mul(t, t, b);
        }
    }
    /* Now e = 0, mask = 1, and r^mask * b^0 = r^mask is the result we want */
    me.set(r, t);
}

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
    uint64_t mask;
    Residue t(me);
    int i = e_nrwords;

    while (i > 0 && e[i - 1] == 0)
        i--;

    if (i == 0) {
        me.set1(r);
        return;
    }
    i--;

    /* Find highest set bit in e[i]. */
    mask = (uint64_t(1) << 63) >> u64arith_clz(e[i]);
    /* t = 1, so t^(mask/2) * b^e = t^mask * b^e  */

    /* Exponentiate */

    me.set(t, b); /* (r*b)^mask * b^(e-mask) = r^mask * b^e */
    mask >>= 1;

    for (; i >= 0; i--) {
        while (mask > 0) {
            me.sqr(t, t);
            if (e[i] & mask)
                me.mul(t, t, b);
            mask >>= 1; /* (r^2)^(mask/2) * b^e = r^mask * b^e */
        }
        mask = uint64_t(1) << 63;
    }
    me.set(r, t);
}

template <typename layer>
void arithxx_details::api<layer>::pow(Residue & r,
                                                   Residue const & b,
                                                   Integer const & e) const
{
    auto const & me = downcast();
    size_t i = e.size_in_words();
    int const bits = layer::Integer::word_bits;
    auto const msb = (typename Integer::value_type)1 << (bits - 1);
    typename Integer::value_type word, mask;
    Residue t(me);

    while (i > 0 && e.getWord(i - 1) == 0)
        i--;

    if (i == 0) {
        me.set1(r);
        return;
    }

    word = e.getWord(i - 1);
    mask = msb >> u64arith_clz(word);

    /* Exponentiate */

    me.set(t, b);
    mask >>= 1;

    for (; i > 0; i--) {
        word = e.getWord(i - 1);
        while (mask > 0) {
            me.sqr(t, t);
            if (word & mask)
                me.mul(t, t, b);
            mask >>= 1;
        }
        mask = msb;
    }
    me.set(r, t);
}

/* Compute r = V_k (b) and rp1 = V_{k+1} (b) if rp1 != NULL
 * r or rp1 can be the same variable as b.
 */
template <typename layer>
void arithxx_details::api<layer>::V(Residue & r, Residue * rp1,
                                                 Residue const & b,
                                                 uint64_t const k) const
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
        unsigned long mask;

        mask = ~(0UL);
        mask -= mask / 2; /* Now the most significant bit of i is set */
        while ((mask & k) == 0)
            mask >>= 1;

        /* Most significant bit of k is 1, do it outside the loop */
        me.set(t0, b);        /* starting value t0 = V_1 (b) = b */
        me.V_dbl(t1, b, two); /* starting value t1 = V_2 (b) */
        mask >>= 1;

        /* If the second most significant bit of k is 0, then we do the
         * iteration manually (to avoid to compute again V_2 (b)) As k >= 3, we
         * know that in this case k has at least 3 bits.
         */
        if (!(k & mask)) /* (t0,t1) <- (V_2 (b), V_3 (b)) */
        {
            me.set(t0, t1);
            me.V_dadd(t1, t1, b, b);
            mask >>= 1;
        }

        for (; mask > 1; mask >>= 1) /* t0 = V_j (b) and t1 = V_{j+1} (b) */
        {
            if (k & mask) /* j -> 2*j+1. Compute V_{2j+1} and V_{2j+2} */
            {
                me.V_dadd(t0, t1, t0, b);
                me.V_dbl(t1, t1, two);
            } else /* j -> 2*j. Compute V_{2j} and V_{2j+1} */
            {
                me.V_dadd(t1, t1, t0, b);
                me.V_dbl(t0, t0, two);
            }
        }

        /* Deal with least significant bit outside the loop */
        if (k & mask) {
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


#endif	/* UTILS_ARITHXX_API_IMPL_HPP_ */
