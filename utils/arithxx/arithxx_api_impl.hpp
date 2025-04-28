#ifndef CADO_UTILS_ARITHXX_API_IMPL_HPP
#define CADO_UTILS_ARITHXX_API_IMPL_HPP

#include <cstdint>
#include <cstddef>

#include <algorithm>

#include "arithxx_common.hpp"
#include "u64arith.h"
#include "cado_math_aux.hpp"
#include "macros.h"

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
    int k = ell.ctz();
    ASSERT_ALWAYS(k);
    ell >>= k;

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

/* Modulus::divn */
/* Division by small integer n.
 * Returns 1 if n is invertible modulo m, 0 if not.
 */

template<typename layer>
template<int n>
bool arithxx_details::api<layer>::divn(Residue & r, Residue const & a) const
{
    using namespace cado_math_aux;
    /*
     * inv_n contains -1/i (mod n) if i is coprime to n, or 0 if i is not
     * coprime to n, for 0 <= i < n c = n^(-1) (mod word base)
     */
    constexpr auto const & inv_n = array_map<n, minus_inverse_mod_n>::value;

    auto const & me = downcast();

    uint64_t const an = a.r.template mod_n<n>();
    uint64_t const mn = me.m.template mod_n<n>();

    if (inv_n[mn] == 0)
        return false;

    Residue t = a;

    /* Make t[1]:t[0] == a+km (mod w^2) with a+km divisible by n */
    /* We want a+km == 0 (mod n), so k = -a*m^{-1} (mod n) */
    uint64_t k = (inv_n[mn] * an) % n;
    ASSERT_ALWAYS((an + k * mn) % n == 0);

    if (k == 0) {
        /* nothing to do */
    } else if (k == 1) {
        t.r += me.m;
        // u64arith_add_2_2(t.r.data(), t.r.data() + 1, me.m[0], me.m[1]);
    } else if (k == 2) {
        t.r += me.m;
        t.r += me.m;
        // u64arith_add_2_2(t.r.data(), t.r.data() + 1, me.m[0], me.m[1]);
        // u64arith_add_2_2(t.r.data(), t.r.data() + 1, me.m[0], me.m[1]);
    } else {
        t.r += me.m * k;
        /*
        u64arith_mul_1_1_2(t.r.data(), t.r.data() + 1, me.m[0], k);
        t.r[1] += me.m[1] * k;
        u64arith_add_2_2(t.r.data(), t.r.data() + 1, a.r[0], a.r[1]);
        */
    }

    // now t is congruent to a multiple of n modulo w^2. We might have
    // encountered a carry in the previous additions, though (the sum can
    // be as large as k times the max modulus).

    constexpr bool has_carry = !multiplier_fits_in_overflow_bits<n, layer>::value;
    typedef at_most<layer::mul_c_cutoff::value> chooser_mul;

    r.r = Integer::template reduce_multiple<n, has_carry, chooser_mul>(t.r);

#ifdef WANT_ASSERT_EXPENSIVE
    {
        uint64_t i;
        me.set(t, r);
        for (i = 1; i < n; i++)
            me.add(t, t, r);
        ASSERT_EXPENSIVE(me.equal(t, a));
    }
#endif

    return true;
}

/* the mod 3 case is a bit special because there ain't many choices for
 * the residues, and in the end we can save a % operation.
 */
template<typename layer>
bool arithxx_details::api<layer>::div3(Residue & r, Residue const & a) const
{
    constexpr int n = 3;

    using namespace cado_math_aux;

    auto const & me = downcast();

    uint64_t const an = a.r.template mod_n<n>();
    uint64_t const mn = me.m.template mod_n<n>();

    if (mn == 0)
        return false;

    Residue t = a;

    if (an != 0) {
        if (an + mn == 3) {
            t.r += me.m;
        } else {
            t.r += me.m;
            t.r += me.m;
        }
    }

    // now t is congruent to a multiple of n modulo w^2. We might have
    // encountered a carry in the previous additions, though.

    constexpr bool has_carry = !multiplier_fits_in_overflow_bits<n, layer>::value;
    typedef at_most<layer::mul_c_cutoff::value> chooser_mul;

    r.r = Integer::template reduce_multiple<n, has_carry, chooser_mul>(t.r);

#ifdef WANT_ASSERT_EXPENSIVE
    {
        uint64_t i;
        me.set(t, r);
        for (i = 1; i < n; i++)
            me.add(t, t, r);
        ASSERT_EXPENSIVE(me.equal(t, a));
    }
#endif

    return true;
}

template <typename layer> bool arithxx_details::api<layer>::div5(Residue & r, Residue const & a) const { return divn<5>(r, a); }
template <typename layer> bool arithxx_details::api<layer>::div7(Residue & r, Residue const & a) const { return divn<7>(r, a); }
template <typename layer> bool arithxx_details::api<layer>::div11(Residue & r, Residue const & a) const { return divn<11>(r, a); }
template <typename layer> bool arithxx_details::api<layer>::div13(Residue & r, Residue const & a) const { return divn<13>(r, a); }

template<typename layer>
int arithxx_details::api<layer>::jacobi(Residue const & a_par) const
{
    auto const & me = downcast();

    /* If we're in Montgomery form, then the modulus is odd and the
     * number of bits is even, so that the Montgomery scaling does not
     * change the squareness, so we can pull a_par.r without converting.
     * If we're not, then we just pull the Integer from the Residue.
     */
    Integer x = a_par.r;
    Integer mm = me.m;

    /* Jacobi(0,1) = 1, Jacobi(0,>1) = 0 */
    if (x == 0)
        return mm == 1 ? 1 : 0;

    ASSERT(x < mm);
    ASSERT(mm & 1);

    unsigned int j, s;

    j = x.ctz();
    x >>= j;
    s = (j << 1) & (mm[0] ^ (mm[0] >> 1));
    /* If we divide by an odd power of 2, and 2 is a QNR, flip sign */
    /* 2 is a QNR (mod mm) iff mm = 3,5 (mod 8)
     * mm = 1 = 001b:   1
     * mm = 3 = 011b:  -1
     * mm = 5 = 101b:  -1
     * mm = 7 = 111b:   1
     * Hence we can store in s the exponent of -1, i.e., s=0 for jacobi()=1
     * and s=1 for jacobi()=-1, and update s ^= (mm>>1) & (mm>>2) & 1.
     * We can do the &1 at the very end.
     * 
     * Small optimization: we store the exponent of -1 in the second bit
     * of s.  The s ^= ((j<<1) & (mm ^ (mm>>1))) still needs 2 shift but
     * one of them can be done with LEA, and f = s ^ (x&mm) needs no
     * shift.
     */

    while (x > 1) {
        /* Here, x < mm, x and mm are odd */

        /* Implicitly swap by reversing roles of x and mm in next loop */
        /* Flip sign if both are 3 (mod 4) */
        s ^= (x[0] & mm[0]);

        /* Make mm<x by subtracting and shifting */
        do {
            mm -= x;    /* Difference is even */
            if (mm == 0)
                break;
            /* Make odd again */
            j = mm.ctz();
            s ^= (j << 1) & (x[0] ^ (x[0] >> 1));
            mm >>= j;
        } while (mm >= x);

        if (mm <= 1) {
            x = mm;
            break;
        }

        /* Flip sign if both are 3 (mod 4) */
        /* Implicitly swap again */
        s ^= (x[0] & mm[0]);

        /* Make x<mm by subtracting and shifting */
        do {
            x -= mm; /* Difference is even */
            if (x == 0)
                break;
            /* Make odd again */
            j = x.ctz();
            s ^= (j << 1) & (mm[0] ^ (mm[0] >> 1));
            x >>= j;
        } while (x >= mm);
    }

    if (x == 0)
        return 0;
    return (s & 2) ? -1 : 1;
}

template<typename layer>
void
arithxx_details::api<layer>::gcd(Integer & r, const Residue & A) const
{
    auto const & me = downcast();

    if (me.is0(A)) {
        r = me.getmod();
        return;
    }

    Integer a = A.r, b = me.m;

    unsigned int s = 0;
    if (layer::even_moduli_allowed::value) {
        const unsigned int ja = a.ctz();
        const unsigned int jb = b.ctz();
        a >>= ja;
        b >>= jb;
        s = std::min(ja, jb);
        if (b < a)
            std::swap(a, b);
    }

    if (s + layer::overflow_bits::value < 2) {
        if (b.high_word() & uint64_t(-(int64_t(1) << 62))) {
            /* make a odd. */
            a >>= a.ctz();
            do {
                ASSERT_EXPENSIVE(a[0] % 2 == 1);
                ASSERT_EXPENSIVE(b[0] % 2 == 1);
                ASSERT_EXPENSIVE(a < b);
                b -= a;
                if (b == 0) {
                    r = a << s;
                    return;
                }
                b >>= b.ctz();
                if (b < a)
                    std::swap(a, b);
            } while (b.high_word() & uint64_t(-(int64_t(1)<<62)));
        }
    }

    /* Each transformation step changes the pair (a,b) into a pair
     * (a',b') such that max(|a'|, |b'|) <= max(|a|, |b|), provided that
     * the inputs A.r and me.m are both < 2^126.
     *
     * Proof:
     * |b±a| <= 2 max(|a|, |b|) < 2^127, and this ± result is followed by
     * a right shift of at least one bit. Since the bound is less than
     * 2^127, then b±a is correctly represented with the signed 2-word
     * representation, and the signed right shift is correct. In
     * particular, the bound is preserved.
     */
    while (a != 0) {
        /* Make a odd */
        a.signed_shift_right(int(a.ctz()));

        /* Try to make the low two bits of b[0] zero */
        ASSERT_EXPENSIVE(a[0] % 2 == 1);
        ASSERT_EXPENSIVE(b[0] % 2 == 1);
        if ((a[0] ^ b[0]) & 2)
            b += a;
        else
            b -= a;

        if (b == 0) {
            r = a.sign_bit() ? -a : a;
            r <<= s;
            return;
        }

        /* Make b odd */
        b.signed_shift_right(int(b.ctz()));
        ASSERT_EXPENSIVE(a[0] % 2 == 1);
        ASSERT_EXPENSIVE(b[0] % 2 == 1);

        if ((a[0] ^ b[0]) & 2)
            a += b;
        else
            a -= b;
    }

    r = b.sign_bit() ? -b : b;
    r <<= s;
}


#endif	/* UTILS_ARITHXX_API_IMPL_HPP_ */
