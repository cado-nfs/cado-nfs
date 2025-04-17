#ifndef UTILS_ARITHXX_API64_IMPL_HPP_
#define UTILS_ARITHXX_API64_IMPL_HPP_

#include <cstdint>

#include "arithxx_common.hpp"
#include "arithxx_api64.hpp"
#include "u64arith.h"
#include "macros.h"

/*  */

// {{{ div{3,5,7,11,13}
template <typename layer>
bool arithxx_details::api64<layer>::div3(Residue & r,
                                                        Residue const & a) const
{
    auto const & me = downcast();
    uint64_t const a3 = a.r % 3;
    uint64_t ml, m3;
    Residue t(me);

    me.assertValid(a);

    ml = me.getmod_u64();
    m3 = ml % 3;
    if (m3 == 0)
        return false;

    me.set(t, a);
    if (a3 != 0) {
        if (a3 + m3 == 3) /* Hence a3 == 1, m3 == 2 or a3 == 2, m3 == 1 */
            t.r += ml;
        else /* a3 == 1, m3 == 1 or a3 == 2, m3 == 2 */
            t.r += 2 * ml;
    }

    /* Now t == a+k*m (mod 2^64) so that a+k*m is divisible by 3.
       (a+k*m)/3 < 2^64, so doing a division (mod 2^64) produces the
       correct result. */

    t.r *= UINT64_C(0xAAAAAAAAAAAAAAAB);

#ifdef WANT_ASSERT_EXPENSIVE
    Residue c(me);
    me.assertValid(t);
    me.add(c, t, t);
    me.add(c, c, t);
    ASSERT_EXPENSIVE(me.equal(c, a));
#endif

    me.set(r, t);

    return true;
}

template <typename layer>
bool arithxx_details::api64<layer>::div5(Residue & r,
                                                        Residue const & a) const
{
    auto const & me = downcast();
    uint64_t ml, m5, k;
    Residue t(me);
    uint64_t const a5 = a.r % 5;
    uint64_t const inv5[5] = {0, 4, 2, 3, 1}; /* inv5[i] = -1/i (mod 5) */

    me.assertValid(a);

    ml = me.getmod_u64();
    m5 = ml % 5;
    if (m5 == 0)
        return false;

    me.set(t, a);
    if (a5 != 0) {
        /* We want a+km == 0 (mod 5), so k = -a*m^{-1} (mod 5) */
        k = (a5 * inv5[m5]) % 5;
        ASSERT_EXPENSIVE((k * m5 + a5) % 5 == 0);
        t.r = a.r + k * ml;
    }

    /* Now t == a+k*m (mod 2^w) so that a+k*m is divisible by 5.
       (a+k*m)/5 < 2^w, so doing a division (mod 2^w) produces the
       correct result. */

    t.r *= UINT64_C(0xcccccccccccccccd); /* 1/5 (mod 2^64) */

#ifdef WANT_ASSERT_EXPENSIVE
    Residue c(me);
    me.assertValid(t);
    me.add(c, t, t);
    me.add(c, c, c);
    me.add(c, c, t);
    ASSERT_EXPENSIVE(me.equal(c, a));
#endif

    me.set(r, t);

    return true;
}

template <typename layer>
bool arithxx_details::api64<layer>::div7(Residue & r,
                                                        Residue const & a) const
{
    auto const & me = downcast();
    uint64_t ml, m7, k;
    Residue t(me);
    uint64_t const a7 = a.r % 7;
    uint64_t const inv7[7] = {0, 6, 3, 2, 5, 4, 1}; /* inv7[i] = -1/i (mod 7) */

    me.assertValid(a);

    ml = me.getmod_u64();
    m7 = ml % 7;
    if (m7 == 0)
        return false;

    me.set(t, a);
    if (a7 != 0) {
        /* We want a+km == 0 (mod 7), so k = -a*m^{-1} (mod 7) */
        k = (a7 * inv7[m7]) % 7;
        ASSERT_EXPENSIVE((k * m7 + a7) % 7 == 0);
        t.r = a.r + k * ml;
    }

    /* Now t == a+k*m (mod 2^w) so that a+k*m is divisible by 7.
       (a+k*m)/7 < 2^w, so doing a division (mod 2^w) produces the
       correct result. */

    t.r *= UINT64_C(0x6db6db6db6db6db7); /* 1/7 (mod 2^64) */

#ifdef WANT_ASSERT_EXPENSIVE
    Residue c(me);
    me.assertValid(t);
    me.add(c, t, t);
    me.add(c, c, c);
    me.add(c, c, c);
    me.sub(c, c, t);
    ASSERT_EXPENSIVE(me.equal(c, a));
#endif

    me.set(r, t);

    return true;
}

template <typename layer>
bool arithxx_details::api64<layer>::div11(
    Residue & r, Residue const & a) const
{
    auto const & me = downcast();
    uint64_t ml, m11, k;
    Residue t(me);
    uint64_t const a11 = a.r % 11;
    /* inv11[i] = -1/i (mod 11) */
    uint64_t const inv11[11] = {0, 10, 5, 7, 8, 2, 9, 3, 4, 6, 1};

    me.assertValid(a);

    ml = me.getmod_u64();
    m11 = ml % 11;
    if (m11 == 0)
        return false;

    me.set(t, a);
    if (a11 != 0) {
        /* We want a+km == 0 (mod 11), so k = -a*m^{-1} (mod 11) */
        k = (a11 * inv11[m11]) % 11;
        ASSERT_EXPENSIVE((k * m11 + a11) % 11 == 0);
        t.r = a.r + k * ml;
    }

    /* Now t == a+k*m (mod 2^w) so that a+k*m is divisible by 11.
       (a+k*m)/11 < 2^w, so doing a division (mod 2^w) produces the
       correct result. */

    t.r *= UINT64_C(0x2e8ba2e8ba2e8ba3); /* 1/11 (mod 2^64) */

#ifdef WANT_ASSERT_EXPENSIVE
    Residue c(me);
    me.assertValid(t);
    me.add(c, t, t); /* c = 2*t */
    me.add(c, c, c); /* c = 4*t */
    me.add(c, c, t); /* c = 5*t */
    me.add(c, c, c); /* c = 10*t */
    me.add(c, c, t); /* c = 11*t */
    ASSERT_EXPENSIVE(me.equal(c, a));
#endif

    me.set(r, t);

    return true;
}

template <typename layer>
bool arithxx_details::api64<layer>::div13(
    Residue & r, Residue const & a) const
{
    auto const & me = downcast();
    uint64_t ml, m13, k;
    Residue t(me);
    uint64_t const a13 = a.r % 13;
    /* inv13[i] = -1/i (mod 13) */
    uint64_t const inv13[13] = {0, 12, 6, 4, 3, 5, 2, 11, 8, 10, 9, 7, 1};

    me.assertValid(a);

    ml = me.getmod_u64();
    m13 = ml % 13;
    if (m13 == 0)
        return false;

    me.set(t, a);
    if (a13 != 0) {
        /* We want a+km == 0 (mod 13), so k = -a*m^{-1} (mod 13) */
        k = (a13 * inv13[m13]) % 13;
        ASSERT_EXPENSIVE((k * m13 + a13) % 13 == 0);
        t.r = a.r + k * ml;
    }

    /* Now t == a+k*m (mod 2^w) so that a+k*m is divisible by 13.
       (a+k*m)/13 < 2^w, so doing a division (mod 2^w) produces the
       correct result. */

    t.r *= UINT64_C(0x4ec4ec4ec4ec4ec5); /* 1/13 (mod 2^64) */

#ifdef WANT_ASSERT_EXPENSIVE
    Residue c(me);
    me.assertValid(t);
    me.add(c, t, t); /* c = 2*t */
    me.add(c, c, t); /* c = 3*t */
    me.add(c, c, c); /* c = 6*t */
    me.add(c, c, c); /* c = 12*t */
    me.add(c, c, t); /* c = 13*t */
    ASSERT_EXPENSIVE(me.equal(c, a));
#endif

    me.set(r, t);

    return true;
}
// }}}

template<typename layer>
void
arithxx_details::api64<layer>::gcd (Integer &g, const Residue &r) const
{
    auto const & me = downcast();
    uint64_t t;

    auto a = uint64_t(r.r); /* This works the same for "a" in plain or Montgomery
                         representation */
    uint64_t b = me.getmod_u64 ();
    /* ASSERT (a < b); Should we require this? */
    ASSERT (b > 0);

    if (a >= b)
        a %= b;

    while (a > 0)
    {
        /* Here 0 < a < b */
        t = b % a;
        b = a;
        a = t;
    }

    g = Integer(b);
}

template <typename layer>
int
arithxx_details::api64<layer>::jacobi (const Residue &a_par) const
{
    auto const & me = downcast();
    unsigned int s, j;

    /* Get residue in Montgomery form directly without converting */
    auto x = uint64_t(a_par.r);
    uint64_t mm = me.getmod_u64 ();
    if (x == 0) {
        return (mm == 1) ? 1 : 0; /* Jacobi(0,1) = 1, Jacobi(0,>1) = 0 */
    }
    ASSERT (x < mm);
    ASSERT(mm % 2 == 1);

    j = u64arith_ctz(x);
    x = x >> j;
    /* If we divide by an odd power of 2, and 2 is a QNR, flip sign */
    /* 2 is a QNR (mod mm) iff mm = 3,5 (mod 8)
       mm = 1 = 001b:   1
       mm = 3 = 011b:  -1
       mm = 5 = 101b:  -1
       mm = 7 = 111b:   1
       Hence we can store in s the exponent of -1, i.e., s=0 for jacobi()=1
       and s=1 for jacobi()=-1, and update s ^= (mm>>1) & (mm>>2) & 1.
       We can do the &1 at the very end.
       In fact, we store the exponent of -1 in the second bit of s.
       The s ^= ((j<<1) & (mm ^ (mm>>1))) still needs 2 shift but one of them can
       be done with LEA, and f = s ^ (x&mm) needs no shift */

    s = ((j<<1) & (mm ^ (mm>>1)));

    while (x > 1) {
        /* Here, x < mm, x and mm are odd */

        /* Implicitly swap by reversing roles of x and mm in next loop */
        /* Flip sign if both are 3 (mod 4) */
        s = s ^ (x&mm);

        /* Make mm smaller by subtracting and shifting */
        do {
            mm -= x; /* Difference is even */
            if (mm == 0)
                break;
            /* Make odd again */
            j = u64arith_ctz(mm);
            s ^= ((j<<1) & (x ^ (x>>1)));
            mm >>= j;
        } while (mm >= x);

        if (mm <= 1) {
            x = mm;
            break;
        }

        /* Flip sign if both are 3 (mod 4) */
        /* Implicitly swap again */
        s = s ^ (x&mm);

        /* Make x<mm by subtracting and shifting */
        do {
            x -= mm; /* Difference is even */
            if (x == 0)
                break;
            /* Make odd again */
            j = u64arith_ctz(x);
            s ^= ((j<<1) & (mm ^ (mm>>1)));
            x >>= j;
        } while (x >= mm);
    }

    if (x == 0)
        return 0;
    return ((s & 2) == 0) ? 1 : -1;
}

#endif	/* UTILS_ARITHXX_API64_IMPL_HPP_ */
