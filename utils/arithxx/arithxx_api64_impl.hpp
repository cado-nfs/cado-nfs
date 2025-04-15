#ifndef UTILS_ARITHXX_API64_IMPL_HPP_
#define UTILS_ARITHXX_API64_IMPL_HPP_

#include <cstddef>
#include <cstdint>

#include "arithxx_common.hpp"
#include "arithxx_api64.hpp"
#include "u64arith.h"
#include "macros.h"

/* {{{ sprp, sprp2, and isprime */

/* Returns 1 if m is a strong probable prime wrt base b, 0 otherwise.
   We assume m is odd. */
template <typename layer>
bool arithxx_details::api64<layer>::sprp(Residue const & b) const
{
    auto const & me = downcast();
    Residue r1(me), minusone(me);
    int i = 0, po2 = 1;
    uint64_t mm1;

    mm1 = me.getmod_u64();

    /* Set mm1 to the odd part of m-1 */
    mm1 = (mm1 - 1) >> 1;
    while (mm1 % 2 == 0) {
        po2++;
        mm1 >>= 1;
    }
    /* Hence, m-1 = mm1 * 2^po2 */

    me.set1(minusone);
    me.neg(minusone, minusone);

    /* Exponentiate */
    me.pow(r1, b, mm1);

    /* Now r1 == b^mm1 (mod m) */
#ifdef PARI
    printf("(Mod(%lu,%lu) ^ %lu) == %lu /* PARI */\n", get_u64(b, m),
           getmod_u64(m), mm1, get_u64(r1, m));
#endif

    i = me.find_minus1(r1, minusone, po2);

    return i;
}

/* Returns 1 if m is a strong probable prime wrt base 2, 0 otherwise.
   Assumes m > 1 and is odd.
 */
template <typename layer>
bool arithxx_details::api64<layer>::sprp2() const
{
    auto const & me = downcast();
    Residue r(me), minusone(me);
    int i = 0, po2 = 1;
    uint64_t mm1;

    mm1 = me.getmod_u64();

    /* If m == 1,7 (mod 8), then 2 is a quadratic residue, and we must find
       -1 with one less squaring. This does not reduce the number of
       pseudo-primes because strong pseudo-primes are also Euler pseudo-primes,
       but makes identifying composites a little faster on average. */
    if (mm1 % 8 == 1 || mm1 % 8 == 7)
        po2--;

    /* Set mm1 to the odd part of m-1 */
    mm1 = (mm1 - 1) >> 1;
    while (mm1 % 2 == 0) {
        po2++;
        mm1 >>= 1;
    }
    /* Hence, m-1 = mm1 * 2^po2 */

    me.set1(minusone);
    me.neg(minusone, minusone);

    /* Exponentiate */
    me.pow2(r, mm1);

    /* Now r == b^mm1 (mod m) */
#ifdef PARI
    printf("(Mod(2,%lu) ^ %lu) == %lu /* PARI */\n", getmod_u64(m), mm1,
           get_u64(r, m));
#endif

    i = me.find_minus1(r, minusone, po2);

    return i;
}

template <typename layer>
bool arithxx_details::api64<layer>::isprime() const
{
    auto const & me = downcast();
    Residue b(me), minusone(me), r1(me);
    uint64_t const n = me.getmod_u64();
    uint64_t mm1;
    int r = 0, po2;

#ifdef PARI
    auto dummy =
        call_dtor([&]() { printf("isprime(%lu) == %d /* PARI */\n", n, r); });
#endif

    if (n == 1)
        return false;

    if (n % 2 == 0) {
        r = (n == 2);
        return r;
    }

    /* Set mm1 to the odd part of m-1 */
    mm1 = n - 1;
    po2 = u64arith_ctz(mm1);
    mm1 >>= po2;

    me.set1(minusone);
    me.neg(minusone, minusone);

    /* Do base 2 SPRP test */
    me.pow2(r1, mm1); /* r = 2^mm1 mod m */
    /* If n is prime and 1 or 7 (mod 8), then 2 is a square (mod n)
       and one less squaring must suffice. This does not strengthen the
       test but saves one squaring for composite input */
    if (n % 8 == 7) {
        if (!me.is1(r1))
            return r;
    } else if (!me.find_minus1(r1, minusone, po2 - ((n % 8 == 1) ? 1 : 0)))
        return r; /* Not prime */

    if (n < 2047) {
        r = 1;
        return r;
    }

    /* Base 3 is poor at identifying composites == 1 (mod 3), but good at
       identifying composites == 2 (mod 3). Thus we use it only for 2 (mod 3) */
    if (n % 3 == 1) {
        me.set1(b);
        me.add(b, b, b);
        me.add(b, b, b);
        me.add(b, b, b);
        me.add(b, b, minusone); /* b = 7 */
        me.pow(r1, b, mm1);     /* r = 7^mm1 mod m */
        if (!me.find_minus1(r1, minusone, po2))
            return r; /* Not prime */

        if (n < 2269093) {
            r = (n != 314821);
            return r;
        }

        /* b is still 7 here */
        me.add(b, b, b);        /* 14 */
        me.sub(b, b, minusone); /* 15 */
        me.add(b, b, b);        /* 30 */
        me.add(b, b, b);        /* 60 */
        me.sub(b, b, minusone); /* 61 */
        me.pow(r1, b, mm1);     /* r = 61^mm1 mod m */
        if (!me.find_minus1(r1, minusone, po2))
            return r; /* Not prime */

        if (n != UINT64_C(4759123141) && n != UINT64_C(8411807377) &&
            n < UINT64_C(11207066041)) {
            r = 1;
            return r;
        }

        me.set1(b);
        me.add(b, b, b);
        me.add(b, b, b);
        me.sub(b, b, minusone); /* b = 5 */
        me.pow(r1, b, mm1);     /* r = 5^mm1 mod m */
        if (!me.find_minus1(r1, minusone, po2))
            return r; /* Not prime */

        /* These are the base 5,7,61 SPSP < 10^13 and n == 1 (mod 3) */
        r = (n != UINT64_C(30926647201) && n != UINT64_C(45821738881) &&
             n != UINT64_C(74359744201) && n != UINT64_C(90528271681) &&
             n != UINT64_C(110330267041) && n != UINT64_C(373303331521) &&
             n != UINT64_C(440478111067) && n != UINT64_C(1436309367751) &&
             n != UINT64_C(1437328758421) && n != UINT64_C(1858903385041) &&
             n != UINT64_C(4897239482521) && n != UINT64_C(5026103290981) &&
             n != UINT64_C(5219055617887) && n != UINT64_C(5660137043641) &&
             n != UINT64_C(6385803726241));
    } else {
        /* Case n % 3 == 0, 2 */

        me.pow3(r1, mm1); /* r = 3^mm1 mod m */
        if (!me.find_minus1(r1, minusone, po2))
            return r; /* Not prime */

        if (n < UINT64_C(102690677) && n != UINT64_C(5173601) &&
            n != UINT64_C(16070429) && n != UINT64_C(54029741)) {
            r = 1;
            return r;
        }

        me.set1(b);
        me.add(b, b, b);
        me.add(b, b, b);
        me.sub(b, b, minusone); /* b = 5 */
        me.pow(r1, b, mm1);     /* r = 5^mm1 mod m */
        if (!me.find_minus1(r1, minusone, po2))
            return r; /* Not prime */

        /* These are the base 3,5 SPSP < 10^13 with n == 2 (mod 3) */
        r = (n != UINT64_C(244970876021) && n != UINT64_C(405439595861) &&
             n != UINT64_C(1566655993781) && n != UINT64_C(3857382025841) &&
             n != UINT64_C(4074652846961) && n != UINT64_C(5783688565841));
    }
    return r;
}

/* }}} */

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
template <typename value_type>
inline void
arithxx_details::api64<layer>::pow2_oneWord(
    value_type mask, const value_type word, Residue &t) const
{
    auto const & me = downcast();
    while (mask > 0) {
        me.sqr (t, t);
        if (word & mask) {
            me.add (t, t, t);
        }
        mask >>= 1;
    }
}

/* Compute r = 2^e. Here, e is a uint64_t */
template <typename layer>
inline void
arithxx_details::api64<layer>::pow2(Residue &r, const uint64_t e) const
{
    auto const & me = downcast();
    uint64_t mask;
    Residue t(me);

    if (e == 0)
    {
        me.set1 (r);
        return;
    }

    mask = (UINT64_C(1) << 63) >> u64arith_clz (e);

    me.set1 (t);
    me.add (t, t, t);
    mask >>= 1;

    me.pow2_oneWord(mask, e, t);

    me.set (r, t);
}


/* Compute r = 3^e. Here, e is a uint64_t */
template <typename layer>
inline void
arithxx_details::api64<layer>::pow3 (
        Residue &r, const uint64_t e) const
{
    auto const & me = downcast();
    uint64_t mask;
    Residue t(me), u(me);

    if (e == 0)
    {
        me.set1 (r);
        return;
    }

    mask = (UINT64_C(1) << 63) >> u64arith_clz (e);

    me.set1 (u);
    me.add (t, u, u);
    me.add (t, t, u);
    mask >>= 1;

    while (mask > 0)
    {
        me.sqr (t, t);
        me.add (u, t, t);
        me.add (u, u, t);
        if (e & mask)
            me.set (t, u);
        mask >>= 1;
    }
    me.set (r, t);
}


/* Computes 2^e (mod m), where e is a multiple precision integer.
   Requires e != 0. The value of 2 in Montgomery representation
   (i.e. 2*2^w (mod m) must be passed. */

template <typename layer>
inline void
arithxx_details::api64<layer>::pow2 (
        Residue &r, const uint64_t *e, const size_t e_nrwords) const
{
    auto const & me = downcast();
    Residue t(me);
    uint64_t mask;
    int i = e_nrwords;

    while (i > 0 && e[i - 1] == 0)
        i--;

    if (i == 0) {
        me.set1(r);
        return;
    }
    i--;

    mask = (UINT64_C(1) << 63) >> u64arith_clz (e[i]);

    me.set1 (t);
    me.add (t, t, t);
    mask >>= 1;

    for ( ; i >= 0; i--) {
        me.pow2_oneWord(mask, e[i], t);
        mask = (UINT64_C(1) << 63);
    }

    me.set (r, t);
}

template <typename layer>
inline void
arithxx_details::api64<layer>::pow2(
        Residue &r, const Integer &e) const
{
    auto const & me = downcast();
    size_t i = e.size_in_words();
    const int bits = Integer::word_bits;
    auto const msb = (typename Integer::value_type) 1 << (bits-1);
    typename Integer::value_type word, mask;

    while (i > 0 && e.getWord(i - 1) == 0)
        i--;

    if (i == 0) {
        me.set1(r);
        return;
    }

    word = e.getWord(i - 1);
    mask = msb >> u64arith_clz (word);

    me.set1 (r);
    me.add (r, r, r);
    mask >>= 1;

    for ( ; i > 0; i--) {
        word = e.getWord(i - 1);
        me.pow2_oneWord(mask, word, r);
        mask = msb;
    }
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
