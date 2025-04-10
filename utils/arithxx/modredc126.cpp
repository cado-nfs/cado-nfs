#include "cado.h" // IWYU pragma: keep

#include <cstdint>
#include <cstdlib> // for abort

#include "macros.h"
#include "modredc126.hpp"
#include "u64arith.h"

#include "arithxx_common.hpp"

/* Only the .cpp source files that emit the non-inline symbols will
 * include this impl header file. So even though it does not look like
 * we're using it, in fact we are!  */
#include "arithxx_api_impl.hpp"      // IWYU pragma: keep

// scan-headers: stop here

/* {{{ implementation of fixed-exponent power npow (utility template function,
 * static) */

/* Simple addition chains for small multipliers, used in the powering
   functions with small bases. r and a may overlap, t must not overlap
   with anything. */
template <>
inline void arithxx_modredc126::Modulus::simple_mul<2>(Residue & r, Residue const & a,
                   Residue & t MAYBE_UNUSED) const
{
    add(r, a, a); /* r = 2*a */
}
template <>
inline void arithxx_modredc126::Modulus::simple_mul<3>(Residue & r, Residue const & a,
                   Residue & t MAYBE_UNUSED) const
{
    add(t, a, a); /* t = 2*a */
    add(r, t, a); /* r = 3*a */
}
template <>
inline void arithxx_modredc126::Modulus::simple_mul<5>(Residue & r, Residue const & a,
                   Residue & t MAYBE_UNUSED) const
{
    add(t, a, a); /* t = 2*a */
    add(t, t, t); /* t = 4*a */
    add(r, r, a); /* r = 5*a */
}
template <>
inline void arithxx_modredc126::Modulus::simple_mul<7>(Residue & r, Residue const & a,
                   Residue & t MAYBE_UNUSED) const
{
    add(t, a, a); /* r = 2*a */
    add(t, t, t); /* r = 4*a */
    add(t, t, t); /* r = 8*a */
    sub(r, t, a); /* r = 7*a */
}

template <int B, typename WordType>
inline void arithxx_modredc126::Modulus::npow_oneWord(
        WordType mask, WordType const word, Residue & t, Residue & u) const
{
    while (mask > 0) {
        sqr(t, t);
        if (word & mask) {
            simple_mul<B>(t, t, u);
        }
        mask >>= 1;
    }
}

/* Compute r = b^e, where b is a small integer, currently b=2,3,5,7 are
   implemented.  Here e is a multiple precision integer
   sum_{i=0}^{e_nrwords-1} e[i] * (machine word base)^i. e_nrwords must be
   minimal, i.e., either e_nrwords == 0 or e[e_nrwords - 1] != 0. */
template <int B>
inline void arithxx_modredc126::Modulus::npow(Residue & r, uint64_t const * e,
                        size_t const e_nrwords) const
{
    Modulus::Residue t(*this), u(*this);
    size_t i = e_nrwords;
    uint64_t mask;

    ASSERT(i == 0 || e[i - 1] != 0);

    if (i == 0) {
        set1(r);
        return;
    }

    set1(t);
    simple_mul<B>(t, t, u); /* t = b */

    mask = (UINT64_C(1) << 63) >> u64arith_clz(e[i - 1]);
    mask >>= 1;

    for (; i > 0; i--) {
        npow_oneWord<B>(mask, e[i - 1], t, u);
        mask = UINT64_C(1) << 63;
    }
    set(r, t);
}

/* Compute r = b^e, where b is a small integer, currently b=2,3,5,7 are
   implemented. Here, e is a uint64_t */
template <int B>
inline void arithxx_modredc126::Modulus::npow(Residue & r, uint64_t const e) const
{
    uint64_t mask;
    Residue t(*this), u(*this);

    if (e == 0) {
        set1(r);
        return;
    }

    set1(t);
    simple_mul<B>(t, t, u); /* t = b */

    mask = (uint64_t(1) << 63) >> u64arith_clz(e);
    ASSERT(e & mask);
    mask >>= 1;

    npow_oneWord<B>(mask, e, t, u);
    set(r, t);
}

template <int B>
inline void arithxx_modredc126::Modulus::npow(Residue & r, Integer const & e) const
{
    if (e.size_in_words() == 2) {
        uint64_t t[2];
        e.get(t, 2);
        npow<B>(r, t, e.size_in_words()); /* r = b^e mod m */
    } else if (e.size_in_words() <= 1) {
        npow<B>(r, e.getWord(0));
    } else {
        abort();
    }
}
// }}}

/* Compute r = 2^e mod m. Here, e is a uint64_t */
template<>
void arithxx_details::api<arithxx_modredc126>::pow2(Residue & r, uint64_t const e) const
{
    downcast().npow<2>(r, e);
}

/* Compute r = 2^e mod m.  Here e is a multiple precision integer
   sum_{i=0}^{e_nrwords-1} e[i] * (machine word base)^i */
template<>
void arithxx_details::api<arithxx_modredc126>::pow2(Residue & r, uint64_t const * e,
                   size_t const e_nrwords) const
{
    downcast().npow<2>(r, e, e_nrwords);
}

template<>
void arithxx_details::api<arithxx_modredc126>::pow2(Residue & r, Integer const & e) const
{
    downcast().npow<2>(r, e);
}

/* {{{ sprp, sprp2, and isprime */

/* Returns 1 if m is a strong probable prime wrt base b, 0 otherwise.
   We assume m is odd. */
template <>
bool arithxx_details::api<
    arithxx_modredc126>::sprp(Residue const & b) const
{
    auto const & me = downcast();
    Residue r(me), minusone(me);
    Integer mm1;
    int po2 = 0;

    me.getmod(mm1);

    if (mm1 == 1)
        return false;

    /* Let m-1 = 2^l * k, k odd. Set mm1 = k, po2 = l */
    mm1 -= 1;
    po2 = mm1.ctz();
    ASSERT_ALWAYS((unsigned int)po2 < Integer::max_bits);
    mm1 >>= po2;

    me.set1(minusone);
    me.neg(minusone, minusone);

    /* Exponentiate */
    me.pow(r, b, mm1);
    return me.find_minus1(r, minusone, po2);
}

/* Returns 1 if m is a strong probable prime wrt base 2, 0 otherwise.
   We assume m is odd. */
template <>
bool arithxx_details::api<arithxx_modredc126>::sprp2() const
{
    auto const & me = downcast();
    Residue r(me), minusone(me);
    int i = 0, po2 = 0;
    Integer mm1;

    me.getmod(mm1);
    if (mm1 == 1)
        return false;

    /* Let m-1 = 2^l * k, k odd. Set mm1 = k, po2 = l */
    --mm1;
    po2 = mm1.ctz();

    ASSERT_ALWAYS((unsigned int)po2 < Integer::max_bits);
    mm1 >>= po2;

    me.set1(minusone);
    me.neg(minusone, minusone);

    /* Exponentiate */
    me.pow2(r, mm1);
    i = find_minus1(r, minusone, po2);

    return i;
}

template <>
bool arithxx_details::api<arithxx_modredc126>::isprime() const
{
    auto const & me = downcast();
    Residue b(me), minusone(me), r1(me);
    Integer n, mm1;
    bool r = false;
    int po2 = 0, i;

    me.getmod(n);

    if (n == 1)
        return false;

    if (n.getWord(0) % 2 == 0)
        return n == 2;

    /* Set mm1 to the odd part of m-1 */
    mm1 = n - 1;
    po2 = mm1.ctz();
    ASSERT_ALWAYS((unsigned int)po2 < Integer::max_bits);
    mm1 >>= po2;

    me.set1(minusone);
    me.neg(minusone, minusone);

    /* Do base 2 SPRP test */
    me.pow2(r1, mm1);

    do {
        /* If n is prime and 1 or 7 (mod 8), then 2 is a square (mod n)
           and one fewer squarings must suffice. This does not strengthen the
           test but saves one squaring for composite input */
        if (n.getWord(0) % 8 == 7) {
            if (!me.is1(r1))
                break;
        } else if (!find_minus1(r1, minusone,
                                po2 - ((n.getWord(0) % 8 == 7) ? 1 : 0)))
            break; /* Not prime */

        /* Base 3 is poor at identifying composites == 1 (mod 3), but
         * good at identifying composites == 2 (mod 3). Thus we use it
         * only for 2 (mod 3) */
        i = n.getWord(0) % 3 + n.getWord(1) % 3;
        if (i == 1 || i == 4) {
            me.npow<7>(r1, mm1); /* r = 7^mm1 mod m */
            if (!find_minus1(r1, minusone, po2))
                break; /* Not prime */

            me.set_reduced(b, 61); /* Use addition chain? */
            pow(r1, b, mm1);    /* r = 61^mm1 mod m */
            if (!find_minus1(r1, minusone, po2))
                break; /* Not prime */

            me.npow<5>(r1, mm1); /* r = 5^mm1 mod m */
            if (!find_minus1(r1, minusone, po2))
                break; /* Not prime */

            /* These are the base 2,5,7,61 SPSP < 10^13 and n == 1 (mod 3) */

            r = n != UINT64_C(30926647201) && n != UINT64_C(45821738881) &&
                n != UINT64_C(74359744201) && n != UINT64_C(90528271681) &&
                n != UINT64_C(110330267041) && n != UINT64_C(373303331521) &&
                n != UINT64_C(440478111067) && n != UINT64_C(1436309367751) &&
                n != UINT64_C(1437328758421) && n != UINT64_C(1858903385041) &&
                n != UINT64_C(4897239482521) && n != UINT64_C(5026103290981) &&
                n != UINT64_C(5219055617887) && n != UINT64_C(5660137043641) &&
                n != UINT64_C(6385803726241);
        } else {
            /* Case n % 3 == 0, 2 */

            me.npow<3>(r1, mm1); /* r = 3^mm1 mod m */
            if (!find_minus1(r1, minusone, po2))
                break; /* Not prime */

            me.npow<5>(r1, mm1); /* r = 5^mm1 mod m */
            if (!find_minus1(r1, minusone, po2))
                break; /* Not prime */

            /* These are the base 2,3,5 SPSP < 10^13 and n == 2 (mod 3) */

            r = n != UINT64_C(244970876021) && n != UINT64_C(405439595861) &&
                n != UINT64_C(1566655993781) && n != UINT64_C(3857382025841) &&
                n != UINT64_C(4074652846961) && n != UINT64_C(5783688565841);
        }
    } while (false);

#if defined(PARI)
    printf("isprime(%lu) == %d /* PARI */\n", n, r);
#endif
    return r;
}
/* }}} */

/* {{{ Modulus::divn */
/* Division by small integer n, where (n-1)*m may overflow the most
   significant word. Returns 1 if n is invertible modulo m, 0 if not.

   w_mod_n is word base (e.g., 2^32 or  2^64) mod n
   inv_n contains -1/i (mod n) if i is coprime to n, or 0 if i is not
   coprime to n, for 0 <= i < n c = n^(-1) (mod word base)
   */

bool arithxx_modredc126::Modulus::divn(Residue & r, Residue const & a, uint64_t const n,
                          uint64_t const w_mod_n, uint64_t const * inv_n,
                          uint64_t const c) const
{
    uint64_t const an = ((a.r[1] % n) * w_mod_n + a.r[0] % n) % n;
    uint64_t const mn = ((m[1] % n) * w_mod_n + m[0] % n) % n;

    Residue t(*this), t2(*this);
    uint64_t k;
#ifdef WANT_ASSERT_EXPENSIVE
    Residue a_backup(*this);
#endif

    if (inv_n[mn] == 0)
        return false;

#ifdef WANT_ASSERT_EXPENSIVE
    set(a_backup, a);
#endif
    set(t, a);

    /* Make t[1]:t[0] == a+km (mod w^2) with a+km divisible by n */
    /* We want a+km == 0 (mod n), so k = -a*m^{-1} (mod n) */
    k = (inv_n[mn] * an) % n;
    ASSERT_EXPENSIVE((an + k * mn) % n == 0);
    u64arith_mul_1_1_2(&(t.r[0]), &(t.r[1]), m[0], k);
    t.r[1] += m[1] * k;
    u64arith_add_2_2(&(t.r[0]), &(t.r[1]), a.r[0], a.r[1]);

    /* We want r = (a+km)/n. */

    /* May overwrite a */
    r.r[0] = t.r[0] * c;

    /* r0 == (a+km)/n (mod w)
       (r1*w + r0) * n = (a+km)
       (r1*w + r0) * n == t (mod w^2)
       r1*w*n == t - n*r0 (mod w^2)
       t - n*r0 == 0 (mod w), thus
       r1*n == (t - n*r0)/w (mod w) */

    u64arith_mul_1_1_2(&(t2.r[0]), &(t2.r[1]), r.r[0], n);
    u64arith_sub_2_2(&(t.r[0]), &(t.r[1]), t2.r[0], t2.r[1]);
    ASSERT_EXPENSIVE(t.r[0] == 0);
    r.r[1] = t.r[1] * c;

#ifdef WANT_ASSERT_EXPENSIVE
    {
        uint64_t i;
        set(t, r);
        for (i = 1; i < n; i++)
            add(t, t, r);
        ASSERT_EXPENSIVE(equal(t, a_backup));
    }
#endif

    return true;
}
/* }}} */

// {{{ div{3,5,7,11,13}
/* Divide residue by 3. Returns 1 if division is possible, 0 otherwise.
   Assumes that a+3m does not overflow */
template <>
bool arithxx_details::api<
    arithxx_modredc126>::div3(Residue & r,
                                                      Residue const & a) const
{
    auto const & me = downcast();
    Residue t(me);
    uint64_t const a3 =
        (a.r[1] % 256 + a.r[1] / 256 + a.r[0] % 256 + a.r[0] / 256) % 3;
    uint64_t const m3 = (me.m[0] % 256 + me.m[0] / 256 + me.m[1] % 256 + me.m[1] / 256) % 3;
#ifdef WANT_ASSERT_EXPENSIVE
    Residue a_backup(*this);
    set(a_backup, a);
#endif

    if (m3 == 0)
        return false;

    me.set(t, a);

    /* Make t[1]:t[0] divisible by 3 */
    if (a3 != 0) {
        if (a3 + m3 == 3) /* Hence a3 == 1, m3 == 2 or a3 == 2, m3 == 1 */
        {
            u64arith_add_2_2(&(t.r[0]), &(t.r[1]), me.m[0], me.m[1]);
        } else /* a3 == 1, m3 == 1 or a3 == 2, m3 == 2 */
        {
            u64arith_add_2_2(&(t.r[0]), &(t.r[1]), me.m[0], me.m[1]);
            u64arith_add_2_2(&(t.r[0]), &(t.r[1]), me.m[0], me.m[1]);
        }

        /* Now t[1]:t[0] is divisible by 3 */
        ASSERT_EXPENSIVE((t.r[0] % 3 + t.r[1] % 3) % 3 == 0);
    }

    /* a = a1 * 2^w + a0, 3|a
       Let a = a' * 3 * 2^w + a'', a'' < 3 * 2^w.
       3 | a'', a'' / 3 < 2^w
       So a / 3 = a' * w + a'' / 3
       a' = trunc(a1 / 3)
       a'' = a0 * 3^{-1} (mod 2^w)
       Hence we get the correct result with one one-word multiplication
       and one one-word truncating division by a small constant.
    */

    r.r[1] = t.r[1] / 3;
    r.r[0] = t.r[0] * UINT64_C(0xaaaaaaaaaaaaaaab); /* 1/3 (mod 2^64) */

#ifdef WANT_ASSERT_EXPENSIVE
    add(t, r, r);
    add(t, t, r);
    ASSERT_EXPENSIVE(equal(a_backup, t));
#endif

    return true;
}

/* Divide residue by 5. Returns 1 if division is possible, 0 otherwise */
template <>
bool arithxx_details::api<
    arithxx_modredc126>::div5(Residue & r,
                                                      Residue const & a) const
{
    auto const & me = downcast();
    /* inv_5[i] = -1/i (mod 5) */
    uint64_t const inv_5[5] = {0, 4, 2, 3, 1};
    uint64_t const c = UINT64_C(0xcccccccccccccccd); /* 1/5 (mod 2^64) */

    return me.divn(r, a, 5, 1, inv_5, c);
}

/* Divide residue by 7. Returns 1 if division is possible, 0 otherwise */
template <>
bool arithxx_details::api<
    arithxx_modredc126>::div7(Residue & r,
                                                      Residue const & a) const
{
    auto const & me = downcast();
    /* inv_7[i] = -1/i (mod 7) */
    uint64_t const inv_7[7] = {0, 6, 3, 2, 5, 4, 1};
    uint64_t const c = UINT64_C(0x6db6db6db6db6db7); /* 1/7 (mod 2^64) */
    return me.divn(r, a, 7, 2, inv_7, c);
}

/* Divide residue by 11. Returns 1 if division is possible, 0 otherwise */
template <>
bool arithxx_details::api<
    arithxx_modredc126>::div11(Residue & r,
                                                       Residue const & a) const
{
    auto const & me = downcast();
    /* inv_11[i] = -1/i (mod 11) */
    uint64_t const inv_11[11] = {0, 10, 5, 7, 8, 2, 9, 3, 4, 6, 1};
    uint64_t const c = UINT64_C(0x2e8ba2e8ba2e8ba3); /* 1/11 (mod 2^64) */
    return me.divn(r, a, 11, 5, inv_11, c);
}

/* Divide residue by 13. Returns 1 if division is possible, 0 otherwise */
template <>
bool arithxx_details::api<
    arithxx_modredc126>::div13(Residue & r,
                                                       Residue const & a) const
{
    auto const & me = downcast();
    /* inv_13[i] = -1/i (mod 13) */
    uint64_t const inv_13[13] = {0, 12, 6, 4, 3, 5, 2, 11, 8, 10, 9, 7, 1};
    uint64_t const c = UINT64_C(0x4ec4ec4ec4ec4ec5); /* 1/13 (mod 2^64) */
    return me.divn(r, a, 13, 3, inv_13, c);
}

// }}}

template<>
void
arithxx_details::api<arithxx_modredc126>::gcd(Integer & r, const Residue & A) const
{
    auto const & me = downcast();
    uint64_t a[2], b[2];
    int sh;

    /* Since we do REDC arithmetic, we must have m odd */
    ASSERT_EXPENSIVE(me.m[0] % 2 != 0);

    if (me.is0(A)) {
        me.getmod(r);
        return;
    }

    a[0] = A.r[0];
    a[1] = A.r[1];
    b[0] = me.m[0];
    b[1] = me.m[1];

    while (a[1] != 0 || a[0] != 0) {
        /* Make a odd */
#if LOOKUP_TRAILING_ZEROS
        do {
            sh = trailing_zeros[(unsigned char)a[0]];
            u64arith_shrd(&(a[0]), a[1], a[0], sh);
            *(int64_t *)&(a[1]) >>= sh;
        } while (sh == 8);
#else
        if (a[0] == 0) /* ctzl does not like zero input */
        {
            a[0] = a[1];
            a[1] = ((int64_t)a[1] < 0L) ? (uint64_t)(-1L) : 0;
        }
        sh = u64arith_ctz(a[0]);
        u64arith_shrd(&(a[0]), a[1], a[0], sh);
        *(int64_t *)&(a[1]) >>= sh;
#endif

        /* Try to make the low two bits of b[0] zero */
        ASSERT_EXPENSIVE(a[0] % 2 == 1);
        ASSERT_EXPENSIVE(b[0] % 2 == 1);
        if ((a[0] ^ b[0]) & 2)
            u64arith_add_2_2(&(b[0]), &(b[1]), a[0], a[1]);
        else
            u64arith_sub_2_2(&(b[0]), &(b[1]), a[0], a[1]);

        if (b[0] == 0 && b[1] == 0) {
            if ((int64_t)a[1] < 0) {
                a[1] = -a[1];
                if (a[0] != 0)
                    a[1]--;
                a[0] = -a[0];
            }
            r = Integer(a[0], a[1]);
            return;
        }

        /* Make b odd */
#if LOOKUP_TRAILING_ZEROS
        do {
            sh = trailing_zeros[(unsigned char)b[0]];
            u64arith_shrd(&(b[0]), b[1], b[0], sh);
            *(int64_t *)&(b[1]) >>= sh;
        } while (sh == 8);
#else
        if (b[0] == 0) /* ctzl does not like zero input */
        {
            b[0] = b[1];
            b[1] = ((int64_t)b[1] < 0) ? (uint64_t)(-1) : 0;
        }
        sh = u64arith_ctz(b[0]);
        u64arith_shrd(&(b[0]), b[1], b[0], sh);
        *(int64_t *)&(b[1]) >>= sh;
#endif
        ASSERT_EXPENSIVE(a[0] % 2 == 1);
        ASSERT_EXPENSIVE(b[0] % 2 == 1);

        if ((a[0] ^ b[0]) & 2)
            u64arith_add_2_2(&(a[0]), &(a[1]), b[0], b[1]);
        else
            u64arith_sub_2_2(&(a[0]), &(a[1]), b[0], b[1]);
    }

    if ((int64_t)b[1] < 0) {
        b[1] = -b[1];
        if (b[0] != 0)
            b[1]--;
        b[0] = -b[0];
    }
    r = Integer(b[0], b[1]);
}

int arithxx_modredc126::Modulus::jacobi(Residue const & a_par) const
{
    Integer a, m, s;
    int t = 1;

    get(a, a_par);
    getmod(m);

    while (a != 0) {
        while ((a & 1) == 0) { /* TODO speedup */
            a >>= 1;
            if ((m & 7) == 3 || (m & 7) == 5)
                t = -t;
        }
        s = a; /* swap a and m */
        a = m;
        m = s;
        if ((a & 3) == 3 && (m & 3) == 3)
            t = -t;

        /* m is odd here */
        if (a >= m)
            a %= m;
    }
    if (m != 1)
        t = 0;

    return t;
}

bool arithxx_modredc126::Modulus::inv(Residue & r, Residue const & A) const
{
    Integer a, b, u, v;
    int t, lsh;
#ifdef WANT_ASSERT_EXPENSIVE
    Residue tmp(*this);

    set(tmp, A);
#endif

    assertValid(A);
    ASSERT_EXPENSIVE(m[0] % 2 != 0);

    if (is0(A))
        return false;

    getmod(b);

    /* Let A = x*2^{2w}, so we want the Montgomery representation of 1/x,
       which is 2^{2w}/x. We start by getting a = x */

    /* We simply set a = x/2^{2w} and t=0. The result before correction
       will be 2^(2w+t)/x so we have to divide by t, which may be >64,
       so we may have to do one or more full and a variable width REDC. */
    /* TODO: If b[1] > 1, we could skip one of the two REDC */
    {
        Residue x(*this);
        set(x, A);
        redc1(x, x);
        get(a, x);
    }
    /* Now a = x/2^w */
    t = -64;

    u = 1;
    v = 0; /* 0 is a valid pointer */

    /* make a odd */
    lsh = a.ctz();
    t += lsh;
    a >>= lsh;

    // Here a and b are odd, and a < b
    do {
        /* Here, a and b are odd, 0 < a < b, u is odd and v is even */
        ASSERT_EXPENSIVE(a < b);
        ASSERT_EXPENSIVE((a & 1) == 1);
        ASSERT_EXPENSIVE((b & 1) == 1);
        ASSERT_EXPENSIVE((u & 1) == 1);
        ASSERT_EXPENSIVE((v & 1) == 0);

        do {
            b -= a;
            v += u;
            ASSERT_EXPENSIVE((b & 1) == 0);

            lsh = b.ctz();
            t += lsh;
            b >>= lsh;
            u <<= lsh;
        } while (a < b); /* ~50% branch taken :( */

        /* Here, a and b are odd, 0 < b =< a, u is even and v is odd */
        ASSERT_EXPENSIVE((a & 1) == 1);
        ASSERT_EXPENSIVE((b & 1) == 1);
        ASSERT_EXPENSIVE((u & 1) == 0);
        ASSERT_EXPENSIVE((v & 1) == 1);

        if (a == b)
            break;
        ASSERT_EXPENSIVE(a > b);

        /* Here, a and b are odd, 0 < b < a, u is even and v is odd */
        do {
            a -= b;
            u += v;

            ASSERT_EXPENSIVE((a & 1) == 0);
            lsh = a.ctz();
            a >>= lsh;
            t += lsh;
            v <<= lsh;
        } while (b < a); /* about 50% branch taken :( */
        /* Here, a and b are odd, 0 < a =< b, u is odd and v is even */
    } while (a != b);

    if (a != 1) /* Non-trivial GCD */
        return false;

    ASSERT_ALWAYS(t >= 0);

    /* Here, the inverse of a is u/2^t mod m. To do the division by t,
       we use a variable-width REDC. We want to add a multiple of m to u
       so that the low t bits of the sum are 0 and we can right-shift by t
       with impunity. */
    for (; t >= 64; t -= 64)
        redc1(u, u);

    if (t > 0) {
        uint64_t s[5], k;
        k = ((u.getWord(0) * invm) &
             ((UINT64_C(1) << t) - 1)); /* tlow <= 2^t-1 */
        u64arith_mul_1_1_2(&(s[0]), &(s[1]), k, m[0]);
        /* s[1]:s[0] <= (2^w-1)*(2^t-1) <= (2^w-1)*(2^(w-1)-1) */
        u64arith_add_2_2(&(s[0]), &(s[1]), u.getWord(0), u.getWord(1));
        /* s[1]:s[0] <= (2^w-1)*(2^(w-1)-1) + (m-1) < 2^(2w) */
        /* s[0] == 0 (mod 2^t) */
        ASSERT_EXPENSIVE((s[0] & ((1UL << t) - 1)) == 0);
        s[2] = 0;
        u64arith_mul_1_1_2(&(s[3]), &(s[4]), k, m[1]);
        u64arith_add_2_2(&(s[1]), &(s[2]), s[3], s[4]);

        /* Now shift s[2]:s[1]:s[0] right by t */
        u64arith_shrd(&(s[0]), s[1], s[0], t);
        u64arith_shrd(&(s[1]), s[2], s[1], t);

        u = Integer(s[0], s[1]);
        // t = 0;
    }

    u.get(r.r, 2);
#ifdef WANT_ASSERT_EXPENSIVE
    mul(tmp, tmp, r);
    ASSERT_EXPENSIVE(is1(tmp));
#endif

    return true;
}

bool arithxx_modredc126::Modulus::batchinv_redc(Residue * r, uint64_t const * a, size_t const n,
                       Integer const * c) const
{
    auto const & me = *this;
    Residue R(me);

    if (n == 0)
        return true;

    r[0] = a[0];

    /* beta' = 2^64, beta = 2^128 */
    for (size_t i = 1; i < n; i++) {
        me.mul_ul(r[i], r[i - 1], a[i]);
        /* r[i] = beta'^{-i} \prod_{0 <= j <= i} a[j] */
    }

    /* Computes R = beta^2/r[n-1] */
    if (!me.inv(R, r[n - 1]))
        return false;
    /* R = beta^2 beta'^{n-1} \prod_{0 <= j < n} a[j]^{-1} */

    if (c != nullptr) {
        Residue t(me);
        t = *c;
        me.mul(R, R, t);
    } else {
        me.redc1(R, R); /* Assume c=1 */
        me.redc1(R, R);
    }
    /* R = beta beta'^{n-1} c \prod_{0 <= j < n} a[j]^{-1} */

    me.redc1(R, R);
    /* R = beta beta'^{n-2} c \prod_{0 <= j < n} a[j]^{-1} */

    for (size_t i = n - 1; i > 0; i--) {
        /* Invariant: R = beta beta'^{i-1} c \prod_{0 <= j <= i} a[j]^{-1} */

        me.mul(r[i], R, r[i - 1]);
        /* r[i] := R * r[i-1] / beta
                = (beta beta'^{i-1} c \prod_{0 <= j <= i} a[j]^{-1}) *
           (1/beta'^{i-1} \prod_{0 <= j <= i-1} a[j]) / beta = c a[i]^{-1} */

        me.mul_ul(R, R, a[i]);
        /* R := R * a[i] / beta'
             = (beta beta'^{i-1} c \prod_{0 <= j <= i} a[j]^{-1}) * a[i] / beta'
             = beta beta'^{i-2} c \prod_{0 <= j < i} a[j]^{-1},
           thus satisfying the invariant for i := i - 1 */
    }
    /* Here have R = beta * beta'^{-1} / a[0]. We need to convert the factor
       beta to a factor of beta', so that the beta' cancel. */
    me.redc1(R, R); /* R := beta * beta'^{-1} / a[0] / beta',
                    with beta = beta'^2, this is 1/a[0] */
    me.set(r[0], R);
    return true;
}

#if 0

Modulus::batch_Q_to_Fp_context_t *
Modulus::batch_Q_to_Fp_init (const Integer &num, const Integer &den) const
{
  batch_Q_to_Fp_context_t *context;
  Integer ratio, remainder;

  context = (batch_Q_to_Fp_context_t *) malloc(sizeof(batch_Q_to_Fp_context_t));
  if (context == NULL)
    return NULL;

  mod_initmod_int(context->m, den);

  /* Compute ratio = floor(num / den), remainder = num % den. We assume that
    ratio fits into uint64_t, and abort if it does not. We need only the
    low word of remainder. */
  remainder = num % den;
  ratio = num - remainder;
  ratio = ratio.divexact(den);
  ASSERT_ALWAYS(ratio.size() == 1);
  ratio.get(&(context->ratio_ul), 1);
  // ASSERT_ALWAYS(remainder.size() == 1);
  remainder.get(&(context->rem_ul), 1);
  if (remainder != 0)
    context->c = den - remainder; /* c = -remainder (mod den) */

  context->den_inv = u64arith_invmod(den.get()[0]);

  return context;
}


void
modredc2ul2_batch_Q_to_Fp_clear (modredc2ul2_batch_Q_to_Fp_context_t * context)
{
  mod_clearmod(context->m);
  mod_intclear(context->c);
  free(context);
}


int
modredc2ul2_batch_Q_to_Fp (uint64_t *r,
                           const modredc2ul2_batch_Q_to_Fp_context_t *context,
                           const uint64_t k, const int neg,
                           const uint64_t *p, const size_t n)
{
  Residue *tr;
  int rc = 1;

  tr = (Residue *) malloc(n * sizeof(Residue));
  for (size_t i = 0; i < n; i++) {
    mod_init_noset0(tr[i], context->m);
  }

  if (!modredc2ul2_batchinv_ul(tr, p, n, context->c, context->m)) {
    rc = 0;
    goto clear_and_exit;
  }

  for (size_t i = 0; i < n; i++) {
    uint64_t t;
    t = ularith_post_process_inverse(mod_intget_ul(tr[i]), p[i],
                                     context->rem_ul, context->den_inv,
                                     context->ratio_ul, k);
    if (neg && t != 0)
      t = p[i] - t;
    r[i] = t;
  }

clear_and_exit:
  for (size_t i = 0; i < n; i++) {
    mod_clear(tr[i], context->m);
  }
  free(tr);
  return rc;
}
#endif

template struct arithxx_details::api<arithxx_modredc126>;
