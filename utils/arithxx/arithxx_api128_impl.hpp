#ifndef CADO_UTILS_ARITHXX_API128_IMPL_HPP
#define CADO_UTILS_ARITHXX_API128_IMPL_HPP

#include <cstdint>

#include <array>

#include "arithxx_common.hpp"
#include "arithxx_api128.hpp"
#include "modint.hpp"
#include "macros.h"
#include "u64arith.h"

/* TODO: clean this up. This implementation does not belong here, since
 * it has traces of redc and possibly of max 126 bits as well.
 */
template<typename layer>
bool arithxx_details::api_bysize<layer, Integer128>::inv(Residue & r, Residue const & A) const
{
    static_assert(layer::uses_montgomery_representation::value,
            "This code assumes that the current layer uses Montgomery representation");

    auto const & me = downcast();
    Integer a, u, v;
    int t;
#ifdef WANT_ASSERT_EXPENSIVE
    Residue tmp(*this);

    set(tmp, A);
#endif

    me.assertValid(A);
    ASSERT_EXPENSIVE(m[0] % 2 != 0);

    if (me.is0(A))
        return false;

    Integer b = me.getmod();

    /* Let A = x*2^{2w}, so we want the Montgomery representation of 1/x,
       which is 2^{2w}/x. We start by getting a = x */

    /* We simply set a = x/2^{2w} and t=0. The result before correction
       will be 2^(2w+t)/x so we have to divide by t, which may be >64,
       so we may have to do one or more full and a variable width REDC. */
    /* TODO: If b[1] > 1, we could skip one of the two REDC */
    {
        Residue x(me);
        me.redc1(x, A);
        a = me.get(x);
    }
    /* Now a = x/2^w */
    t = -64;

    u = 1;
    v = 0; /* 0 is a valid pointer */

    /* make a odd */
    auto lsh = int(a.ctz());
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

            lsh = int(b.ctz());
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
            lsh = int(a.ctz());
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
        me.redc1(u, u);

    if (t > 0) {
        uint64_t s[5], k;
        k = ((u[0] * me.invm) & ((uint64_t(1) << t) - 1)); /* tlow <= 2^t-1 */
        u64arith_mul_1_1_2(&s[0], &s[1], k, me.m[0]);
        /* s[1]:s[0] <= (2^w-1)*(2^t-1) <= (2^w-1)*(2^(w-1)-1) */
        u64arith_add_2_2(&s[0], &s[1], u[0], u[1]);
        /* s[1]:s[0] <= (2^w-1)*(2^(w-1)-1) + (m-1) < 2^(2w) */
        /* s[0] == 0 (mod 2^t) */
        ASSERT_EXPENSIVE((s[0] & ((1UL << t) - 1)) == 0);
        s[2] = 0;
        u64arith_mul_1_1_2(&(s[3]), &(s[4]), k, me.m[1]);
        u64arith_add_2_2(&(s[1]), &(s[2]), s[3], s[4]);

        /* Now shift s[2]:s[1]:s[0] right by t */
        u64arith_shrd(&(s[0]), s[1], s[0], t);
        u64arith_shrd(&(s[1]), s[2], s[1], t);

        u = Integer(std::array<uint64_t, 2> {s[0], s[1]});
        // t = 0;
    }

    r.r = u.get();
#ifdef WANT_ASSERT_EXPENSIVE
    mul(tmp, tmp, r);
    ASSERT_EXPENSIVE(is1(tmp));
#endif

    return true;
}

#endif	/* UTILS_ARITHXX_API128_IMPL_HPP_ */
