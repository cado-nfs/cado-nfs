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
    Integer a, u;
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

    /* The loop works on plain words, as the old arith code did: with
     * Integer128 objects, gcc spills to memory in this loop, and it is
     * then 10% slower than the old code.
     */
    uint64_t a0 = a[0], a1 = a[1], b0 = b[0], b1 = b[1];
    uint64_t u0 = 1, u1 = 0, v0 = 0, v1 = 0;
    /* x < y. The best form depends on the compiler: clang is 20% faster
     * with plain C than with the asm subtraction with borrow, and gcc
     * (13 to 15) is 1 to 6% faster with the latter. */
    auto lt = [](uint64_t x0, uint64_t x1, uint64_t y0, uint64_t y1) {
#ifdef __clang__
        return x1 < y1 || (x1 == y1 && x0 < y0);
#else
        return u64arith_sub_2_2_cy(&x0, &x1, y0, y1) != 0;
#endif
    };

    /* make a odd */
    if (a0 == 0) {
        a0 = a1;
        a1 = 0;
        t += 64;
    }
    auto lsh = int(u64arith_ctz(a0));
    t += lsh;
    u64arith_shr_2(&a0, &a1, lsh);

    // Here a and b are odd, and a < b
    do {
        /* Here, a and b are odd, 0 < a < b, u is odd and v is even */
        ASSERT_EXPENSIVE(lt(a0, a1, b0, b1));
        ASSERT_EXPENSIVE(a0 & b0 & u0 & 1);
        ASSERT_EXPENSIVE((v0 & 1) == 0);

        do {
            u64arith_sub_2_2(&b0, &b1, a0, a1);
            u64arith_add_2_2(&v0, &v1, u0, u1);
            /* a zero low word is a word shift, so that all shifts
             * below are by less than 64 bits */
            if (b0 == 0) {
                b0 = b1;
                b1 = 0;
                ASSERT_EXPENSIVE(u1 == 0);
                u1 = u0;
                u0 = 0;
                t += 64;
            }
            lsh = int(u64arith_ctz(b0));
            t += lsh;
            u64arith_shr_2(&b0, &b1, lsh);
            u64arith_shl_2(&u0, &u1, lsh);
        } while (lt(a0, a1, b0, b1)); /* ~50% branch taken :( */

        /* Here, a and b are odd, 0 < b =< a, u is even and v is odd */
        ASSERT_EXPENSIVE(a0 & b0 & v0 & 1);
        ASSERT_EXPENSIVE((u0 & 1) == 0);

        if (a0 == b0 && a1 == b1)
            break;

        /* Here, a and b are odd, 0 < b < a, u is even and v is odd */
        do {
            u64arith_sub_2_2(&a0, &a1, b0, b1);
            u64arith_add_2_2(&u0, &u1, v0, v1);
            if (a0 == 0) {
                a0 = a1;
                a1 = 0;
                ASSERT_EXPENSIVE(v1 == 0);
                v1 = v0;
                v0 = 0;
                t += 64;
            }
            lsh = int(u64arith_ctz(a0));
            t += lsh;
            u64arith_shr_2(&a0, &a1, lsh);
            u64arith_shl_2(&v0, &v1, lsh);
        } while (lt(b0, b1, a0, a1)); /* about 50% branch taken :( */
        /* Here, a and b are odd, 0 < a =< b, u is odd and v is even */
    } while (a0 != b0 || a1 != b1);

    a = Integer(a0, a1);
    u = Integer(u0, u1);

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
