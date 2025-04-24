#include "cado.h" // IWYU pragma: keep

#include <cstddef>
#include <cstdint>

#include <vector>

#include "arithxx_common.hpp"
#include "macros.h"
#include "modredc64.hpp"
#include "u64arith.h"

/* Only the .cpp source files that emit the non-inline symbols will
 * include this impl header file. So even though it does not look like
 * we're using it, in fact we are!  */
#include "arithxx_api_impl.hpp"      // IWYU pragma: keep
#include "arithxx_api64_impl.hpp"  // IWYU pragma: keep
#include "arithxx_redc_impl.hpp"  // IWYU pragma: keep
#include "arithxx_redc64_impl.hpp"  // IWYU pragma: keep
#include "arithxx_batch_Q_to_Fp_impl.hpp"  // IWYU pragma: keep

// scan-headers: stop here

template<>
bool arithxx_details::api<arithxx_modredc64>::inv(Residue & r, const Residue & A) const
{
    auto const & me = downcast();
    Integer x = me.m;
    uint64_t u, v;
    int t, lsh;

    ASSERT(A.r < x);
    ASSERT(x & 1);

    if (A.r == 0)
        return false;

    /* Let A = a*2^w, so we want the Montgomery representation of 1/a,
       which is 2^w/a. We start by getting y = a */
    Integer y = me.frommontgomery(A.r);

    /* We simply set y = a/2^w and t=0. The result before
       correction will be 2^(w+t)/a so we have to divide by t, which
       may be >64, so we may have to do a full and a variable width REDC. */
    y = me.frommontgomery(y);

    /* Now y = a/2^w */
    t = 0;

    u = 1;
    v = 0;

    // make y odd
    lsh = int(y.ctz());
    y >>= lsh;
    t += lsh;
    /* v <<= lsh; ??? v is 0 here */

    // Here y and x are odd, and y < x
    do {
        /* Here, y and x are odd, 0 < y < x, u is odd and v is even */
        do {
            x -= y;
            v += u;
            if (x == 0)
                break;
            lsh = int(x.ctz());
            ASSERT_EXPENSIVE(lsh > 0);
            x >>= lsh;
            t += lsh;
            u <<= lsh;
        } while (x > y); /* ~50% branch taken :( */
        /* Here, y and x are odd, 0 < x =< y, u is even and v is odd */

        /* x is the one that got reduced, test if we're done */

        if (x <= 1)
            break;

        /* Here, y and x are odd, 0 < x < y, u is even and v is odd */
        do {
            y -= x;
            u += v;
            if (y == 0)
                break;
            lsh = int(y.ctz());
            ASSERT_EXPENSIVE(lsh > 0);
            y >>= lsh;
            t += lsh;
            v <<= lsh;
        } while (x < y); /* about 50% branch taken :( */
        /* Here, y and x are odd, 0 < y =< x, u is odd and v is even */
        /* y is the one that got reduced, test if we're done */
    } while (y > 1);

    if ((x & y) == 0) /* Non-trivial GCD */
        return false;

    if (y != 1) {
        /* So x is the one that reached 1.
           We maintained ya == u2^t (mod m) and xa = -v2^t (mod m).
           So 1/a = -v2^t.
         */
        u = me.m[0] - v;
        /* Now 1/a = u2^t */
    }

    /* Here, u = 2^w * 2^t / a. We want 2^w / a. */

    /* Here, the inverse of y is u/2^t mod x. To do the division by t,
       we use a variable-width REDC. We want to add a multiple of m to u
       so that the low t bits of the sum are 0 and we can right-shift by t. */
    if (t >= 64) {
        uint64_t tlow, thigh;
        tlow = u * me.invm; /* tlow <= 2^w-1 */
        u64arith_mul_1_1_2(&tlow, &thigh, tlow, me.m[0]); /* thigh:tlow <= (2^w-1)*m */
        u = thigh + ((u != 0) ? 1 : 0);
        /* thigh:tlow + u < (2^w-1)*m + m < 2^w*m. No correction necesary */
        t -= 64;
    }

    ASSERT(t < 64);
    if (t > 0) {
        uint64_t tlow, thigh;
        /* Necessarily t < 64, so the shift is ok */
        /* Doing a left shift first and then a full REDC needs a modular
           addition at the end due to larger summands and thus is probably
           slower */
        tlow = ((u * me.invm) & (((uint64_t)1 << t) - 1)); /* tlow <= 2^t-1 */
        u64arith_mul_1_1_2(&tlow, &thigh, tlow, me.m[0]); /* thigh:tlow <= m*(2^t-1) */
        u64arith_add_1_2(&tlow, &thigh,
                         u); /* thigh:tlow <= m*2^t-1 (since u<m) */
        /* Now the low t bits of tlow are 0 */
        ASSERT_EXPENSIVE((tlow & (((uint64_t)1 << t) - 1)) == 0);
        u64arith_shrd(&tlow, thigh, tlow, t);
        u = tlow;
        ASSERT_EXPENSIVE((thigh >> t) == 0 && u < m[0]);
    }

    r.r = u;
    return true;
}

/* same as inv, but for classical representation (not Montgomery). */
template<>
bool arithxx_details::api<arithxx_modredc64>::intinv(Integer & r, Integer const & A) const
{
    auto const & me = downcast();
    uint64_t x = me.m[0], y, u, v;
    int t, lsh;

    ASSERT(x & 1);

    if (A == 0)
        return false;

    y = uint64_t(A);

    /* We don't expect it's going to happen normally, except maybe for tests */
    if (y >= me.m[0])
        y %= me.m[0];

    t = 0;

    u = 1;
    v = 0;

    // make y odd
    lsh = u64arith_ctz(y);
    y >>= lsh;
    t += lsh;

    do {
        /* Here, x and y are odd, 0 < y < x, u is odd and v is even */
        ASSERT_EXPENSIVE(y % 2 == 1);
        ASSERT_EXPENSIVE(u % 2 == 1);
        ASSERT_EXPENSIVE(v % 2 == 0);
        ASSERT_EXPENSIVE(0 < y);
        ASSERT_EXPENSIVE(y < x);
        do {
            ASSERT_EXPENSIVE(x % 2 == 1);
            x -= y;
            v += u;
            lsh = u64arith_ctz(x);
            ASSERT_EXPENSIVE(lsh > 0);
            x >>= lsh;
            t += lsh;
            u <<= lsh;
        } while (x > y); /* about 50% branch taken :( */

        /* x is the one that got reduced, test if we're done */
        /* Here, x and y are odd, 0 < x <= y, u is even and v is odd */
        ASSERT_EXPENSIVE(0 < x);
        ASSERT_EXPENSIVE(x <= y);
        ASSERT_EXPENSIVE(x % 2 == 1);
        ASSERT_EXPENSIVE(y % 2 == 1);
        ASSERT_EXPENSIVE(u % 2 == 0);
        ASSERT_EXPENSIVE(v % 2 == 1);

        if (x == y)
            break;

        /* Here, x and y are odd, 0 < x < y, u is even and v is odd */
        do {
            ASSERT_EXPENSIVE(y % 2 == 1);
            y -= x;
            u += v;
            lsh = u64arith_ctz(y);
            ASSERT_EXPENSIVE(lsh > 0);
            y >>= lsh;
            t += lsh;
            v <<= lsh;
        } while (x < y); /* about 50% branch taken :( */
        /* Here, x and y are odd, 0 < y <= x, u is odd and v is even */

        /* y is the one that got reduced, test if we're done */
    } while (x != y);

    /* when we exit, x == y */
    ASSERT_EXPENSIVE(x == y);

    if (x > 1) /* Non-trivial GCD */
        return false;

    if (u % 2 == 0) {
        /* We exited the loop after reducing x */
        /* We maintained ya == u2^t (mod m) and xa = -v2^t (mod m).
           So 1/a = -v2^t. */
        u = me.m[0] - v;
        /* Now 1/a = u2^t */
    }

    /* Here, u = 2^w * 2^t / a. We want 2^w / a. */

    /* Here, the inverse of y is u/2^t mod x. To do the division by t,
       we use a variable-width REDC. We want to add a multiple of m to u
       so that the low t bits of the sum are 0 and we can right-shift by t. */
    if (t >= 64) {
        uint64_t tlow, thigh;
        tlow = u * me.invm; /* tlow <= 2^w-1 */
        u64arith_mul_1_1_2(&tlow, &thigh, tlow, me.m[0]); /* thigh:tlow <= (2^w-1)*m */
        u = thigh + ((u != 0) ? 1 : 0);
        /* thigh:tlow + u < (2^w-1)*m + m < 2^w*m. No correction necesary */
        t -= 64;
    }

    ASSERT(t < 64);
    if (t > 0) {
        uint64_t tlow, thigh;
        /* Necessarily t < 64, so the shift is ok */
        /* Doing a left shift first and then a full REDC needs a modular
           addition at the end due to larger summands and thus is probably
           slower */
        tlow = ((u * me.invm) & (((uint64_t)1 << t) - 1)); /* tlow <= 2^t-1 */
        u64arith_mul_1_1_2(&tlow, &thigh, tlow, me.m[0]); /* thigh:tlow <= m*(2^t-1) */
        u64arith_add_1_2(&tlow, &thigh, u); /* thigh:tlow <= m*2^t-1 (since u<m) */
        /* Now the low t bits of tlow are 0 */
        ASSERT_EXPENSIVE((tlow & (((uint64_t)1 << t) - 1)) == 0);
        u64arith_shrd(&tlow, thigh, tlow, t);
        u = tlow;
        ASSERT_EXPENSIVE((thigh >> t) == 0 && u < m[0]);
    }

    r = u;
    return true;
}

#if 0
/* Let v = lo + 2^64 * hi and
   subtrahend = subtrahend_lo + subtrahend_hi * 2^64.
   Return 1 if (v - subtrahend) / divisor is a non-negative integer less than
   2^64, and 0 otherwise */
MAYBE_UNUSED static inline int check_divisible(uint64_t const lo,
                                               uint64_t const hi,
                                               uint64_t const subtrahend_lo,
                                               uint64_t const subtrahend_hi,
                                               uint64_t const divisor)
{
    /* Test for subtraction underflow */
    if (u64arith_gt_2_2(subtrahend_lo, subtrahend_hi, lo, hi))
        return 0;
    uint64_t diff_lo = lo, diff_hi = hi;
    u64arith_sub_2_2(&diff_lo, &diff_hi, subtrahend_lo, subtrahend_hi);
    /* Test for division overflow */
    if (diff_hi >= divisor)
        return 0;
    uint64_t q, r;
    u64arith_divqr_2_1_1(&q, &r, diff_lo, diff_hi, divisor);
    return r == 0;
}
#endif

template<>
arithxx_details::batch_Q_to_Fp_context<arithxx_modredc64>::batch_Q_to_Fp_context(Integer const & num, Integer const & den)
        : remainder(num[0] % den[0])
        , quotient(num[0] / den[0])
        , D(den)
    { }

template struct arithxx_details::batch_Q_to_Fp_context<arithxx_modredc64>;
template struct arithxx_details::redc<arithxx_modredc64>;
template struct arithxx_details::redc64<arithxx_modredc64>;
template struct arithxx_details::api<arithxx_modredc64>;
template struct arithxx_details::api_bysize<arithxx_modredc64>;
