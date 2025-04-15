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

// scan-headers: stop here

template<>
bool arithxx_details::api<arithxx_modredc64>::inv(Residue & r, const Residue & A) const
{
    auto const & me = downcast();
    uint64_t x = me.m, y, u, v;
    int t, lsh;

    ASSERT(A.r < x);
    ASSERT(x & 1);

    if (A.r == 0)
        return false;

    /* Let A = a*2^w, so we want the Montgomery representation of 1/a,
       which is 2^w/a. We start by getting y = a */
    y = me.get_u64(A);

    /* We simply set y = a/2^w and t=0. The result before
       correction will be 2^(w+t)/a so we have to divide by t, which
       may be >64, so we may have to do a full and a variable width REDC. */
    me.frommontgomery(y, y);
    /* Now y = a/2^w */
    t = 0;

    u = 1;
    v = 0;

    // make y odd
    lsh = u64arith_ctz(y);
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
            lsh = u64arith_ctz(x);
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
            lsh = u64arith_ctz(y);
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
        u = me.m - v;
        /* Now 1/a = u2^t */
    }

    /* Here, u = 2^w * 2^t / a. We want 2^w / a. */

    /* Here, the inverse of y is u/2^t mod x. To do the division by t,
       we use a variable-width REDC. We want to add a multiple of m to u
       so that the low t bits of the sum are 0 and we can right-shift by t. */
    if (t >= 64) {
        uint64_t tlow, thigh;
        tlow = u * me.invm; /* tlow <= 2^w-1 */
        u64arith_mul_1_1_2(&tlow, &thigh, tlow, me.m); /* thigh:tlow <= (2^w-1)*m */
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
        u64arith_mul_1_1_2(&tlow, &thigh, tlow, me.m); /* thigh:tlow <= m*(2^t-1) */
        u64arith_add_1_2(&tlow, &thigh,
                         u); /* thigh:tlow <= m*2^t-1 (since u<m) */
        /* Now the low t bits of tlow are 0 */
        ASSERT_EXPENSIVE((tlow & (((uint64_t)1 << t) - 1)) == 0);
        u64arith_shrd(&tlow, thigh, tlow, t);
        u = tlow;
        ASSERT_EXPENSIVE((thigh >> t) == 0 && u < m);
    }

    r.r = u;
    return true;
}

/* same as inv, but for classical representation (not Montgomery). */
template<>
bool arithxx_details::api<arithxx_modredc64>::intinv(Integer & r, Integer const & A) const
{
    auto const & me = downcast();
    uint64_t x = me.m, y, u, v;
    int t, lsh;

    ASSERT(x & 1);

    if (A == 0)
        return false;

    y = uint64_t(A);

    /* We don't expect it's going to happen normally, xcept maybe for tests */
    if (y >= me.m)
        y %= me.m;

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
        u = me.m - v;
        /* Now 1/a = u2^t */
    }

    /* Here, u = 2^w * 2^t / a. We want 2^w / a. */

    /* Here, the inverse of y is u/2^t mod x. To do the division by t,
       we use a variable-width REDC. We want to add a multiple of m to u
       so that the low t bits of the sum are 0 and we can right-shift by t. */
    if (t >= 64) {
        uint64_t tlow, thigh;
        tlow = u * me.invm; /* tlow <= 2^w-1 */
        u64arith_mul_1_1_2(&tlow, &thigh, tlow, me.m); /* thigh:tlow <= (2^w-1)*m */
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
        u64arith_mul_1_1_2(&tlow, &thigh, tlow, me.m); /* thigh:tlow <= m*(2^t-1) */
        u64arith_add_1_2(&tlow, &thigh, u); /* thigh:tlow <= m*2^t-1 (since u<m) */
        /* Now the low t bits of tlow are 0 */
        ASSERT_EXPENSIVE((tlow & (((uint64_t)1 << t) - 1)) == 0);
        u64arith_shrd(&tlow, thigh, tlow, t);
        u = tlow;
        ASSERT_EXPENSIVE((thigh >> t) == 0 && u < m);
    }

    r = u;
    return true;
}

/* Compute r[i] = c * a[i]^(-1) mod m, for 0 <= i < n, where a[i] are
   non-negative integers and r[i] are integers with 0 <= r[i] < m.  a[i]
   need not be reduced modulo m. r_ul and a_ul must be non-overlapping.
   If any inverse does not exists, returns 0 and contents of r are
   undefined, otherwise returns 1. */

std::vector<arithxx_modredc64::Integer>
arithxx_modredc64::Modulus::batchinv_redc(std::vector<uint64_t> const & a, Integer c) const
{
    /* We simply don't convert c to or from Montgomery representation.
     * Strangely enough, it all turns out well. */

    if (a.empty())
        return {};

    std::vector<Integer> r;
    r.reserve(a.size());
    /* The a[i]'s need not be reduced. When we multiply them by something
     * (by 1 for a[0], for example), we get a reduced representative */
    Residue R = one;
    for (auto const & x : a) {
        mul_u64_u64(R.r[0], R, x);
        ASSERT_ALWAYS(R.r < m);
        r.push_back(R.r);
    }

    /* r[i] is a reduced representative of a[0]*...*a[i]. It
     * happens to be the Montgomery representative of
     * a'_0*...*a'_i where a'_i = a[i]/beta.
     */
    int const rc = inv(R, R);
    if (rc == 0)
        return {};

    /* R is the Montgomery representative of [a'_0*...*a'_{n-1}]^-1
     * c is the Montgomery representative of c/beta
     */
    mul_u64_u64(R.r[0], R, uint64_t(c));
    frommontgomery(R.r[0], R.r[0]);
    /* R is now [a'_0*...*a'_{n-1}]^-1*c/beta, a.k.a the Montgomery
     * representative of [a'_0*...*a'_{n-1}]^-1*c/beta^2 */

    for (size_t i = a.size() - 1; i > 0; i--) {
        mul_u64_u64(r[i][0], R, r[i - 1][0]);
        /* r[i] is the Montgomery representative of a'_i^-1*c/beta^2
         * i.e.
         * r[i] = a'_i^-1*c/beta == (a_i / beta)^-1 * c/beta = c/a_i
         */
        mul_u64_u64(R.r[0], R, a[i]);
    }
    r[0] = R.r;

    return r;
}

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

arithxx_modredc64::Modulus::batch_Q_to_Fp_context::batch_Q_to_Fp_context(
        Integer const & num, Integer const & den)
    : remainder(num[0] % den[0])
    , quotient(num[0] / den[0])
    , D(den)
{
}

std::vector<uint64_t> arithxx_modredc64::Modulus::batch_Q_to_Fp_context::operator()(std::vector<uint64_t> const & p, int const k) const
{
    /* We use -rem (mod den) here. batchinv_ul() does not
       mandate its c parameter to be fully reduced, which occurs here in the
       case of rem == 0. */
    auto r = D.batchinv_redc(p, Integer(D.m - remainder[0]));
    if (r.empty())
        return {};

    std::vector<uint64_t> ri(p.size());

    for (size_t i = 0; i < p.size(); i++)
        ri[i] = u64arith_post_process_inverse(r[i][0], p[i],
                remainder[0], -D.invm, quotient[0], k);

    return ri;
}

/* For each 0 <= i < n, compute r[i] = num/(den*2^k) mod p[i].
   den must be odd. If k > 0, then all p[i] must be odd.
   The memory pointed to be r and p must be non-overlapping.
   Returns 1 if successful. If any modular inverse does not exist,
   returns 0 and the contents of r are undefined. */
std::vector<uint64_t> arithxx_modredc64::Modulus::batch_Q_to_Fp(Integer const & num,
                            Integer const & den, int const k,
                            std::vector<uint64_t> const & p)
{
    return batch_Q_to_Fp_context(num, den)(p, k);
}

template struct arithxx_details::api<arithxx_modredc64>;
template struct arithxx_details::api64<arithxx_modredc64>;
