#ifndef UTILS_ARITHXX_API64_IMPL_HPP_
#define UTILS_ARITHXX_API64_IMPL_HPP_

#include <cstdint>

#include "arithxx_common.hpp"
#include "arithxx_api64.hpp"
#include "u64arith.h"
#include "macros.h"

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
