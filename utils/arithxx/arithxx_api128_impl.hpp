#ifndef UTILS_ARITHXX_ARITHXX_API128_IMPL_HPP_
#define UTILS_ARITHXX_ARITHXX_API128_IMPL_HPP_

#include <cstdint>

#include "arithxx_common.hpp"
#include "arithxx_api128.hpp"
#include "u64arith.h"
#include "macros.h"

template <typename layer>
bool arithxx_details::api128<layer>::div3(Residue & r,
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
            u64arith_add_2_2(t.r.data(), t.r.data() + 1, me.m[0], me.m[1]);
        } else /* a3 == 1, m3 == 1 or a3 == 2, m3 == 2 */
        {
            u64arith_add_2_2(t.r.data(), t.r.data() + 1, me.m[0], me.m[1]);
            u64arith_add_2_2(t.r.data(), t.r.data() + 1, me.m[0], me.m[1]);
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
    me.add(t, r, r);
    me.add(t, t, r);
    ASSERT_EXPENSIVE(me.equal(a_backup, t));
#endif

    return true;
}


#endif	/* UTILS_ARITHXX_ARITHXX_API128_IMPL_HPP_ */
