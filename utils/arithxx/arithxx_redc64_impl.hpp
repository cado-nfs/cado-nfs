#ifndef CADO_UTILS_ARITHXX_REDC64_IMPL_HPP
#define CADO_UTILS_ARITHXX_REDC64_IMPL_HPP

#include <cstdint>

#include <vector>

#include "arithxx_redc64.hpp"

/* Compute r[i] = c * a[i]^(-1) mod m, for 0 <= i < n, where a[i] are
   non-negative integers and r[i] are integers with 0 <= r[i] < m.  a[i]
   need not be reduced modulo m. r_ul and a_ul must be non-overlapping.
   If any inverse does not exists, returns 0 and contents of r are
   undefined, otherwise returns 1. */

template<typename layer>
auto
arithxx_details::redc64<layer>::batchinv_redc(std::vector<uint64_t> const & a, Integer const & c) const
-> std::vector<Integer> 
{
    auto const & me = downcast();

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
        me.mul_u64_u64(R.r[0], R, x);
        ASSERT_ALWAYS(R.r < me.m[0]);
        r.push_back(R.r);
    }

    /* r[i] is a reduced representative of a[0]*...*a[i]. It
     * happens to be the Montgomery representative of
     * a'_0*...*a'_i where a'_i = a[i]/beta.
     */
    int const rc = me.inv(R, R);
    if (rc == 0)
        return {};

    /* R is the Montgomery representative of [a'_0*...*a'_{n-1}]^-1
     * c is the Montgomery representative of c/beta
     */
    me.mul_u64_u64(R.r[0], R, uint64_t(c));
    me.frommontgomery(R, R);
    /* R is now [a'_0*...*a'_{n-1}]^-1*c/beta, a.k.a the Montgomery
     * representative of [a'_0*...*a'_{n-1}]^-1*c/beta^2 */

    for (size_t i = a.size() - 1; i > 0; i--) {
        me.mul_u64_u64(r[i][0], R, r[i - 1][0]);
        /* r[i] is the Montgomery representative of a'_i^-1*c/beta^2
         * i.e.
         * r[i] = a'_i^-1*c/beta == (a_i / beta)^-1 * c/beta = c/a_i
         */
        me.mul_u64_u64(R.r[0], R, a[i]);
    }
    r[0] = R.r;

    return r;
}



#endif	/* UTILS_ARITHXX_REDC64_IMPL_HPP_ */
