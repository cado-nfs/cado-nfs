#ifndef UTILS_ARITHXX_API64_IMPL_HPP_
#define UTILS_ARITHXX_API64_IMPL_HPP_

#include <cstdint>

#include "arithxx_common.hpp"
#include "arithxx_api64.hpp"
#include "u64arith.h"
#include "macros.h"

template<typename layer>
void
arithxx_details::api_bysize<layer, Integer64>::gcd (Integer &g, const Residue &r) const
{
    auto const & me = downcast();
    uint64_t t;

    auto a = uint64_t(r.r); /* This works the same for "a" in plain or Montgomery
                         representation */
    uint64_t b = me.getmod()[0];
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

#endif	/* UTILS_ARITHXX_API64_IMPL_HPP_ */
