#ifndef CADO_UTILS_ARITHXX_API64_HPP
#define CADO_UTILS_ARITHXX_API64_HPP

#include <type_traits>

#include "arithxx_api.hpp"
#include "modint.hpp"

namespace arithxx_details {
template <typename layer>
struct api_bysize<layer, Integer64>
: public api<layer>
{
    using typename api<layer>::Modulus;
    using typename api<layer>::Residue;
    using typename api<layer>::Integer;
    using api<layer>::downcast;

    using api<layer>::api;

    static_assert(std::is_same<Integer, Integer64>::value,
            "api64 assumes that our underlying integer type is Integer64");

    bool sprp2_is_enough() const { return downcast().m < 2047; }

    /* a type that is 64-bit based needs reduction when we create a
     * residue from an immediate
     */
    using api<layer>::set;
    void set(Residue & r, uint64_t const s) const {
        auto const & me = downcast();
        r.r = s % me.m[0];
    }

};
} /* namespace arithxx_details */

#endif	/* UTILS_ARITHXX_API64_HPP_ */
