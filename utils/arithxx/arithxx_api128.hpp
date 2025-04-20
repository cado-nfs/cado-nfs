#ifndef UTILS_ARITHXX_ARITHXX_API128_HPP_
#define UTILS_ARITHXX_ARITHXX_API128_HPP_

#include <cstdint>
#include <type_traits>

#include "arithxx_api.hpp"
#include "modint.hpp"

namespace arithxx_details {
template <typename layer>
struct api128
: public api<layer>
{
    using typename api<layer>::Modulus;
    using typename api<layer>::Residue;
    using typename api<layer>::Integer;
    using api<layer>::downcast;

    static_assert(std::is_same<Integer, Integer128>::value,
            "api128 assumes that our underlying integer type is Integer128");

    // not yet
    // void gcd (Integer &g, const Residue &r) const;
    // int jacobi(Residue const &) const;
};
}



#endif	/* UTILS_ARITHXX_ARITHXX_API128_HPP_ */
