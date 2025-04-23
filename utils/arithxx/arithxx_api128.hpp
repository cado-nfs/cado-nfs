#ifndef CADO_UTILS_ARITHXX_API128_HPP
#define CADO_UTILS_ARITHXX_API128_HPP


#include <type_traits>

#include "arithxx_api.hpp"
#include "modint.hpp"

namespace arithxx_details {
template <typename layer>
struct api_bysize<layer, Integer128>
: public api<layer>
{
    using typename api<layer>::Modulus;
    using typename api<layer>::Residue;
    using typename api<layer>::Integer;
    using api<layer>::downcast;

    using api<layer>::api;

    static_assert(std::is_same<Integer, Integer128>::value,
            "api128 assumes that our underlying integer type is Integer128");

    bool inv(Residue &, Residue const &) const;

};

}

#endif	/* UTILS_ARITHXX_API128_HPP_ */
