#ifndef UTILS_ARITHXX_COMMON_HPP_
#define UTILS_ARITHXX_COMMON_HPP_

#include <climits>

#include <type_traits>

#include "arithxx_api.hpp"       // IWYU pragma: export
#include "arithxx_api64.hpp"     // IWYU pragma: export
#include "arithxx_api128.hpp"     // IWYU pragma: export

namespace arithxx_details {
    template <typename layer>
        struct Residue_base {
            typedef typename layer::Modulus Modulus;
            typedef typename layer::Residue Residue;
            typedef typename layer::Integer Integer;
        };

    template<int n, typename layer>
        struct multiplier_fits_in_overflow_bits
        : public std::integral_constant<bool,
            layer::overflow_bits::value == INT_MAX ||
            (n >> layer::overflow_bits::value == 0)> {};
}

#endif /* UTILS_ARITHXX_COMMON_HPP_ */
