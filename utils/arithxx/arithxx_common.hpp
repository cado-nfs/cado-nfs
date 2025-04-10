#ifndef UTILS_ARITHXX_COMMON_HPP_
#define UTILS_ARITHXX_COMMON_HPP_

#include "arithxx_api.hpp"       // IWYU pragma: export
#include "arithxx_api64.hpp"     // IWYU pragma: export

namespace arithxx_details {
    template <typename layer>
        struct Residue_base {
            typedef typename layer::Modulus Modulus;
            typedef typename layer::Residue Residue;
            typedef typename layer::Integer Integer;
        };
}

#endif /* UTILS_ARITHXX_COMMON_HPP_ */
