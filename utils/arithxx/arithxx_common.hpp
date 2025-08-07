#ifndef CADO_UTILS_ARITHXX_COMMON_HPP
#define CADO_UTILS_ARITHXX_COMMON_HPP

#include <limits>
#include <type_traits>

#include "arithxx_api.hpp"       // IWYU pragma: export
#include "arithxx_api64.hpp"     // IWYU pragma: export
#include "arithxx_api128.hpp"    // IWYU pragma: export

namespace arithxx_details {
    template <typename layer>
        struct Residue_base {
            typedef typename layer::Modulus Modulus;
            typedef typename layer::Residue Residue;
            typedef typename layer::Integer Integer;

            /* the default implementation of residues consists in having
             * an Integer (a priori fixed-size), and force the ctor to
             * contain a reference to the modulus, even though we don't
             * actually use it. As a matter of fact, it's not entirely
             * clear that we can ensure that all constructions of
             * Residue types will have a complete Modulus object
             * available!
             */
            Integer r;

            explicit Residue_base(Modulus const &)
                : r(0)
            { }
            Residue_base(Modulus const &, Residue const & s)
                : r(s.r)
            { }

            Residue_base() = delete;
        };

    template<int n, typename layer>
        struct multiplier_fits_in_overflow_bits
        : public std::integral_constant<bool,
            layer::overflow_bits::value >= std::numeric_limits<int>::digits ||
            (n >> layer::overflow_bits::value == 0)> {};

    template<typename layer>
        struct redc;

    template<typename layer>
        struct batch_Q_to_Fp_context;
} /* namespace arithxx_details */

#endif /* UTILS_ARITHXX_COMMON_HPP_ */
