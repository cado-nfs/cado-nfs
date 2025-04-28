#ifndef CADO_UTILS_NUMBERTHEORY_NUMBER_FIELD_ORDER_ELEMENT_HPP
#define CADO_UTILS_NUMBERTHEORY_NUMBER_FIELD_ORDER_ELEMENT_HPP

#include "numbertheory/numbertheory_fwd_types.hpp"
#include "numbertheory/number_field_order.hpp"
#include "mpz_mat.h"

class number_field_order_element {
    friend class number_field;
    friend class number_field_prime_ideal;
    friend struct fmt::formatter<number_field_order_element>;
    number_field_order const & O;
    cxx_mpz_mat coefficients;
    public:
    number_field_order_element(number_field_order const & O, cxx_mpz_mat const & e) : O(O), coefficients(e) {}
    number_field_order const & order() const { return O; }
    number_field_order_element(number_field_order const & O, number_field_element const & e);
    number_field_order_element operator*(number_field_order_element const &) const;
    cxx_mpz_mat multiplication_matrix() const;

    number_field_order_element(number_field_order_element const & a)
        : O(a.O)
        , coefficients(a.coefficients)
    {}

    number_field_order_element(number_field_order_element && a)
        : O(a.O)
        , coefficients(a.coefficients)
    {}

    number_field_order_element& operator=(number_field_order_element const & a)
    {
        ASSERT_ALWAYS(&O == &a.O);
        coefficients = a.coefficients;
        return *this;
    }
    number_field_order_element& operator=(number_field_order_element && a)
    {
        ASSERT_ALWAYS(&O == &a.O);
        coefficients = std::move(a.coefficients);
        return *this;
    }
};

namespace fmt {
    template <>
    struct formatter<number_field_order_element>
        : formatter<string_view>
        , fmt_helper_sagemath<number_field_order_element>
    {
        static constexpr const decltype(custom_format) custom_format_default = SAGEMATH;
        using fmt_helper_sagemath::parse;
        auto format(number_field_order_element const & e, format_context& ctx) const
            -> format_context::iterator;
    };
}
inline std::ostream& operator<<(std::ostream& os, number_field_order_element const & e) { return os << fmt::format("{}", e); }





#endif	/* UTILS_NUMBERTHEORY_NUMBER_FIELD_ORDER_ELEMENT_HPP_ */
