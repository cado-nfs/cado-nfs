#ifndef CADO_UTILS_NUMBERTHEORY_NUMBER_FIELD_PRIME_IDEAL_HPP
#define CADO_UTILS_NUMBERTHEORY_NUMBER_FIELD_PRIME_IDEAL_HPP

#include "numbertheory/numbertheory_fwd_types.hpp"
#include "numbertheory/number_field_fractional_ideal.hpp"
#include "numbertheory/number_field_order_element.hpp"

class number_field_prime_ideal : private number_field_fractional_ideal {
    number_field_order_element valuation_helper;
    cxx_mpz p;
    int e;      /* ramification index */
    public:
        using number_field_fractional_ideal::two_element;
        using number_field_fractional_ideal::order;
        using number_field_fractional_ideal::number_field;

    operator two_element() const;
    int inertia_degree() const;
    inline int ramification_index() const { return e; };

    private:
    /* This ctor computes the helper. The tricky thing is that we may not
     * want to expose it, except to friends...
     */
    number_field_prime_ideal(number_field_fractional_ideal const & I, cxx_mpz const & p, int e);
    int valuation(number_field_fractional_ideal const & I) const;
    friend class number_field_order;
    friend class number_field_fractional_ideal;
};

namespace fmt {
    template <>
    struct formatter<number_field_prime_ideal>
        : formatter<string_view>
        , fmt_helper_sagemath<number_field_prime_ideal>
    {
        using fmt_helper_sagemath::parse;
        static constexpr const decltype(custom_format) custom_format_default = SAGEMATH;
        auto format(number_field_prime_ideal const & e, format_context& ctx) const -> format_context::iterator;
    };
}
inline std::ostream& operator<<(std::ostream& os, number_field_prime_ideal const & e) { return os << fmt::format("{}", e); }


#endif	/* UTILS_NUMBERTHEORY_NUMBER_FIELD_PRIME_IDEAL_HPP_ */
