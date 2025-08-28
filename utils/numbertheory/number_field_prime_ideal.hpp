#ifndef CADO_UTILS_NUMBERTHEORY_NUMBER_FIELD_PRIME_IDEAL_HPP
#define CADO_UTILS_NUMBERTHEORY_NUMBER_FIELD_PRIME_IDEAL_HPP

#include <gmp.h>
#include "fmt/base.h"
#include "fmt/format.h"

#include "numbertheory/numbertheory_fwd_types.hpp"
#include "numbertheory/number_field_fractional_ideal.hpp"
#include "numbertheory/number_field_order_element.hpp"
#include "fmt_helper_sagemath.hpp"
#include "cxx_mpz.hpp"

class number_field_prime_ideal : private number_field_fractional_ideal {
    number_field_order_element valuation_helper;
    cxx_mpz p;
    int e;      /* ramification index */
    public:
        using number_field_fractional_ideal::two_element;
        using number_field_fractional_ideal::order;
        using number_field_fractional_ideal::number_field;

    explicit operator two_element() const;
    int inertia_degree() const;
    int ramification_index() const { return e; };

    // The valuation helper a of I is such that (a/p)*I is in O, yet a is
    // not in p*O.
    // For an element u whose I-valuation is v, we thus have that u*(a/p)^v
    // maps to a non-zero element in O/pO.
    number_field_order_element const & get_valuation_helper() const {
        return valuation_helper;
    }
    
    /* Compute the valuation of I at the prime ideal *this */
    int valuation(number_field_fractional_ideal const & I) const;

    /* This ctor computes the helper.
     *
     * We expose it because otherwise we have no way to create a prime
     * ideal from external data. Wheher it's a good idea to offer this
     * possiblity isn't totally clear.
     *
     * We could also explore alternatives that include some checks that
     * we don't do at the moment (starting with the determinant of I
     * being some power of p)
     */
    number_field_prime_ideal(number_field_fractional_ideal I, cxx_mpz p, int e);

    friend class number_field_order;
    friend class number_field_fractional_ideal;

    int cmp(number_field_prime_ideal const & I) const {
        int r;
        r = mpz_cmp(p, I.p);
        if (r) return r;
        r = e - I.e;
        if (r) return r;
        r = number_field_fractional_ideal::cmp(I);
        return r;
    }
    bool operator<(number_field_prime_ideal const & I) const { return cmp(I) < 0; }
    bool operator==(number_field_prime_ideal const & I) const { return cmp(I) == 0; }
};

namespace fmt {
    template <>
    struct formatter<number_field_prime_ideal>
        : fmt_helper_sagemath<number_field_prime_ideal>
    {
        static constexpr const decltype(custom_format) custom_format_default = SAGEMATH;
        auto format(number_field_prime_ideal const & e, format_context& ctx) const -> format_context::iterator;
    };
}
inline std::ostream& operator<<(std::ostream& os, number_field_prime_ideal const & e) { return os << fmt::format("{}", e); }


#endif	/* UTILS_NUMBERTHEORY_NUMBER_FIELD_PRIME_IDEAL_HPP_ */
