#ifndef UTILS_NUMBERTHEORY_NUMBER_FIELD_FRACTIONAL_IDEAL_HPP_
#define UTILS_NUMBERTHEORY_NUMBER_FIELD_FRACTIONAL_IDEAL_HPP_

#include <utility>
#include "fmt/format.h"
#include "numbertheory/numbertheory_fwd_types.hpp"
#include "numbertheory/number_field_order.hpp"
#include "cxx_mpz.hpp"
#include "mpz_mat.h"


class number_field_fractional_ideal {
    friend struct fmt::formatter<number_field_fractional_ideal>;
    friend class number_field_order;
    friend class number_field_prime_ideal;
    number_field_order const & O;
    /* the basis matrix is in HNF form, with respect to the basis of the
     * order */
    cxx_mpz_mat ideal_basis_matrix;
    cxx_mpz denominator;
    public:
    /* oddly enough, we have no code for 2-element except for prime
     * ideals. Should we do it just by CRT? It's an option, but we don't
     * have much need for it, really.
     */
    typedef std::pair<cxx_mpz, number_field_order_element> two_element;
    private:
    mutable std::unique_ptr<two_element> cached_two_element;
    number_field_fractional_ideal(number_field_order const & O,
            cxx_mpz_mat const& I,
            cxx_mpz const & d = 1)
        : O(O)
        , ideal_basis_matrix(I)
        , denominator(d)
    {}
    public:

    /* return the order this ideal is defined to be a fractional ideal
     * of. Note that the actual endomorphism ring of the ideal may of
     * course be bigger!
     */
    inline number_field_order const & order() const { return O; }

    /* return the parent number field */
    inline class number_field const & number_field() const { return O.number_field(); }

    // unimplemented (for now)
    // operator two_element() const;

    /* return the fractional O-ideal generated by the given generators
     */
    number_field_fractional_ideal(number_field_order const & O, std::vector<number_field_element> const & gens);

    /* return the valuation at the given prime ideal
     */
    int valuation(number_field_prime_ideal const & fkp) const;

    number_field_fractional_ideal(number_field_fractional_ideal const & a)
        : O(a.O)
        , ideal_basis_matrix(a.ideal_basis_matrix)
        , denominator(a.denominator)
    {}

    number_field_fractional_ideal(number_field_fractional_ideal && a)
        : O(a.O)
        , ideal_basis_matrix(a.ideal_basis_matrix)
        , denominator(a.denominator)
    {}

    number_field_fractional_ideal& operator=(number_field_fractional_ideal const & a)
    {
        ASSERT_ALWAYS(&O == &a.O);
        ideal_basis_matrix = a.ideal_basis_matrix;
        denominator = a.denominator;
        return *this;
    }
    number_field_fractional_ideal& operator=(number_field_fractional_ideal && a)
    {
        ASSERT_ALWAYS(&O == &a.O);
        ideal_basis_matrix = std::move(a.ideal_basis_matrix);
        denominator = std::move(a.denominator);
        return *this;
    }
};

namespace fmt {
    template <> struct formatter<number_field_fractional_ideal> : formatter<string_view>{
        auto format(number_field_fractional_ideal const & e, format_context& ctx) const -> format_context::iterator;
    };
}
inline std::ostream& operator<<(std::ostream& os, number_field_fractional_ideal const & e) { return os << fmt::format("{}", e); }
#endif	/* UTILS_NUMBERTHEORY_NUMBER_FIELD_FRACTIONAL_IDEAL_HPP_ */
