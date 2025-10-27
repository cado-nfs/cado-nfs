#include "cado.h" // IWYU pragma: keep

#include <climits>

#include <memory>
#include <utility>

#include "fmt/base.h"

#include "cxx_mpz.hpp"
#include "mpz_mat.h"
#include "gmp_aux.h"
#include "numbertheory/number_field.hpp"
#include "numbertheory/number_field_element.hpp"        // IWYU pragma: keep
#include "numbertheory/number_field_fractional_ideal.hpp"
#include "numbertheory/number_field_order.hpp"
#include "fmt_helper_sagemath.hpp"
#include "numbertheory/number_field_prime_ideal.hpp"
#include "numbertheory/numbertheory_internals.hpp"

number_field_prime_ideal::operator two_element() const
{
    number_field_order const & O = order();
    if (!cached_two_element) {
        auto T = numbertheory_internals::prime_ideal_two_element(O.basis_matrix, O.number_field().defining_polynomial(), O.multiplication_table, ideal_basis_matrix);
        cached_two_element = std::make_unique<two_element>(T.first, O(T.second));
    }
    return *cached_two_element;
}

int number_field_prime_ideal::valuation(number_field_fractional_ideal const & I) const
{
    if (&order() != &I.order())
        throw number_field_inconsistency();

    number_field_order const & O(order());
    cxx_mpz_mat const & M(O.multiplication_table);
    cxx_mpz_mat const & a(valuation_helper.coefficients);

    int const v = numbertheory_internals::valuation_of_ideal_at_prime_ideal(M, I.ideal_basis_matrix, a, p);
    if (v == INT_MAX) return v;

    int const w = mpz_p_valuation(I.denominator, p);
    return v - w * e;

}

number_field_prime_ideal::number_field_prime_ideal(number_field_fractional_ideal I, cxx_mpz p, int e)
    : number_field_fractional_ideal(std::move(I))
    , valuation_helper(
            order(),
            numbertheory_internals::valuation_helper_for_ideal(
                order().multiplication_table,
                ideal_basis_matrix,
                p))
    , p(std::move(p))
    , e(e)
{
}

int number_field_prime_ideal::inertia_degree() const
{
    return numbertheory_internals::prime_ideal_inertia_degree(ideal_basis_matrix);
}

namespace fmt {
    auto formatter<number_field_prime_ideal>::format(number_field_prime_ideal const & I, format_context& ctx) const -> format_context::iterator
    {
        number_field_prime_ideal::two_element uv(I);
        number_field_order const & O(I.order());
        number_field const & K(O.number_field());

        if (custom_format == SAGEMATH || custom_format == TEXT) {
            fmt::format_to(ctx.out(), "{}.fractional_ideal([{}, {}])", O.name, uv.first, K(uv.second));
        } else if (custom_format == MACHINE) {
            fmt::format_to(ctx.out(), "{} {:M}",
                    uv.first,
                    K(uv.second));
        }
        return ctx.out();
    }
}
