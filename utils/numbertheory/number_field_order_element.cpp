#include "cado.h" // IWYU pragma: keep

#include "mpz_mat.h"
#include "cxx_mpz.hpp"
#include "numbertheory/number_field_order_element.hpp"
#include "numbertheory/number_field_element.hpp"
#include "numbertheory/number_field_order.hpp"
#include "numbertheory/number_field.hpp"
#include "numbertheory/numbertheory_internals.hpp"

number_field_order_element::number_field_order_element(number_field_order const & O, number_field_element const & e)
    : O(O)
{
    cxx_mpq_mat c;
    mpq_mat_mul(c, e.coefficients, O.inv_basis_matrix);
    cxx_mpz denom;
    mpq_mat_numden(coefficients, denom, c);
    if (denom != 1) {
        throw element_not_integral();
    }
}

number_field_order_element number_field_order_element::operator*(number_field_order_element const & a) const
{
    if (&order() != &a.order())
        throw number_field_inconsistency();

    return { order(),
            numbertheory_internals::multiply_elements_in_order(order().multiplication_table,
                coefficients,
                a.coefficients) };
}

cxx_mpz_mat number_field_order_element::multiplication_matrix() const
{
    number_field_order const & O(order());
    /* return the matrix of the multiplication by e in O */
    int const n = O.number_field().degree();
    cxx_mpz_mat eM;
    mpz_mat_mul(eM, coefficients, O.multiplication_table);
    cxx_mpz_mat ret(n, n);
    for(int i = 0 ; i < n ; i++)
        mpz_mat_submat_swap(ret, i, 0, eM, 0, n*i, 1, n);
    return ret;
}


namespace fmt {
    auto formatter<number_field_order_element>::format(number_field_order_element const & e, format_context& ctx) const -> format_context::iterator
    {
        if (custom_format == SAGEMATH || custom_format == TEXT) {
            fmt::format_to(ctx.out(), "{}([", e.order().name);
            for(int i = 0 ; i < e.order().number_field().degree() ; i++) {
                if (i) fmt::format_to(ctx.out(), ", ");
                fmt::format_to(ctx.out(), "{}", * (cxx_mpz const *) e.coefficients(0,i));
            }
            fmt::format_to(ctx.out(), "])");
        } else if (custom_format == MACHINE) {
            for(int i = 0 ; i < e.order().number_field().degree() ; i++) {
                if (i) fmt::format_to(ctx.out(), " ");
                fmt::format_to(ctx.out(), "{}", * (cxx_mpz const *) e.coefficients(0,i));
            }
        }
        return ctx.out();
    }
}
