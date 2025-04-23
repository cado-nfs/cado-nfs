#include "cado.h" // IWYU pragma: keep

#include <cstddef>

#include <vector>
#include <tuple>

#include "fmt/base.h"

#include "mpz_mat.h"
#include "numbertheory/number_field.hpp"
#include "numbertheory/number_field_fractional_ideal.hpp"
#include "numbertheory/number_field_prime_ideal.hpp"
#include "numbertheory/number_field_element.hpp"
#include "numbertheory/number_field_order.hpp"
#include "numbertheory/numbertheory_internals.hpp"
#include "runtime_numeric_cast.hpp"

int number_field_fractional_ideal::valuation(number_field_prime_ideal const & fkp) const
{
    return fkp.valuation(*this);
}


number_field_fractional_ideal::number_field_fractional_ideal(number_field_order const & O, std::vector<number_field_element> const & gens)
    : O(O)
{
    int const n = O.number_field().degree();
    cxx_mpq_mat G(runtime_numeric_cast<int>(gens.size()), n);
    for(size_t i = 0 ; i < gens.size() ; i++)
        mpq_mat_submat_set(G, i, 0, gens[i].coefficients, 0, 0, 1, n);
    std::tie(ideal_basis_matrix, denominator) = numbertheory_internals::generate_ideal(O.basis_matrix, O.multiplication_table, G);
}

namespace fmt {
    auto formatter<number_field_fractional_ideal>::format(number_field_fractional_ideal const & I, format_context& ctx) const -> format_context::iterator
    {
        number_field_order const & O(I.order());
        int const n = O.number_field().degree();
        number_field const & K(O.number_field());

        if (I.denominator != 1) fmt::format_to(ctx.out(), "(");
        fmt::format_to(ctx.out(), "{}.fractional_ideal([", O.name);
        cxx_mpz_mat a(1, n);
        for(int i = 0 ; i < n ; i++) {
            if (i) fmt::format_to(ctx.out(), ", ");
            mpz_mat_submat_set(a, 0, 0, I.ideal_basis_matrix, i, 0, 1, n);
            fmt::format_to(ctx.out(), "{}", K(O(a)));
        }
        fmt::format_to(ctx.out(), "])");
        if (I.denominator != 1)
            if (I.denominator != 1) fmt::format_to(ctx.out(), "/{})", I.denominator);
        return ctx.out();
    }
}
