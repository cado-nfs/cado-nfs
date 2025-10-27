#include "cado.h" // IWYU pragma: keep

#include <utility>
#include <stdexcept>
#include <vector>

#include "fmt/base.h"
#include "fmt/format.h"

#include "cxx_mpz.hpp"
#include "mpz_mat.h"
#include "gmp_aux.h"
#include "numbertheory/number_field.hpp"
#include "numbertheory/number_field_order.hpp"
#include "numbertheory/number_field_order_element.hpp"
#include "numbertheory/number_field_element.hpp"
#include "numbertheory/number_field_fractional_ideal.hpp"
#include "numbertheory/number_field_prime_ideal.hpp"
#include "numbertheory/numbertheory_internals.hpp"
#include "fmt_helper_sagemath.hpp"

number_field_order::number_field_order(class number_field const & K, cxx_mpq_mat mat)
    : K(K)
    , name("$")
    , basis_matrix(std::move(mat))
{
    multiplication_table = numbertheory_internals::multiplication_table_of_order(basis_matrix, K.defining_polynomial());
    mpq_mat_inv(inv_basis_matrix, basis_matrix);
}

std::vector<std::pair<number_field_prime_ideal, int>>
    number_field_order::factor(cxx_mpz const & p, cxx_gmp_randstate & state) const
{
    auto R = numbertheory_internals::factorization_of_prime(basis_matrix,
            number_field().defining_polynomial(),
            p,
            state);

    std::vector<std::pair<number_field_prime_ideal, int>> ret;
    for(auto const & bd : R) {
        number_field_fractional_ideal const I(*this, bd.first);
        number_field_prime_ideal const fkp(I, p, bd.second);
        ret.emplace_back(fkp, bd.second);
    }
    std::sort(ret.begin(), ret.end());
    return ret;
}

std::vector<number_field_prime_ideal>
    number_field_order::factor_radical(cxx_mpz const & p, cxx_gmp_randstate & state) const
{
    std::vector<number_field_prime_ideal> ret;
    for(auto const & F : factor(p, state)) {
        ret.emplace_back(F.first);
    }
    return ret;
}

std::vector<std::pair<number_field_prime_ideal, int>> number_field_order::factor(cxx_mpz const & p) const {
    cxx_gmp_randstate state;
    return factor(p, state);
}

std::vector<number_field_prime_ideal> number_field_order::factor_radical(cxx_mpz const & p) const {
    cxx_gmp_randstate state;
    return factor_radical(p, state);
}

number_field_fractional_ideal number_field_order::fractional_ideal(std::vector<number_field_element> const & gens) const
{
    return { *this, gens };
}

number_field_fractional_ideal number_field_order::p_radical(cxx_mpz const& p) const
{
    return { *this, numbertheory_internals::p_radical_of_order(multiplication_table, p)};
}


number_field_element number_field_order::operator[](int i) const
{
    class number_field const & K(number_field());
    int const n = K.degree();
    cxx_mpq_mat a(1, n);
    mpq_mat_submat_set(a, 0, 0, basis_matrix, i, 0, 1, n);
    return { K, a };
}

std::vector<number_field_element> number_field_order::basis() const
{
    std::vector<number_field_element> ret;
    int const n = K.degree();
    ret.reserve(n);
    for(int i = 0 ; i < n ; i++)
        ret.push_back((*this)[i]);
    return ret;
}


number_field_order_element number_field_order::operator()(cxx_mpz_mat const & a) const
{
    return { *this, a };
}

namespace fmt {
    auto formatter<number_field_order>::format(number_field_order const & O, format_context& ctx) const -> format_context::iterator {
        if (custom_format == SAGEMATH) {
            int const n = O.number_field().degree();
            std::string s = "[";
            for(int i = 0 ; i < n ; i++) {
                if (i) s += ", ";
                s += fmt::format("{}", O[i]);
            }
            s += "]";
            return fmt::format_to(ctx.out(), "{}.order({})", O.number_field().name, s);
        } else if (custom_format == TEXT) {
            return fmt::format_to(ctx.out(), "Order {} in {}", O.name, O.number_field());
        } else {
            throw std::runtime_error("bad format");
        }
    }
}

number_field_order number_field_order::p_maximal_order(cxx_mpz const& p) const
{
    cxx_mpz_poly g;
    mpz_poly_to_monic(g, K.defining_polynomial());
    
    cxx_mpq_mat B = K.basis_matrix_from_f_to_monic(basis_matrix);

    cxx_mpq_mat D = numbertheory_internals::p_maximal_order(B, g, p);

    D = K.basis_matrix_from_monic_to_f(D);

    // Put D into HNF.
    cxx_mpz_mat Dz;
    cxx_mpz den;
    mpq_mat_numden(Dz, den, D);
    mpz_mat_hermite_form_rev(Dz);
    mpq_mat_set_mpz_mat_denom(D, Dz, den);

    return { K, D };
}

cxx_mpz number_field_order::index(number_field_element const & a) const
{
    cxx_mpq_mat c;
    mpq_mat_mul(c, a.coefficients, inv_basis_matrix);
    cxx_mpz denom;
    mpq_mat_numden(nullptr, denom, c);
    return denom;
}

static bool sl_equivalent_matrices(cxx_mpq_mat const& M, cxx_mpq_mat const& A, cxx_mpz const& p)/*{{{*/
{
    /* This is over SL_n(Z_p) */
    if (M->m != A->m) return false;
    if (M->n != A->n) return false;
    cxx_mpq_mat Mi;
    mpq_mat_inv(Mi, M);
    cxx_mpq_mat AMi;
    mpq_mat_mul(AMi, A, Mi);
    /* check that the p-valuation is zero */
    for(unsigned int i = 0 ; i < AMi->m ; i++) {
        for(unsigned int j = 0 ; j < AMi->n ; j++) {
            mpq_srcptr mij = mpq_mat_entry_const(AMi, i, j);
            if (mpz_divisible_p(mpq_denref(mij), p)) return false;
        }
    }
    return true;
}/*}}}*/

bool number_field_order::equal_mod(number_field_order const & O, cxx_mpz const & p) const
{
    return sl_equivalent_matrices(basis_matrix, O.basis_matrix, p);
}


#if 0
bool sl_equivalent_matrices(cxx_mpq_mat const& M, cxx_mpq_mat const& A)/*{{{*/
{
    /* unimplemented for the moment, since unneeded.  We would compute
     * A*M^-1, and see whether we have a denominator. */
    if (M->m != A->m) return false;
    if (M->n != A->n) return false;
    cxx_mpq_mat Mi;
    mpq_mat_inv(Mi, M);
    cxx_mpq_mat AMi;
    mpq_mat_mul(AMi, A, Mi);
    for(unsigned int i = 0 ; i < AMi->m ; i++) {
        for(unsigned int j = 0 ; j < AMi->n ; j++) {
            mpq_srcptr mij = mpq_mat_entry_const(AMi, i, j);
            if (mpz_cmp_ui(mpq_denref(mij), 1) != 0) return false;
        }
    }
    return true;
}/*}}}*/
#endif

