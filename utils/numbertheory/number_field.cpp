#include "cado.h" // IWYU pragma: keep

#include <memory>

#include <gmp.h>
#include "fmt/base.h"
#include "fmt/format.h"

#include "cxx_mpz.hpp"
#include "macros.h"
#include "misc.h"
#include "mpz_mat.h"
#include "mpz_poly.h"
#include "numbertheory/fmt_helpers.hpp"
#include "numbertheory/number_field.hpp"
#include "numbertheory/number_field_element.hpp"
#include "numbertheory/number_field_order.hpp"
#include "numbertheory/number_field_order_element.hpp"
#include "numbertheory/numbertheory_internals.hpp"

number_field::number_field(cxx_mpz_poly const & f)
    : f(f)
    , name("$")
    , varname("$.1")
{
    mpz_poly_to_monic(f_hat, f);
}

void number_field::bless(std::string const & name, std::string const & varname)
{
    this->name = name;
    this->varname = varname;
}

void number_field::bless(std::string const & name)
{
    this->name = name;
    this->varname = fmt::format("{}.1", name);
}

void number_field::bless(const char * name, const char * varname)
{
    this->name = name;
    if (varname)
        this->varname = varname;
    else
        this->varname = fmt::format("{}.1", name);
}

number_field_element number_field::gen() const {
    return { *this, "x" };
}

number_field_order number_field::order(number_field_element const & a) const
{
    /* Here we probably want to make sure that a is integral */
    ASSERT_ALWAYS(0);
    number_field_element x = (*this)(1);
    int const n = degree();
    cxx_mpq_mat B(n, n);
    for(int i = 0; i < n; i++) {
        mpq_set_ui(B(0,i), i == 0, 1);
    }
    for(int i = 1; i < n; i++) {
        x = x * a;
        mpq_mat_submat_set(B, i, 0, x.coefficients, 0, 0, 1, n);
    }
    return { *this, std::move(B) };
}

number_field_order number_field::equation_order() const
{
    cxx_mpz x;
    int const n = degree();
    mpz_set_ui(x, 1);
    cxx_mpq_mat B(n, n);
    mpq_mat_set_ui(B, 1);
    for(int i = 0; i < n; i++) {
        mpq_set_z(B(i,i), x);
        mpz_mul(x, x, mpz_poly_lc(f));
    }
    return { *this, std::move(B) };
}

number_field_order number_field::p_maximal_order(cxx_mpz const & p) const
{
    return { *this, numbertheory_internals::p_maximal_order(defining_polynomial(), p) };
}

static cxx_mpq_mat companion_matrix(cxx_mpz_poly const & f)
{
    int const n = f.degree();
    cxx_mpq_mat M(n, n);
    for(int i = 0 ; i < n-1 ; i++)
        mpq_set_ui(M(i, i+1), 1, 1);
    for(int j = 0 ; j < n ; j++) {
        mpz_neg(mpq_numref(M(n-1, j)), mpz_poly_coeff_const(f, j));
        mpz_set(mpq_denref(M(n-1, j)), mpz_poly_lc(f));
        mpq_canonicalize(M(n-1, j));
    }
    return M;
}

cxx_mpq_mat number_field::trace_matrix() const
{
    if (cached_trace_matrix == nullptr) {
        cxx_mpz_poly const & f(defining_polynomial());
        int const n = degree();
        cxx_mpq_mat C = companion_matrix(f);
        cxx_mpq_mat M = C;
        cxx_mpq_mat T(1, n);
        mpq_set_ui(T(0, 0), n, 1);
        if (n == 1)
            return T;
        mpq_set(T(0, 1), C(n-1, n-1));
        for(int i = 2 ; i < n ; i++) {
            mpq_mat_mul(M, M, C);
            mpq_mat_trace(T(0, i), M);
        }
        cached_trace_matrix = std::unique_ptr<cxx_mpq_mat>(new cxx_mpq_mat(T));
    }
    return *cached_trace_matrix;
}


number_field_order const& number_field::maximal_order(unsigned long prime_limit) const
{
    if (cached_maximal_order == nullptr) {

        number_field_order O = equation_order();

        cxx_mpz disc;
        mpz_poly_discriminant(disc, f);
        mpz_mul(disc, disc, mpz_poly_lc(f));

        /* We're not urged to use ecm here */
        for(auto const & pe : trial_division(disc, prime_limit, disc)) {
            fmt::print("{} {}\n", pe.first, O);
        }

        ASSERT_ALWAYS(0);

        cached_maximal_order = std::unique_ptr<number_field_order>(new number_field_order(O));
    }

    return *cached_maximal_order;
}

number_field_element number_field::operator()(cxx_mpz_poly const & a, cxx_mpz const & d) const
{
    return { *this, a, d };
}
number_field_element number_field::operator()(cxx_mpq_mat const & a) const
{
    return { *this, a };
}
number_field_element number_field::operator()(cxx_mpz_mat const & a, cxx_mpz const & d) const
{
    return { *this, a, d };
}

number_field_element number_field::operator()(number_field_order_element const & e) const
{
    ASSERT_ALWAYS(this == &e.order().number_field());
    cxx_mpq_mat c = e.coefficients;
    mpq_mat_mul(c, c, e.order().basis_matrix);
    return { *this, c };
}


auto fmt::formatter<number_field>::format(number_field const & K, format_context& ctx) const -> format_context::iterator
{
    if (custom_format == TEXT) {
        fmt::format_to(ctx.out(), "Number Field {} in variable {} defined by {}", K.name, K.varname, K.defining_polynomial());
    } else if (custom_format == SAGEMATH) {
        fmt::format_to(ctx.out(), "NumberField({}, name=(\"{}\",))", K.defining_polynomial(), K.varname);
    } else if (custom_format == MACHINE) {
        fmt::format_to(ctx.out(), "{}", K.defining_polynomial());
    }
    return ctx.out();
}
