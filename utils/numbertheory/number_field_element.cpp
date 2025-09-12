#include "cado.h" // IWYU pragma: keep

#include <utility>

#include <gmp.h>
#include "fmt/base.h"

#include "cxx_mpz.hpp"
#include "numbertheory/number_field_element.hpp"
#include "numbertheory/number_field.hpp"
#include "fmt_helper_sagemath.hpp"
#include "mpz_mat_accessors.h"
#include "mpz_poly.h"

number_field_element::number_field_element(class number_field const & K, cxx_mpz_poly const & a, cxx_mpz const & d)
    : K(K)
{
    cxx_mpz_poly r;
    cxx_mpz_poly const & f(K.defining_polynomial());
    cxx_mpz denom;
    int const n = f.degree();

    if (a.degree() >= n) {
        cxx_mpz c;
        mpz_poly_pseudo_remainder(r, a, f);
        mpz_pow_ui(denom, mpz_poly_lc(f), a.degree() - n + 1);
        mpz_poly_content(c, r);
        mpz_gcd(c, c, denom);
        mpz_poly_divexact_mpz(r, r, c);
        mpz_divexact(denom, denom, c);
        mpz_mul(denom, denom, d);
    } else {
        r = a;
        denom = d;
    }

    coefficients = cxx_mpq_mat(1, n);
    for(int i = 0 ; i < n ; i++) {
        mpq_ptr ci = coefficients(0, i);
        mpz_set(mpq_numref(ci), mpz_poly_coeff_const(r, i));
        mpz_set(mpq_denref(ci), denom);
        mpq_canonicalize(ci);
    }
}

number_field_element::number_field_element(class number_field const & K, cxx_mpz_mat const & a, cxx_mpz const & d)
    : K(K)
    , coefficients(1, K.degree())
{
    int const n = K.degree();
    for(int i = 0 ; i < n ; i++) {
        mpq_ptr ci = coefficients(0, i);
        mpz_set(mpq_numref(ci), a(0, i));
        mpz_set(mpq_denref(ci), d);
        mpq_canonicalize(ci);
    }
}

std::pair<cxx_mpz_poly, cxx_mpz> number_field_element::as_polynomial() const
{
    std::pair<cxx_mpz_poly, cxx_mpz> ret;
    mpq_mat_row_to_poly(ret.first, ret.second, coefficients, 0);
    return ret;
}

/* return the matrix of the multiplication by e in K */
cxx_mpq_mat number_field_element::multiplication_matrix() const
{
    class number_field const & K(number_field());
    cxx_mpz_poly const & f(K.defining_polynomial());
    int const n = K.degree();
    cxx_mpq_mat M(n, n);
    cxx_mpz denom;
    cxx_mpz_poly a;
    std::tie(a, denom) = as_polynomial();
    mpq_poly_to_mat_row(M, 0, a, denom);
    for(int i = 1 ; i < n ; i++) {
        mpz_poly_mul_xi(a, a, 1);
        if (a.degree() == n) {
            mpz_poly_pseudo_remainder(a, a, f);
            mpz_mul(denom, denom, mpz_poly_lc(f));
        }
        mpq_poly_to_mat_row(M, i, a, denom);
    }
    return M;
}


cxx_mpq number_field_element::trace() const
{
    cxx_mpq_mat const & T(number_field().trace_matrix());
    int const n = number_field().degree();
    cxx_mpq t = 0, c;
    for(int i = 0 ; i < n ; i++) {
        mpq_mul(c, coefficients(0, i), T(0, i));
        mpq_add(t, t, c);
    }
    return t;
}

number_field_element number_field_element::operator*(number_field_element const & o) const
{
    class number_field const & K(number_field());
    cxx_mpz_poly const & f(K.defining_polynomial());
    int const n = K.degree();
    cxx_mpz_poly a,b,c;
    cxx_mpz ad,bd,cd;
    std::tie(a, ad) = as_polynomial();
    std::tie(b, bd) = o.as_polynomial();
    mpz_poly_mul(c, a, b);
    mpz_mul(cd, ad, bd);
    int const k = c.degree() - n;
    if (k >= 0) {
        mpz_poly_pseudo_remainder(c, c, f);
        cxx_mpz lck;
        mpz_pow_ui(lck, mpz_poly_lc(f), k+1);
        mpz_mul(cd, cd, lck);
    }
    number_field_element ret(K, c);
    for(int i = 0 ; i < n ; i++) {
        mpq_ptr ri = ret.coefficients(0, i);
        mpz_mul(mpq_denref(ri), mpq_denref(ri), cd);
        mpq_canonicalize(ri);
    }
    return ret;
}

number_field_element number_field_element::operator*(cxx_mpz const & m) const
{
    cxx_mpq_mat cm;
    mpq_mat_mul_mpz(cm, coefficients, m);
    return { K, cm };
}

number_field_element number_field_element::operator/(cxx_mpz const & m) const
{
    cxx_mpq_mat cm;
    mpq_mat_div_mpz(cm, coefficients, m);
    return { K, cm };
}

auto fmt::formatter<number_field_element>::format(number_field_element const & e, format_context& ctx) const -> format_context::iterator
{
    auto y = e.as_polynomial();
    if (custom_format == SAGEMATH || custom_format == TEXT) {
        if (y.second == 1) {
            fmt::format_to(ctx.out(), "{}", y.first.print_poly(e.number_field().varname));
        } else {
            fmt::format_to(ctx.out(), "({})/{}", y.first.print_poly(e.number_field().varname), y.second);
        }
    } else if (custom_format == MACHINE) {
        cxx_mpz_mat mz;
        cxx_mpz d;
        mpq_mat_numden(mz, d, e.coefficients);
        fmt::format_to(ctx.out(), "{}", d);
        for(int i = 0 ; i < e.number_field().degree() ; i++)
            fmt::format_to(ctx.out(), " {}", cxx_mpz(mz(0,i)));
    }
    return ctx.out();
}
