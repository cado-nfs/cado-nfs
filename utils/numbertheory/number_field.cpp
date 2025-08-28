#include "cado.h" // IWYU pragma: keep

#include <memory>
#include <utility>
#include <stdexcept>

#include <gmp.h>
#include "fmt/base.h"
#include "fmt/format.h"

#include "cxx_mpz.hpp"
#include "macros.h"
#include "misc.h"
#include "mpz_mat.h"
#include "mpz_poly.h"

#include "fmt_helper_sagemath.hpp"
#include "numbertheory/number_field.hpp"
#include "numbertheory/number_field_element.hpp"
#include "numbertheory/number_field_order.hpp"
#include "numbertheory/number_field_order_element.hpp"

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
    bless(name, fmt::format("{}.1", name));
}

void number_field::bless(const char * name, const char * varname)
{
    bless(std::string(name), varname ? std::string(varname) : fmt::format("{}.1", name));
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
    cxx_mpq_mat B(degree(), degree());
    mpq_mat_set_ui(B, 1);
    return { *this, basis_matrix_from_monic_to_f(B) };
}

number_field_order number_field::p_maximal_order(cxx_mpz const & p) const
{
    return equation_order().p_maximal_order(p);
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
        cached_trace_matrix = std::make_unique<cxx_mpq_mat>(T);
    }
    return *cached_trace_matrix;
}

template<typename iterator>
static number_field_order maximize_recursively(number_field_order && O, cxx_mpz const & disc, iterator begin, iterator end)
{
    if (begin == end)
        return std::move(O);

    auto p = begin->first;

    if (mpz_p_valuation(disc, p) == 1)
        return maximize_recursively(std::move(O), disc, ++begin, end);
    else
        return maximize_recursively(O.p_maximal_order(p), disc, ++begin, end);
}

/* the maximal order in itself isn't necessarily something very useful.
 * And anyway we only compute an approximation of it, given that we don't
 * expect that we'll factor the discriminant completely.
 *
 * on top of that, we have an implementation difficulty caused by the
 * fact that assigning to a number_field_order object isn't supported.
 * Beyond the recursive kludge above, alternatives include changing
 * references to shared_ptr's after all, or make p_maximal_order an
 * in-place operation.
 */
number_field_order const& number_field::maximal_order(unsigned long prime_limit) const
{
    if (cached_maximal_order == nullptr) {
        cxx_mpz disc, cofac;
        mpz_poly_discriminant(disc, f);
        mpz_mul(disc, disc, mpz_poly_lc(f));

        /* We're not urged to use ecm here */
        auto d_fac = trial_division(disc, prime_limit, cofac);

        cached_maximal_order = std::make_unique<number_field_order>(
                    maximize_recursively(equation_order(),
                        disc,
                        d_fac.begin(),
                        d_fac.end()));

    }

    return *cached_maximal_order;
}

number_field_element number_field::operator()(cxx_mpz_poly const & a, cxx_mpz const & d) const
{
    return { *this, a, d };
}
number_field_element number_field::operator()(cxx_mpq_mat const & a) const
{
    if (a.ncols() != static_cast<unsigned int>(degree()))
        throw std::out_of_range("wrong number of coefficients");
    return { *this, a };
}
number_field_element number_field::operator()(cxx_mpz_mat const & a, cxx_mpz const & d) const
{
    if (a.ncols() != static_cast<unsigned int>(degree()))
        throw std::out_of_range("wrong number of coefficients");
    return { *this, a, d };
}

number_field_element number_field::operator()(number_field_order_element const & e) const
{
    ASSERT_ALWAYS(this == &e.order().number_field());
    cxx_mpq_mat c = e.coefficients;
    mpq_mat_mul(c, c, e.order().basis_matrix);
    return { *this, c };
}

cxx_mpq_mat number_field::basis_matrix_from_f_to_monic(cxx_mpq_mat const & B) const
{
    /* given the basis of some Z-lattice in K that is expressed with
     * respect to the polynomial basis defined by f, return a basis of
     * the same Z-lattice, but with respect to the polynomial basis
     * defined by make_monic(f)
     *
     * For example, if K=Q(alpha) with alpha a root of f = ell*x^2-1, and
     * B is [a,b,c,d] representing the Z-lattice with generating elements
     * a+b*alpha and c+d*alpha, then since g = x^2-ell has the root
     * alpha_hat = ell*alpha, we return the matrix [a, b/ell, c, d/ell]
     */
    if (mpz_poly_is_monic(defining_polynomial()))
        return B;

    unsigned int const n = degree();

    cxx_mpq_mat C = B;
    cxx_mpz x;
    mpz_set_ui(x, 1);
    for(unsigned int j = 0 ; j < n ; j++) {
        for(unsigned int i = 0; i < n; i++) {
            mpq_ptr dij = mpq_mat_entry(C, i, j);
            mpz_mul(mpq_denref(dij), mpq_denref(dij), x);
            mpq_canonicalize(dij);
        }
        mpz_mul(x, x, mpz_poly_lc(defining_polynomial()));
    }
    return C;
}

cxx_mpq_mat number_field::basis_matrix_from_monic_to_f(cxx_mpq_mat const & B) const
{
    /* does the converse of basis_matrix_from_f_to_monic */
    if (mpz_poly_is_monic(defining_polynomial()))
        return B;

    unsigned int const n = degree();

    cxx_mpq_mat C = B;
    cxx_mpz x;
    mpz_set_ui(x, 1);
    for(unsigned int j = 0 ; j < n ; j++) {
        for(unsigned int i = 0; i < n; i++) {
            mpq_ptr dij = mpq_mat_entry(C, i, j);
            mpz_mul(mpq_numref(dij), mpq_numref(dij), x);
            mpq_canonicalize(dij);
        }
        mpz_mul(x, x, mpz_poly_lc(defining_polynomial()));
    }
    return C;
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

std::pair<unsigned int, unsigned int> number_field::signature() const
{
    if (!cached_signature) {
        int const r1 = mpz_poly_number_of_real_roots(defining_polynomial());
        ASSERT_ALWAYS((degree() - r1) % 2 == 0);
        int const r2 = (degree() - r1) / 2;
        cached_signature = std::make_unique<std::pair<unsigned int, unsigned int>>(r1, r2);
    }
    return *cached_signature;
}
