#include "cado.h" // IWYU pragma: keep

#include <sstream>
#include <string>

#include <gmp.h>

#include "cado_expression_parser.hpp"
#include "gmp_aux.h"
#include "cado_mp_conversions.hpp"

/* This compilation unit is a complement to gmp_aux.c ; some of the
 * utilities there are actually implemented by c++ code
 */
long double
mpz_get_ld (mpz_srcptr z)
{
    return cado_math_aux::mpz_get<long double>(z);
}

struct mpz_parser_traits {
    static constexpr const int accept_literals = 0;
    using type = cxx_mpz;
    using number_type = cxx_mpz;
    static void add(cxx_mpz & c, cxx_mpz const & a, cxx_mpz const & b) {
        mpz_add(c, a, b);
    }
    static void sub(cxx_mpz & c, cxx_mpz const & a, cxx_mpz const & b) {
        mpz_sub(c, a, b);
    }
    static void neg(cxx_mpz & c, cxx_mpz const & a) {
        mpz_neg(c, a);
    }
    static void mul(cxx_mpz & c, cxx_mpz const & a, cxx_mpz const & b) {
        mpz_mul(c, a, b);
    }
    static void pow(cxx_mpz & c, cxx_mpz const & a, unsigned long e) {
        mpz_pow_ui(c, a, e);
    }
    static void swap(cxx_mpz & a, cxx_mpz & b) {
        mpz_swap(a, b);
    }
    static void set(cxx_mpz & a, cxx_mpz const & z) {
        mpz_set(a, z);
    }
    static void set_literal_power(cxx_mpz &, std::string const&, unsigned long) {
        // never called. we could do some gymnastics to statically elide
        // this call, but that does not seem to be worth it.
    }
};

using integer_parser = cado_expression_parser<mpz_parser_traits>;

cxx_mpz mpz_from_expression(std::string const & value)
{
    std::istringstream is(value);
    integer_parser P;
    P.tokenize(is);
    return P.parse();
}

int mpz_set_from_expression(mpz_ptr f, const char * value)
{
    try {
        mpz_set(f, mpz_from_expression(value));
    } catch (cado_expression_parser_details::token_error const & p) {
        return 0;
    } catch (cado_expression_parser_details::parse_error const & p) {
        return 0;
    }
    return 1;
}

bool cado::params::parser<mpz_t>::operator()(std::string const & s, mpz_t value) const
{
    cxx_mpz z;
    const bool r = parse(s, z);
    if (r)
        mpz_set(value, z);
    return r;
}
