#define __STDCPP_MATH_SPEC_FUNCS__ 201003L
#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1       /* for expint() */
#include "cado.h" // IWYU pragma: keep
#include <cmath>
#include <sstream>
#include "cxx_mpz.hpp"
#include "cado_expression_parser.hpp"

#include "misc.h"

double nprimes_interval(double p0, double p1)
{
#ifdef HAVE_STDCPP_MATH_SPEC_FUNCS
    return std::expint(log(p1)) - std::expint(log(p0));
#else
    /* that can't be sooo wrong... */
    double l0 = log(p0);
    double l1 = log(p1);
    double s1 = p1*(1/l1+1/pow(l1,2)+2/pow(l1,3)+6/pow(l1,4));
    double s0 = p0*(1/l0+1/pow(l0,2)+2/pow(l0,3)+6/pow(l0,4));
    return s1 - s0;
#endif
}

struct mpz_parser_traits {
    static constexpr const int accept_literals = 0;
    typedef cxx_mpz type;
    void add(cxx_mpz & c, cxx_mpz const & a, cxx_mpz const & b) {
        mpz_add(c, a, b);
    }
    void sub(cxx_mpz & c, cxx_mpz const & a, cxx_mpz const & b) {
        mpz_sub(c, a, b);
    }
    void mul(cxx_mpz & c, cxx_mpz const & a, cxx_mpz const & b) {
        mpz_mul(c, a, b);
    }
    void pow_ui(cxx_mpz & c, cxx_mpz const & a, unsigned long e) {
        mpz_pow_ui(c, a, e);
    }
    void swap(cxx_mpz & a, cxx_mpz & b) {
        mpz_swap(a, b);
    }
    void set_mpz(cxx_mpz & a, cxx_mpz const & z) {
        mpz_set(a, z);
    }
    void set_literal_power(cxx_mpz &, char, unsigned long) {
        // never called. we could do some gymnastics to statically elide
        // this call, but that does not seem to be worth it.
    }
};

typedef cado_expression_parser<mpz_parser_traits> integer_parser;

int mpz_set_from_expression(mpz_ptr f, const char * value)
{
    std::istringstream is(value);

    integer_parser P;
    try {
        P.tokenize(is);
        cxx_mpz tmp = P.parse();
        mpz_set(f, tmp);
    } catch (integer_parser::token_error const & p) {
        return 0;
    } catch (integer_parser::parse_error const & p) {
        return 0;
    }
    return 1;
}
