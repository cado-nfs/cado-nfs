#include "cado.h" // IWYU pragma: keep

#include <ios>
#include <istream>

#include "cado_expression_parser.hpp"
#include "cxx_mpfr.hpp"

std::istream & operator>>(std::istream & is, cxx_mpfr::input_with_precision xp)
{
    using cado_expression_parser_details::number_literal;
    using traits = cado_expression_parser_details::number_traits<cxx_mpfr>;
    number_literal N;
    if (is >> N) {
        try {
            xp.x = traits::from_number_literal(N, xp.p);
        } catch (cado_expression_parser_details::parse_error const &) {
            is.setstate(std::ios::failbit);
        }
    }
    return is;
}
