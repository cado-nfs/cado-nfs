#ifndef UTILS_NUMBER_LITERAL_HPP_
#define UTILS_NUMBER_LITERAL_HPP_

#include <cctype>

#include <string>
#include <istream>
#include <ios>

#include "cado_parsing_base.hpp"

/*
 * https://en.cppreference.com/w/cpp/language/floating_literal
 *
 * we also take inspiration from
 *
 * https://en.cppreference.com/w/cpp/numeric/complex/operator%2522%2522i.html
 */

namespace cado {

struct number_literal {
    std::string full;
    bool is_hex = false;
    bool has_point = false;
    bool has_exponent = false;
    bool is_imaginary = false;
    // long exponent = false;
    // std::string::size_type begin_integral = 0,
    // std::string::size_type begin_fractional = 0, end_fractional = 0;
    // std::string::size_type begin_exponent = 0, end_exponent = 0;
    private:
    std::string::size_type begin_exponent = 0;
    std::string::size_type end_integral = 0;
    public:
    number_literal strip_i() const {
        if (!is_imaginary) return *this;
        number_literal N = *this;
        N.is_imaginary = false;
        if (N.full.back() != 'i') throw cado::token_error();
        N.full.erase(N.full.size()-1);
        if (!N.full.empty() && std::isspace(N.full.back()))
            N.full.erase(N.full.size()-1);
        return N;
    }
    std::string integral_part() const {
        return full.substr(0, end_integral);
    }
    std::string exponent_part() const {
        return full.substr(begin_exponent);
    }
    bool isdigit(int c) const {
        return std::isdigit(c) ||
            (is_hex && ((c >= 'a' && c <= 'f') || (c >= 'A' && c <= 'F')));
    }
    /* I don't think we want to recognize the f and l suffix, but it
     * could certainly be done
     */
    long as_exponent() const {
        if (is_hex || has_exponent || has_point)
            throw token_error();
        return std::stol(exponent_part());
    }
    public:
    /* recognize a number literal, return true if successful, and is is
     * stopped at the first character that is not part of the recognized
     * literal. [is] may or may not be in eof state on return, it does not
     * affect the boolean return value of this function.
     */
    static bool recognize(number_literal & N, std::istream& is);
};

inline std::istream& operator>>(std::istream& is, number_literal & N)
{
    if (!number_literal::recognize(N, is))
        is.setstate(std::ios::failbit);
    return is;
}

} /* namespace cado */

#endif	/* UTILS_NUMBER_LITERAL_HPP_ */
