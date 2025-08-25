#ifndef UTILS_NUMBER_LITERAL_HPP_
#define UTILS_NUMBER_LITERAL_HPP_

#include "cado_config.h"        // IWYU pragma: keep

#include <string>
#include <istream>
#include <ios>
#include <cctype>

/*
 * https://en.cppreference.com/w/cpp/language/floating_literal
 */

namespace cado {

struct number_literal {
    struct parse_error: public std::exception {
        const char * what() const noexcept override { return "parse error"; }
    };
    struct token_error: public std::exception {
        const char * what() const noexcept override { return "token error"; }
    };
    std::string full;
    bool is_hex = false;
    bool has_point = false;
    bool has_exponent = false;
    // long exponent = false;
    // std::string::size_type begin_integral = 0,
    // std::string::size_type begin_fractional = 0, end_fractional = 0;
    // std::string::size_type begin_exponent = 0, end_exponent = 0;
    private:
    std::string::size_type begin_exponent = 0;
    std::string::size_type end_integral = 0;
    public:
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

    bool swallow_hits_end(int & c, std::istream & is) {
        full += static_cast<char>(c);
        is.get();
        c = is.peek();
        return is.eof();
    }
    public:
    static bool recognize(number_literal & N, std::istream& is)// {{{
    {
        N = number_literal {};
        int c = is.peek();
        bool no_mantissa = true;

        if (c == '-')
            if (N.swallow_hits_end(c, is))
                throw token_error();

        if (c == 'n' || c == 'N') {
            if (N.swallow_hits_end(c, is))
                throw token_error();
            if (c != 'a' && c != 'A')
                throw token_error();
            if (N.swallow_hits_end(c, is))
                throw token_error();
            if (c != 'n' && c != 'N') 
                throw token_error();
            N.swallow_hits_end(c, is);
            return true;
        }

        if (c == 'i') {
            if (N.swallow_hits_end(c, is))
                throw token_error();
            if (c != 'n')
                throw token_error();
            if (N.swallow_hits_end(c, is))
                throw token_error();
            if (c != 'f')
                throw token_error();
            N.swallow_hits_end(c, is);
            return true;
        }

        // std::string * ps = &N.integral;
        if (!std::isdigit(c)) {
            if (c == '.') {
                N.end_integral = N.full.size();
                N.has_point = true;
                if (N.swallow_hits_end(c, is)) throw token_error();
                // *ps = "0";
                // ps = &N.fractional;
            } else {
                return false;
            }
        } else if (c == '0') {
            if (N.swallow_hits_end(c, is)) return true;
            if (std::tolower(c) == 'x') {
                N.is_hex = true;
                if (N.swallow_hits_end(c, is)) throw token_error();
            } else {
                no_mantissa = false;
            }
        }

        for( ; N.isdigit(c) || (!N.has_point && c == '.') ; ) {
            if (c == '.') N.has_point = true;
            no_mantissa = false;
            if (!N.has_point)
                N.end_integral = N.full.size() + 1;
            if (N.swallow_hits_end(c, is)) return true;
        }
        if (no_mantissa)
            throw token_error();

        if (std::tolower(c) == (N.is_hex ? 'p' : 'e')) {
            N.has_exponent = true;
            N.begin_exponent = N.full.size() + 1;
            if (N.swallow_hits_end(c, is)) return true;
            if (c == '+' || c == '-')
                if (N.swallow_hits_end(c, is)) return true;
            for( ; std::isdigit(c) ; )
                if (N.swallow_hits_end(c, is)) return true;
        }
        return true;
    }// }}}
};

inline std::istream& operator>>(std::istream& is, number_literal & N)
{
    if (!number_literal::recognize(N, is))
        is.setstate(std::ios::failbit);
    return is;
}

} /* namespace cado */

#endif	/* UTILS_NUMBER_LITERAL_HPP_ */
