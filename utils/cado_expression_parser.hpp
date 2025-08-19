#ifndef CADO_EXPRESSION_PARSER_HPP
#define CADO_EXPRESSION_PARSER_HPP

#include "cado_config.h"

#include <cstddef>

#include <algorithm>
#include <cctype>
#include <exception>
#include <istream>
#include <ios>
#include <string>
#include <vector>
#include <utility>
#include <functional>

#include <gmp.h>

#include "cxx_mpz.hpp"
// #include "cado_math_aux.hpp"

#ifdef HAVE_MPFR
#include <mpfr.h>
#include "cxx_mpfr.hpp"
#endif

/* This structure has only two public functions: tokenize and parse */

/* A typical use can be found in utils/mpz_poly.cpp */

namespace cado_expression_parser_details {
    struct parse_error: public std::exception {
        const char * what() const noexcept override { return "parse error"; }
    };
    struct token_error: public std::exception {
        const char * what() const noexcept override { return "token error"; }
    };

    /*
     * https://en.cppreference.com/w/cpp/language/floating_literal
     */

    struct number_literal {
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

    struct cado_expression_parser_base {
        // typedef cado_expression_parser_details::token_error token_error;
        // typedef cado_expression_parser_details::parse_error parse_error;
        protected:
        enum expression_token {
            LEFT_PAREN,
            RIGHT_PAREN,
            PLUS,
            MINUS,
            TIMES,
            POWER,
            POSITIVE_NUMBER,
            LITERAL
        };
        std::vector<expression_token> tokens;
        std::vector<std::string> literals;
        std::vector<cado_expression_parser_details::number_literal> numbers;

        typename std::vector<expression_token>::const_iterator ctok;
        std::vector<std::string>::const_iterator clit;
        std::vector<cado_expression_parser_details::number_literal>::const_iterator cnumber;

        void next() { if (ctok != tokens.end()) ctok++; }
        bool test(expression_token s) { return ctok != tokens.end() && *ctok == s; }

        bool accept(expression_token s) {
            if (!test(s)) return false;
            next();
            return true;
        }

        int expect(expression_token s) {
            if (accept(s))
                return 1;
            throw cado_expression_parser_details::parse_error();
            return 0;
        }

        long exponent() {
            /* a^b^c parenthesizes as a^(b^c)... */
            using namespace cado_expression_parser_details;
            if (accept(POWER)) {
                long sign = 1;
                if (accept(MINUS)) sign = -1;
                expect(POSITIVE_NUMBER);
                const long e = (*cnumber++).as_exponent();
                const long ne = exponent();
                /* now we want to return sign * e ^ ne */
                if (ne == 1) {
                    return sign*e;
                } else if (ne < 0) {
                    throw parse_error();
                } else if (ne == 0) {
                    return sign;
                } else {
                    cxx_mpz ze = e;
                    mpz_pow_ui(ze, ze, ne);
                    mpz_mul_si(ze, ze, sign);
                    if (!mpz_fits_slong_p(ze)) throw parse_error();
                    return mpz_get_si(ze);
                }
            } else {
                return 1;
            }
        }

        bool tokenize(std::istream& is, int accept_literals) {
            using namespace cado_expression_parser_details;
            tokens.clear();
            literals.clear();
            numbers.clear();
            for( ; !is.eof() ; ) {
                int c;
                number_literal z;
                for(;;is.get()) {
                    c = is.peek();
                    if (is.eof() || !isspace(c)) break;
                }
                /* c is the next non-whitespace character */
                if (is.eof()) { break;
                } else if (c == '+') { is.get(); tokens.push_back(PLUS);
                } else if (c == '-') { is.get(); tokens.push_back(MINUS);
                } else if (c == '*') { is.get(); tokens.push_back(TIMES);
                } else if (c == '(') { is.get(); tokens.push_back(LEFT_PAREN);
                } else if (c == ')') { is.get(); tokens.push_back(RIGHT_PAREN);
                } else if (c == '^') { is.get(); tokens.push_back(POWER);
                } else if (number_literal::recognize(z, is)) {
                    numbers.push_back(z);
                    tokens.push_back(POSITIVE_NUMBER);
                } else if (accept_literals && isalpha(c)) {
                    is.get();
                    std::string lit(1, (char) c);
                    /* we should be able to collect literals such as x0 */
                    for( ;; lit += (char) c, is.get()) {
                        c = is.peek();
                        if (is.eof() || !isalnum(c)) break;
                    }
                    literals.push_back(lit);
                    tokens.push_back(LITERAL);
                } else {
                    throw token_error();
                }
            };
            ctok = tokens.begin();
            clit = literals.begin();
            cnumber = numbers.begin();
            return true;
        }

    };
}

/* TODO: I wonder whether this could be expanded to do the parsing of the
 * fantastic galois actions that we have here and there... */

template<typename T>
struct cado_expression_parser : public T, public cado_expression_parser_details::cado_expression_parser_base {
    using typename T::type;
    using typename T::number_type;

    template<typename... Args>
    explicit cado_expression_parser(Args&& ...args)
    /* we need braced initialization, otherwise the compiler (at least
     * clang++15 on macs) will look for an explicitly defined ctor
     */
    : T { std::forward<Args>(args)... }
    {}

private:
    template<typename... Args>
    type parse_factor(Args... args) {
        using cado_expression_parser_details::parse_error;
        type p(args...);

        if (T::accept_literals && accept(LITERAL)) {
            T::set_literal_power(p, *clit++, exponent());
        } else if (accept(MINUS)) {
            p = parse_factor(args...);
            T::neg(p, p);
        } else if (accept(POSITIVE_NUMBER)) {
            T::set(p, *cnumber++);
        } else if (accept(LEFT_PAREN)) {
            p = parse_expression(args...);
            expect(RIGHT_PAREN);
        } else {
            throw parse_error();
        }
        long e = exponent();
        /* At this point it would possibly be doable to accept negative
         * (integer!) exponents. real exponents such as ^2^2^-2 would
         * require expanding to fuller symbolic parsing.
         */
        if (e < 0)
            throw parse_error();
        if (e != 1)
            T::pow_ui(p, p, e);
        return p;
    }

    template<typename... Args>
    type parse_term(Args... args) {
        type p = parse_factor(args...);
        for( ; accept(TIMES) ; )
            T::mul(p, p, parse_factor(args...));
        return p;
    }

    template<typename... Args>
    type parse_expression(Args... args) {
        type p(args...);
        if (test(PLUS)) {
            next();
            p = parse_term(args...);
        } else if (!test(MINUS)) {
            p = parse_term(args...);
        }

        for(;;) {
            if (test(PLUS)) {
                next();
                T::add(p, p, parse_term(args...));
            } else if (test(MINUS)) {
                next();
                T::sub(p, p, parse_term(args...));
            } else
                break;
        }
        return p;
    }
public:
    template<typename... Args>
    type parse(Args&&... args) {
        using cado_expression_parser_details::parse_error;
        {
            /* count literals, see if we have the right number (at most) */
            auto lcopy = literals;
            std::sort(lcopy.begin(), lcopy.end());
            std::string c;
            int nlit = 0;
            for(auto const & l : lcopy) {
                if (l != c) nlit++;
                c = l;
                if (nlit > T::accept_literals)
                    throw parse_error();
            }
        }
        type p = parse_expression(std::forward<Args>(args)...);
        if (ctok != tokens.end())
            throw parse_error();
        return p;
    }
    bool tokenize(std::istream& is) {
        return cado_expression_parser_details::cado_expression_parser_base::tokenize(is, T::accept_literals);
    }
};


#endif	/* CADO_EXPRESSION_PARSER_HPP_ */
