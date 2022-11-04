#ifndef CADO_EXPRESSION_PARSER_HPP_
#define CADO_EXPRESSION_PARSER_HPP_

#include <exception>
#include <vector>
#include <algorithm>
#include "cxx_mpz.hpp"

/* This structure has only two public functions: tokenize and parse */

/* A typical use can be found in utils/mpz_poly.cpp */

/* TODO: I wonder whether this could be expanded to do the parsing of the
 * fantastic galois actions that we have here and there... */

template<typename T>
struct cado_expression_parser : public T {
    using typename T::type;

    struct parse_error: public std::exception {
        const char * what() const noexcept override { return "parse error"; }
    };
    struct token_error: public std::exception {
        const char * what() const noexcept override { return "token error"; }
    };

private:

    enum expression_token {
        LEFT_PAREN,
        RIGHT_PAREN,
        PLUS,
        MINUS,
        TIMES,
        POWER,
        POSITIVE_INTEGER,
        LITERAL
    };
    std::vector<expression_token> tokens;
    std::vector<char> literals;
    std::vector<cxx_mpz> integers;

    typename std::vector<expression_token>::const_iterator ctok;
    std::vector<char>::const_iterator clit;
    std::vector<cxx_mpz>::const_iterator cint;

    inline void next() { if (ctok != tokens.end()) ctok++; }
    inline bool test(expression_token s) { return ctok != tokens.end() && *ctok == s; }

    bool accept(expression_token s) {
        if (!test(s)) return false;
        next();
        return true;
    }

    int expect(expression_token s) {
        if (accept(s))
            return 1;
        throw parse_error();
        return 0;
    }

    unsigned long exponent() {
        if (accept(POWER)) {
            expect(POSITIVE_INTEGER);
            cxx_mpz e = *cint++;
            while (accept(POWER)) {
                expect(POSITIVE_INTEGER);
                if (!mpz_fits_ulong_p(*cint)) throw parse_error();
                unsigned long ei = mpz_get_ui(*cint++);
                mpz_pow_ui(e, e, ei);
            }
            if (!mpz_fits_ulong_p(e)) throw parse_error();
            return mpz_get_ui(e);
        } else {
            return 1;
        }
    }

    type parse_factor() {
        type p;
        if (T::accept_literals && accept(LITERAL)) {
            T::set_literal_power(p, *clit++, exponent());
        } else if (accept(POSITIVE_INTEGER)) {
            T::set_mpz(p, *cint++);
        } else if (accept(LEFT_PAREN)) {
            p = parse_expression();
            expect(RIGHT_PAREN);
        } else {
            throw parse_error();
        }
        unsigned long e = exponent();
        if (e != 1)
            T::pow_ui(p, p, e);
        return p;
    }

    type parse_term() {
        type p = parse_factor();
        for( ; accept(TIMES) ; )
            T::mul(p, p, parse_factor());
        return p;
    }

    type parse_expression() {
        type p;
        if (test(PLUS)) {
            next();
            p = parse_term();
        } else if (!test(MINUS)) {
            p = parse_term();
        }

        for(;;) {
            if (test(PLUS)) {
                next();
                T::add(p, p, parse_term());
            } else if (test(MINUS)) {
                next();
                T::sub(p, p, parse_term());
            } else
                break;
        }
        return p;
    }
public:
    type parse() {
        /* count literals, see if we have the right number (at most) */
        std::sort(literals.begin(), literals.end());
        char c = '\0';
        int nlit = 0;
        for(auto const & l : literals) {
            if (l != c) nlit++;
            c = l;
            if (nlit > T::accept_literals)
                throw parse_error();
        }
        type p = parse_expression();
        if (ctok != tokens.end())
            throw parse_error();
        return p;
    }

    bool tokenize(std::istream& is) {
        tokens.clear();
        literals.clear();
        integers.clear();
        for( ; !is.eof() ; ) {
            int c;
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
            } else if (isdigit(c)) {
                /* gmp's mpz parser really wants only an mpz, nothing
                 * else. We have to collect digits first.
                 */
                std::string s;
                for(;!is.eof() && isdigit(c);is.get(), c=is.peek()) {
                    s += c;
                }
                cxx_mpz z;
                mpz_set_str(z, s.c_str(), 0);

                integers.push_back(z);
                tokens.push_back(POSITIVE_INTEGER);
            } else if (T::accept_literals && isalpha(c)) {
                is.get();
                literals.push_back(c);
                tokens.push_back(LITERAL);
            } else {
                throw token_error();
            }
        };
        ctok = tokens.begin();
        clit = literals.begin();
        cint = integers.begin();
        return true;
    }
};


#endif	/* CADO_EXPRESSION_PARSER_HPP_ */
