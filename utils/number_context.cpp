#include "cado.h"       // IWYU pragma: keep

#include <cctype>

#include <string>
#include <complex>
#include <vector>
#include <istream>
#include <sstream>
#include <stdexcept>

#ifdef HAVE_MPC
#include "cxx_mpc.hpp"
#endif

#include "number_literal.hpp"
#include "number_context.hpp"
#include "cado_parsing_base.hpp"

using cado::parse_error;
using cado::token_error;
using cado::number_context;
using cado::number_literal;

/* parsing of complex numbers is a bit annoying. They're written as
 * algebraic expressions, really, so in a sense we probably want to
 * fire up our expression parser -- although at the same time it
 * feels a bit overkill.
 *
 * It's also remarkable that the output of formatting a complex
 * number (fmt produces: "(2+3i)" when there's a real part) is not a
 * number literal at all.
 *
 * As boring as it may seem, we must also include provision for the
 * case of "inf i" and "nan i", which are both printed by fmt.
 *
 * The approach we take is:
 *  - we separate the calls operator()(std::string) and
 *  operator()(number_literal).
 *  - we tokenize mildly for the former, allowing a restricted set of
 *  expressions that yield a complex number, the regexp being (MINUS ?
 *  POSITIVE_REAL_OR_TOTALLY_IMAGINARY_NUMBER | \( MINUS ? POSITIVE_REAL_NUMBER [PLUS | MINUS]
 *  POSITIVE_TOTALLLY_IMAGINARY_NUMBER \))
 *  - the call operator()(number_literal) only recognizes real
 *  numbers or totally imaginary numbers. The
 *  number_literal::recognize tokenizer 
 *  terminating i in the definition of POSITIVE_NUMBER, *BUT* 2+3i is
 *  still not a number literal. It's an expression, which can only be
 *  parsed if it's in the form (2+3i), using operator()(std::string)
 */

struct complex_tokenizer {
    enum token {
        LEFT_PAREN,
        RIGHT_PAREN,
        PLUS,
        MINUS,
        POSITIVE_NUMBER
    };
    std::vector<token> tokens;
    std::vector<number_literal> numbers;

    std::vector<token>::const_iterator ctok;
    std::vector<number_literal>::const_iterator cnumber;

    void reset() {
        ctok = tokens.begin();
        cnumber = numbers.begin();
    }
    void next() { if (ctok != tokens.end()) ctok++; }
    bool test(token s) { return ctok != tokens.end() && *ctok == s; }

    bool accept(token s) {
        if (!test(s)) return false;
        next();
        return true;
    }

    int expect(token s) {
        if (accept(s))
            return 1;
        throw parse_error();
        return 0;
    }

    bool tokenize(std::istream& is) {
        tokens.clear();
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
            } else if (c == '(') { is.get(); tokens.push_back(LEFT_PAREN);
            } else if (c == ')') { is.get(); tokens.push_back(RIGHT_PAREN);
            } else if (number_literal::recognize(z, is)) {
                numbers.push_back(z);
                tokens.push_back(POSITIVE_NUMBER);
            } else {
                throw token_error();
            }
        };
        reset();
        return true;
    }
    template<typename T>
        auto parse(number_context<T> const & ctx) -> decltype(ctx(number_literal())) {
            number_literal h;
            decltype(ctx(h)) z;
            if (accept(MINUS)) {
                expect(POSITIVE_NUMBER);
                z = -ctx(*cnumber++);
            } else if (accept(POSITIVE_NUMBER)) {
                /* This will also recognize an imaginary number */
                z = ctx(*cnumber++);
            } else if (accept(LEFT_PAREN)) {
                /* This *MUST* be the sum (or difference) of two
                 * numbers, the first one real, and the second one
                 * imaginary. Any other form raises a parse error */
                bool neg = false;
                if (accept(MINUS))
                    neg = true;
                expect(POSITIVE_NUMBER);
                h = *cnumber++;
                z = ctx(h);
                if (neg) z = -z;
                if (h.is_imaginary)
                    throw parse_error();
                if (accept(PLUS)) {
                    neg = false;
                } else if (accept(MINUS)) {
                    neg = true;
                } else {
                    throw parse_error();
                }
                expect(POSITIVE_NUMBER);
                h = *cnumber++;
                if (!h.is_imaginary)
                    throw parse_error();
                if (neg)
                    z += ctx(h);
                else
                    z -= ctx(h);
                expect(RIGHT_PAREN);
            }
            if (ctok != tokens.end())
                throw parse_error();
            return z;
        }
};


template<typename T>
std::complex<T> 
number_context<std::complex<T>>::operator()(std::string const & s) const
{
    complex_tokenizer P;
    std::istringstream is(s);
    P.tokenize(is);
    is >> std::ws;
    if (!is.eof())
        throw std::invalid_argument(s);
    return P.parse(*this);
}

#ifdef HAVE_MPC
template<>
cxx_mpc number_context<cxx_mpc>::operator()(std::string const & s) const
{
    complex_tokenizer P;
    std::istringstream is(s);
    P.tokenize(is);
    is >> std::ws;
    if (!is.eof())
        throw std::invalid_argument(s);
    return P.parse(*this);
}
#endif

namespace cado {

CADO_NUMBER_CONTEXT_DEFINE_EXPLICIT_INSTANTIATION(std::complex<float>)
CADO_NUMBER_CONTEXT_DEFINE_EXPLICIT_INSTANTIATION(std::complex<double>)
CADO_NUMBER_CONTEXT_DEFINE_EXPLICIT_INSTANTIATION(std::complex<long double>)

} /* namespace cado */
