#include "cado.h" // IWYU pragma: keep

#include <istream>
#include <string>
#include <cctype>

#include "cado_parsing_base.hpp"
#include "number_literal.hpp"
#include "macros.h"

static void swallow_not_terminal(int & c, std::string & dst, std::istream & is)
{
    dst += static_cast<char>(c);
    is.get();
    c = is.peek();
    if (is.eof())
        throw cado::token_error();
}
static bool swallow_hits_end(int & c, std::string & dst, std::istream & is)
{
    dst += static_cast<char>(c);
    is.get();
    c = is.peek();
    return is.eof();
}

static bool recognize_constant_string(int & c, cado::number_literal & N,
                                      std::istream & is, std::string const & S)
{
    ASSERT_ALWAYS(!S.empty());
    if (c != S[0])
        return false;
    for (auto s: S) {
        c = is.get();
        if (is.eof() || c != s)
            throw cado::token_error();
    }
    c = is.peek();
    N.full += S;
    return true;
}

/* recognize a number literal, return true if successful, and is is
 * stopped at the first character that is not part of the recognized
 * literal. is may or may not be in eof state on return, it does not
 * affect the boolean return value of this function.
 */
bool cado::number_literal::recognize(number_literal & N,
                                     std::istream & is) // {{{
{
    N = number_literal {};
    int c = is.peek();
    bool no_mantissa = true;

    if (c == '-')
        swallow_not_terminal(c, N.full, is);

    if (recognize_constant_string(c, N, is, "nan") ||
        recognize_constant_string(c, N, is, "NaN") ||
        recognize_constant_string(c, N, is, "inf")) {
        /* "nan i" or "inf i", both with one preceding space, are
         * valid literals */
        if (recognize_constant_string(c, N, is, " i"))
            N.is_imaginary = true;
        return true;
    }

    if (!std::isdigit(c)) {
        if (c == '.') {
            N.end_integral = N.full.size();
            N.has_point = true;
            swallow_not_terminal(c, N.full, is);
            // *ps = "0";
            // ps = &N.fractional;
        } else {
            return false;
        }
    } else if (c == '0') {
        if (swallow_hits_end(c, N.full, is))
            return true;
        if (std::tolower(c) == 'x') {
            N.is_hex = true;
            swallow_not_terminal(c, N.full, is);
        } else {
            no_mantissa = false;
        }
    }

    for (; N.isdigit(c) || (!N.has_point && c == '.');) {
        if (c == '.')
            N.has_point = true;
        no_mantissa = false;
        if (!N.has_point)
            N.end_integral = N.full.size() + 1;
        if (swallow_hits_end(c, N.full, is))
            return true;
    }
    if (no_mantissa)
        throw token_error();

    if (std::tolower(c) == (N.is_hex ? 'p' : 'e')) {
        N.has_exponent = true;
        N.begin_exponent = N.full.size() + 1;
        if (swallow_hits_end(c, N.full, is))
            return true;
        if (c == '+' || c == '-')
            if (swallow_hits_end(c, N.full, is))
                return true;
        for (; std::isdigit(c);)
            if (swallow_hits_end(c, N.full, is))
                return true;
    }

    /* a suffix i must come without a preceding space */
    if (c == 'i') {
        N.is_imaginary = true;
        swallow_hits_end(c, N.full, is);
    }
    return true;
} // }}}
