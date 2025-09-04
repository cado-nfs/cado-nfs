#include "cado.h" // IWYU pragma: keep

#include <cstddef>

#include <string>
#include <vector>
#include <stdexcept>
#include <sstream>

#include <gmp.h>
#include "fmt/base.h"

#include "macros.h"
#include "cxx_mpz.hpp"
#include "number_literal.hpp"
#include "number_context.hpp"
#include "cado_parsing_base.hpp"

using cado::token_error;
using cado::parse_error;

/* this only tests our number (floating point and integer) literal parser
 *
 * testing of the actual parsers is done in passing, since we use it in
 * many tests that touch on polynomials of various kinds.
 */
int main()
{
    // NOLINTBEGIN(modernize-use-designated-initializers)
    struct test_case {
        std::string input;
        size_t tail_size;
        std::string parse_flags;
        double d;
        cxx_mpz z;
    };

    const std::vector<test_case> test_cases
    {
        /* these are all valid floating-point literals, but they do not
         * parse completely as integers. So we should get /Z as a result
         */
        { "0x000p+1",     0, "/Z", 0,                   0 },
        { "0x000p-1",     0, "/Z", 0,                   0 },
        { "0x0001p-1",    0, "/Z", 0.5,                 0 },
        { "0x0001p-2",    0, "/Z", 0.25,                0 },
        { "0x1p-2",       0, "/Z", 0.25,                0 },
        { "0x19p1",       0, "/Z", 50,                  0 },
        { "0x19p-1",      0, "/Z", 12.5,                0 },
        { "0x1.9p-1",     0, "/Z", 0.78125,             0 },
        { "3e9",          0, "/Z", 3000000000,          0 },
        { ".2",           0, "/Z", 0.20000000000000001, 0 },
        { "0x.2p1",       0, "/Z", 0.25,                0 },
        { "0x1.921fap+1", 0, "/Z", 3.1415901184082031,  0 },
        { "0x1.921fAp+1", 0, "/Z", 3.1415901184082031,  0 },
        { "0001.2",       0, "/Z", 1.2,                 0 },
        { "0001.2e0",     0, "/Z", 1.2,                 0 },
        // currently our code does not understand this as an integer.
        { "0x1p256",      0, "/Z", 1.157920892373162e+77,0 },
        { ".0",           0, "/Z", 0,                   0 },
        { "0.",           0, "/Z", 0,                   0 },
        { "0x1p2",        0, "/Z", 4,                   0 },
        { "0x2p0",        0, "/Z", 2,                   0 },
        { "0x2p1",        0, "/Z", 4,                   0 },
        { "0x.1p2",       0, "/Z", 0.25,                0 },
        { "0e+1",         0, "/Z",  0,                  0 },

        /* same for these, except that they do have an unrecognized
         * suffix that we choose not to parse */
        { "0x1.ffp-1f",   1, "/Z", .998046875,          0 },
        { "0u.2p1",       5, "/",  0,                   0 },
        { "0p.2p1",       5, "/",  0,                   0 },
        // in the cases below it's because we have a p exponent suffix
        // without a 0x marker
        { ".2p1",         2, "/Z", 0.2,                 0 },
        { "0001.2p",      1, "/Z", 1.2,                 0 },
        { "0001.2p0",     2, "/Z", 1.2,                 0 },

        /* This one parsers incompletely and thus has a tail that is not
         * recognized. Beyond that, 1 is of course a recognized integer
         * or floating-point number
         */
        { "1x1p-2",     5, "/",  1,                     1 },

        /* Here, the hexadecimal marker is too far down, so only the head
         * zeros are recognized
         */
        { "000x1.2p",   5, "/",  0,                     0 },

        /* Here the tokenization is successful but the parse as a
         * floating point number does not work, because p alone is not
         * good (note that this is among the corner cases that gdb
         * recognizes, though).
         * 
         * XXX these could be token_error, really
         */
        { "0x1.2p",     0, "/DZ", 0,                    0 },
        { "0x1.2pe",    1, "/DZ", 0,                    0 },
        { "0x2p",       0, "/DZ", 0,                    0 },
        { "0x2.p+",     0, "/DZ", 0,                    0 },
        { "0x2.p-",     0, "/DZ", 0,                    0 },
        { "0e",         0, "/DZ", 0,                    0 },
        { "0x.",        0, "/DZ", 0,                    0 },

        /* These are perfectly valid integers! */
        { "0xe2",       0, "/",   226,                  226 },
        { "0x1e2",      0, "/",   482,                  482 },

        /* These are all token errors because tokenization notices that
         * there's no mantissa, which is obviously bad.
         */
        { ".",          0, "t",   0,                    0 },
        { "..",         0, "t",   0,                    0 },
        { "0x",         0, "t",   0,                    0 },
        { "0xp2",       0, "t",   0,                    0 },
    };
    // NOLINTEND(modernize-use-designated-initializers)

    for(auto s : test_cases) {
        std::istringstream is(s.input);
        cado::number_literal N;
        std::string exc;
        double d = 0;
        cxx_mpz z = 0;
        std::string tail;

        using cado::number_literal;
        using cado::number_context;

        try {
            const bool b = number_literal::recognize(N, is);
            is >> tail;

            /* we have a / if the recognition of a number literal was
             * successful. Note that this only means tokenization,
             * really. Depending on the context in which this is going to
             * be interpreted, we can still fail to parse (e.g. if a
             * floating-point literal is parsed as an integer)
             */
            exc += '/';

            ASSERT_ALWAYS(b == (s.input.size() > tail.size()));
            try {
                d = number_context<double>()(N);
            } catch (parse_error const & e) {
                exc += 'D';
            }
            try {
                z = number_context<cxx_mpz>()(N);
            } catch (parse_error const & e) {
                exc += 'Z';
            }
        } catch (token_error const & e) {
            exc += 't';
        } catch (parse_error const & e) {
            exc += 'p';
        } catch (std::runtime_error const & e) {
            exc += 'R';
        }
        const bool good =
            exc == s.parse_flags
            && tail.size() == s.tail_size
            && d == s.d
            && z == s.z;
        fmt::print("{}\t{}\t[{}]\t{}\t{}\t{}\n", s.input, exc, tail, d, z, good);
        ASSERT_ALWAYS(good);
    }
}
