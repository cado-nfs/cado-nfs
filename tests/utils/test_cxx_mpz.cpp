#include "cado.h" // IWYU pragma: keep

#include <string>

#include <gmp.h>
#include "fmt/base.h"
#include "fmt/format.h"

#include "cxx_mpz.hpp"

/* we've had this cxx_mpz type around for years, and yet no proper
 * unit testing is in place. At the very least, the specifics of fmt
 * printing deserve their own test.
 */

static void check(long z, fmt::format_string<cxx_mpz> const & format, std::string const & fcopy, std::string const & exp)
{
    cxx_mpz zz = z;
    std::string gotz = fmt::format(fmt::runtime(format), z);
    std::string gotzz = fmt::format(fmt::runtime(format), zz);
    if (gotz == exp && gotzz == exp)
        return;
    throw cado::error("format error while printing {} with format {}:"
            " long -> {} ; cxx_mpz -> {} ; expected {}",
            z, fcopy, gotz, gotzz, exp);
}

#define checkcheck(a,b,c) check(a,b,b,c)

int main()
{
    {
        const long z = 1234;
        checkcheck(z, "{}", "1234");
        checkcheck(z, "{:x}", "4d2");
        checkcheck(z, "{:#x}", "0x4d2");
        checkcheck(z, "{:X}", "4D2");
        checkcheck(z, "{:#X}", "0X4D2");
        checkcheck(-z, "{:#X}", "-0X4D2");
        checkcheck(-z, "{:#o}", "-02322");
    }
}
