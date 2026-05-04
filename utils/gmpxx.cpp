#include "cado.h" // IWYU pragma: keep

#include <cstdlib>

#include <istream>
#include <ostream>
#include <string>
#include <memory>
#include <ios>

#include <gmp.h>
#include "gmpxx.hpp"

#include "utils_cxx.hpp"

using namespace std;

/* In case we don't have gmpxx, we have to provide the C++ I/O functions
 * by ourselves...
 *
 * We make no effort to do this very accurately, though.
 */

namespace {
inline int getbase(ostream const& o)
{
    switch(o.flags() & std::ios_base::basefield) {
        case std::ios::hex:
            return 16;
        case std::ios::dec:
            return 10;
        case std::ios::oct:
            return 8;
        default:
            return 10;
    }
}
} // namespace

ostream& operator<<(ostream& os, mpz_srcptr x)
{
    int b = getbase(os);
    const unique_ptr<char[], free_delete<char>> str(mpz_get_str(nullptr, b, x));
    if ((os.flags() & std::ios::showbase) && (b == 8 || b == 16)) {
        if (str[0] == '-')
            os << "-0" << (b == 8 ? 'b' : 'x') << str.get() + 1;
        else
            os << "0" << (b == 8 ? 'b' : 'x') << str.get();
    } else {
        os << str.get();
    }
    return os;
}

ostream& operator<<(ostream& os, mpq_srcptr x)
{
    const unique_ptr<char, free_delete<char>> str(mpq_get_str(nullptr, getbase(os), x));
    os << str.get();
    return os;
}

istream& operator>>(istream& is, mpz_ptr x)
{
    string s;
    is >> s;
    mpz_set_str(x, s.c_str(), 0);
    return is;
}


istream& operator>>(istream& is, mpq_ptr x)
{
    string s;
    is >> s;
    mpq_set_str(x, s.c_str(), 0);
    return is;
}

