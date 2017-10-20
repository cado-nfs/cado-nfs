#include "cado.h"
#include <gmp.h>
#include <istream>
#include <ostream>
#include <string>
#include <stdlib.h>
#include "gmpxx.hpp"

using namespace std;

/* In case we don't have gmpxx, we have to provide the C++ I/O functions
 * by ourselves...
 *
 * We make no effort to do this very accurately, though.
 */

static inline int getbase(ostream const& o)
{
    std::ios_base::fmtflags ff;
    ff = o.flags();
    ff &= std::ios_base::basefield;
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

ostream& operator<<(ostream& os, mpz_srcptr x)
{
    char * str;
    os << (str = mpz_get_str(NULL, getbase(os), x));
    free(str);
    return os;
}

ostream& operator<<(ostream& os, mpq_srcptr x)
{
    char * str;
    os << (str = mpq_get_str(NULL, getbase(os), x));
    free(str);
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

