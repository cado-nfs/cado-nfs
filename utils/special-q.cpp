#include "cado.h" // IWYU pragma: keep

#include <ostream>     // for operator<<, ostream, basic_ostream, basic_ostr...
#include <istream>     // for operator>>, istream, ...
#include <sstream>     // for istringstream
#include <string>
#include <ios>

#include <gmp.h>

#include "cxx_mpz.hpp"
#include "getprime.h"  // for getprime_mt, prime_info_clear, prime_info_init
#include "gmp_aux.h"
#include "special-q.hpp"
#include "macros.h"

void special_q::find_prime_factors()
{
    prime_factors.clear();

    /* It sometimes happens that we're only dealing with the basis, not
     * the special-q itself. This is a sign of the fact that we should
     * detach q and the basis, really.
     */
    if (mpz_cmp_ui(p, 0) == 0)
        return;

    /* We really do not want composites to be considered as primes... */
    if (mpz_probab_prime_p(p, 10)) {
        /* This is rubbish if p does not fit in 64 bits */
        if (mpz_fits_uint64_p(p))
            prime_factors.push_back(mpz_get_uint64(p));
        else
            prime_factors.push_back(0);

        return;
    }

    // Need to pre-compute the prime factors of the special-q in order to
    // skip them while sieving.
    prime_info pi;
    prime_info_init (pi);
    unsigned long f = 2;

    cxx_mpz B;
    mpz_sqrt(B, p);
    unsigned long const bound = mpz_get_ui(B);

    cxx_mpz pp;
    mpz_init_set(pp, p);
    while ((f <= bound) && (mpz_cmp_ui(pp, 1) > 0)) {
        if (mpz_divisible_ui_p(pp, f)) {
            mpz_divexact_ui(pp, pp, f);
            // Powers are not allowed in special-q
            ASSERT_ALWAYS(!mpz_divisible_ui_p(pp, f));
            prime_factors.push_back(f);
            if (mpz_probab_prime_p(pp, 10)) {
                ASSERT_ALWAYS(mpz_fits_ulong_p(pp));
                prime_factors.push_back(mpz_get_ui(pp));
                mpz_set_ui(pp, 1);
            }
        }
        f = getprime_mt(pi);
    }
    prime_info_clear (pi);

    ASSERT_ALWAYS(mpz_cmp_ui(pp, 1) == 0);
}

/* This format is also parsed by read_sq_comment in dupsup.cpp ! */
std::ostream& operator<<(std::ostream& os, special_q const & doing)
{
    os << "side-" << doing.side << " q=" << doing.p;
    if (!doing.is_prime()) {
        char c = '=';
        for(auto const & p : doing.prime_factors) {
            os << c << p;
            c = '*';
        }
    }
    os << "; rho=" << doing.r;
    return os;
}

std::istream& operator>>(std::istream& is, special_q & doing)
{
    doing = special_q();
    is >> std::ws >> expect("side-") >> doing.side;
    is >> std::ws >> expect("q=");
    std::string token;
    /* read until space-or-semicolon, because old-style q printing did
     * not have the semicolon...
     */
    for(char c ; is.get(c) && c != ';' && c != ' ' ; token.push_back(c));
    for(char & c : token)
        if (c == '=' || c == '*') c = ' ';
    std::istringstream iss(token);
    iss >> doing.p;
    for( ; iss.peek() != EOF ; ) {
        uint64_t g;
        iss >> g;
        if (iss.fail()) {
            is.setstate(std::ios::failbit);
            return is;
        }
        doing.prime_factors.push_back(g);
    }
    if (doing.prime_factors.empty())
        doing.find_prime_factors();
    is >> std::ws >> expect("r");
    if (!is) return is;
    /* support old-style q printing */
    int const c = is.peek();
    if (c == 'h')
        is >> expect("ho");
    is >> expect("=");
    if (!is) return is;
    std::getline(is, token, ';');
    iss.clear();
    iss.str(token);
    if (!(iss >> doing.r)) {
        is.setstate(std::ios::failbit);
    }
    return is;
}

