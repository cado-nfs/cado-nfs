#include "cado.h" // IWYU pragma: keep
#include <ostream>     // for operator<<, ostream, basic_ostream, basic_ostr...
#include <gmp.h>
#include "getprime.h"  // for getprime_mt, prime_info_clear, prime_info_init
#include "gmp_aux.h"
#include "las-todo-entry.hpp"
#include "macros.h"

void las_todo_entry::find_prime_factors()
{
    prime_factors.clear();

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
    unsigned long bound = mpz_get_ui(B);

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
std::ostream& operator<<(std::ostream& os, las_todo_entry const & doing)
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

