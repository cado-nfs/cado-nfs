#include "cado.h" // IWYU pragma: keep

#include <cstdlib>
#include <cstdio>

#include <utility>
#include <map>
#include <mutex>

#include <gmp.h>

#include "powers_of_p.hpp"
#include "cxx_mpz.hpp"

mpz_srcptr power_lookup_table::operator()(int i)
{
    const std::lock_guard<std::mutex> dummy(mx);
    mpz_srcptr res = inside(i);
    return res;
}

mpz_srcptr power_lookup_table::operator()(int i) const
{
    const std::lock_guard<std::mutex> dummy(mx);
    auto const px = m.find(i);
    if (px == m.end()) {
        fprintf(stderr, "Fatal error: we would have expected p^%d to have been computed already\n", i);
        abort();
    }
    mpz_srcptr res = z[px->second];
    return res;
}

mpz_srcptr power_lookup_table::inside(int i)
{
    auto const px = m.find(i);
    if (px != m.end()) {
        mpz_srcptr res = z[px->second];
        return res;
    }
    cxx_mpz q;
    // XXX valgrind says that this sometimes leaks. I don't understand
    // why.  Perhaps it's obvious.
    if (i == 0) {
        mpz_set_ui(q, 1);
    } else if (i == 1) {
        mpz_init_set_ui(q, p);
    } else if (i & 1) {
        /* we do it in such a way that the newton lift encounters exactly
         * this sequence of primes -- and we apologize for the division
         * by p, which could quite probably be saved (store the whole
         * chain in a ladder manner) */
        mpz_srcptr ph = inside(i - i/2);
        mpz_divexact_ui(q,ph,p);
        mpz_mul(q,q,ph);
    } else {
        mpz_srcptr pl = inside(i/2);
        mpz_mul(q,pl,pl);
    }
    z.push_back(std::move(q));
    m[i] = z.size()-1;
    return z.back();
}

