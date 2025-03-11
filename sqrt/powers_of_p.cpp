
/* This file uses the standard C++ functions to provide a look-up table
 * for primes of a power p. It must be given an opaque structure, which
 * in fact is of type power_lookup_table
 */
#include "cado.h" // IWYU pragma: keep

#include <cstdlib>
#include <cstdio>

#include <utility>                 // pair
#include <map>
#include <mutex>

#include <gmp.h>

#include "powers_of_p.h"
#include "utils_cxx.hpp"

using namespace std;

struct power_lookup_table {
    mutable std::mutex mx;
    unsigned long p;
    typedef map<int, int> m_t;
    mpz_ptr * z = nullptr;
    unsigned int alloc = 0;
    unsigned int nz = 0;
    m_t m;
    int extra_power_swapstore(mpz_ptr w) {
        if (nz == alloc) {
            alloc += alloc ? alloc / 4 : 16;
            checked_realloc(z, alloc);
        }
        z[nz] = (mpz_ptr) malloc(sizeof(mpz_t));
        mpz_init(z[nz]);
        mpz_swap(z[nz], w);
        nz++;
        return nz-1;
    }
    power_lookup_table(unsigned long p) : p(p) { }
    ~power_lookup_table() {
        for(unsigned int i = 0 ; i < nz ; i++) {
            mpz_clear(z[i]);
            free(z[i]);
        }
        free(z);
        z = nullptr;
    }
    mpz_srcptr operator()(int i);
    mpz_srcptr operator()(int i) const;
    private: mpz_srcptr inside(int i);
};

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
    mpz_t q;
    // XXX valgrind says that this sometimes leaks. I don't understand
    // why.  Perhaps it's obvious.
    mpz_init(q);
    if (i == 0) {
        mpz_set_ui(q, 1);
    } else if (i == 1) {
        mpz_init_set_ui(q, p);
    } else if (i & 1) {
        /* we do it in such a way that the newton lift encounter exactly
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
    int const r = extra_power_swapstore(q);
    m[i]=r;
    mpz_clear(q);  // has been swapped with the other one.
    mpz_srcptr res = z[r];
    return res;
}

// C entry functions

void * power_lookup_table_init(unsigned long p)
{
    return new power_lookup_table(p);
}

void power_lookup_table_clear(void * t)
{
    auto * pt = static_cast<power_lookup_table *>(t);
    delete pt;
}

mpz_srcptr power_lookup(void * t, int i)
{
    auto * pt = static_cast<power_lookup_table *>(t);
    return (*pt)(i);
}
mpz_srcptr power_lookup_const(const void * t, int i)
{
    auto const * pt = static_cast<const power_lookup_table *>(t);
    return (*pt)(i);
}
