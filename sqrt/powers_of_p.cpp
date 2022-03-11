
/* This file uses the standard C++ functions to provide a look-up table
 * for primes of a power p. It must be given an opaque structure, which
 * in fact is of type power_lookup_table
 */
#include "cado.h" // IWYU pragma: keep
#include <cstddef>      /* see https://gcc.gnu.org/gcc-4.9/porting_to.html */
#include <cstdlib>
#include <cstdio>
#include <utility>                 // pair
#include <map>
#include <gmp.h>
#include <pthread.h>
#include "powers_of_p.h"

using namespace std;

struct power_lookup_table {
    mutable pthread_mutex_t mx[1];
    unsigned long p;
    typedef map<int, int> m_t;
    mpz_ptr * z;
    unsigned int alloc;
    unsigned int nz;
    m_t m;
    int extra_power_swapstore(mpz_ptr w) {
        if (nz == alloc) {
            alloc += alloc ? alloc / 4 : 16;
            z = (mpz_ptr *) realloc(z, alloc * sizeof(mpz_t));
        }
        z[nz] = (mpz_ptr) malloc(sizeof(mpz_t));
        mpz_init(z[nz]);
        mpz_swap(z[nz], w);
        nz++;
        return nz-1;
    }
    power_lookup_table(unsigned long p) : p(p), z(NULL), alloc(0), nz(0) {
        pthread_mutex_init(mx, NULL);
    }
    ~power_lookup_table() {
        pthread_mutex_destroy(mx);
        for(unsigned int i = 0 ; i < nz ; i++) {
            mpz_clear(z[i]);
            free(z[i]);
        }
        free(z);
        z = NULL;
    }
    mpz_srcptr operator()(int i);
    mpz_srcptr operator()(int i) const;
    private: mpz_srcptr inside(int i);
};

mpz_srcptr power_lookup_table::operator()(int i)
{
    pthread_mutex_lock(mx);
    mpz_srcptr res = inside(i);
    pthread_mutex_unlock(mx);
    return res;
}

mpz_srcptr power_lookup_table::operator()(int i) const
{
    pthread_mutex_lock(mx);
    m_t::const_iterator px = m.find(i);
    if (px == m.end()) {
        pthread_mutex_unlock(mx);
        fprintf(stderr, "Fatal error: we would have expected p^%d to have been computed already\n", i);
        abort();
    }
    mpz_srcptr res = z[px->second];
    pthread_mutex_unlock(mx);
    return res;
}

mpz_srcptr power_lookup_table::inside(int i)
{
    m_t::const_iterator px = m.find(i);
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
    int r = extra_power_swapstore(q);
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
    power_lookup_table * pt = static_cast<power_lookup_table *>(t);
    return delete pt;
}

mpz_srcptr power_lookup(void * t, int i)
{
    power_lookup_table * pt = static_cast<power_lookup_table *>(t);
    return (*pt)(i);
}
mpz_srcptr power_lookup_const(const void * t, int i)
{
    const power_lookup_table * pt = static_cast<const power_lookup_table *>(t);
    return (*pt)(i);
}
