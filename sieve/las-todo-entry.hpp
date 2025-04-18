#ifndef CADO_LAS_TODO_ENTRY_HPP
#define CADO_LAS_TODO_ENTRY_HPP

#include <cstdint>      // for uint64_t
#include <iosfwd>       // for ostream
#include <vector>       // for vector
#include <utility>

#include <gmp.h>        // for mpz_cmp, mpz_cmp_ui, mpz_mul_ui, mpz_set_ui
#include "fmt/ostream.h"
#include "fmt/base.h"

#include "macros.h"     // for ASSERT_ALWAYS

#include "cxx_mpz.hpp"  // for cxx_mpz


struct las_todo_entry {
    cxx_mpz p; /* this is the 'special-q', despite the 'p' name... */
    /* even for a rational side, the field below is used, since
     * it is needed for the initialization of the q-lattice. All callers
     * of las_todo_push must therefore make sure that a proper argument
     * is provided.
     */
    cxx_mpz r;
    int side;
    /* The array of prime_factors is valid unless p is a prime above 64
     * bits.
     */
    std::vector<uint64_t> prime_factors;
    bool is_prime() const { return prime_factors.size() == 1; }

    /* some fields which are specific to the descent */
    int depth;
    int iteration;      /* number of times we failed on this prime */

    las_todo_entry() : side(0), depth(0), iteration(0) {}

    bool operator<(las_todo_entry const & o) const {
        int c = mpz_cmp(p, o.p);
        if (c) return c < 0;
        c = mpz_cmp(r, o.r);
        if (c) return c < 0;
        return false;
    }

    /* Empty p, r is used for "closing brace" */
    las_todo_entry(const int side, const int depth) : side(side), depth(depth), iteration(0) { }

    las_todo_entry(cxx_mpz p, cxx_mpz r, const int side, const int depth = 0, const int iteration = 0)
        : p(std::move(p))
        , r(std::move(r))
        , side(side)
        , depth(depth)
        , iteration(iteration)
    {
        find_prime_factors();
    }

    las_todo_entry(cxx_mpz p, cxx_mpz r, const int side, std::vector<uint64_t> & prime_factors, const int depth = 0, const int iteration = 0)
        : p(std::move(p))
        , r(std::move(r))
        , side(side)
        , prime_factors(prime_factors)
        , depth(depth)
        , iteration(iteration)
    {
        /* make sure that the factorization provided is correct */
        cxx_mpz t = 1;
        for(auto x : prime_factors)
            mpz_mul_ui(t, t, x);
        ASSERT_ALWAYS(mpz_cmp(this->p, t) == 0);
    }

    // Given a *prime* ell, check whether ell is coprime to current
    // entry.
    bool is_coprime_to(unsigned long ell) const {
        if (is_prime()) {
            return (mpz_cmp_ui(p, ell) != 0);
        } else {
            for (auto p : prime_factors) {
                if (p == ell)
                    return false;
            }
            return true;
        }
    }
private:
    friend std::istream& operator>>(std::istream&, las_todo_entry &);
    void find_prime_factors();
};

std::ostream& operator<<(std::ostream&, las_todo_entry const &);

/* this parsing routine is used exclusively in fake_rels.cpp, and it goes
 * to some level of complication in order to parse what operator<<
 * prints. At least we have symmetry, sure. But overall, I don't think we
 * should expose this kind of parsing globally
 */
std::istream& operator>>(std::istream&, las_todo_entry &);

namespace fmt {
    template <> struct formatter<las_todo_entry>: ostream_formatter {};
}

#endif	/* CADO_LAS_TODO_ENTRY_HPP */
