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
    cxx_mpz p = 0; /* this is the 'special-q', despite the 'p' name... */
    /* even for a rational side, the field below is used, since
     * it is needed for the initialization of the q-lattice. All callers
     * of las_todo_push must therefore make sure that a proper argument
     * is provided.
     */
    cxx_mpz r = 0;
    int side = 0;

    bool operator<(las_todo_entry const & o) const {
        int c;
        c = mpz_cmp(p, o.p); if (c) return c < 0;
        c = side - o.side;   if (c) return c < 0;
        c = mpz_cmp(r, o.r); if (c) return c < 0;
        return false;
    }

    las_todo_entry() = default;

    explicit operator bool() const { return p != 0; }
    bool operator!() const { return p == 0; }

    /********************************************************************/
    /* what comes below is only contextual information -- really, a
     * todo entry must be thought of as (p,r,side)
     */

    /* The array of prime_factors is valid unless p is a prime above 64
     * bits.
     */
    std::vector<uint64_t> prime_factors;
    bool is_prime() const { return prime_factors.size() == 1; }

    /* some fields which are specific to the descent */
    int depth = 0;
    int iteration = 0;      /* number of times we failed on this prime */

    private:
    las_todo_entry(const int side, const int depth)
        : side(side)
        , depth(depth)
    { }
    public:

    static las_todo_entry closing_brace(int depth)
    {
        return { -1, depth };
    }
    bool is_closing_brace() const { return side < 0; }

    las_todo_entry(cxx_mpz p, cxx_mpz r, const int side)
        : p(std::move(p))
        , r(std::move(r))
        , side(side)
    {
        find_prime_factors();
    }

    las_todo_entry(cxx_mpz p, cxx_mpz r, const int side,
            const int depth, const int iteration)
        : p(std::move(p))
        , r(std::move(r))
        , side(side)
        , depth(depth)
        , iteration(iteration)
    {
        find_prime_factors();
    }

    /* it's only used in finishbatch */
    las_todo_entry(cxx_mpz p, cxx_mpz r, const int side,
            std::vector<uint64_t> const & prime_factors)
        : p(std::move(p))
        , r(std::move(r))
        , side(side)
        , prime_factors(prime_factors)
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
