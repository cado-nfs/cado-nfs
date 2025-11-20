#ifndef CADO_UTILS_SPECIAL_Q_HPP
#define CADO_UTILS_SPECIAL_Q_HPP

#include <cstdint>

#include <iosfwd>
#include <ostream>
#include <string>
#include <vector>
#include <utility>

#include <gmp.h>
#include "fmt/ostream.h"
#include "fmt/format.h"
#include "fmt/base.h"

#include "macros.h"
#include "relation.hpp"

#include "cxx_mpz.hpp"


struct special_q;

std::ostream& operator<<(std::ostream&, special_q const &);

/* this parsing routine is used exclusively in fake_rels.cpp and
 * dupsup.cpp, and it goes to some level of complication in order to
 * parse what operator<< prints. At least we have symmetry, sure. But
 * overall, I don't think we should expose this kind of parsing globally
 */
std::istream& operator>>(std::istream&, special_q &);

namespace fmt {
    template <> struct formatter<special_q>: ostream_formatter {};
}


struct special_q {
    cxx_mpz p = 0; /* this is the 'special-q', despite the 'p' name... */
    /* even for a rational side, the field below is used, since
     * it is needed for the initialization of the q-lattice. All callers
     * of las_todo_push must therefore make sure that a proper argument
     * is provided.
     */
    cxx_mpz r = 0;
    int side = 0;

    bool operator<(special_q const & o) const {
        int c;
        c = mpz_cmp(p, o.p); if (c) return c < 0;
        c = side - o.side;   if (c) return c < 0;
        c = mpz_cmp(r, o.r); if (c) return c < 0;
        return false;
    }
    bool operator!=(special_q const & o) const {
        return (*this) < o || o < (*this);
    }
    bool operator==(special_q const & o) const {
        return !((*this) != o);
    }
    std::string shortname() const {
        return fmt::format("{}@{}", mpz_sizeinbase(p, 2), side);
    }
    std::string fullname() const {
        return fmt::format("{}", *this);
    }


    special_q() = default;

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

    public:

    special_q(cxx_mpz p, cxx_mpz r, const int side)
        : p(std::move(p))
        , r(std::move(r))
        , side(side)
    {
        find_prime_factors();
    }

    explicit special_q(int side, relation::pr const & pr)
        : p(pr.p)
        , r(pr.r)
        , side(side)
    {
        find_prime_factors();
    }

    /* it's only used in finishbatch */
    special_q(cxx_mpz p, cxx_mpz r, const int side,
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
    friend std::istream& operator>>(std::istream&, special_q &);
    void find_prime_factors();
};

#endif	/* CADO_UTILS_SPECIAL_Q_HPP */
