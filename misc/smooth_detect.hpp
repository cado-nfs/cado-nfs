#ifndef CADO_SMOOTH_DETECT_HPP
#define CADO_SMOOTH_DETECT_HPP

#include <climits>

#include <gmp.h>

#include "cxx_mpz.hpp"
#include "macros.h"

/*
 * Smoothness detector: take a list of pairs of integers and a smoothness
 * bound, test them for smoothness and returns one which is smooth.
 * It does not claim exhaustivity: if a candidate does not look promising
 * after a few ECMs, it is just skipped.
 *
 * The list of candidates is not given as a table, but a function must be
 * passed that returns the next candidate (computed or read from a
 * stream, whatever). A void * parameter is passed to this function, so
 * that it can have a context, and update it after each call.
 *
 * The next_cand() function must create candidates with something like:
 *   descent_init_candidate C(u0, v0, e);
 * where u0 and v0 are the two integers that must be checked for
 * simultaneous smoothness.
 *
 * An example of usage is given in descent_init_Fp.cpp .
 */

// Type and basic functions for candidates
struct descent_init_candidate
{
    cxx_mpz u0, v0;        // original integers to be tested for smoothness
    cxx_mpz u, v;          // unfactored parts of u, v, must be composite or 1
    unsigned int lpu = 0, lpv = 0; // bitsize of largest primes on both sides
    double effort = 0;     // sum of the B1 already tried on these numbers
    unsigned long e;       // this is a reconstruction of target^e

    descent_init_candidate() = default;
    descent_init_candidate(descent_init_candidate const&) = default;
    descent_init_candidate& operator=(descent_init_candidate const&) = default;
    descent_init_candidate(descent_init_candidate&&) = default;
    descent_init_candidate& operator=(descent_init_candidate&&) = default;

    void check_prime()
    {
        if (mpz_probab_prime_p(u, 1)) {
            lpu = MAX(lpu, mpz_sizeinbase(u, 2));
            u = 1;
        }
        if (mpz_probab_prime_p(v, 1)) {
            lpv = MAX(lpv, mpz_sizeinbase(v, 2));
            mpz_set_ui(v, 1);
        }
    }

    descent_init_candidate(cxx_mpz const& u0,
                           cxx_mpz const& v0,
                           cxx_mpz const& u,
                           cxx_mpz const& v,
                           unsigned int lpu,
                           unsigned int lpv,
                           unsigned long e)
      : u0(u0)
      , v0(v0)
      , u(u)
      , v(v)
      , lpu(lpu)
      , lpv(lpv)
      , effort(0.0)
      , e(e)
    {
        check_prime();
    }

    descent_init_candidate(cxx_mpz const& u0,
                           cxx_mpz const& v0,
                           unsigned long e)
      : descent_init_candidate(u0, v0, u0, v0, 0, 0, e)
    {
    }

    bool is_factored() const
    {
        return (mpz_cmp_ui(u, 1) == 0) && (mpz_cmp_ui(v, 1) == 0);
    }

    int cost() const
    {
        if (is_factored())
            return INT_MAX;
        else
            return MAX(mpz_sizeinbase(u, 2), mpz_sizeinbase(v, 2));
    }

    unsigned int maxlp() const { return MAX(lpu, lpv); }

    auto operator<=>(descent_init_candidate const& c) const
    {
        return cost() <=> c.cost();
    }
    bool operator==(descent_init_candidate const& c) const
    {
        return cost() == c.cost();
    }
    bool is_probably_not_smooth(unsigned int bound) const;
};

std::ostream&
operator<<(std::ostream&, descent_init_candidate const&);

// Type for tuning parameters for smooth_detect.
//   min_effort: the effort at the start (effort = sum of the B1 already tried)
//   max_effort: when this value is reached, don't increase effort anymore
//   max_pool_size: number of candidates to keep in mind at the same time
//   verbose: 0 = no verbose, 1 verbose
//   minB1: first B1 to try, in the ECM chain.
// Default values: {2000, +inf, 10, 1, 100.0}.
struct smooth_detect_params
{
    double min_effort;
    double max_effort;
    unsigned int max_pool_size;
    int verbose;
    double minB1;
};

// The main exported function. Smooth candidate is put in C.
// Last argument is for changing default strategy. NULL can be passed.
descent_init_candidate
smooth_detect(int (*next_cand)(descent_init_candidate&, const void*),
              const void* param_next_cand,
              unsigned long bound,
              smooth_detect_params const& param);

#endif /* CADO_SMOOTH_DETECT_HPP */
