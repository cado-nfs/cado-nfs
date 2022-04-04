#ifndef FACUL_STRATEGIES_HPP_
#define FACUL_STRATEGIES_HPP_

#include <cstdio>      // for FILE
#include <array>
#include <vector>
#include <map>

#include "facul_method.hpp"

struct cxx_mpz;

extern int nb_curves (const unsigned int lpb, const unsigned int mfb);

/* All prime factors in the input number must be > fb. A factor of the 
   input number is assumed to be prime if it is < fb^2.
   The input number is taken to be not smooth if it has a 
   prime factor > 2^lpb. */

/* used in testbench -- it should _probably_ be used in factor_one as
 * well.
 */
struct facul_strategy_oneside {
  unsigned long B;
  unsigned int lpb;  /* Large prime bound 2^lpb */
  double BB;         /* The factor base bound squared.
                      * We assume that primes <= fbb have already been
                      * removed, thus any factor <= BB is assumed prime
                      * without further test.
                      * Note that BB is _only_ used for that, and that
                      * setting BB=0 is an effective way to consider all
                      * numbers as potential composites.
                      */
  double BBB;        /* The factor base bound cubed. */
  unsigned int mfb;
  std::vector<facul_method> methods;  /* List of methods to try */

  /* B lpb mfb ncurves verbose */
  facul_strategy_oneside(
          unsigned long B, unsigned int lpb, unsigned int mfb,
          int ncurves,
          int verbose);
  facul_strategy_oneside(
          unsigned long B, unsigned int lpb, unsigned int mfb,
          std::vector<facul_method::parameters> const & mps,
          int verbose);
  facul_strategy_oneside() = default;

  /* This only returns the parameters! The methods themselves must be
   * instantiated afterwards */
  static std::vector<facul_method::parameters> default_strategy (int n);
};
int facul (std::vector<cxx_mpz> &, cxx_mpz const &, facul_strategy_oneside const &);


struct facul_strategies_base {
    std::vector<unsigned long> B; /* The factor base bounds */
    std::vector<unsigned int> lpb; /* Large prime bounds (in bits) */
    std::vector<double> BB; /* The factor base bounds squared.
				 We assume that primes <= fbb have
				 already been removed, thus any
				 factor <= assume_prime_thresh is
				 assumed prime without further
				 test. */
    std::vector<double> BBB; /* The factor base bounds cubed. */
    std::vector<unsigned int> mfb; /* The cofactor bounds (bits) */
    facul_strategies_base(
            std::vector<unsigned long> const & lim,
            std::vector<unsigned int> const & lpb,
            std::vector<unsigned int> const & mfb,
            bool perfectly_sieved);
};

struct facul_strategies : public facul_strategies_base {
    std::vector<facul_method_side> const & operator()(unsigned int r, unsigned int a) const;

private:
    /* Optimization for facul_make_strategies () */
    std::map<facul_method::parameters, facul_method> precomputed_methods;

    /* now in reality, this cache could even be shared across many siever
     * configs, right?
     */
    std::map<
        std::vector<facul_method::parameters_with_side>,
        std::vector<facul_method_side>
    > precomputed_strategies;

    /* two quick strategies as a default (when no strategy file is
     * provided) */
    std::vector<std::vector<facul_method_side>> uniform_strategy;

public:

    typedef std::map<
        std::array<unsigned int, 2>,
        std::vector<facul_method::parameters_with_side>> strategy_file;

    facul_strategies(
            std::vector<unsigned long> const & lim,
            std::vector<unsigned int> const & lpb,
            std::vector<unsigned int> const & mfb,
            std::vector<int> ncurves,
            bool,
            const int);

    facul_strategies(
            std::vector<unsigned long> const & lim,
            std::vector<unsigned int> const & lpb,
            std::vector<unsigned int> const & mfb,
            bool,
            FILE *,
            const int);

    facul_strategies(
            std::vector<unsigned long> const & lim,
            std::vector<unsigned int> const & lpb,
            std::vector<unsigned int> const & mfb,
            bool,
            strategy_file const &,
            const int);

    void print(FILE *) const;

    /* used (internally) in ctors and (externally) in tests */
    void precompute_method(facul_method::parameters const & mp, int verbose);

private:
    /* This is the backing store of operator() */
    std::vector<const std::vector<facul_method_side> *> direct_access;
    std::vector<facul_method_side> const * & direct_access_get(unsigned int r, unsigned int a);
    std::vector<facul_method_side> const * const & direct_access_get(unsigned int r, unsigned int a) const;
};


#endif	/* FACUL_STRATEGIES_HPP_ */
