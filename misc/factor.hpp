#ifndef CADO_FACTOR_HPP
#define CADO_FACTOR_HPP

#include <list>
#include <set>
#include <utility>
#include <vector>

#include "ecm.h"

#include "cxx_mpz.hpp"
#include "fmt/format.h"

/*
 * Class with methods to fully factor (or at least try to) the absolute value of
 * integer.
 */
class fully_factor {
    cxx_mpz N; /* input integer that we want to factor */
    unsigned long B; /* all primes <= B have been trial divided */
    cxx_mpz BB; /* square of the trial division bound */

    ecm_params params;

    std::set<std::pair<cxx_mpz, int>> primes;

    /* Invariant: prod(p^e for p,e in primes) * cofactor == N */
    cxx_mpz cofactor;
    /* List of **coprime** composites dividing cofactor.
     * Note that we can have prod(composites) != cofactor if one the composite
     * has multiplicity >= 2.
     */
    std::list<cxx_mpz> composites;

    /* contains B1 and #curves for targeting a factor of 5*i+30 bits for the
     * ith element.
     */
    static const std::vector<std::pair<unsigned int, unsigned int>> ECM_DATA;

    /* Set B1 and ncurves according to ECM_DATA to look for factors of n. */
    static void get_B1_ncurves(unsigned int & B1, unsigned int & ncurves,
                               cxx_mpz const & n);

    void print_progress(const char * header) const;

    bool isprime(cxx_mpz const & n) const;

    /* write n as r^e, with e >= 1 as large as possible.
     * r and n cannot be a ref to the same variable.
     */
    static void write_as_power(cxx_mpz & r, int & e, cxx_mpz const & n);

    /* Assumes
     *  - p is prime, no check is performed
     *  - p divides cofactor
     *  - p is not a reference to cofactor
     * /!\ Does not take care of removing p from composites if necessary.
     */
    void push_prime(cxx_mpz const & p);

    /* return iterator to the next composite (either just ++it or
     * composites.erase(it).
     */
    std::list<cxx_mpz>::iterator ecmlib_wrapper(std::list<cxx_mpz>::iterator it,
                                                unsigned int B1, int method,
                                                gmp_randstate_t randgen);

    void apply_method_to_all_composites(int method,
                                        std::vector<unsigned int> const & B1s,
                                        char const * method_str,
                                        gmp_randstate_t randgen);

public:
    int ecmlib_verbose = 0;
    int isprime_niter = 15; /* number of iterations for primality testing */

    explicit fully_factor(cxx_mpz n);

    /* return iterator to the next composite (either just ++it or
     * composites.erase(it).
     */
    std::list<cxx_mpz>::iterator split_composite(std::list<cxx_mpz>::iterator,
                                                 cxx_mpz const &);

    /* trial division for all primes p < B */
    void trial_division_up_to(unsigned long B);

    /* compute gcd between all hints and composites */
    void use_hints(std::vector<cxx_mpz> const & hints);

    /* Run PM1 with the given B1 value on all composites */
    void do_PM1(unsigned int B1, gmp_randstate_t randgen)
    {
        apply_method_to_all_composites(ECM_PM1, {B1}, "PM1", randgen);
    }

    /* Run PP1 with the given B1 value on all composites */
    void do_PP1(unsigned int B1, gmp_randstate_t randgen)
    {
        apply_method_to_all_composites(ECM_PP1, {B1}, "PP1", randgen);
    }

    /* Do niter pass over all composites; B1 and number of curves for ECM are
     * chosen based on the size of the composites and the data from ECM_DATA.
     */
    void do_ECM_based_on_composites_size(unsigned int niter,
                                         gmp_randstate_t randgen);

    /* default method to factor. It performs:
     *  - trial division
     *  - use hints
     *  - P-1
     *  - P+1
     *  - first round of ECM with one curve and increasing values of B1
     *  - multiple rounds of ECM with B1 and ncurves taken from ECM_DATA.
     */
    void factor_using_default_strategy(unsigned long trialdiv_bound,
                                       std::vector<cxx_mpz> const & hints,
                                       unsigned long PM1_B1,
                                       unsigned long PP1_B1,
                                       unsigned int nECM_rounds,
                                       gmp_randstate_t randgen);

    ~fully_factor()
    {
        ecm_clear(params);
    }

    std::set<std::pair<cxx_mpz, int>> const & prime_factors() const
    {
        return primes;
    }

    bool is_complete () const
    {
        return composites.empty(); /* should be equivalent to cofactor == 1U */
    }

    friend std::ostream & operator<<(std::ostream & o, fully_factor const & f);
};

namespace fmt {
    template <> struct formatter<fully_factor>: ostream_formatter {};
}

#endif /* CADO_FACTOR_HPP */
