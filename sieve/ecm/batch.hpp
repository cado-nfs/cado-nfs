#ifndef SIEVE_ECM_BATCH_HPP
#define SIEVE_ECM_BATCH_HPP

#include <cstdint>
#include <cstdio>

#include <list>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <gmp.h>

#include "macros.h"
#include "cado_poly.h"
#include "mpz_poly.h"
#include "cxx_mpz.hpp"
#include "special-q.hpp"

struct relation; // IWYU pragma: keep


struct cofac_candidate {
    int64_t a = 0;
    uint64_t b = 0;
    std::vector<cxx_mpz> cofactor;
    cofac_candidate() = default;
    cofac_candidate(
            int64_t a, uint64_t b,
            std::vector<cxx_mpz> & cofactor)
        : a(a)
        , b(b)
        , cofactor(std::move(cofactor))
    { }
};

/*
 * These functions add to to the double& extra_time argument the
 * cpu time (RUSAGE_THREAD, seconds_thread()) spent in openmp helper
 * threads, NOT counting the time spent in the main thread.
 */
size_t find_smooth (
        std::list<cofac_candidate> & l,
        std::vector<cxx_mpz> const & batchP,
        std::vector<unsigned int> const & batchlpb,
        std::vector<unsigned int> const & lpb,
        std::vector<unsigned int> const & batchmfb,
        FILE *out,
        int nthreads MAYBE_UNUSED, double &);

size_t
find_smooth (std::list<std::pair<special_q, std::list<cofac_candidate>>> & l,
        std::vector<cxx_mpz> const & batchP,
        std::vector<unsigned int> const & batchlpb,
        std::vector<unsigned int> const & lpb,
        std::vector<unsigned int> const & batchmfb,
        FILE *out,
        int nthreads MAYBE_UNUSED, double & extra_time);

std::list<relation> factor (
        std::list<cofac_candidate> const &,
        cxx_cado_poly const&,
        special_q const &,
        std::vector<unsigned int> const & batchlpb,
        std::vector<unsigned int> const & lpb,
        int max_ncurves,
        FILE* output,
        int loose,
        double& extra_time,
        int);
void create_batch_file (std::string const &, cxx_mpz &, unsigned long, unsigned long,
                        cxx_mpz_poly const &, FILE*, int, double &);

#endif /* SIEVE_ECM_BATCH_HPP */
