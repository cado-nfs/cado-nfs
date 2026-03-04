#ifndef CADO_LAS_DIVIDE_PRIMES_HPP
#define CADO_LAS_DIVIDE_PRIMES_HPP

#include <algorithm>  // for max
#include <cstdint>    // for uint64_t, int64_t
#include <cstdio>     // for FILE
#include <vector>     // for vector
#include "fb.hpp"     // for fb_factorbase
class bucket_array_complete;
class bucket_primes_t;
struct cxx_mpz;
struct trialdiv_data;


typedef std::vector<uint64_t> factor_list_t;
extern void display_bucket_prime_stats();

extern int factor_list_fprint(FILE *f, factor_list_t const & fl);

extern void
divide_known_primes (std::vector<uint64_t> & fl, cxx_mpz & norm, const unsigned int N, unsigned int x,
           const bool handle_2, bucket_primes_t *primes,
           bucket_array_complete *purged,
	   trialdiv_data const & td,
           int64_t a, uint64_t b,
           fb_factorbase::slicing const & fbs,
           bool very_verbose);

#endif	/* CADO_LAS_DIVIDE_PRIMES_HPP */
