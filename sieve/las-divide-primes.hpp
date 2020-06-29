#ifndef LAS_DIVIDE_PRIMES_HPP_
#define LAS_DIVIDE_PRIMES_HPP_

#include <cstdio>
#include <cstdint>
#include <vector>
#include "cxx_mpz.hpp"
#include "bucket.hpp"
#include "trialdiv.hpp"
#include "fb.hpp"

typedef std::vector<uint64_t> factor_list_t;
extern void display_bucket_prime_stats();

extern int factor_list_fprint(FILE *f, factor_list_t const & fl);

extern void
divide_known_primes (std::vector<uint64_t> & fl, cxx_mpz & norm, const unsigned int N, unsigned int x,
           const bool handle_2, bucket_primes_t *primes,
           bucket_array_complete *purged,
	   trialdiv_data const & td,
           int64_t a, uint64_t b,
           fb_factorbase::slicing const & fbs);

#endif	/* LAS_DIVIDE_PRIMES_HPP_ */
