#ifndef LAS_FILL_IN_BUCKETS_HPP_
#define LAS_FILL_IN_BUCKETS_HPP_

#include "las-types.hpp"
#include "fb-types.h"
#include "las-plattice.hpp"
#include "las-threads.hpp"
#include "tdict.hpp"

// This one is used for keeping information of middle primes.
typedef std::vector<plattices_vector_t *> precomp_plattice_t [2][FB_MAX_PARTS];

// This one is for remembering the FK basis in sublat mode, between two
// different congruences of (i,j) mod m.
// For simplicity, we remember them only for the toplevel.
typedef std::vector<plattices_dense_vector_t *> precomp_plattice_dense_t;

template <int LEVEL>
void
downsort_tree(timetree_t&,
        uint32_t bucket_index,
        uint32_t first_region0_index,
        thread_workspaces &ws,
        thread_pool &pool,
        sieve_info& si,
        precomp_plattice_t precomp_plattice);
void fill_in_buckets(timetree_t&, thread_pool &, thread_workspaces &, sieve_info &, int side);

#endif
