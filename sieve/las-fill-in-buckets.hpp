#ifndef LAS_FILL_IN_BUCKETS_HPP_
#define LAS_FILL_IN_BUCKETS_HPP_

#include <array>
#include "las-types.hpp"
#include "fb-types.h"
#include "las-plattice.hpp"
#include "las-threads.hpp"
#include "tdict.hpp"

// This one is used for keeping information of middle primes.
struct precomp_plattice_t {
    typedef std::vector<plattices_vector_t> vec_type;
    std::array<std::array<vec_type, FB_MAX_PARTS>, 2> v;
    precomp_plattice_t(precomp_plattice_t const&) = delete;
    precomp_plattice_t() = default;
    void push(int side, int level, plattices_vector_t&& x) {
        v[side][level].push_back(std::move(x));
    }
    ~precomp_plattice_t() = default;
    std::vector<plattices_vector_t> & operator()(int side, int level) {
        return v[side][level];
    }
    std::vector<plattices_vector_t> const & operator()(int side, int level) const {
        return v[side][level];
    }
};


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
        precomp_plattice_t & precomp_plattice);

void fill_in_buckets(timetree_t&, thread_pool &, thread_workspaces &, sieve_info &, int side);

void fill_in_buckets_prepare_precomp_plattice(
        thread_pool &pool,
        int side,
        int level,
        sieve_info const & si MAYBE_UNUSED,
        precomp_plattice_t & precomp_plattice);

#endif
