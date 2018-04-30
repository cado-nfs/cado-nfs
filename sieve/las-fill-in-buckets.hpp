#ifndef LAS_FILL_IN_BUCKETS_HPP_
#define LAS_FILL_IN_BUCKETS_HPP_

#include <array>
#include "fb-types.h"
#include "las-plattice.hpp"
#include "tdict.hpp"
#include "threadpool.hpp"

#include "las-forwardtypes.hpp"

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
typedef std::vector<plattices_dense_vector_t> precomp_plattice_dense_t;

template <int LEVEL>
void
downsort_tree(
        nfs_work &ws,
        nfs_aux &aux,
        thread_pool &pool,
        uint32_t bucket_index,
        uint32_t first_region0_index,
        sieve_info& si,
        precomp_plattice_t & precomp_plattice,
        where_am_I & w);

void fill_in_buckets(
        nfs_work &ws,
        nfs_aux &aux,
        thread_pool &pool,
        sieve_info& si,
        int side,
        where_am_I & w);

void fill_in_buckets_prepare_precomp_plattice(
        nfs_aux &aux,
        thread_pool &pool,
        int side,
        int level,
        sieve_info const & si,
        precomp_plattice_t & precomp_plattice);

#endif
