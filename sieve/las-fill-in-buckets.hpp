#ifndef LAS_FILL_IN_BUCKETS_HPP_
#define LAS_FILL_IN_BUCKETS_HPP_

#include <array>
#include <memory>
#include "las-config.h" // FB_MAX_PARTS
#include "fb-types.h"
#include "las-plattice.hpp"
#include "tdict.hpp"
#include "threadpool.hpp"
#include "las-threads-work-data.hpp"
#include "las-forwardtypes.hpp"
#include "multityped_array.hpp"

// This one is used for keeping information of middle primes.
template<int LEVEL>
struct precomp_plattice_t {
    static const int level = LEVEL;
    typedef precomp_plattice_t type;    /* for multityped_array */
    typedef std::vector<plattices_vector_t<LEVEL>> vec_type;
    std::array<vec_type, 2> v;
    precomp_plattice_t(precomp_plattice_t<LEVEL> const&) = delete;
    precomp_plattice_t() = default;
    void push(int side, vec_type&& x) {
        std::swap(v[side], x);
    }
    ~precomp_plattice_t() = default;
    /* This allows us to access the contents with range-based for loops */
    vec_type & operator()(int side) {
        return v[side];
    }
    vec_type const & operator()(int side) const {
        return v[side];
    }
};


template <int LEVEL>
void
downsort_tree(
        nfs_work &ws,
        std::shared_ptr<nfs_work_cofac> wc_p,
        std::shared_ptr<nfs_aux> aux_p,
        thread_pool &pool,
        uint32_t bucket_index,
        uint32_t first_region0_index,
        multityped_array<precomp_plattice_t, 1, FB_MAX_PARTS> & precomp_plattice,
        where_am_I & w);

void fill_in_buckets_toplevel(
        nfs_work &ws,
        nfs_aux &aux,
        thread_pool &pool,
        int side,
        where_am_I & w);

void fill_in_buckets_prepare_plattices(
        nfs_work & ws,
        thread_pool &pool,
        int side,
        multityped_array<precomp_plattice_t, 1, FB_MAX_PARTS> & precomp_plattice);
#endif
