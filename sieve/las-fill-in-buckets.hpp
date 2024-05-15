#ifndef LAS_FILL_IN_BUCKETS_HPP_
#define LAS_FILL_IN_BUCKETS_HPP_

#include <cstdint>
#include <vector>
#include <array>
#include <memory>
#include "las-config.h" // FB_MAX_PARTS
class nfs_aux; // IWYU pragma: keep
class nfs_work;
class nfs_work_cofac;
class thread_pool;
struct where_am_I;
class plattices_vector_t;
template <template <int> class F, int n0, int n1> struct multityped_array;

// This one is used for keeping information of middle primes.
template<int LEVEL>
struct precomp_plattice_t {
    static const int level = LEVEL;
    typedef precomp_plattice_t type;    /* for multityped_array */
    typedef std::vector<plattices_vector_t> vec_type;
    std::vector<vec_type> v;
    precomp_plattice_t(precomp_plattice_t<LEVEL> const&) = delete;
    precomp_plattice_t(int nsides) : v(nsides) {}
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
