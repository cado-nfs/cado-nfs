#ifndef CADO_LAS_FILL_IN_BUCKETS_HPP
#define CADO_LAS_FILL_IN_BUCKETS_HPP

#include <cstdint>

#include <vector>
#include <memory>

#include "las-config.hpp"
#include "las-plattice.hpp"
#include "multityped_array.hpp"

class nfs_aux; // IWYU pragma: keep
class nfs_work;
class nfs_work_cofac;
class thread_pool;
struct where_am_I;
class plattices_vector_t;

// This one is used for keeping information of middle primes.
template<int LEVEL>
struct precomp_plattice_t :
    public std::vector<plattices_vector_t>
{
    static const int level = LEVEL;
    using type = precomp_plattice_t;    /* for multityped_array */
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
        std::vector<cado::multityped_array<precomp_plattice_t, 1, FB_MAX_PARTS - 1>> & precomp_plattice,
        where_am_I & w);

void fill_in_buckets_toplevel_multiplex(
        nfs_work &ws,
        nfs_aux &aux,
        thread_pool &pool,
        int side,
        where_am_I & w);

/* This prepares the p-lattices that will be used several times, which
 * means in particular: not at the toplevel!
 *
 * For this reason, it is enough to slice the precomp_plattice object,
 * and not look at the topmost level.
 */
void fill_in_buckets_prepare_plattices(
        nfs_work & ws,
        thread_pool &pool,
        int side,
        cado::multityped_array<precomp_plattice_t, 1, FB_MAX_PARTS - 1> & precomp_plattice);
#endif
