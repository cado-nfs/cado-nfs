#ifndef CADO_LAS_FILL_IN_BUCKETS_HPP
#define CADO_LAS_FILL_IN_BUCKETS_HPP

#include <cstdint>

#include <vector>
#include <memory>

#include "las-config.hpp"
#include "multityped_array.hpp"
#ifdef SIQS_SIEVE
#include "siqs-largesieve.hpp"
#else
#include "las-plattice.hpp"
#endif
#include "sieve-methods.hpp"

/* Ideally this #ifdef should only be in las.ccp and not here */
#ifdef SIQS_SIEVE
    using ALGO = SIQS;
#else
    using ALGO = NFS;
#endif

class nfs_aux; // IWYU pragma: keep
class nfs_work;
class nfs_work_cofac;
class thread_pool;
struct where_am_I;

class plattices_vector_t : public std::vector<ALGO::largesieve>
{
    /* The index here is the global index, across all fb parts */
    slice_index_t index = -1;
    double weight = 0;

  public:
    plattices_vector_t() = default;
    plattices_vector_t(slice_index_t index, double weight)
        : index(index)
        , weight(weight)
    {
    }
    /* returns a global index */
    slice_index_t get_index() const { return index; };
    slice_index_t get_weight() const { return weight; };
};


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
        ALGO::special_q_data const & Q,
        thread_pool &pool,
        uint32_t bucket_index,
        uint32_t first_region0_index,
        std::vector<cado::multityped_array<precomp_plattice_t, 1, FB_MAX_PARTS - 1>> & precomp_plattice,
        where_am_I & w);

void fill_in_buckets_toplevel_multiplex(
        nfs_work &ws,
        nfs_aux &aux,
        ALGO::special_q_data const & Q,
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
        ALGO::special_q_data const & Q,
        thread_pool &pool,
        int side,
        cado::multityped_array<precomp_plattice_t, 1, FB_MAX_PARTS - 1> & precomp_plattice);
#endif
