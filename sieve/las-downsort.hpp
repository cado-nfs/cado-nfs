#ifndef SIEVE_LAS_DOWNSORT_HPP_
#define SIEVE_LAS_DOWNSORT_HPP_

#include <memory>
#include <vector>

#include "las-config.hpp"
#include "las-fill-in-buckets.hpp"
#include "multityped_array.hpp"
#ifdef SIQS_SIEVE
#include "siqs-largesieve.hpp"
#endif
#include "sieve-methods.hpp"

/* Ideally this #ifdef should only be in las.ccp and not here */
#ifdef SIQS_SIEVE
    using ALGO = SIQS;
#else
    using ALGO = NFS;
#endif

class nfs_aux;
class nfs_work;
class nfs_work_cofac;
class thread_pool;
struct where_am_I;

extern void downsort_toplevel(
        nfs_work &ws,
        std::shared_ptr<nfs_work_cofac> const & wc_p,
        std::shared_ptr<nfs_aux> const & aux_p,
        ALGO::special_q_data const & Q,
        thread_pool &pool,
        std::vector<cado::multityped_array<precomp_plattice_t, 1, FB_MAX_PARTS - 1>> & precomp_plattices,
        where_am_I & w);

#endif	/* SIEVE_LAS_DOWNSORT_HPP_ */
