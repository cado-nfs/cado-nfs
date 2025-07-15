#ifndef CADO_LAS_PROCESS_BUCKET_REGION_HPP
#define CADO_LAS_PROCESS_BUCKET_REGION_HPP

#include <memory>

#include "las-where-am-i-proxy.hpp"
#include "las-auxiliary-data.hpp"

class nfs_work;
class nfs_work_cofac;
class thread_pool;
class worker_thread;


/* process_many_bucket_regions is found in las.cpp, currently */

/* {{{ process_one_bucket_region */

struct process_bucket_region_spawn {
    nfs_work & ws;
    std::shared_ptr<nfs_work_cofac> wc_p;
    std::shared_ptr<nfs_aux> aux_p;

    where_am_I w_saved;
    
    /* These two indices are set from within process_many_bucket_regions,
     * prior to spawning all threads.
     *
     * first_region0_index is the bucket index of the first region for
     * which we filled the buckets.
     *
     * already done is, relative to first_region0_index, the index of the
     * first region for which the small sieve start position are
     * available in ssdpos_many.
     *
     * The i-th process_bucket_region task thus handles the bucket region
     * of index first_region0_index + already_done + i
     *
     * the two fields below must be set BY HAND before operator() is
     * called. See the .cpp file.
     */
    int first_region0_index = 0;
    int already_done = 0;

    process_bucket_region_spawn(
            nfs_work & ws,
            std::shared_ptr<nfs_work_cofac> wc_p,
            std::shared_ptr<nfs_aux> aux_p,
            where_am_I const & w)
    : ws(ws)
    , wc_p(std::move(wc_p))
    , aux_p(std::move(aux_p))
    , w_saved(w)
    {}

    void operator()(worker_thread * worker, int id);
};

/*}}}*/

extern void process_many_bucket_regions(nfs_work & ws, std::shared_ptr<nfs_work_cofac> wc_p, std::shared_ptr<nfs_aux> aux_p, thread_pool & pool, int first_region0_index, where_am_I & w);

#endif	/* CADO_LAS_PROCESS_BUCKET_REGION_HPP */
