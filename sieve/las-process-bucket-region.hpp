#ifndef LAS_PROCESS_BUCKET_REGION_HPP_
#define LAS_PROCESS_BUCKET_REGION_HPP_

#include "threadpool.hpp"
#include "las-threads.hpp"
#include "las-auxiliary-data.hpp"
#include <stdint.h>

struct process_bucket_region_parameters: public task_parameters {
    /* aux is just stats and so on. This will become a shared pointer
     * someday */
    nfs_work & ws;
    nfs_aux & aux;
    sieve_info & si;
    where_am_I w;

    /* I _think_ that first_region0_index is really the "smallest N" (N
     * as in trace_Nx, that is, the bucket number) which is being
     * processed here. Of course when si.toplevel==1, we have only one
     * level of buckets, so we don't need that */
    uint32_t first_region0_index=0;
    process_bucket_region_parameters(nfs_work & ws, nfs_aux & aux, sieve_info & si, where_am_I const& w)
        : ws(ws), aux(aux), si(si), w(w)
    {}
};

/* This is in las.cpp */
extern task_result * process_bucket_region(const worker_thread * worker, task_parameters * _param);

#endif	/* LAS_PROCESS_BUCKET_REGION_HPP_ */
