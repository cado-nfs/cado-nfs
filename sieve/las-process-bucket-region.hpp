#ifndef LAS_PROCESS_BUCKET_REGION_HPP_
#define LAS_PROCESS_BUCKET_REGION_HPP_

#include <stdint.h>
#include <memory>
#include <array>
#include "threadpool.hpp"
#include "las-threads-work-data.hpp"
#include "las-auxiliary-data.hpp"
#include "las-report-stats.hpp"

#if 0
struct process_bucket_region_parameters: public task_parameters {
    /* ws is used in the synchronous phase, so we don't need to retain
     * ownership via a shared pointer. For wc and aux, we do, because
     * they're used late */
    nfs_work & ws;
    std::shared_ptr<nfs_work_cofac> wc_p;
    std::shared_ptr<nfs_aux> aux_p;
    sieve_info & si;
    where_am_I w;

    /* For good-old process_bucket_region, which actually does a whole
     * congruence class of bucket regions, this is the smallest N that is
     * considered. For the one-bucket-region-at-a-time approach, this is
     * simply the bucket index itself.
     *
     * Note that in the former (good-old) case, when si.toplevel==1, we
     * have only one level of buckets, so we don't need this field */
    uint32_t first_region0_index=0;

    process_bucket_region_parameters(nfs_work & ws, std::shared_ptr<nfs_work_cofac> wc_p, std::shared_ptr<nfs_aux> aux_p, sieve_info & si, where_am_I const& w)
        : ws(ws), wc_p(wc_p), aux_p(aux_p), si(si), w(w)
    {}
};

/* This is in las.cpp */
extern task_result * process_bucket_regions_congruence_class(worker_thread * worker, task_parameters * _param, int);

#endif

extern void process_many_bucket_regions(nfs_work & ws, std::shared_ptr<nfs_work_cofac> wc_p, std::shared_ptr<nfs_aux> aux_p, thread_pool & pool, int first_region0_index, sieve_info & si, where_am_I const & w);

/* {{{ process_one_bucket_region */

struct process_bucket_region_spawn {
    nfs_work & ws;
    std::shared_ptr<nfs_work_cofac> wc_p;
    std::shared_ptr<nfs_aux> aux_p;
    sieve_info & si;
    where_am_I w_saved;
    int first_region0_index;

    process_bucket_region_spawn(
            nfs_work & ws,
            std::shared_ptr<nfs_work_cofac> wc_p,
            std::shared_ptr<nfs_aux> aux_p,
            sieve_info & si,
            where_am_I w,
            int first_region0_index
            ) :
        ws(ws), wc_p(wc_p), aux_p(aux_p), si(si), w_saved(w), first_region0_index(first_region0_index) {}

    void operator()(worker_thread * worker, int id);
};

struct process_bucket_region_run : public process_bucket_region_spawn {
    worker_thread * worker;
    nfs_aux::thread_data & taux;
    nfs_work::thread_data & tws;
    timetree_t & timer;
    int bucket_relative_index;
    timetree_t::accounting_activate dummy;
    las_report& rep;
    unsigned char * S[2];
    /* We will have this point to the thread's where_am_I data member.
     * (within nfs_aux::th). However it might be just as easy to let this
     * field be defined here, and drop the latter.
     */
    where_am_I & w;

    /* A note on SS versus S[side]
     *
     * SS is temp data. It's only used here, and it could well be defined
     * here only. We declare it at the thread_data level to avoid
     * constant malloc()/free().
     *
     * S[side] is where we compute the norm initialization. Some
     * tolerance is subtracted from these lognorms to account for
     * accepted cofactors.
     *
     * SS is the bucket region where we apply the buckets, and also later
     * where we do the small sieve.
     *
     * as long as SS[x] >= S[side][x], we are good.
     */

    unsigned char *SS;
    
    struct side_data {/*{{{*/
        bucket_array_complete purged;   /* for purge_buckets */
        bucket_primes_t primes;         /* for resieving */
        side_data() :
            purged(bucket_array_complete(BUCKET_REGION)),
            primes(bucket_primes_t(BUCKET_REGION))
        {}
    };/*}}}*/

    std::array<side_data, 2> sides;

    process_bucket_region_run(process_bucket_region_spawn const & p, worker_thread * worker, int id);

    /* will be passed as results of functions
    std::vector<uint32_t> survivors;
    std::vector<bucket_update_t<1, shorthint_t>::br_index_t> survivors2;
     * */

    /* most probably useless, I guess
    int N;
    int cpt;
    int copr;
    */

    void init_norms(int side);
    void apply_buckets(int side);
    void small_sieve(int side);
    void SminusS(int side);
    typedef std::vector<uint32_t> surv1_t;
    typedef std::vector<bucket_update_t<1, shorthint_t>::br_index_t> surv2_t;
    surv1_t search_survivors();
    surv2_t convert_survivors(surv1_t&& survivors);
    void purge_buckets(int side);
    void resieve(int side);
    void cofactoring_sync (surv2_t & survivors2);
    void operator()();
};


/*}}}*/

#endif	/* LAS_PROCESS_BUCKET_REGION_HPP_ */
