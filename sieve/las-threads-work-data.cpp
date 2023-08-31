#include "cado.h" // IWYU pragma: keep

#include <algorithm>
#include <string.h>        // for size_t, NULL, memset
#include <ostream>         // for operator<<, basic_ostream, ostringstream
#include <string>          // for char_traits, basic_string, string
#include <type_traits>     // for is_same

#include "las-threads-work-data.hpp"

#include "macros.h"        // for ASSERT_ALWAYS, iceildiv

#include "ecm/batch.hpp"       // for cofac_list
#include "bucket.hpp"      // for bucket_array_t, emptyhint_t (ptr only)
#include "las-info.hpp"    // for las_info
#include "las-memory.hpp"  // for las_memory_accessor
#include "threadpool.hpp"  // for thread_pool
#include "verbose.h"
class nfs_aux; // IWYU pragma: keep

nfs_work::thread_data::thread_data(thread_data && o) : ws(o.ws)
{
    for(int side = 0 ; side < 2 ; side++) {
        sides[side].bucket_region = o.sides[side].bucket_region;
        o.sides[side].bucket_region = NULL;
    }
    SS = o.SS;
    o.SS = NULL;
}

#if 0
/* See comment in las-threads-work-data.hpp */
nfs_work::thread_data::thread_data(thread_data const & o)
    : ws(o.ws)
{
    for(int side = 0 ; side < 2 ; side++) {
        ASSERT_ALWAYS(o.sides[side].bucket_region == NULL);
        sides[side].bucket_region = NULL;
    }
    ASSERT_ALWAYS(o.SS == NULL);
    SS = NULL;
#if 0
    /* We do not need MEMSET_MIN here. However we're making life easier
     * for the memory allocator if we allocate and free always the same
     * size.
     */
    for(int side = 0 ; side < 2 ; side++) {
        if (!sides[side].bucket_region)
            sides[side].bucket_region = ws.local_memory.alloc_bucket_region();
        memcpy(sides[side].bucket_region, o.sides[side].bucket_region, BUCKET_REGION);
    }
    SS = ws.local_memory.alloc_bucket_region();
    memcpy(SS, o.SS, BUCKET_REGION);
#endif
}
#endif

nfs_work::thread_data::thread_data(nfs_work & ws)
    : ws(ws)
{
    for(int side = 0 ; side < 2 ; side++) {
        sides[side].bucket_region = NULL;
        if (ws.sides[side].no_fb()) continue;
        sides[side].bucket_region = ws.local_memory.alloc_bucket_region();
    }
    /* Allocate memory for the intermediate sum (only one for both sides) */
    SS = ws.local_memory.alloc_bucket_region();
}

nfs_work::thread_data::~thread_data()
{
    for(int side = 0 ; side < 2 ; side++) {
        if (ws.sides[side].no_fb()) continue;
        ws.local_memory.free_bucket_region(sides[side].bucket_region);
        sides[side].bucket_region = NULL;
    }
    ws.local_memory.free_bucket_region(SS);
    SS = NULL;
}

/* only used once a siever_config has been attached to the parent structure */
void nfs_work::thread_data::allocate_bucket_regions()
{
    for(int side = 0 ; side < 2 ; side++) {
        if (ws.sides[side].no_fb()) {
            if (sides[side].bucket_region)
                ws.local_memory.free_bucket_region(sides[side].bucket_region);
        } else {
            if (!sides[side].bucket_region)
                sides[side].bucket_region = ws.local_memory.alloc_bucket_region();
        }
    }
    /* Allocate memory for the intermediate sum (only one for both sides) */
    if (!SS)
        SS = ws.local_memory.alloc_bucket_region();
}

void nfs_work::zeroinit_defaults()
{
    memset(nb_buckets, 0, sizeof(nb_buckets));
    toplevel = 0;
}

nfs_work::nfs_work(las_info & _las)
    : nfs_work(_las, NUMBER_OF_BAS_FOR_THREADS(_las.number_of_threads_per_subjob()))
{
}
nfs_work::nfs_work(las_info & _las, int nr_workspaces)
    : las(_las),
    local_memory(_las.local_memory_accessor()),
    nr_workspaces(nr_workspaces),
    sides {{ {nr_workspaces}, {nr_workspaces} }}
{
    zeroinit_defaults();
    // we cannot do this because thread_data has no copy ctor (on
    // purpose)
    // th(_las.number_of_threads_per_subjob(), thread_data(*this))
    th.reserve(_las.number_of_threads_per_subjob());
    for(int x = _las.number_of_threads_per_subjob() ; x-- ; )
        th.emplace_back(*this);
    sides[0].dumpfile.open(las.dump_filename, Q.doing, 0);
    sides[1].dumpfile.open(las.dump_filename, Q.doing, 1);
}

nfs_work_cofac::nfs_work_cofac(las_info & las, nfs_work const & ws) :
    las(las),
    sc(ws.conf),
    doing(ws.Q.doing)
{
    strategies = las.get_strategies(sc);
}

/* Prepare to work on sieving a special-q as described by _si.
   This implies allocating all the memory we need for bucket arrays,
   sieve regions, etc. */
void nfs_work::allocate_buckets(nfs_aux & aux, thread_pool & pool)
{
    /* We always allocate the max number of buckets (i.e., as if we were
     * using the max value for J), even if we use a smaller J due to a
     * poor q-lattice basis */ 

    verbose_output_print(0, 2, "# Reserving buckets with a multiplier of %s\n",
            bk_multiplier.print_all().c_str());

    bool do_resieve = conf.sides[0].lim && conf.sides[1].lim;

    for (int side = 0; side < 2; side++) {
        side_data & wss(sides[side]);
        if (wss.no_fb()) continue;
        wss.group.allocate_buckets(
                local_memory,
                nb_buckets,
                bk_multiplier,
                wss.fbs->stats.weight,
                conf.logI,
                aux, pool, do_resieve);
    }
    pool.drain_queue(2);
}

/* only used once a siever_config has been attached to the structure */
void nfs_work::allocate_bucket_regions() {
    for(auto & T : th)
        T.allocate_bucket_regions();
}

template <int LEVEL, typename HINT>
double
nfs_work::buckets_max_full()
{
    /* find the most full bucket across all buckets in the bucket array */
    double maxfull_ratio = 0;
    int maxfull_side = -1;
    unsigned int maxfull_index = 0;
    size_t maxfull_updates = 0;
    size_t maxfull_room = 0;
    typedef bucket_array_t<LEVEL, HINT> BA_t;
    for(int side = 0 ; side < 2 ; side++) {
        side_data & wss(sides[side]);
        for (auto const & BA : wss.bucket_arrays<LEVEL, HINT>()) {
            unsigned int index;
            const double ratio = BA.max_full(&index);
            if (ratio > maxfull_ratio) {
                maxfull_ratio = ratio;
                maxfull_side = side;
                maxfull_index = index;
                maxfull_updates = BA.nb_of_updates(index);
                maxfull_room = BA.room_allocated_for_updates(index);
            }
        }
    }
    if (maxfull_ratio > 1) {
        int side = maxfull_side;
        side_data & wss(sides[side]);
        auto const & BAs = wss.bucket_arrays<LEVEL, HINT>();
        std::ostringstream os;
        os << "bucket " << maxfull_index << " on side " << maxfull_side << ":";
        size_t m = 0;
        for (auto const & BA : wss.bucket_arrays<LEVEL, HINT>()) {
            if (BA.nb_of_updates(maxfull_index) >= m)
                m = BA.nb_of_updates(maxfull_index);
        }
        for (auto const & BA : wss.bucket_arrays<LEVEL, HINT>()) {
            size_t z = BA.nb_of_updates(maxfull_index);
            os << " " << z;
            if (z == m) os << "*";
        }
        os << " /" << BAs[0].room_allocated_for_updates(maxfull_index);

        auto k = bkmult_specifier::getkey<typename BA_t::update_t>();
        verbose_output_print(0, 1, "# Error: %s buckets are full, worst is %s\n",

                bkmult_specifier::printkey(k).c_str(),
                os.str().c_str());

        throw buckets_are_full(
                k,
                maxfull_index,
                maxfull_updates,
                maxfull_room);
    }
    return maxfull_ratio;
}
double nfs_work::check_buckets_max_full()
{
    double mf0 = 0, mf;
    mf = buckets_max_full<3, shorthint_t>(); if (mf > mf0) mf0 = mf;
    mf = buckets_max_full<2, shorthint_t>(); if (mf > mf0) mf0 = mf;
    mf = buckets_max_full<2, longhint_t>();  if (mf > mf0) mf0 = mf;
    mf = buckets_max_full<1, shorthint_t>(); if (mf > mf0) mf0 = mf;
    mf = buckets_max_full<1, longhint_t>();  if (mf > mf0) mf0 = mf;
    mf = buckets_max_full<3, emptyhint_t>(); if (mf > mf0) mf0 = mf;
    mf = buckets_max_full<2, emptyhint_t>(); if (mf > mf0) mf0 = mf;
    mf = buckets_max_full<2, logphint_t>();  if (mf > mf0) mf0 = mf;
    mf = buckets_max_full<1, emptyhint_t>(); if (mf > mf0) mf0 = mf;
    mf = buckets_max_full<1, logphint_t>();  if (mf > mf0) mf0 = mf;
    return mf0;
}

template <typename HINT>
double nfs_work::check_buckets_max_full(int level)
{
    switch(level) {
        case 3:
            if (std::is_same<HINT, shorthint_t>::value)
                return buckets_max_full<3, shorthint_t>();
            else if (std::is_same<HINT, emptyhint_t>::value)
                return buckets_max_full<3, emptyhint_t>();
            else
                ASSERT_ALWAYS(0);
            /* shouldn't be needed, but g++-12 likes to complain for nothing */
            no_break();
        case 2: return buckets_max_full<2, HINT>();
        case 1: return buckets_max_full<1, HINT>();
    }
    return 0;
}

template double nfs_work::buckets_max_full<1, shorthint_t>();
template double nfs_work::buckets_max_full<2, shorthint_t>();
template double nfs_work::buckets_max_full<3, shorthint_t>();
template double nfs_work::buckets_max_full<1, longhint_t>();
template double nfs_work::buckets_max_full<2, longhint_t>();
template double nfs_work::buckets_max_full<1, emptyhint_t>();
template double nfs_work::buckets_max_full<2, emptyhint_t>();
template double nfs_work::buckets_max_full<3, emptyhint_t>();
template double nfs_work::buckets_max_full<1, logphint_t>();
template double nfs_work::buckets_max_full<2, logphint_t>();
template double nfs_work::check_buckets_max_full<shorthint_t>(int);
template double nfs_work::check_buckets_max_full<longhint_t>(int);
template double nfs_work::check_buckets_max_full<emptyhint_t>(int);
template double nfs_work::check_buckets_max_full<logphint_t>(int);

template <int LEVEL, typename HINT>
void
nfs_work::side_data::reset_all_pointers() {
    group.get<LEVEL, HINT>().reset_all_pointers();
}

template void nfs_work::side_data::reset_all_pointers<1, shorthint_t>();
template void nfs_work::side_data::reset_all_pointers<2, shorthint_t>();
template void nfs_work::side_data::reset_all_pointers<3, shorthint_t>();
template void nfs_work::side_data::reset_all_pointers<1, emptyhint_t>();
template void nfs_work::side_data::reset_all_pointers<2, emptyhint_t>();
template void nfs_work::side_data::reset_all_pointers<3, emptyhint_t>();

void nfs_work::compute_toplevel_and_buckets()
{
    // Now that fb have been initialized, we can set the toplevel.
    toplevel = -1;
    for(int side = 0 ; side < 2 ; side++) {
        side_data & wss(sides[side]);
        if (wss.no_fb()) continue;

        int level = wss.fbs->get_toplevel();
        if (level > toplevel) toplevel = level;
    }
    ASSERT_ALWAYS(toplevel >= 1);

    /* update number of buckets at toplevel */
    size_t (&BRS)[FB_MAX_PARTS] = BUCKET_REGIONS;

    for(int i = 0 ; i < FB_MAX_PARTS ; ++i) nb_buckets[i] = 0;

    nb_buckets[toplevel] = iceildiv(1UL << conf.logA, BRS[toplevel]);

    // maybe there is only 1 bucket at toplevel and less than 256 at
    // toplevel-1, due to a tiny J.
    if (toplevel > 1) {
        if (nb_buckets[toplevel] == 1) {
            nb_buckets[toplevel-1] = iceildiv(1UL << conf.logA, BRS[toplevel - 1]);
            // we forbid skipping two levels.
            ASSERT_ALWAYS(nb_buckets[toplevel-1] != 1);
        } else {
            nb_buckets[toplevel-1] = BRS[toplevel]/BRS[toplevel-1];
        }
    }
}

void nfs_work::prepare_for_new_q(las_info & las0) {
    ASSERT_ALWAYS(&las == &las0);
    /* The config on which we're running is now decided. In order to
     * select the factor base to use, we also need the log scale */
    for(int side = 0 ; side < 2 ; side++) {
        nfs_work::side_data & wss(sides[side]);

        /* Even when we have no factor base, we do the lognorm setup
         * because we need it for the norm computation. It's admittedly a
         * bit ridiculous, given that only initializing lognorm_base
         * would be sufficient. The cost should be negligible, though.
         */
        wss.lognorms = lognorm_smart(conf, las.cpoly, side, Q, conf.logI, J);

        if (las.no_fb(side)) {
            wss.fbs = NULL;
            continue;
        }

        wss.fbK = conf.instantiate_thresholds(side);
        wss.fbK.scale = wss.lognorms.scale;
        wss.fbK.nb_threads = las.number_of_threads_per_subjob();

        /* Now possibly trigger the creation of a new slicing. There's a
         * design decision of whether we want the slicing replicated on
         * all sub-jobs, on sub-jobs with the same memory binding only,
         * or in a loose, unbound fashion. I think the second option is
         * better.
         */
        wss.fbs = las0.get_factorbase_slicing(side, wss.fbK);
        wss.td = las0.get_trialdiv_data(side, wss.fbK, wss.fbs);
    }
    bk_multiplier = las0.get_bk_multiplier();
    compute_toplevel_and_buckets();

    jd = las0.get_j_divisibility_helper(J);
    us = las0.get_unsieve_data(conf);

    /* we may now allocate the bucket regions for all threads. Those are
     * quite unsignificant of course, but in cases where we have gobs of
     * threads, we'd rather not allocate RAM for memory we won't actually
     * need because of the absence of factor base on some side.
     */
    allocate_bucket_regions();
}

/* Yes, it's quite unfortunate that we add so much red tape.
 */
struct helper_functor_precomp_plattice_dense_clear {
    template<typename T>
        void operator()(T & x) { x.clear(); }
};
void nfs_work::side_data::precomp_plattice_dense_clear()
{
    helper_functor_precomp_plattice_dense_clear H;
    multityped_array_foreach(H, precomp_plattice_dense);
}
