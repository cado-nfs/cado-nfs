#include "cado.h" // IWYU pragma: keep

#include <cstring>

#include <algorithm>
#include <ostream>
#include <string>
#include <sstream>
#include <utility>

#include "bucket.hpp"
#include "ecm/batch.hpp"
#include "las-config.hpp"
#include "las-bkmult.hpp"
#include "las-info.hpp"
#include "las-memory.hpp"
#include "las-threads-work-data.hpp"
#include "multityped_array.hpp"
#include "macros.h"
#include "threadpool.hpp"
#include "verbose.h"

class nfs_aux; // IWYU pragma: keep

nfs_work::thread_data::thread_data(thread_data && o) noexcept
    : ws(o.ws),
    sides(std::move(o.sides)),
    SS(o.SS)
{
    for(auto & s: o.sides)
        s.bucket_region = nullptr;
    o.SS = nullptr;
}

#if 0
/* See comment in las-threads-work-data.hpp */
nfs_work::thread_data::thread_data(thread_data const & o)
    : ws(o.ws)
{
    for(unsigned int side = 0 ; side < sides.size() ; side++) {
        ASSERT_ALWAYS(o.sides[side].bucket_region == nullptr);
        sides[side].bucket_region = nullptr;
    }
    ASSERT_ALWAYS(o.SS == nullptr);
    SS = nullptr;
#if 0
    /* We do not need MEMSET_MIN here. However we're making life easier
     * for the memory allocator if we allocate and free always the same
     * size.
     */
    for(unsigned int side = 0 ; side < sides.size() ; side++) {
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
    : ws(ws),
    sides(ws.sides.size())
{
    for(unsigned int side = 0 ; side < sides.size() ; side++) {
        sides[side].bucket_region = nullptr;
        if (ws.sides[side].no_fb()) continue;
        sides[side].bucket_region = ws.local_memory.alloc_bucket_region();
    }
    /* Allocate memory for the intermediate sum (only one for both sides) */
    SS = ws.local_memory.alloc_bucket_region();
}

nfs_work::thread_data::~thread_data()
{
    for(unsigned int side = 0 ; side < sides.size() ; side++) {
        if (ws.sides[side].no_fb()) continue;
        ws.local_memory.free_bucket_region(sides[side].bucket_region);
        sides[side].bucket_region = nullptr;
    }
    ws.local_memory.free_bucket_region(SS);
    SS = nullptr;
}

/* only used once a siever_config has been attached to the parent structure */
void nfs_work::thread_data::allocate_bucket_regions()
{
    for(unsigned int side = 0 ; side < sides.size() ; side++) {
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

nfs_work::nfs_work(las_info & _las, sieve_method auto tag)
    : nfs_work(_las, NUMBER_OF_BAS_FOR_THREADS(_las.number_of_threads_per_subjob()), tag)
{
}

template nfs_work::nfs_work(las_info &, NFS);
template nfs_work::nfs_work(las_info &, SIQS);

nfs_work::nfs_work(las_info & _las, int nr_workspaces, sieve_method auto tag)
    : las(_las),
    local_memory(_las.local_memory_accessor()),
    nr_workspaces(nr_workspaces),
    sides {{ {nr_workspaces, tag}, {nr_workspaces, tag} }}
{
    zeroinit_defaults();
    // we cannot do this because thread_data has no copy ctor (on
    // purpose)
    // th(_las.number_of_threads_per_subjob(), thread_data(*this))
    th.reserve(_las.number_of_threads_per_subjob());
    for(int x = _las.number_of_threads_per_subjob() ; x-- ; )
        th.emplace_back(*this);
    for (size_t i = 0; i < sides.size(); ++i) {
      ASSERT_ALWAYS(!las.dump_filename);
      /* Well, I don't think that doing is set to anything non-default
       * at this point: Q is only initialized by choose_sieve_area.
       */
      sides[i].dumpfile.open(las.dump_filename, doing, i);
    }
}

nfs_work_cofac::nfs_work_cofac(las_info & las, nfs_work const & ws) :
    las(las),
    sc(ws.conf),
    doing(ws.doing),
    strategies(las.get_strategies(sc))
{
}

/* Prepare to work on sieving a special-q as described by _si.
   This implies allocating all the memory we need for bucket arrays,
   sieve regions, etc. */
void nfs_work::allocate_buckets(nfs_aux & aux, thread_pool & pool)
{
    /* We always allocate the max number of buckets (i.e., as if we were
     * using the max value for J), even if we use a smaller J due to a
     * poor q-lattice basis */ 

    verbose_fmt_print(0, 2,
            "# Reserving buckets with a multiplier of {}\n",
            bk_multiplier.print_all());

    bool const do_resieve = conf.needs_resieving();

    for (auto & wss : sides) {
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
nfs_work::buckets_max_full() const
{
    /* find the most full bucket across all buckets in the bucket array */
    double maxfull_ratio = 0;
    int maxfull_side = -1;
    unsigned int maxfull_index = 0;
    size_t maxfull_updates = 0;
    size_t maxfull_room = 0;
    using BA_t = bucket_array_t<LEVEL, HINT>;
    for(unsigned int side = 0 ; side < sides.size() ; side++) {
        side_data  const& wss(sides[side]);
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
        int const side = maxfull_side;
        side_data  const& wss(sides[side]);
        auto const & BAs = wss.bucket_arrays<LEVEL, HINT>();
        std::ostringstream os;
        os << "bucket " << maxfull_index << " on side " << maxfull_side << ":";
        size_t m = 0;
        for (auto const & BA : wss.bucket_arrays<LEVEL, HINT>()) {
            if (BA.nb_of_updates(maxfull_index) >= m)
                m = BA.nb_of_updates(maxfull_index);
        }
        for (auto const & BA : wss.bucket_arrays<LEVEL, HINT>()) {
            size_t const z = BA.nb_of_updates(maxfull_index);
            os << " " << z;
            if (z == m) os << "*";
        }
        os << " /" << BAs[0].room_allocated_for_updates(maxfull_index);

        auto k = bkmult_specifier::getkey<typename BA_t::update_t>();
        verbose_fmt_print(0, 1,
                "# Error: {} buckets are full, worst is {}\n",
                bkmult_specifier::printkey(k),
                os.str());

        throw buckets_are_full(
                k,
                maxfull_index,
                maxfull_updates,
                maxfull_room);
    }
    return maxfull_ratio;
}
double nfs_work::check_buckets_max_full() const
{
    double mf = 0;

    mf = std::max(mf, buckets_max_full<1, shorthint_t>());
    mf = std::max(mf, buckets_max_full<1, emptyhint_t>());

#if MAX_TOPLEVEL >= 2
    mf = std::max(mf, buckets_max_full<1, longhint_t>());
    mf = std::max(mf, buckets_max_full<1, logphint_t>());
    mf = std::max(mf, buckets_max_full<2, shorthint_t>());
    mf = std::max(mf, buckets_max_full<2, emptyhint_t>());
#endif

#if MAX_TOPLEVEL >= 3
    mf = std::max(mf, buckets_max_full<2, longhint_t>());
    mf = std::max(mf, buckets_max_full<2, logphint_t>());
    mf = std::max(mf, buckets_max_full<3, shorthint_t>());
    mf = std::max(mf, buckets_max_full<3, emptyhint_t>());
#endif

    static_assert(MAX_TOPLEVEL == 3);

    return mf;
}

double nfs_work::check_buckets_max_full_toplevel(int level) const
{

    double mf = 0;

    if (level == 1) {
        mf = std::max(mf, buckets_max_full<1, shorthint_t>());
        mf = std::max(mf, buckets_max_full<1, emptyhint_t>());
        return mf;
    }

#if MAX_TOPLEVEL >= 2
    if (level == 2) {
        mf = std::max(mf, buckets_max_full<2, shorthint_t>());
        mf = std::max(mf, buckets_max_full<2, emptyhint_t>());
        return mf;
    }
#endif

#if MAX_TOPLEVEL >= 3
    if (level == 3) {
        mf = std::max(mf, buckets_max_full<3, shorthint_t>());
        mf = std::max(mf, buckets_max_full<3, emptyhint_t>());
        return mf;
    }
#endif
    static_assert(MAX_TOPLEVEL == 3);

    ASSERT_ALWAYS(0);

}

#if 0
template <typename HINT>
double nfs_work::check_buckets_max_full(int level)
    requires (HINT::allowed_at_toplevel)
{
    static_assert(MAX_TOPLEVEL == 3);

#if MAX_TOPLEVEL >= 3
    if (level == 3)
        return buckets_max_full<3, HINT>();
#endif
#if MAX_TOPLEVEL >= 2
    if (level == 2)
        return buckets_max_full<2, HINT>();
#endif
    if (level == 1)
        return buckets_max_full<1, HINT>();
    ASSERT_ALWAYS(0);
}
template <typename HINT>
double nfs_work::check_buckets_max_full(int level)
    requires (!(HINT::allowed_at_toplevel))
{
    static_assert(MAX_TOPLEVEL == 3);

#if MAX_TOPLEVEL >= 3
    if (level == 2)
        return buckets_max_full<2, HINT>();
#endif

#if MAX_TOPLEVEL >= 2
    if (level == 1)
        return buckets_max_full<1, HINT>();
#endif

    ASSERT_ALWAYS(0);
}
#endif

template double nfs_work::buckets_max_full<1, shorthint_t>() const;
template double nfs_work::buckets_max_full<1, emptyhint_t>() const;

#if MAX_TOPLEVEL >= 2
template double nfs_work::buckets_max_full<2, shorthint_t>() const;
template double nfs_work::buckets_max_full<2, emptyhint_t>() const;
template double nfs_work::buckets_max_full<1, longhint_t>() const;
template double nfs_work::buckets_max_full<1, logphint_t>() const;
#endif

#if MAX_TOPLEVEL >= 3
template double nfs_work::buckets_max_full<3, shorthint_t>() const;
template double nfs_work::buckets_max_full<3, emptyhint_t>() const;
template double nfs_work::buckets_max_full<2, longhint_t>() const;
template double nfs_work::buckets_max_full<2, logphint_t>() const;
#endif

static_assert(MAX_TOPLEVEL == 3);

#if 0
template double nfs_work::check_buckets_max_full<shorthint_t>(int) const;
template double nfs_work::check_buckets_max_full<emptyhint_t>(int) const;
template double nfs_work::check_buckets_max_full<longhint_t>(int) const;
template double nfs_work::check_buckets_max_full<logphint_t>(int) const;
#endif

template <int LEVEL, typename HINT>
void
nfs_work::side_data::reset_all_pointers() {
    group.get<LEVEL, HINT>().reset_all_pointers();
}

template void nfs_work::side_data::reset_all_pointers<1, shorthint_t>();
template void nfs_work::side_data::reset_all_pointers<1, emptyhint_t>();
#if MAX_TOPLEVEL >= 2
template void nfs_work::side_data::reset_all_pointers<2, shorthint_t>();
template void nfs_work::side_data::reset_all_pointers<2, emptyhint_t>();
#endif
#if MAX_TOPLEVEL >= 3
template void nfs_work::side_data::reset_all_pointers<3, shorthint_t>();
template void nfs_work::side_data::reset_all_pointers<3, emptyhint_t>();
#endif
static_assert(MAX_TOPLEVEL == 3);

void nfs_work::compute_toplevel_and_buckets()
{
    // Now that fb have been initialized, we can set the toplevel.
    // XXX TODO: we should decouple the toplevel and the sieving side,
    // really.
    toplevel = -1;
    for(unsigned int side = 0 ; side < sides.size() ; side++) {
        side_data  const& wss(sides[side]);
        if (wss.no_fb()) continue;

        toplevel = std::max(toplevel, wss.fbs->get_toplevel());
    }
    ASSERT_ALWAYS(toplevel >= 1 && toplevel <= MAX_TOPLEVEL);

    /* update number of buckets at toplevel */
    size_t  const(&BRS)[FB_MAX_PARTS] = BUCKET_REGIONS;

    std::fill_n(nb_buckets, FB_MAX_PARTS, 0);

    size_t above = 1UL << conf.logA;

    for(int level = toplevel ; level >= 1 ; level--) {
        // For level < toplevel, we typically have BRS[level]/BRS[level-1]
        // buckets. However we tolerate the case where we have only one
        // toplevel bucket. It is a waste of resources, however. Using a
        // less intricate sieving process would be better.
        nb_buckets[level] = iceildiv(above, BRS[level]);
        if (nb_buckets[level] > 1)
            above = BRS[level];
        else
            fmt::print(stderr, "# WARNING: only 1 bucket at level {}. Consider dropping bkthresh{} setting.\n", level, level-1);

        /* Having this suboptimal situation for two consecutive levels
         * really makes no sense.
         */
        ASSERT_ALWAYS(level == toplevel || nb_buckets[level] != 1);
    }
}

template<sieve_method Algo>
void nfs_work::prepare_for_new_q(las_info & las0, special_q_task * task, typename Algo::special_q_data const & Q)
{
    ASSERT_ALWAYS(&las == &las0);
    this->task = task;

    /* The config on which we're running is now decided. In order to
     * select the factor base to use, we also need the log scale */
    for(int side = 0 ; side < las.cpoly->nb_polys ; side++) {
        nfs_work::side_data & wss(sides[side]);

        /* Even when we have no factor base, we do the lognorm setup
         * because we need it for the norm computation. It's admittedly a
         * bit ridiculous, given that only initializing lognorm_base
         * would be sufficient. The cost should be negligible, though.
         */
        wss.lognorms = lognorm_smart(conf, las.cpoly, side, Q, conf.logI, J);

        if (las.no_fb(side)) {
            wss.fbs = nullptr;
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

template void nfs_work::prepare_for_new_q<NFS>(las_info & las0, special_q_task * task, NFS::special_q_data const & Q);
template void nfs_work::prepare_for_new_q<SIQS>(las_info & las0, special_q_task * task, SIQS::special_q_data const & Q);

void nfs_work::side_data::precomp_plattice_dense_clear()
{
    precomp_plattice_dense.foreach([](auto & x) { x.clear(); });
}
