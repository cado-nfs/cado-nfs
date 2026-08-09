#include "cado.h" // IWYU pragma: keep

/* This compilation units reacts to TRACK_CODE_PATH and uses macros
 * such as WHERE_AM_I_UPDATE.
 * This compilation unit _must_ produce different object files depending
 * on the value of TRACK_CODE_PATH.
 * The WHERE_AM_I_UPDATE macro itself is defined in las-where-am-i.hpp
 */

#include <cstddef>
#include <cstdint>

#include <algorithm>
#include <array>
#include <functional>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

#include "bucket.hpp"
#include "bucket-push-update.hpp"       // IWYU pragma: keep
#include "chronograms.hpp"
#include "fb-types.hpp"
#include "fb.hpp"
#include "las-auxiliary-data.hpp"
#include "las-bkmult.hpp"
#include "las-config.hpp"
#include "las-fill-in-buckets.hpp"
#include "las-globals.hpp"
#include "las-process-bucket-region.hpp"
#include "las-report-stats.hpp"
#include "las-siever-config.hpp"
#include "las-threads-work-data.hpp"
#include "las-where-am-i-proxy.hpp"
#include "las-where-am-i.hpp"
#include "macros.h"
#include "multityped_array.hpp"
#include "tdict.hpp"
#include "threadpool.hpp"
#include "verbose.hpp"

/***************************************************************************/
/********        Main bucket sieving functions                    **********/

/* {{{ Big question: shall we enable bucket-sieving for powers ?
 *
 * There are several difficulties, in fact. One rationale that yields a
 * straight "no" answer is that such primes make very little difference
 * to the smooth part, so we'd better skip them anyway.
 *
 * But it's not the hardest thing.
 *
 * For the small sieve, we create the small_sieve_data from the factor
 * base entries, and we compute the logp accordingly, per entry.
 *
 * For the bucket-sieve, we use the fact that the factor base is sorted
 * in increasing log(p) order, and we create slices with ranges of primes
 * that have the same round(scale*log(p)).
 *
 * Currently, the factor base is sorted by q=p^k. A power that makes the
 * p-valuation go from p^k0 to p^k1 contributes
 * round(k1*log(p))-round(k0*log(p)). Therefore, sorting by q does not
 * mean that log(p)'s are sorted, and we're in trouble because when we
 * take powers aboard in a slice, their log(p) value is not correctly
 * represented.
 *
 * Previously, we had the behaviour of setting powlim to bucket_thresh-1,
 * effectively preventing powers from appearing in the bucket-sieve.
 *
 * Now powlim is a factor base parameter, and bucket_thresh comes later,
 * so such a default does not work.
 *
 * The strategy we take here is that *if* we see powers here (and we know
 * that will happen only for the fairly rare general entries), then we do
 * something special:
 *  - either we say that we skip over this entry
 *  - or we defer to apply_buckets the expensive computation of a proper
 *    logp value.
 *
 * Currently we do the former. The latter would be viable since only a
 * small fraction of the apply_one_bucket time is devoted to dealing with
 * general entries, so we could imagine having a branch in there for
 * dealing with them. But that would be quite painful. Furthermore, it
 * would then be mandatory to split the entries with same q, same p, but
 * different k0,k1 pairs (we do encounter these), so that the hint would
 * still open the possibility to infer the value of log(p).
 *
 *
 * Note that it would not be possible to sidestep the issue by sorting
 * the vectors of entries by (k1-k0)*log(p) (which would make a
 * difference only for general entries anyway). This is because even
 * sorting by increasing (k1-k0)*log(p) does not guarantee that
 * round(s*k1*log(p))-round(s*k0*log(p)) increases. (counter-example:
 * s=1, k1*log(p)=0.51, k0*log(p)=0.49 diff=0.02 round-round=1
 *      k1*log(p)=1.49, k0*log(p)=0.51 diff=0.98 round-round=0
 * )
 * }}} */
template <class FB_ENTRY_TYPE>
static inline bool discard_power_for_bucket_sieving(FB_ENTRY_TYPE const &)
{
    /* the entry is not a general entry, therefore k is a const thing
     * equal to 1.
     */
    return false;
}
#ifndef BUCKET_SIEVE_POWERS
template <>
inline bool
discard_power_for_bucket_sieving<fb_entry_general>(fb_entry_general const & e)
{
    return e.k > 1;
}
#endif

/***********************************************************************/
/* multithreaded processing of make_lattice_bases (a.k.a
 * precomp_plattices)
 */


#ifdef SIQS_SIEVE
#include "siqs-fill-in-buckets.inl"
#else
#include "las-fill-in-buckets.inl"
#endif

void fill_in_buckets_prepare_plattices(
        nfs_work & ws,
        ALGO::special_q_data const & Q,
        thread_pool & pool,
        int side,
        cado::multityped_array<precomp_plattice_t, 1, FB_MAX_PARTS-1> & precomp_plattice)
{
    /* this will *not* do anything for level==ws.toplevel, by design */
    precomp_plattice.foreach([&](auto & V) {
        /* T is precomp_plattice_t<n> for some level n */
        using T = std::remove_reference_t<decltype(V)>;
        if (T::level >= ws.toplevel)
            return;

        nfs_work::side_data const & wss(ws.sides[side]);
        fb_factorbase::slicing::part const & P = wss.fbs->get_part(T::level);
        V.clear();
        V.resize(P.nslices());
        P.slices.foreach([&](auto const & sl) {
                for(auto const & s : sl) {
                    using E = std::remove_reference_t<decltype(s)>::entry_t;
                    // Forward function pointer & arguments directly!
                    pool.add_task(thread_pool::QUEUE_GENERIC, 0.0,
                            make_lattice_bases<T::level, E>,
                            side, std::ref(ws), std::ref(Q), std::ref(V), std::cref(s));
                }
        });
    });
}

/* {{{ */


// For internal levels, the fill-in is not exactly the same as for
// top-level, since the plattices have already been precomputed.
template <int LEVEL, typename TARGET_HINT>
static void
fill_in_buckets_one_slice_internal(worker_thread * worker,
        nfs_work & ws,
        nfs_aux & aux,
        ALGO::special_q_data const & Q,
        int side,
        fb_slice_interface *,
        plattices_vector_t * plattices_vector,
        plattices_dense_vector_t *,
        uint32_t const first_region0_index)
{
    static_assert(!TARGET_HINT::is_long_v);

    /* Import some contextual stuff */
    int const id = worker->rank();
    nfs_aux::thread_data & taux(aux.th[id]);
    timetree_t & timer(aux.get_timer(worker));
    ENTER_THREAD_TIMER(timer);
    where_am_I & w(taux.w);
    nfs_work::side_data & wss(ws.sides[side]);

    MARK_TIMER_FOR_SIDE(timer, side);

    // we're declaring the timer here, but really the work happens below
    // in fill_in_buckets_lowlevel. We happen to have access to
    // param->side here, so we use it to provide a nicer timing report.
    CHILD_TIMER(timer,
            fmt::format("fill_in_buckets_one_slice_internal<{}>", LEVEL));

    WHERE_AM_I_UPDATE(w, side, side);
    WHERE_AM_I_UPDATE(w, i, plattices_vector->get_index());
    WHERE_AM_I_UPDATE(w, N, first_region0_index);

    try {
        auto acquired = wss.reserve_BA<LEVEL, TARGET_HINT>();
        auto tt = worker->trace(chronograms::FIB(
                    side,
                    LEVEL,
                    wss.rank_BA(acquired.access()),
                    plattices_vector->get_index()));
            
        /* Get an unused bucket array that we can write to */
        /* clearly, reserve_BA() possibly throws. As it turns out,
         * fill_in_buckets_lowlevel<> does not, at least currently. One
         * could imagine that it could throw, so let's wrap it too.
         */
        fill_in_buckets_lowlevel<LEVEL, TARGET_HINT>(
                acquired.access(),
                ws, Q, *plattices_vector,
                first_region0_index, w);
    } catch (buckets_are_full & e) {
        e.side = side;
        throw e;
    }
}

// At top level.
// We need to interleave the root transforms and the FK walk,
// otherwise, we spend all the time waiting for memory.
// Hence the ugly de-templatization.
// At some point, the code should be re-organized, I'm afraid.
template <int LEVEL, class FB_ENTRY_TYPE, typename TARGET_HINT>
static void
fill_in_buckets_toplevel_wrapper(worker_thread * worker,
        int side,
        nfs_work & ws,
        nfs_aux & aux,
        ALGO::special_q_data const & Q,
        plattices_dense_vector_t * plattices_dense_vector,
        fb_slice<FB_ENTRY_TYPE> const & slice)
{
    static_assert(!TARGET_HINT::is_long_v);

    /* Import some contextual stuff */
    int const id = worker->rank();
    nfs_aux::thread_data & taux(aux.th[id]);
    timetree_t & timer(aux.get_timer(worker));
    where_am_I & w(taux.w);
    nfs_work::side_data & wss(ws.sides[side]);

    ENTER_THREAD_TIMER(timer);
    MARK_TIMER_FOR_SIDE(timer, side);

#ifndef DISABLE_TIMINGS
    /* This is one of the places where helgrind is likely to complain. We
     * use thread-safe statics. Helgrind can't cope with it,
     * unfortunately. So the error is a false positive.
     *
     * https://sourceforge.net/p/valgrind/mailman/message/32434015/
     */
    timetree_t::accounting_child const local_timer_sentry(timer,
                                                          tdict_slot_for_fibt);
#endif

    WHERE_AM_I_UPDATE(w, side, side);
    WHERE_AM_I_UPDATE(w, i, slice.get_index());
    WHERE_AM_I_UPDATE(w, N, 0);

    try {
        /* Get an unused bucket array that we can write to */
        // bucket_array_t<LEVEL, TARGET_HINT> &;
        auto acquired = wss.reserve_BA<LEVEL, TARGET_HINT>();
        auto tt = worker->trace(chronograms::FIB(
                    side,
                    LEVEL,
                    wss.rank_BA(acquired.access()),
                    slice.get_index()));

        fill_in_buckets_toplevel<LEVEL, FB_ENTRY_TYPE, TARGET_HINT>(
                acquired.access(),
                ws, slice, Q, plattices_dense_vector, w);
        return;
    } catch (buckets_are_full & e) {
        e.side = side;
        throw e;
    }
}
/* same for sublat */
template <int LEVEL, class FB_ENTRY_TYPE, typename TARGET_HINT>
static void
fill_in_buckets_toplevel_sublat_wrapper(worker_thread * worker,
        int side,
        nfs_work & ws,
        nfs_aux & aux,
        ALGO::special_q_data const & Q,
        plattices_dense_vector_t * plattices_dense_vector,
        fb_slice<FB_ENTRY_TYPE> const & slice
        )
{
    static_assert(!TARGET_HINT::is_long_v);

    /* Import some contextual stuff */
    int const id = worker->rank();
    nfs_aux::thread_data & taux(aux.th[id]);
    timetree_t & timer(aux.get_timer(worker));
    where_am_I & w(taux.w);
    nfs_work::side_data & wss(ws.sides[side]);

    ENTER_THREAD_TIMER(timer);
    MARK_TIMER_FOR_SIDE(timer, side);

#ifndef DISABLE_TIMINGS
    /* This is one of the places where helgrind is likely to complain. We
     * use thread-safe statics. Helgrind can't cope with it,
     * unfortunately. So the error is a false positive.
     *
     * https://sourceforge.net/p/valgrind/mailman/message/32434015/
     */
    timetree_t::accounting_child const local_timer_sentry(timer,
                                                          tdict_slot_for_fibt);
#endif

    WHERE_AM_I_UPDATE(w, side, side);
    WHERE_AM_I_UPDATE(w, i, slice.get_index());
    WHERE_AM_I_UPDATE(w, N, 0);

    try {
        /* Get an unused bucket array that we can write to */
        auto acquired = wss.reserve_BA<LEVEL, TARGET_HINT>();
        auto tt = worker->trace(chronograms::FIB(
                    side,
                    LEVEL,
                    wss.rank_BA(acquired.access()),
                    slice.get_index()));
        fill_in_buckets_toplevel_sublat<LEVEL, FB_ENTRY_TYPE>(
                acquired.access(),
                ws, Q,
                plattices_dense_vector,
                slice, w);
    } catch (buckets_are_full & e) {
        e.side = side;
        throw e;
    }
}

// Static helper function outside loop to avoid lambda closure bloat inside foreach_slice
template <int LEVEL, class FB_ENTRY_TYPE, typename TARGET_HINT>
static void run_fill_in_buckets_toplevel(worker_thread* worker, int side, nfs_work & ws,
                                         nfs_aux & aux, ALGO::special_q_data const & Q,
                                         plattices_dense_vector_t * pre,
                                         fb_slice<FB_ENTRY_TYPE> const & s,
                                         std::shared_ptr<where_am_I> w_copy)
{
    int const id = worker->rank();
    aux.th[id].w = *w_copy;
    if (pre) {
        fill_in_buckets_toplevel_sublat_wrapper<LEVEL, FB_ENTRY_TYPE, TARGET_HINT>(
            worker, side, ws, aux, Q, pre, s);
    } else {
        fill_in_buckets_toplevel_wrapper<LEVEL, FB_ENTRY_TYPE, TARGET_HINT>(
            worker, side, ws, aux, Q, nullptr, s);
    }
}

template <int LEVEL, typename TARGET_HINT>
static void fill_in_buckets_one_side(nfs_work & ws, nfs_aux & aux,
                                     ALGO::special_q_data const & Q,
                                     thread_pool & pool, int const side,
                                     where_am_I & w)
{
    timetree_t & timer(aux.rt.timer);
    nfs_work::side_data & wss(ws.sides[side]);

    BOOKKEEPING_TIMER(timer);

    auto const & BA_ins = wss.bucket_arrays<LEVEL, TARGET_HINT>();

    verbose_fmt_print(0, 3,
            "# Filling the side-{} {}{} buckets ({} groups of {} buckets)\n",
            side,
            LEVEL, TARGET_HINT::rtti[0],
            BA_ins.size(), BA_ins[0].n_bucket);

    fb_factorbase::slicing::part const & P = wss.fbs->get_part(LEVEL);

    typename precomp_plattice_dense_t<LEVEL>::type * Vpre = nullptr;

    if (Q.sublat.m) {
        auto & Vpre_ref(wss.precomp_plattice_dense.get<LEVEL>());
        if (Q.sublat.i0 == 0 && Q.sublat.j0 == 1) {
            Vpre_ref = typename precomp_plattice_dense_t<LEVEL>::type(P.nslices());
        }
        ASSERT_ALWAYS(Vpre_ref.size() == P.nslices());
        Vpre = &Vpre_ref;
    }

    size_t pushed = 0;
    P.foreach_slice([&](auto const & s) {
        auto w_copy = std::make_shared<where_am_I>(w);
        slice_index_t const idx = s.get_index();
        ASSERT_ALWAYS(P.first_slice_index + pushed == idx);
        plattices_dense_vector_t * pre = Vpre ? &((*Vpre)[idx]) : nullptr;
        using entry_t = std::decay_t<decltype(s)>::entry_t;

        pool.add_task(thread_pool::QUEUE_GENERIC, s.get_weight(),
            run_fill_in_buckets_toplevel<LEVEL, entry_t, TARGET_HINT>,
            side, std::ref(ws), std::ref(aux), std::ref(Q), pre, std::cref(s), w_copy);

        pushed++;
    });
}

/* This is a compile-time loop over the possible values from 1 to level,
 * and 0 errors out. */
template<int level, hint_type hint_t>
struct fib1s_caller_s : public fib1s_caller_s<level-1, hint_t> {
    template<typename... Args>
    void operator()(nfs_work & ws, Args&& ...args) const {
        if (ws.toplevel == level)
            fill_in_buckets_one_side<level, hint_t>(ws, std::forward<Args>(args)...);
        else
            fib1s_caller_s<level-1, hint_t>::operator()(ws, std::forward<Args>(args)...);
    }
};
template<hint_type hint_t>
struct fib1s_caller_s<0, hint_t> {
    template<typename... Args>
    void operator()(nfs_work &, Args&& ...) const {
        ASSERT_ALWAYS(0);
    }
};

template<int level, hint_type hint_t, typename... Args>
inline void fib_one_side(nfs_work & ws, Args&& ...args)
{
    fib1s_caller_s<level, hint_t>()(ws, std::forward<Args>(args)...);
}

void fill_in_buckets_toplevel_multiplex(nfs_work & ws, nfs_aux & aux,
        ALGO::special_q_data const & Q, thread_pool & pool, int side,
        where_am_I & w)
{
    // per se, we're not doing anything here.
    // CHILD_TIMER(timer, __func__);
    if (ws.conf.needs_resieving()) {
        fib_one_side<MAX_TOPLEVEL, shorthint_t>(ws, aux, Q, pool, side, w);
    } else {
        fib_one_side<MAX_TOPLEVEL, emptyhint_t>(ws, aux, Q, pool, side, w);
    }
}

/* }}} */

/* multithreaded implementation of the downsort procedure. It becomes a
 * bottleneck sooner than one might think.
 *
 */

/* This is auxiliary only. We downsort stuff that we already downsorted.
 * So it applies only if LEVEL+1 is itself not the toplevel.
 * For this reason, we must have a specific instantiation that reduces
 * this to a no-op if LEVEL+1>=3, because there's no longhint_t for level
 * 3 presently.
 */
template <int LEVEL, bool WITH_HINTS>
static void downsort_aux(nfs_work & ws,
                         nfs_aux & aux, thread_pool & pool, task_group & tg, int side,
                         uint32_t bucket_index, where_am_I & w)
{
    static_assert(LEVEL <= MAX_TOPLEVEL - 1);

    using my_longhint_t = hints_proxy<WITH_HINTS>::l;

    nfs_work::side_data & wss(ws.sides[side]);

    auto const & BA_ins = wss.bucket_arrays<LEVEL + 1, my_longhint_t>();
    auto & BA_outs = wss.bucket_arrays<LEVEL, my_longhint_t>();
    ASSERT_ALWAYS(BA_ins.size() == BA_outs.size());

    verbose_fmt_print(0, 3,
            "# Downsorting the side-{} {}{} buckets ({} groups of {} buckets"
            ", taking only bucket {}/{})"
            " to {}{} buckets ({} groups of {} buckets)\n",
            side,
            LEVEL + 1, my_longhint_t::rtti[0],
            BA_ins.size(), BA_ins[0].n_bucket,
            bucket_index, BA_ins[0].n_bucket,
            LEVEL, my_longhint_t::rtti[0],
            BA_outs.size(), BA_outs[0].n_bucket);


    // What comes from already downsorted data above:
    for (auto const & BA_in: BA_ins) {
        pool.add_task(
            tg,
            [&aux, &ws, side, bucket_index, &BA_in](worker_thread * worker, where_am_I && w) {
                nfs_work::side_data & wss(ws.sides[side]);
                fb_factorbase::slicing const & fbs(*wss.fbs);
                int const id = worker->rank();
                nfs_aux::thread_data & taux(aux.th[id]);
                timetree_t & timer(aux.get_timer(worker));
                ENTER_THREAD_TIMER(timer);
                MARK_TIMER_FOR_SIDE(timer, side);
                taux.w = std::move(w);
                CHILD_TIMER(timer, fmt::format("downsort<{}>", LEVEL));
                auto tt = worker->trace(chronograms::DS(side,LEVEL,bucket_index));
                downsort<LEVEL + 1>(fbs,
                        wss.acquire_BA<LEVEL, my_longhint_t>(wss.rank_BA(BA_in)),
                        BA_in, bucket_index, taux.w);
            }, where_am_I(w));
    }
}

#if MAX_TOPLEVEL == 2
template <>
void downsort_aux<1, false>(nfs_work &, nfs_aux &,
                     thread_pool &, task_group &, int, uint32_t, where_am_I &)
{
}
template <>
void downsort_aux<1, true>(nfs_work &, nfs_aux &,
                     thread_pool &, task_group &, int, uint32_t, where_am_I &)
{
}
#endif
#if MAX_TOPLEVEL == 3
template <>
void downsort_aux<2, false>(nfs_work &, nfs_aux &,
                     thread_pool &, task_group &, int, uint32_t, where_am_I &)
{
}
template <>
void downsort_aux<2, true>(nfs_work &, nfs_aux &,
                     thread_pool &, task_group &, int, uint32_t, where_am_I &)
{
}
#endif
static_assert(MAX_TOPLEVEL == 3);

// first_region0_index is a way to remember where we are in the tree.
// The depth-first is a way to process all the the regions of level 0 in
// increasing order of j-value.
// first_region0_index * nb_lines_per_region0 therefore gives the j-line
// where we are. This is what is called N by WHERE_AM_I and friends.

template <int LEVEL, bool WITH_HINTS>
static void downsort_tree_inner(
    nfs_work & ws,
    std::shared_ptr<nfs_work_cofac> wc_p,
    std::shared_ptr<nfs_aux> aux_p,
    ALGO::special_q_data const & Q,
    thread_pool & pool,
    uint32_t bucket_index, /* for the current level ! */
    uint32_t first_region0_index,
    std::vector<cado::multityped_array<precomp_plattice_t, 1, FB_MAX_PARTS - 1>> & precomp_plattices,
    where_am_I & w)
{
    /* LEVEL is not the toplevel here, so we must have the following: */
    static_assert(LEVEL <= MAX_TOPLEVEL - 1);

    int const nsides = ws.las.cpoly.nsides();
    nfs_aux & aux(*aux_p);
    timetree_t & timer(aux.rt.timer);

    using my_longhint_t = hints_proxy<WITH_HINTS>::l;
    using my_shorthint_t = hints_proxy<WITH_HINTS>::s;

    CHILD_TIMER(timer, fmt::format("downsort_tree<{}>", LEVEL));
    TIMER_CATEGORY(timer, sieving_mixed());
    ASSERT_ALWAYS(LEVEL > 0);

    WHERE_AM_I_UPDATE(w, N, first_region0_index);

    std::vector<task_group> ds_tgs(nsides);
    std::vector<task_group> ds_aux_tgs(nsides);
    std::vector<task_group> fib_tgs(nsides);
    std::vector<task_group> sss_tgs(nsides);

    for (int side = 0; side < nsides; ++side) {
        nfs_work::side_data & wss(ws.sides[side]);
        if (wss.no_fb())
            continue;

        WHERE_AM_I_UPDATE(w, side, side);
        TIMER_CATEGORY(timer, sieving(side));

        auto & ds_tg = ds_tgs[side];
        auto & ds_aux_tg = ds_aux_tgs[side];
        auto & fib_tg = fib_tgs[side];

        /* FIRST: Downsort what is coming from the level above, for this
         * bucket index */
        // All these BA are global stuff; see reservation_group.
        // We reserve those where we write, and access the ones for
        // reading without reserving. We require that things at level
        // above are finished before entering here.

        {
            auto const & BA_ins = wss.bucket_arrays<LEVEL + 1, my_shorthint_t>();
            auto & BA_outs = wss.bucket_arrays<LEVEL, my_longhint_t>();
            /* otherwise the code here can't work */
            /* see also the comment above in downsort_aux */
            ASSERT_ALWAYS(BA_ins.size() == BA_outs.size());

            verbose_fmt_print(0, 3,
                    "# Downsorting the side-{} {}{} buckets ({} groups of {} buckets"
                    ", taking only bucket {}/{})"
                    " to {}{} buckets ({} groups of {} buckets)\n",
                    side,
                    LEVEL + 1, my_shorthint_t::rtti[0],
                    BA_ins.size(), BA_ins[0].n_bucket,
                    bucket_index, BA_ins[0].n_bucket,
                    LEVEL, my_longhint_t::rtti[0],
                    BA_outs.size(), BA_outs[0].n_bucket);

            ASSERT_ALWAYS(BA_ins.size() == BA_outs.size());

            wss.reset_all_pointers<LEVEL, my_longhint_t>();

            for (auto const & BA_in: BA_ins) {
                pool.add_task(
                    ds_tg,
                    [&aux, &ws, side, &BA_in, bucket_index](worker_thread * worker, where_am_I && w) {
                        nfs_work::side_data & wss(ws.sides[side]);
                        fb_factorbase::slicing const & fbs(*wss.fbs);
                        int const id = worker->rank();
                        nfs_aux::thread_data & taux(aux.th[id]);
                        timetree_t & timer(aux.get_timer(worker));
                        taux.w = std::move(w);
                        ENTER_THREAD_TIMER(timer);
                        MARK_TIMER_FOR_SIDE(timer, side);
                        CHILD_TIMER(timer, fmt::format("downsort<{}>", LEVEL));
                        auto tt = worker->trace(chronograms::DS(side, LEVEL, bucket_index));

                        downsort<LEVEL + 1>(fbs,
                                wss.acquire_BA<LEVEL, my_longhint_t>(wss.rank_BA(BA_in)),
                                BA_in, bucket_index,
                                taux.w);
                    },
                    where_am_I(w));
            }
        }

        /* Once ds_tg completes for this side, enqueue all FIB slice tasks into fib_tg */
        ds_tg.on_complete(
            [&ws, &aux, &pool, &precomp_plattices, &ds_aux_tg, &fib_tg, side, first_region0_index, &Q, bucket_index, w_val = where_am_I(w)]() mutable {
                /* SECOND: fill in buckets at this level, for this region. */
                /* We do so only once ds_tg completes for this side */
                auto do_fib = [&ws, &aux, &pool, &precomp_plattices, &fib_tg, side, first_region0_index, &Q]() {
                    nfs_work::side_data & wss(ws.sides[side]);
                    wss.reset_all_pointers<LEVEL, my_shorthint_t>();

                    auto & BA_outs = wss.bucket_arrays<LEVEL, my_shorthint_t>();
                    auto & lattices = precomp_plattices[side].get<LEVEL>();

                    verbose_fmt_print(0, 3,
                            "# Filling the side-{} {}{} buckets ({} groups of {} buckets)"
                            " using {} precomputed lattices\n",
                            side,
                            LEVEL, my_shorthint_t::rtti[0],
                            BA_outs.size(), BA_outs[0].n_bucket,
                            lattices.size());
                    if (!lattices.empty()) {
                        verbose_fmt_print(0, 3,
                                "#   lattices go from slice {} ({} primes) to slice {} ({} primes)\n",
                                lattices.front().get_index(), lattices.front().size(),
                                lattices.back().get_index(), lattices.back().size()
                            );
                    }

                    for (auto & it: lattices) {
                        pool.add_task(fib_tg, thread_pool::QUEUE_GENERIC, it.get_weight(),
                                fill_in_buckets_one_slice_internal<LEVEL, my_shorthint_t>,
                                std::ref(ws), std::ref(aux), std::ref(Q), side, nullptr, &it, nullptr, first_region0_index);
                    }
                };

                if (LEVEL < ws.toplevel - 1) {
                    downsort_aux<LEVEL, WITH_HINTS>(ws, aux, pool, ds_aux_tg, side, bucket_index, w_val);
                    ds_aux_tg.on_complete(do_fib);
                } else {
                    do_fib();
                }
            });
    }

    if (LEVEL == 1) {
        /* Prepare for PBR: we need to precompute the small sieve positions
         * for all the small sieved primes.
         *
         * For ws.toplevel==1, we don't reach here, of course, and the
         * corresponding initialization is done with identical code in
         * las.cpp
         */
        ASSERT(ws.toplevel > 1);
        for (int side = 0; side < nsides; side++) {
            nfs_work::side_data const & wss(ws.sides[side]);
            if (wss.no_fb())
                continue;

            auto & sss_tg = sss_tgs[side];
            wss.ssd->small_sieve_prepare_many_start_positions(
                    pool, &sss_tg,
                    first_region0_index,
                    std::min(SMALL_SIEVE_START_POSITIONS_MAX_ADVANCE, ws.nb_buckets[1]),
                    ws.conf.logI, Q.sublat);
            sss_tg.on_complete([&wss]() {
                wss.ssd->small_sieve_activate_many_start_positions();
            });
        }
    }

    for (int side = 0; side < nsides; side++) {
        nfs_work::side_data const & wss(ws.sides[side]);
        if (wss.no_fb())
            continue;
        ds_tgs[side].wait();
        ds_aux_tgs[side].wait();
        fib_tgs[side].wait();
        if (LEVEL == 1)
            sss_tgs[side].wait();
    }

    /* RECURSE */
    if (LEVEL > 1) {
        size_t const(&BRS)[FB_MAX_PARTS] = BUCKET_REGIONS;
        verbose_fmt_print(0, 3,
                "# recursively downsort level-{} buckets ({} buckets)"
                " to level {} (+ fill {}{} buckets). Target bucket indices: {}..{}\n",
                LEVEL, ws.nb_buckets[LEVEL], LEVEL-1,
                LEVEL - 1, my_shorthint_t::rtti[0],
                first_region0_index,
                first_region0_index + ws.nb_buckets[LEVEL] * BRS[LEVEL] / BRS[1]);

        for (int i = 0; i < ws.nb_buckets[LEVEL]; ++i) {
            /* This is quite suspicious. Shouldn't we do BRS[LEVEL] /
             * BRS[LEVEL - 1] instead?
             */
            uint32_t const N = first_region0_index + i * (BRS[LEVEL] / BRS[1]);
            downsort_tree<LEVEL - 1>(ws, wc_p, aux_p, Q, pool, i, N,
                                     precomp_plattices, w);
        }
    } else {
        /* Now fill_in_buckets has completed for all levels. Time to check
         * that we had no overflow, and move on to process_bucket_region.
         */

        ws.check_buckets_max_full();

        if (print_slice_statistics)
            /* It may make sense to print a diagnosis of the fill ratio
             * of the different buckets. It's quite verbose, though */
            ws.slice_statistics(LEVEL);

        auto exc = pool.get_exceptions<buckets_are_full>(thread_pool::QUEUE_GENERIC);
        if (!exc.empty()) {
            throw *std::ranges::max_element(exc);
        }

        // it seems difficult to compute the max target bucket index, in
        // fact. Well of course it should be ws.nb_buckets[1], but just
        // based on the input that we have, it's less obvious.
        // size_t const(&BRS)[FB_MAX_PARTS] = BUCKET_REGIONS;
        verbose_fmt_print(0, 3,
                "# calling process_bucket_region"
                " on regions of indices {}..\n",
                first_region0_index);

        /* PROCESS THE REGIONS AT LEVEL 0 */
        process_many_bucket_regions(ws, wc_p, aux_p, Q, pool, first_region0_index, w);

        /* We need that, because the next downsort_tree call in the loop
         * above (for LEVEL>1) will reset the pointers while filling the 1l
         * buckets -- and we read the 1l buckets from PBR.
         */
        if (ws.toplevel > 1)
            pool.drain_queue(thread_pool::QUEUE_GENERIC);
    }
}

template <int LEVEL>
void downsort_tree(
    nfs_work & ws,
    std::shared_ptr<nfs_work_cofac> wc_p,
    std::shared_ptr<nfs_aux> aux_p,
    ALGO::special_q_data const & Q,
    thread_pool & pool,
    uint32_t bucket_index, /* for the current level ! */
    uint32_t first_region0_index,
    std::vector<cado::multityped_array<precomp_plattice_t, 1, FB_MAX_PARTS - 1>> & precomp_plattices,
    where_am_I & w)
{
    if (ws.conf.needs_resieving()) {
        downsort_tree_inner<LEVEL, true>(ws, wc_p, aux_p, Q, pool, bucket_index,
                                         first_region0_index, precomp_plattices,
                                         w);
    } else {
        downsort_tree_inner<LEVEL, false>(ws, wc_p, aux_p, Q, pool, bucket_index,
                                          first_region0_index, precomp_plattices,
                                          w);
    }
}
/* Instances to be compiled */

// A fake level 0, to avoid infinite loop during compilation.
template <>
void downsort_tree<0>(nfs_work &,
                      std::shared_ptr<nfs_work_cofac>,
                      std::shared_ptr<nfs_aux>,
                      ALGO::special_q_data const &,
                      thread_pool &, uint32_t,
                      uint32_t,
                      std::vector<cado::multityped_array<precomp_plattice_t, 1, FB_MAX_PARTS - 1>> &,
                      where_am_I &)
{
    ASSERT_ALWAYS(0);
}

// Now the exported instances

template void downsort_tree<1>(
    nfs_work &, std::shared_ptr<nfs_work_cofac>, std::shared_ptr<nfs_aux> aux_p,
    ALGO::special_q_data const &, thread_pool &, uint32_t, uint32_t,
    std::vector<cado::multityped_array<precomp_plattice_t, 1, FB_MAX_PARTS - 1>> &, where_am_I &);


#if MAX_TOPLEVEL >= 3
template void
downsort_tree<2>(nfs_work &, std::shared_ptr<nfs_work_cofac>,
                 std::shared_ptr<nfs_aux>, ALGO::special_q_data const & Q,
                 thread_pool &, uint32_t, uint32_t,
                 std::vector<cado::multityped_array<precomp_plattice_t, 1, FB_MAX_PARTS - 1>> &,
                 where_am_I &);
#endif
static_assert(MAX_TOPLEVEL == 3);
