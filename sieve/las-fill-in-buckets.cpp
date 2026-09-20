#include "cado.h" // IWYU pragma: keep

/* This compilation units reacts to TRACK_CODE_PATH and uses macros
 * such as WHERE_AM_I_UPDATE.
 * This compilation unit _must_ produce different object files depending
 * on the value of TRACK_CODE_PATH.
 * The WHERE_AM_I_UPDATE macro itself is defined in las-where-am-i.hpp
 */

#include <cstddef>

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
#include "las-process-bucket-region.hpp"
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

// At top level.
// We need to interleave the root transforms and the FK walk,
// otherwise, we spend all the time waiting for memory.
// Hence the ugly de-templatization.
// At some point, the code should be re-organized, I'm afraid.
template <int LEVEL, class FB_ENTRY_TYPE, hint_type TARGET_HINT>
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
template <int LEVEL, class FB_ENTRY_TYPE, hint_type TARGET_HINT>
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
template <int LEVEL, class FB_ENTRY_TYPE, hint_type TARGET_HINT>
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

template <int LEVEL, hint_type TARGET_HINT>
static void fill_in_buckets_one_side(nfs_work & ws, nfs_aux & aux,
                                     ALGO::special_q_data const & Q,
                                     thread_pool & pool, int const side,
                                     where_am_I & w)
{
    timetree_t & timer(aux.rt.timer);
    nfs_work::side_data & wss(ws.sides[side]);

    // Early exit if this side does not reach LEVEL + 1
    if (wss.fbs->get_toplevel() < LEVEL)
        return;

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

void fill_in_buckets_toplevel_entry(nfs_work & ws, nfs_aux & aux,
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
