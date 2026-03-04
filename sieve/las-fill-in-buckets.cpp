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
#include <memory>
#include <utility>
#include <vector>
#include <type_traits>

#include "bucket-push-update.hpp"
#include "bucket.hpp"
#include "fb-types.hpp"
#include "fb.hpp"
#include "las-auxiliary-data.hpp"
#include "las-bkmult.hpp"
#include "las-config.hpp"
#include "las-fill-in-buckets.hpp"
#include "las-process-bucket-region.hpp"
#include "las-qlattice.hpp"
#include "las-report-stats.hpp"
#include "las-siever-config.hpp"
#include "las-smallsieve.hpp"
#include "las-threads-work-data.hpp"
#include "las-where-am-i-proxy.hpp"
#include "las-where-am-i.hpp"
#include "macros.h"
#include "multityped_array.hpp"
#include "tdict.hpp"
#include "threadpool.hpp"
#include "utils_cxx.hpp"
#include "verbose.h"

/* is this in the std library or not ? */
template <typename T> static inline T const & const_ref(T & x)
{
    return x;
}

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
 *
 * This creates one control object per slice, with storage ownership of
 * the control object transfered to the called function. Because we
 * depend on the slice, the type of the object is parameterized by the
 * slice type.
 *
 * We may elect to make the "model" a shared_ptr.
 */

template <int LEVEL>
struct make_lattice_bases_parameters_base : public task_parameters {
    int side;
    nfs_work & ws;
    ALGO::special_q_data const & Q;
    precomp_plattice_t<LEVEL> & V;
    make_lattice_bases_parameters_base(int side, nfs_work & ws,
            ALGO::special_q_data const & Q,
            precomp_plattice_t<LEVEL> & V)
        : side(side)
        , ws(ws)
        , Q(Q)
        , V(V)
    {
    }
};
template <int LEVEL, class FB_ENTRY_TYPE>
struct make_lattice_bases_parameters
    : public make_lattice_bases_parameters_base<LEVEL> {
    using super = make_lattice_bases_parameters_base<LEVEL>;
    fb_slice<FB_ENTRY_TYPE> const & slice;
    make_lattice_bases_parameters(super const & model,
                                  fb_slice<FB_ENTRY_TYPE> const & slice)
        : super(model)
        , slice(slice)
    {
    }
};

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
    precomp_plattice.foreach([&](auto & precomp_plattice) {
        /* T is precomp_plattice_t<n> for some level n */
        using T = std::remove_reference_t<decltype(precomp_plattice)>;
        if (T::level >= ws.toplevel)
            return;

        nfs_work::side_data const & wss(ws.sides[side]);
        fb_factorbase::slicing::part const & P = wss.fbs->get_part(T::level);
        /* We pre-assign the result, so that all threads can write to it
         * comfortably.
         *
         * It would be nice to have a way to notify that all threads here are
         * done with their job.
         */
        precomp_plattice.assign(P.nslices(), plattices_vector_t());
        make_lattice_bases_parameters_base<T::level> const model {side, ws, Q, precomp_plattice};
        P.slices.foreach([&](auto const & sl) {
                for(auto const & s : sl) {
                    using E = std::remove_reference_t<decltype(s)>::entry_t;
                    auto param = new make_lattice_bases_parameters<T::level, E>(model, s);
                    task_function_t f = make_lattice_bases<T::level, E>;
                    pool.add_task(f, param, 0);
                }
        });
    });
}

/* {{{ */
template <int LEVEL> class fill_in_buckets_parameters : public task_parameters
{
  public:
    nfs_work & ws;
    nfs_aux & aux;
    ALGO::special_q_data const & Q;
    int const side;
    fb_slice_interface const * slice;
    plattices_vector_t * plattices_vector; // content changed during fill-in
    plattices_dense_vector_t * plattices_dense_vector; // for sublat
    uint32_t const first_region0_index;
    where_am_I w;

    fill_in_buckets_parameters(nfs_work & _ws, nfs_aux & aux,
                               ALGO::special_q_data const & Q, int const _side,
                               fb_slice_interface const * _slice,
                               plattices_vector_t * _platt,
                               plattices_dense_vector_t * _dplatt,
                               uint32_t const _reg0, where_am_I const & w)
        : ws(_ws)
        , aux(aux)
        , Q(Q)
        , side(_side)
        , slice(_slice)
        , plattices_vector(_platt)
        , plattices_dense_vector(_dplatt)
        , first_region0_index(_reg0)
        , w(w)
    {
    }
};

/* short of a better solution. I know some exist, but it seems way
 * overkill to me.
 *
 * This needs constexpr, though... So maybe I could use a more powerful
 * C++11 trick after all.
 */
#define PREPARE_TEMPLATE_INST_NAMES(F, suffix)                                 \
    template <int> struct CADO_CONCATENATE(F, _name) {                         \
    };                                                                         \
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 0);                                  \
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 1);                                  \
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 2);                                  \
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 3);                                  \
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 4);                                  \
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 5);                                  \
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 6);                                  \
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 7);                                  \
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 8);                                  \
    PREPARE_TEMPLATE_INST_NAME(F, suffix, 9)

#define PREPARE_TEMPLATE_INST_NAME(F, suffix, k)                               \
    template <> struct CADO_CONCATENATE(F, _name)<k> {                         \
        static constexpr const char * value = #F "<" #k ">" suffix;            \
    }

/* By tweaking the "" argument, it is possible to have these names
 * embody a suffix like " on side ", so that it's possible tu run
 * parametric timer slots.
 */
PREPARE_TEMPLATE_INST_NAMES(fill_in_buckets_one_slice_internal, "");
PREPARE_TEMPLATE_INST_NAMES(downsort, "");
PREPARE_TEMPLATE_INST_NAMES(downsort_tree, " (dispatcher only)");

#define TEMPLATE_INST_NAME(x, y) CADO_CONCATENATE(x, _name)<y>::value

// For internal levels, the fill-in is not exactly the same as for
// top-level, since the plattices have already been precomputed.
template <int LEVEL, typename TARGET_HINT>
static task_result *
fill_in_buckets_one_slice_internal(worker_thread * worker,
                                   task_parameters * _param, int)
{
    auto * param = static_cast<fill_in_buckets_parameters<LEVEL> *>(_param);

    /* Import some contextual stuff */
    int const id = worker->rank();
    nfs_aux::thread_data & taux(param->aux.th[id]);
    timetree_t & timer(param->aux.get_timer(worker));
    ENTER_THREAD_TIMER(timer);
    nfs_work & ws(param->ws);
    ALGO::special_q_data const & Q(param->Q);
    where_am_I & w(taux.w);
    int const side = param->side;
    nfs_work::side_data & wss(ws.sides[side]);

    MARK_TIMER_FOR_SIDE(timer, side);

    // we're declaring the timer here, but really the work happens below
    // in fill_in_buckets_lowlevel. We happen to have access to
    // param->side here, so we use it to provide a nicer timing report.
    CHILD_TIMER(timer,
                TEMPLATE_INST_NAME(fill_in_buckets_one_slice_internal, LEVEL));

    w = param->w;
    WHERE_AM_I_UPDATE(w, side, param->side);
    WHERE_AM_I_UPDATE(w, i, param->plattices_vector->get_index());
    WHERE_AM_I_UPDATE(w, N, param->first_region0_index);

    try {
        /* Get an unused bucket array that we can write to */
        /* clearly, reserve_BA() possibly throws. As it turns out,
         * fill_in_buckets_lowlevel<> does not, at least currently. One
         * could imagine that it could throw, so let's wrap it too.
         */
        auto & BA = wss.reserve_BA<LEVEL, TARGET_HINT>(-1);

        /* Fill the buckets */
        try {
            fill_in_buckets_lowlevel<LEVEL>(BA, ws, Q, *param->plattices_vector,
                                            param->first_region0_index, w);
        } catch (buckets_are_full & e) {
            wss.release_BA(BA);
            throw e;
        }
        /* Release bucket array again */
        wss.release_BA(BA);
    } catch (buckets_are_full & e) {
        delete param;
        throw e;
    }
    delete param;
    return new task_result;
}

// At top level.
// We need to interleave the root transforms and the FK walk,
// otherwise, we spend all the time waiting for memory.
// Hence the ugly de-templatization.
// At some point, the code should be re-organized, I'm afraid.
template <int LEVEL, class FB_ENTRY_TYPE, typename TARGET_HINT>
static task_result *
fill_in_buckets_toplevel_wrapper(worker_thread * worker MAYBE_UNUSED,
                                 task_parameters * _param, int)
{
    auto * param = static_cast<fill_in_buckets_parameters<LEVEL> *>(_param);

    /* Import some contextual stuff */
    int const id = worker->rank();
    nfs_aux::thread_data & taux(param->aux.th[id]);
    timetree_t & timer(param->aux.get_timer(worker));
    nfs_work & ws(param->ws);
    ALGO::special_q_data const & Q(param->Q);
    where_am_I & w(taux.w);
    int const side = param->side;
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

    WHERE_AM_I_UPDATE(w, side, param->side);
    WHERE_AM_I_UPDATE(w, i, param->slice->get_index());
    WHERE_AM_I_UPDATE(w, N, 0);

    try {
        /* Get an unused bucket array that we can write to */
        bucket_array_t<LEVEL, TARGET_HINT> & BA =
            wss.reserve_BA<LEVEL, TARGET_HINT>(-1);
        ASSERT(param->slice);
        auto const * sl = dynamic_cast<fb_slice<FB_ENTRY_TYPE> const *>(param->slice);
        ASSERT_ALWAYS(sl != NULL);
        fill_in_buckets_toplevel<LEVEL, FB_ENTRY_TYPE, TARGET_HINT>(
            BA, ws, *sl, Q, param->plattices_dense_vector, w);
        /* Release bucket array again */
        wss.release_BA(BA);
        delete param;
        return new task_result;
    } catch (buckets_are_full const & e) {
        delete param;
        throw e;
    }
}
/* same for sublat */
template <int LEVEL, class FB_ENTRY_TYPE, typename TARGET_HINT>
static task_result *
fill_in_buckets_toplevel_sublat_wrapper(worker_thread * worker,
                                        task_parameters * _param, int)
{
    auto * param = static_cast<fill_in_buckets_parameters<LEVEL> *>(_param);

    /* Import some contextual stuff */
    int const id = worker->rank();
    nfs_aux::thread_data & taux(param->aux.th[id]);
    timetree_t & timer(param->aux.get_timer(worker));
    nfs_work & ws(param->ws);
    ALGO::special_q_data const & Q(param->Q);
    where_am_I & w(taux.w);
    int const side = param->side;
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

    WHERE_AM_I_UPDATE(w, side, param->side);
    WHERE_AM_I_UPDATE(w, i, param->slice->get_index());
    WHERE_AM_I_UPDATE(w, N, 0);

    try {
        /* Get an unused bucket array that we can write to */
        bucket_array_t<LEVEL, TARGET_HINT> & BA =
            wss.reserve_BA<LEVEL, TARGET_HINT>(-1);
        ASSERT(param->slice);
        fill_in_buckets_toplevel_sublat<LEVEL, FB_ENTRY_TYPE>(
            BA, ws, Q,
            *dynamic_cast<fb_slice<FB_ENTRY_TYPE> const *>(param->slice),
            param->plattices_dense_vector, w);
        /* Release bucket array again */
        wss.release_BA(BA);
        delete param;
        return new task_result;
    } catch (buckets_are_full const & e) {
        delete param;
        throw e;
    }
}

/* Whether or not we want fill_in_buckets_one_slice to be templatized
 * both for LEVEL and n is not clear. At some point, we're doing code
 * bloat for almost nothing.
 *
 * Now given the code below, it's easy enough to arrange so that we go
 * back to the virtual base fb_slice_interface.
 */
template <int LEVEL, typename TARGET_HINT> struct push_slice_to_task_list {
    thread_pool & pool;
    fill_in_buckets_parameters<LEVEL> model;
    push_slice_to_task_list(thread_pool & pool,
                            fill_in_buckets_parameters<LEVEL> const & m)
        : pool(pool)
        , model(m)
    {
    }
    size_t pushed = 0;
    template <typename T> void operator()(T const & s)
    {
        auto * param = new fill_in_buckets_parameters<LEVEL>(model);
        param->slice = &s;
        using entry_t = typename T::entry_t;
        task_function_t f =
            fill_in_buckets_toplevel_wrapper<LEVEL, entry_t, TARGET_HINT>;
        pool.add_task(f, param, 0, 0, s.get_weight());
        pushed++;
    }
};
template <int LEVEL, typename TARGET_HINT>
struct push_slice_to_task_list_saving_precomp {
    thread_pool & pool;
    fb_factorbase::slicing::part const & P;
    fill_in_buckets_parameters<LEVEL> model;
    /* precomp_plattice_dense_t == std::vector<plattices_dense_vector_t> */
    typename precomp_plattice_dense_t<LEVEL>::type & Vpre;
    size_t pushed = 0;
    push_slice_to_task_list_saving_precomp(
        thread_pool & pool, fb_factorbase::slicing::part const & P,
        fill_in_buckets_parameters<LEVEL> const & m,
        typename precomp_plattice_dense_t<LEVEL>::type & Vpre)
        : pool(pool)
        , P(P)
        , model(m)
        , Vpre(Vpre)
    {
    }
    template <typename T> void operator()(T const & s)
    {
        /* we're pushing the global index, relative to all fb parts */
        slice_index_t const idx = s.get_index();
        ASSERT_ALWAYS((size_t)idx == pushed);

        plattices_dense_vector_t & pre(Vpre[idx]);

        auto * param = new fill_in_buckets_parameters<LEVEL>(model);
        param->slice = &s;
        param->plattices_dense_vector = &pre;

        using entry_t = typename T::entry_t;
        task_function_t f =
            fill_in_buckets_toplevel_sublat_wrapper<LEVEL, entry_t,
                                                    TARGET_HINT>;
        pool.add_task(f, param, 0, 0, s.get_weight());
        pushed++;
    }
};

template <int LEVEL, typename TARGET_HINT>
static void fill_in_buckets_one_side(nfs_work & ws, nfs_aux & aux,
                                     ALGO::special_q_data const & Q,
                                     thread_pool & pool, int const side,
                                     where_am_I & w)
{
    timetree_t & timer(aux.rt.timer);
    nfs_work::side_data & wss(ws.sides[side]);

    /* We're just pushing tasks, here. */
    BOOKKEEPING_TIMER(timer);

    fill_in_buckets_parameters<LEVEL> const model(ws, aux, Q, side, NULL, NULL,
                                                  NULL, 0, w);

    auto const & BA_ins = wss.bucket_arrays<LEVEL, TARGET_HINT>();

    verbose_fmt_print(0, 3,
            "# Filling the side-{} {}{} buckets ({} groups of {} buckets)\n",
            side,
            LEVEL, TARGET_HINT::rtti[0],
            BA_ins.size(), BA_ins[0].n_bucket);

    /* We'd like to also display info on the slices we're about to run
     * FIB on, but the multityped_array interface won't let us do this
     * easily.
     */

    if (!Q.sublat.m) {
        /* This creates a task meant to call
         * fill_in_buckets_toplevel_wrapper */
        push_slice_to_task_list<LEVEL, TARGET_HINT> F(pool, model);
        wss.fbs->get_part(LEVEL).foreach_slice(F);
    } else {
        /* This creates a task meant to call
         * fill_in_buckets_toplevel_sublat_wrapper */
        auto & Vpre(wss.precomp_plattice_dense.get<LEVEL>());
        fb_factorbase::slicing::part const & P = wss.fbs->get_part(LEVEL);

        /* This way we can spare the need to expose the copy contructor
         * of the container's value_type */
        if (Q.sublat.i0 == 0 && Q.sublat.j0 == 1) {
            /* first sublat */
            Vpre = typename precomp_plattice_dense_t<LEVEL>::type(P.nslices());
        }

        ASSERT_ALWAYS(Vpre.size() == P.nslices());
        push_slice_to_task_list_saving_precomp<LEVEL, TARGET_HINT> F(
            pool, P, model, Vpre);
        P.foreach_slice(F);
    }
}


/* This is a compile-time loop over the possible values from 1 to level,
 * and 0 errors out. */
template<int level, typename hint_t>
struct fib1s_caller_s : public fib1s_caller_s<level-1, hint_t> {
    template<typename... Args>
    void operator()(nfs_work & ws, Args&& ...args) const {
        if (ws.toplevel == level)
            fill_in_buckets_one_side<level, hint_t>(ws, std::forward<Args>(args)...);
        else
            fib1s_caller_s<level-1, hint_t>::operator()(ws, std::forward<Args>(args)...);
    }
};
template<typename hint_t>
struct fib1s_caller_s<0, hint_t> {
    template<typename... Args>
    void operator()(nfs_work &, Args&& ...) const {
        ASSERT_ALWAYS(0);
    }
};

template<int level, typename hint_t, typename... Args>
inline void fib_one_side(nfs_work & ws, Args&& ...args)
{
    fib1s_caller_s<level, hint_t>()(ws, std::forward<Args>(args)...);
}

void fill_in_buckets_toplevel_multiplex(nfs_work & ws, nfs_aux & aux,
        ALGO::special_q_data const & Q, thread_pool & pool, int side, where_am_I & w)
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
static void downsort_aux(fb_factorbase::slicing const & fbs, nfs_work & ws,
                         nfs_aux & aux, thread_pool & pool, int side,
                         uint32_t bucket_index, where_am_I & w)
{
    static_assert(LEVEL <= MAX_TOPLEVEL - 1);

    using my_longhint_t = hints_proxy<WITH_HINTS>::l;

    nfs_work::side_data & wss(ws.sides[side]);

    auto const & BA_ins = wss.bucket_arrays<LEVEL + 1, my_longhint_t>();
    auto & BA_outs = wss.bucket_arrays<LEVEL, my_longhint_t>();

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
        pool.add_task_lambda(
            [&, side, w](worker_thread * worker, int bucket_index) {
                nfs_aux::thread_data & taux(aux.th[worker->rank()]);
                timetree_t & timer(aux.get_timer(worker));
                ENTER_THREAD_TIMER(timer);
                MARK_TIMER_FOR_SIDE(timer, side);
                taux.w = w;
                CHILD_TIMER(timer, TEMPLATE_INST_NAME(downsort, LEVEL));
                auto & BA_out(
                    wss.reserve_BA<LEVEL, my_longhint_t>(wss.rank_BA(BA_in)));
                downsort<LEVEL + 1>(fbs, BA_out, BA_in, bucket_index, taux.w);
                wss.template release_BA<LEVEL, my_longhint_t>(BA_out);
            },
            bucket_index, 0);
    }
}

#if MAX_TOPLEVEL == 2
template <>
void downsort_aux<1, false>(fb_factorbase::slicing const &, nfs_work &, nfs_aux &,
                     thread_pool &, int, uint32_t, where_am_I &)
{
}
template <>
void downsort_aux<1, true>(fb_factorbase::slicing const &, nfs_work &, nfs_aux &,
                     thread_pool &, int, uint32_t, where_am_I &)
{
}
#endif
#if MAX_TOPLEVEL == 3
template <>
void downsort_aux<2, false>(fb_factorbase::slicing const &, nfs_work &, nfs_aux &,
                     thread_pool &, int, uint32_t, where_am_I &)
{
}
template <>
void downsort_aux<2, true>(fb_factorbase::slicing const &, nfs_work &, nfs_aux &,
                     thread_pool &, int, uint32_t, where_am_I &)
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

    int const nsides = ws.las.cpoly->nb_polys;
    nfs_aux & aux(*aux_p);
    timetree_t & timer(aux.rt.timer);

    using my_longhint_t = hints_proxy<WITH_HINTS>::l;
    using my_shorthint_t = hints_proxy<WITH_HINTS>::s;

    CHILD_TIMER(timer, TEMPLATE_INST_NAME(downsort_tree, LEVEL));
    TIMER_CATEGORY(timer, sieving_mixed());
    ASSERT_ALWAYS(LEVEL > 0);

    WHERE_AM_I_UPDATE(w, N, first_region0_index);

    for (int side = 0; side < nsides; ++side) {
        nfs_work::side_data & wss(ws.sides[side]);
        if (wss.no_fb())
            continue;

        WHERE_AM_I_UPDATE(w, side, side);
        TIMER_CATEGORY(timer, sieving(side));
        /* FIRST: Downsort what is coming from the level above, for this
         * bucket index */
        // All these BA are global stuff; see reservation_group.
        // We reserve those where we write, and access the ones for
        // reading without reserving. We require that things at level
        // above are finished before entering here.

        {
            /* This is the "dictionary" that maps slice indices to actual fb
             * entries. We rarely need it, except when downsorting short entries
             * in the case where we've eliminated the hint
             */
            fb_factorbase::slicing const & fbs(*wss.fbs);

            auto const & BA_ins = wss.bucket_arrays<LEVEL + 1, my_shorthint_t>();
            auto & BA_outs = wss.bucket_arrays<LEVEL, my_longhint_t>();

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

            /* otherwise the code here can't work */
            ASSERT_ALWAYS(BA_ins.size() == BA_outs.size());

            /* We create one output array for each input array, and we
             * process them in parallel. There would be various ways to
             * achieve that.
             */
            for (auto const & BA_in: BA_ins) {
                pool.add_task_lambda(
                    [&, side, w](worker_thread * worker, int bucket_index) {
                        nfs_aux::thread_data & taux(aux.th[worker->rank()]);
                        timetree_t & timer(aux.get_timer(worker));
                        taux.w = w;
                        ENTER_THREAD_TIMER(timer);
                        MARK_TIMER_FOR_SIDE(timer, side);
                        CHILD_TIMER(timer, TEMPLATE_INST_NAME(downsort, LEVEL));
                        auto & BA_out(wss.reserve_BA<LEVEL, my_longhint_t>(
                            wss.rank_BA(BA_in)));
                        BA_out.reset_pointers();
                        downsort<LEVEL + 1>(fbs, BA_out, BA_in, bucket_index,
                                            taux.w);
                        //wss.template release_BA<LEVEL, my_longhint_t>(BA_out);
                        wss.release_BA(BA_out);
                    },
                    bucket_index, 0);
            }
            // What comes from already downsorted data above. We put this in
            // an external function because we need the code to be elided or
            // LEVEL >= 2.
            if (LEVEL < ws.toplevel - 1) {
                pool.drain_queue(0);
                downsort_aux<LEVEL, WITH_HINTS>(fbs, ws, aux, pool, side, bucket_index, w);
            }
        }

        /* There might be a performance hit here, and honestly I'm not
         * 100% sure it's useful. The F9_sieve_3_levels test wants it,
         * and apparently really wants it here.
         */
        pool.drain_queue(0);

        {
            /* SECOND: fill in buckets at this level, for this region. */
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
                pool.add_task(
                        fill_in_buckets_one_slice_internal<LEVEL, my_shorthint_t>,
                        new fill_in_buckets_parameters<LEVEL> {
                        ws, aux, Q, side, (fb_slice_interface *)NULL, &it, NULL,
                        first_region0_index, const_ref(w)},
                        0, 0, it.get_weight());
            }
        }
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
            pool.add_task_lambda(
                [=, &ws, &aux](worker_thread * worker, int) {
                    timetree_t & timer(aux.get_timer(worker));
                    ENTER_THREAD_TIMER(timer);
                    MARK_TIMER_FOR_SIDE(timer, side);
                    SIBLING_TIMER(timer, "prepare small sieve");
                    nfs_work::side_data & wss(ws.sides[side]);
                    // if (wss.no_fb()) return;
                    SIBLING_TIMER(timer, "small sieve start positions");
                    /* When we're doing 2-level sieving, there is probably
                     * no real point in doing ssdpos initialization in
                     * several passes.
                     */
                    wss.ssd->small_sieve_prepare_many_start_positions(
                        first_region0_index,
                        std::min(SMALL_SIEVE_START_POSITIONS_MAX_ADVANCE,
                                 ws.nb_buckets[1]),
                        ws.conf.logI, Q.sublat);
                    wss.ssd->small_sieve_activate_many_start_positions();
                },
                0);
        }

        pool.drain_queue(0);
        /* Now fill_in_buckets has completed for all levels. Time to check
         * that we had no overflow, and move on to process_bucket_region.
         */

        ws.check_buckets_max_full();
        auto exc = pool.get_exceptions<buckets_are_full>(0);
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
        // first_region0_index + ws.nb_buckets[ws.toplevel] * BRS[ws.toplevel] / BRS[1]);

        /* PROCESS THE REGIONS AT LEVEL 0 */
        process_many_bucket_regions(ws, wc_p, aux_p, Q, pool, first_region0_index, w);

        /* We need that, because the next downsort_tree call in the loop
         * above (for LEVEL>1) will reset the pointers while filling the 1l
         * buckets -- and we read the 1l buckets from PBR.
         */
        if (ws.toplevel > 1)
            pool.drain_queue(0);
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
