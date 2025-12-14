#include "cado.h" // IWYU pragma: keep

/* This compilation units reacts to TRACK_CODE_PATH and uses macros
 * such as WHERE_AM_I_UPDATE.
 * This compilation unit _must_ produce different object files depending
 * on the value of TRACK_CODE_PATH.
 * The WHERE_AM_I_UPDATE macro itself is defined in las-where-am-i.hpp
 */

#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <array>
#include <condition_variable>
#include <fstream>
#include <functional>
#include <list>
#include <map>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <dirent.h>

#include <gmp.h>
#include "fmt/base.h"
#include "fmt/format.h"

#include "bucket.hpp"
#include "cado-sighandlers.h"
#include "cado_poly.h"
#include "cxx_mpz.hpp"
#include "ecm/batch.hpp"
#include "ecm/facul_strategies_stats.hpp"
#include "fb-types.hpp"
#include "fb.hpp"
#include "gmp_aux.h"
#include "json.hpp"
#include "las-auxiliary-data.hpp"
#include "las-bkmult.hpp"
#include "las-choose-sieve-area.hpp"
#include "las-cofactor.hpp"
#include "las-config.hpp"
#include "las-divide-primes.hpp"
#include "las-dlog-base.hpp"
#include "las-duplicate.hpp"
#include "las-fill-in-buckets.hpp"
#include "las-globals.hpp"
#include "las-info.hpp"
#include "las-multiobj-globals.hpp"
#include "las-norms.hpp"
#include "las-output.hpp"
#include "las-parallel.hpp"
#include "las-plattice.hpp"
#include "las-process-bucket-region.hpp"
#include "las-qlattice.hpp"
#include "las-report-stats.hpp"
#include "las-side-config.hpp"
#include "las-sieve-shared-data.hpp"
#include "las-siever-config.hpp"
#include "las-smallsieve.hpp"
#include "las-special-q-task-collection.hpp"
#include "las-special-q-task.hpp"
#include "las-threads-work-data.hpp"
#include "las-todo-list.hpp"
#include "las-where-am-i-proxy.hpp"
#include "las-where-am-i.hpp"
#include "macros.h"
#include "memusage.h"
#include "misc.h"
#include "mpz_poly.h"
#include "multityped_array.hpp"
#include "params.h"
#include "relation.hpp"
#include "special-q.hpp"
#include "tdict.hpp"
#include "threadpool.hpp"
#include "timing.h"
#include "utils_cxx.hpp"
#include "verbose.h"

/*************************** main program ************************************/

static void configure_aliases(cxx_param_list & pl)
{
    las_info::configure_aliases(pl);
    param_list_configure_alias(pl, "log-bucket-region", "B");
    param_list_configure_alias(pl, "log-bucket-region-step", "Bi");
    las_output::configure_aliases(pl);
    tdict::configure_aliases(pl);
}

static void configure_switches(cxx_param_list & pl)
{
    las_info::configure_switches(pl);
    las_output::configure_switches(pl);
    tdict::configure_switches(pl);

    param_list_configure_switch(pl, "-allow-largesq", &allow_largesq);
    if (dlp_descent)
        param_list_configure_switch(pl, "-recursive-descent", &recursive_descent);

    param_list_configure_switch(pl, "-prepend-relation-time", &prepend_relation_time);
    param_list_configure_switch(pl, "-sync", &sync_at_special_q);
    param_list_configure_switch(pl, "-sync-thread-pool", &sync_thread_pool);
    param_list_configure_switch(pl, "-never-discard", &never_discard);
    param_list_configure_switch(pl, "-production", &las_production_mode);
}

static void declare_usage(cxx_param_list & pl)/*{{{*/
{
    param_list_usage_header(pl,
            "In the names and in the descriptions of the parameters, below there are often\n"
            "aliases corresponding to the convention that 0 is the rational side and 1\n"
            "is the algebraic side. If the two sides are algebraic, then the word\n"
            "'rational' just means the side number 0. Note also that for a rational\n"
            "side, the factor base is recomputed on the fly (or cached), and there is\n"
            "no need to provide a fb0 parameter.\n"
            );

    las_info::declare_usage(pl);
    las_parallel_desc::declare_usage(pl);
    las_todo_list::declare_usage(pl);
    las_output::declare_usage(pl);
    tdict::declare_usage(pl);

    param_list_decl_usage(pl, "trialdiv-first-side", "begin trial division on this side");
    param_list_decl_usage(pl, "allow-largesq", "allows large special-q, e.g. for a DL descent");

    param_list_decl_usage(pl, "sublat", "modulus for sublattice sieving");

    param_list_decl_usage(pl, "log-bucket-region", "set bucket region to 2^x");
    param_list_decl_usage(pl, "log-bucket-region-step", "set the number of level-(n-1) buckets inside a level-n bucket to 2^x");

    siever_config::declare_usage(pl);

    param_list_decl_usage(pl, "exit-early", "once a relation has been found, go to next special-q (value==1), or exit (value==2)");
    param_list_decl_usage(pl, "file-cofact", "provide file with strategies for the cofactorization step");
    param_list_decl_usage(pl, "prepend-relation-time", "prefix all relation produced with time offset since beginning of special-q processing");
    param_list_decl_usage(pl, "sync", "synchronize all threads at each special-q");
    param_list_decl_usage(pl, "sync-thread-pool", "synchronize the thread pool (implies -t 1 !!)");
    where_am_I::decl_usage(pl);
    if (dlp_descent) {
        param_list_decl_usage(pl, "recursive-descent", "descend primes recursively");
        param_list_decl_usage(pl, "grace-time-ratio", "Fraction of the estimated further descent time which should be spent processing the current special-q, to find a possibly better relation");

    }
    /* given that this option is dangerous, we used to enable it only for
     * las_descent
     */
    param_list_decl_usage(pl, "never-discard", "Disable the discarding process for special-q's. This is dangerous. See bug #15617");

    param_list_decl_usage(pl, "production", "Sort of an opposite to -v. Disable all diagnostics except the cheap or critical ones. See #21688 and #21825.");
    verbose_decl_usage(pl);
}/*}}}*/

struct round_me {
    slice_index_t initial = 1;
    slice_index_t increase = 1;
    template<typename T>
    static
    round_me from(T const &)
    {
        return { .initial = T::initial, .increase = T::increase };
    }
    /*
    round_me() = default;
    round_me(round_me const &) = default;
    round_me& operator=(round_me const &) = default;
    round_me(round_me &&) = default;
    round_me& operator=(round_me &&) = default;
    ~round_me() = default;
    */
    slice_index_t operator()(slice_index_t y) const {
        return std::max(initial, increase * iceildiv(y, increase));
    }
};

/* Our fetching of the siever_config fields is definitely wrong here. We
 * should only use logA, logI, and the siever thresholds.
 *
 * Note that we are sowewhat limited in our queries of the las_info type
 * here, as it has not been *completely* initialized yet. Well, almost,
 * but the shared structure cache has not loaded the factor base yet.
 */
static size_t expected_memory_usage_per_binding_zone(siever_config const & sc,/*{{{*/
        las_info const & las,
        int print)
{
    int const hush = print ? 0 : 3;
    cado_poly_srcptr cpoly = las.cpoly;

    /*
    int logImin = (1+sc.logA)/2;
    if (las.adjust_strategy == 2) {
        logImin -= ADJUST_STRATEGY2_MAX_SQUEEZE;
    sc.instantiate_thresholds(0);
    sc.instantiate_thresholds(1);
    */

    /* do the estimate based on the typical config stuff provided.
     * This is most often going to give a reasonable rough idea anyway.
     */
    size_t memory = 0;

    for(int side = 0 ; side < las.cpoly->nb_polys ; side++) {
        if (!sc.sides[side].lim) continue;
        double const p1 = sc.sides[side].lim;
        double const p0 = 2;
        /* in theory this should depend on the galois group and so on.
         * Here we're counting only with respect to a full symmetric
         * Galois group.
         * 
         * The average number of roots modulo primes is the average
         * number of fixed points of a permutation, and that is 1. If we
         * average over permutations with at least one fixed point, then
         * we have n! / (n! - D_n), and D_n/n! = \sum_{0\leq k\leq
         * n}(-1)^k/k! (which tends to 1/e, so the average number of
         * roots whenever there's at least one root tends to
         * 1/(1-1/e) = 1.5819767...)
         */
        int const d = cpoly->pols[side]->deg;
        double ideals_per_prime = 1;
        double fac=1;
        for(int k = 1 ; k <= d ; k++) {
            fac *= -k;
            ideals_per_prime += 1/fac;
        }
        ideals_per_prime = 1/(1-ideals_per_prime);
        size_t const nideals = nprimes_interval(p0, p1);
        size_t const nprimes = nideals / ideals_per_prime;
        size_t more = 0;
        /* we have nideals/ideals_per_prime prime numbers, totalling
         * nideals roots.
         * Per prime, we have:
         *      fbprime_t
         *      redc_invp_t
         *      double  (for the weight_cdf table).
         */
        more += sizeof(fbprime_t) * nprimes;
        more += sizeof(redc_invp_t) * nprimes;
        more += sizeof(double) * nprimes;
        /* Per root we have:
         *      fbroot_t
         */
        more += sizeof(fbroot_t) * nideals;

        verbose_fmt_print(0, 3 + hush,
                "# side {}, lim={}, {} fb primes"
                " (d={}, {} roots per p if G=S_d): {}\n",
                side,
                sc.sides[side].lim,
                nideals,
                d, ideals_per_prime,
                more >> 20);
        memory += more;

        /* Maybe also count the cost of storing a few slicings ? */
    }
    return memory;
}/*}}}*/
/* This does not count the footprint per binding zone */
static size_t expected_memory_usage_per_subjob(siever_config const & sc,/*{{{*/
        las_info const & las,
        int nthreads,
        int print)
{
    int const hush = print ? 0 : 3;
    bkmult_specifier const bkmult = las.get_bk_multiplier();

    /* FIXME: I think that this code misses the case of sublat. */

    /* sc.instantiate_thresholds() depends on sc.logI */
    std::vector<fb_factorbase::key_type> K;

    const bool do_resieve = sc.needs_resieving();
    K.reserve(sc.sides.size());
    for(int side = 0 ; side < (int) sc.sides.size() ; side++)
        K.emplace_back(sc.instantiate_thresholds(side));

    size_t memory = 0;
    size_t more;

#if 0 && defined(__linux__) && defined(HAVE_GLIBC) && defined(__x86_64__)
    /* count threads. Each costs 8M+4k for the stack, 64MB for the
     * private heap. However, this is only virtual address space.
     * Therefore it's not really clear how we should count it. Maybe "not
     * at all" is a good start, in fact.
     *
     * 64MB is actually transiently 128MB, then 64MB.
     */
    if (0) {
        verbose_fmt_print(0, 3 + hush, "# {} threads: {}\n",
                nthreads,
                (more = nthreads * 0x4801000) >> 20);
        memory += more;
    }
#endif

    /*
    // toplevel is computed by fb_factorbase::slicing::slicing, based on
    // thresholds in fbK
    int toplevel = -1;
    for(int side = 0 ; side < las.cpoly->nb_polys ; side++) {
        int m;
        for(m = 0 ; m < FB_MAX_PARTS && K[side].thresholds[m] < sc.sides[side].lim; ++m);
        toplevel = std::max(m, toplevel);
    }

    if (toplevel == 0) toplevel++;

    ASSERT_ALWAYS(1 <= toplevel && toplevel <= MAX_TOPLEVEL);
    */

    int const nba = NUMBER_OF_BAS_FOR_THREADS(nthreads);

    std::array<double, MAX_TOPLEVEL + 1> ms, ss;
    std::array<round_me, MAX_TOPLEVEL + 1> rs;
    std::array<double, MAX_TOPLEVEL> ml, sl;
    std::array<round_me, MAX_TOPLEVEL> rl;

#if 0
    double m1s, s1s;
    round_me round1s;

#if MAX_TOPLEVEL >= 2
    double m2s, s2s;
    double m1l, s1l;
    round_me round2s, round1l;
#endif
#if MAX_TOPLEVEL >= 3
    double m3s, s3s;
    double m2l, s2l;
    round_me round3s, round2l;
#endif
    static_assert(MAX_TOPLEVEL == 3);
#endif



    if (do_resieve) {
        using T1s = bucket_update_t<1, shorthint_t>;
        ss[1] = sizeof(T1s);
        ms[1] = bkmult.get<T1s>();
        rs[1] = round_me::from(bucket_slice_alloc_defaults<1, shorthint_t>());

#if MAX_TOPLEVEL >= 2
        using T2s = bucket_update_t<2, shorthint_t>;
        using T1l = bucket_update_t<1, longhint_t>;
        ss[2] = sizeof(T2s);
        ms[2] = bkmult.get<T2s>();
        rs[2] = round_me::from(bucket_slice_alloc_defaults<2, shorthint_t>());
        sl[1] = sizeof(T1l);
        ml[1] = bkmult.get<T1l>();
        rl[1] = round_me::from(bucket_slice_alloc_defaults<1, longhint_t>());
#endif
#if MAX_TOPLEVEL >= 3
        using T3s = bucket_update_t<3, shorthint_t>;
        using T2l = bucket_update_t<2, longhint_t>;
        ss[3] = sizeof(T3s);
        ms[3] = bkmult.get<T3s>();
        rs[3] = round_me::from(bucket_slice_alloc_defaults<3, shorthint_t>());
        sl[2] = sizeof(T2l);
        ml[2] = bkmult.get<T2l>();
        rl[2] = round_me::from(bucket_slice_alloc_defaults<2, longhint_t>());
#endif
        static_assert(MAX_TOPLEVEL == 3);
    } else {
        using T1s = bucket_update_t<1, emptyhint_t>;
        ss[1] = sizeof(T1s);
        ms[1] = bkmult.get<T1s>();
        rs[1] = round_me::from(bucket_slice_alloc_defaults<1, emptyhint_t>());
#if MAX_TOPLEVEL >= 2
        using T2s = bucket_update_t<2, emptyhint_t>;
        using T1l = bucket_update_t<1, logphint_t>;
        ss[2] = sizeof(T2s);
        ms[2] = bkmult.get<T2s>();
        rs[2] = round_me::from(bucket_slice_alloc_defaults<2, emptyhint_t>());
        sl[1] = sizeof(T1l);
        ml[1] = bkmult.get<T1l>();
        rl[1] = round_me::from(bucket_slice_alloc_defaults<1, logphint_t>());
#endif      
#if MAX_TOPLEVEL >= 3
        using T3s = bucket_update_t<3, emptyhint_t>;
        using T2l = bucket_update_t<2, logphint_t>;
        ss[3] = sizeof(T3s);
        ms[3] = bkmult.get<T3s>();
        rs[3] = round_me::from(bucket_slice_alloc_defaults<3, emptyhint_t>());
        sl[2] = sizeof(T2l);
        ml[2] = bkmult.get<T2l>();
        rl[2] = round_me::from(bucket_slice_alloc_defaults<2, logphint_t>());
#endif
    }

    // very large factor base primes, between bkthresh1 and lim = we
    // compute all bucket updates in one go. --> those updates are
    // bucket_update_t<2, shorthint_t>, that is, an XSIZE2 position
    // (should be 24 bits, is actually 32) and a short hint.
    //
    // For each big bucket region (level-2), we then transform these
    // updates into updates for the lower-level buckets. We thus
    // create bucket_update_t<1, longhint_t>'s with the downsort<>
    // function. The long hint is because we have the full fb_slice
    // index. the position in such a bucket update is shorter, only
    // XSIZE1.

    // For moderately large factor base primes (between bkthresh and
    // bkthresh1), we precompute the FK lattices (at some cost), and we
    // fill the buckets locally with short hints (and short positions:
    // bucket_update_t<1, shorthint_t>)


    for(int side = 0 ; side < las.cpoly->nb_polys ; side++) {
        if (!sc.sides[side].lim) continue;

        int toplevel = INT_MAX;

        /* Iterate through the factor base, looking for the level at
         * which a given set of primes is processed by a fill_in_buckets
         * call, as opposed to a downsort<> call.
         */
        for(int fib_level = MAX_TOPLEVEL ; fib_level >= 1 ; fib_level--) {
            if (K[side].thresholds[fib_level - 1] >= sc.sides[side].lim)
                continue;

            bool use_precomputed_lattices = false;

            if (toplevel == INT_MAX) {
                toplevel = fib_level;
            } else {
                use_precomputed_lattices = true;
            }

            double const p1 = K[side].thresholds[fib_level];
            double const p0 = K[side].thresholds[fib_level - 1];
            ASSERT_ALWAYS(p0 <= p1);

            size_t const nprimes = nprimes_interval(p0, p1);
            double const w = (std::log(std::log(p1)) - std::log(std::log(p0)));

            /* we duplicate code that is found in allocate_memory. TODO:
             * refactor that */
            size_t const nreg = (fib_level == toplevel) ?
                (1UL << (sc.logA - LOG_BUCKET_REGIONS[fib_level]))
                : (1 << (LOG_BUCKET_REGIONS[fib_level + 1] - LOG_BUCKET_REGIONS[fib_level]));
            size_t nup_per_reg = 0.25 * w * BUCKET_REGIONS[fib_level] / nba;
            /* assume LOG_BUCKET_REGIONS[fib_level] > logI */
            nup_per_reg *= 3;
            nup_per_reg += NB_DEVIATIONS_BUCKET_REGIONS * sqrt(nup_per_reg);
            size_t const nupdates = nup_per_reg * nreg * nba;
            {
                verbose_fmt_print(0, 3 + hush,
                        "# level {}, side {}:"
                        " {} primes, room for {} {}-updates [{}s] in {} arrays: {}\n",
                        fib_level, side, nprimes, nupdates,
                        fib_level, fib_level,
                        nba,
                        size_disp(more = ms[fib_level] * nupdates * ss[fib_level]));
                memory += more;
            }
            if (use_precomputed_lattices) {
                verbose_fmt_print(0, 3 + hush,
                        "# level {}, side {}: {} primes => precomp_plattices: {}\n",
                        fib_level,
                        side, nprimes,
                        size_disp(more = nprimes * sizeof(plattice_enumerator)));
            }
            {
                /* Count the slice_start pointers as well. We need to know
                 * how many slices will be processed in each bucket
                 * array. A rough rule of thumb probably works.
                 */
                size_t const nslices_estim = iceildiv(nprimes >> 16, nba);
                size_t const nslices_alloc = rs[fib_level](nslices_estim);
                std::string comment;
                size_t const waste = (nslices_alloc - nslices_estim) * nba * nreg * sizeof(void*);
                if (waste > (100<<20))
                    comment = fmt::format(" [note: using coarse-grain value of "
                            "{} slices instead; {}% waste ({} MB) !]",
                            nslices_alloc,
                            100.0*double_ratio(nslices_alloc-nslices_estim, nslices_alloc),
                            (waste>>20));
                verbose_fmt_print(0, 3 + hush,
                        "# level {}, side {}:"
                        " expect {} slices per array,"
                        " {} pointers each, in {} arrays:"
                        " {}{}\n",
                        fib_level,
                        side,
                        nslices_estim,
                        nreg,
                        nba,
                        size_disp(more = nba * nreg * MAX(nslices_estim, nslices_alloc) * sizeof(void*)),
                        comment);
                memory += more;
            }
            for(int level = fib_level - 1 ; level >= 1 ; level--) {
                {
                    // how many downsorted updates are alive at a given point in
                    // time ?
                    size_t const nupdates_D = nupdates >> 8;
                    verbose_fmt_print(0, 3 + hush,
                            "# level {}, side {}:"
                            " {} downsorted {}-updates [{}l]: {}\n",
                            level,
                            side, nupdates_D,
                            level, level,
                            size_disp(more = ml[level] * nupdates_D * sl[level]));
                    memory += more;
                }
                {
                    size_t const nslices_estim = 1;
                    size_t const nslices_alloc = rl[level](nslices_estim);
                    size_t const nreg = 1 << (LOG_BUCKET_REGIONS[level + 1] - LOG_BUCKET_REGIONS[level]);
                    std::string comment;
                    size_t const waste = (nslices_alloc - nslices_estim) * nba * nreg * sizeof(void*);
                    if (waste > (100<<20))
                        comment = fmt::format(" [note: using coarse-grain value of "
                                "{} slices instead; {}% waste ({} MB) !]",
                                nslices_alloc,
                                100.0*double_ratio(nslices_alloc-nslices_estim, nslices_alloc),
                                (waste>>20));
                    verbose_fmt_print(0, 3 + hush,
                            "# level {}, side {}:"
                            " expect {} slices per array,"
                            " {} pointers each, in {} arrays: {}{}\n",
                            level,
                            side,
                            nslices_estim,
                            nreg,
                            nba,
                            size_disp(more = nba * nreg * MAX(nslices_estim, nslices_alloc) * sizeof(void*)),
                            comment);
                    memory += more;
                }
            }
        }
    }

    return memory;
}/*}}}*/
static size_t expected_memory_usage_per_subjob_worst_logI(siever_config const & sc0, las_info const & las, int nthreads, int print)/*{{{*/
{
    /* We're not getting number_of_threads_per_subjob() from las, because
     * this function gets called before the las_parallel_desc layer of
     * the las_info structure is set.
     */
    int const hush = print ? 0 : 3;
    /* do the estimate based on the typical config stuff provided.
     * This is most often going to give a reasonable rough idea anyway.
     */
    siever_config sc = sc0;
    /* See for which I we'll have the most expensive setting */
    int logImin, logImax;
    if (sc0.adjust_strategy == 2) {
        logImin = (1+sc.logA)/2 - ADJUST_STRATEGY2_MAX_SQUEEZE;
        logImax = (1+sc.logA)/2 - ADJUST_STRATEGY2_MIN_SQUEEZE;
    } else {
        logImin = logImax = (1+sc.logA)/2;
    }
    size_t max_memory = 0;
    int logI_max_memory = 0;
    for(int logI = logImin ; logI <= logImax ; logI++) {
        sc.logI = logI;
        verbose_fmt_print(0, 3 + hush,
                "# Expected memory usage per subjob for logI={} [{} threads]:\n",
                sc.logI, nthreads);

        size_t const memory = expected_memory_usage_per_subjob(sc, las, nthreads, print);

        verbose_fmt_print(0, 2 + hush,
                "# Expected memory usage per subjob for logI={}: {}\n",
                sc.logI, size_disp(memory));

        if (memory > max_memory) {
            logI_max_memory = sc.logI;
            max_memory = memory;
        }
    }
    if (logImin != logImax || main_output->verbose < 2 + hush)
        verbose_fmt_print(0, 0 + hush,
                "# Expected memory use per subjob (max reached for logI={}):"
                " {}\n",
                logI_max_memory, size_disp(max_memory));
    return max_memory;
}/*}}}*/

static size_t expected_memory_usage(siever_config const & sc,/*{{{*/
        las_info & las,
        int print,
        size_t base_memory = 0)
{
    /* Contrary to the previous functions which are used very early on in
     * order to decide on the parallel setting, this function can safely
     * assume that this parallel setting is now complete.
     */
    int const hush = print ? 0 : 3;

    size_t const fb_memory = expected_memory_usage_per_binding_zone(sc, las, print);
    verbose_fmt_print(0, 0 + hush,
            "# Expected memory usage per binding zone for the factor base: {}\n",
            size_disp(fb_memory));

    size_t const subjob_memory = expected_memory_usage_per_subjob_worst_logI(sc, las, las.number_of_threads_per_subjob(), print);

    size_t memory;
    memory = subjob_memory;
    memory *= las.number_of_subjobs_per_memory_binding_zone();
    memory += fb_memory;
    memory *= las.number_of_memory_binding_zones();
    memory += base_memory;

    verbose_fmt_print(0, 0 + hush,
            "# Expected memory use for {} binding zone(s) and {} {}-threaded jobs per zone, counting {} MB of base footprint: {}\n",
            las.number_of_memory_binding_zones(),
            las.number_of_subjobs_per_memory_binding_zone(),
            las.number_of_threads_per_subjob(),     /* per subjob, always */
            base_memory >> 20,
            size_disp(memory));

    return memory;

}/*}}}*/

static void check_whether_q_above_large_prime_bound(siever_config const & conf, special_q const & doing)/*{{{*/
{
    /* Check whether q has a factor larger than the large prime bound.
     * This can create some problems, for instance in characters.
     * By default, this is not allowed, but the parameter
     * -allow-largesq is a by-pass to this test.
     */
    if (allow_largesq) return;

    if (doing.is_prime()) {
        /* Note: a prime special-q will have a 1-item prime_factors list,
         * but this item will be zero if the special-q does not fit in
         * 64 bits */
        if (mpz_sizeinbase(doing.p, 2) > conf.sides[doing.side].lpb) {
            fmt::print(stderr, "ERROR: The special q ({} bits) is larger than"
                    " the large prime bound on side {} ({} bits).\n",
                    (int) mpz_sizeinbase(doing.p, 2),
                    doing.side,
                    conf.sides[doing.side].lpb);
            fmt::print(stderr, "       You can disable this check with "
                    "the -allow-largesq argument,\n");
            fmt::print(stderr, "       It is for instance useful for the "
                    "descent.\n");
            fmt::print(stderr, "       Use tasks.sieve.allow_largesq=true.\n");
            exit(EXIT_FAILURE);
        }
    } else for (auto const & f: doing.prime_factors) {
        if ((unsigned int) nbits(f) > conf.sides[doing.side].lpb) {
            fmt::print(stderr, "ERROR: The special q ({} bits) has a factor {}"
                    " ({} bits) which is larger than"
                    " the large prime bound on side {} ({} bits).\n",
                    (int) mpz_sizeinbase(doing.p, 2),
                    f, nbits(f),
                    doing.side,
                    conf.sides[doing.side].lpb);
            fmt::print(stderr, "       You can disable this check with "
                    "the -allow-largesq argument,\n");
            fmt::print(stderr, "       It is for instance useful for the "
                    "descent.\n");
            fmt::print(stderr, "       Use tasks.sieve.allow_largesq=true.\n");
            exit(EXIT_FAILURE);
        }
    }
}
/*}}}*/

static void check_whether_special_q_is_root(cado_poly_srcptr cpoly, special_q const & doing)/*{{{*/
{
    cxx_mpz const & p(doing.p);
    cxx_mpz const & r(doing.r);
    ASSERT_ALWAYS(mpz_poly_is_root(cpoly->pols[doing.side], r, p));
}
/*}}}*/
static void per_special_q_banner(special_q const & doing)
{
    // arrange so that we don't have the same header line as the one
    // which prints the q-lattice basis
    verbose_fmt_print(0, 2, "#\n");
    verbose_fmt_print(0, 1, "# Now sieving {}\n", doing);
}

/* This is the core of the sieving routine. We do fill-in-buckets,
 * downsort, apply-buckets, lognorm computation, small sieve computation,
 * and survivor search and detection, all from here.
 */
static void do_one_special_q_sublat(nfs_work & ws, std::shared_ptr<nfs_work_cofac> const & wc_p, std::shared_ptr<nfs_aux> const & aux_p, thread_pool & pool)/*{{{*/
{
    int const nsides = ws.sides.size();

    nfs_aux & aux(*aux_p);
    timetree_t & timer_special_q(aux.rt.timer);
    where_am_I & w(aux.w);

    /* essentially update the fij polynomials and the max log bounds */
    if (main_output->verbose >= 2) {
        verbose_output_start_batch();
        for (int side = 0; side < nsides; ++side) {
            verbose_fmt_print (0, 1, "# f_{}'(x) = ", side);
            mpz_poly_fprintf(main_output->output, ws.sides[side].lognorms.fij);
        }
        verbose_output_end_batch();
    }

    where_am_I::begin_special_q(ws);

    /* TODO: is there a way to share this in sublat mode ? */

    std::vector<cado::multityped_array<precomp_plattice_t, 1, FB_MAX_PARTS - 1>> precomp_plattices(nsides);

    {
        {
            /* allocate_bucket_regions is probably ridiculously cheap in
             * comparison to allocate_buckets */
            // ws.allocate_bucket_regions();
            ws.allocate_buckets(*aux_p, pool);
        }

        for(int side = 0 ; side < nsides ; side++) {
            nfs_work::side_data  const& wss(ws.sides[side]);
            if (wss.no_fb()) continue;

            fill_in_buckets_toplevel_multiplex(ws, aux, pool, side, w);

            fill_in_buckets_prepare_plattices(ws, pool, side, precomp_plattices[side]);

        }

        /*
         * Mixing-and-matching threads here with the fill-in-buckets threads
         * might lead to unbalance.
         */
        BOOKKEEPING_TIMER(timer_special_q);

        for(int side = 0 ; side < nsides ; side++) {
            nfs_work::side_data  const& wss(ws.sides[side]);
            if (wss.no_fb()) continue;
            pool.add_task_lambda([&ws,aux_p,side](worker_thread * worker,int){
                    timetree_t & timer(aux_p->get_timer(worker));
                    ENTER_THREAD_TIMER(timer);
                    MARK_TIMER_FOR_SIDE(timer, side);

                    SIBLING_TIMER(timer, "prepare small sieve");

                    nfs_work::side_data & wss(ws.sides[side]);
                    // if (wss.no_fb()) return;

                    wss.ssd->small_sieve_init(
                            wss.fbs->small_sieve_entries.resieved,
                            wss.fbs->small_sieve_entries.rest,
                            ws.conf.logI,
                            side,
                            wss.fbK,
                            ws.Q,
                            wss.lognorms.scale);

                    wss.ssd->small_sieve_info("small sieve", side);

                    if (ws.toplevel == 1) {
                        /* when ws.toplevel > 1, this start_many call is done
                         * several times.
                         */
                        SIBLING_TIMER(timer, "small sieve start positions ");
                        wss.ssd->small_sieve_prepare_many_start_positions(
                                0,
                                std::min(SMALL_SIEVE_START_POSITIONS_MAX_ADVANCE, ws.nb_buckets[1]),
                                ws.conf.logI, ws.Q.sublat);
                        wss.ssd->small_sieve_activate_many_start_positions();
                    }
            },0);
        }

        /* Note: we haven't done any downsorting yet ! */

        pool.drain_queue(0);

        ws.check_buckets_max_full_toplevel(ws.toplevel);

        auto exc = pool.get_exceptions<buckets_are_full>(0);
        if (!exc.empty())
            throw *std::ranges::max_element(exc);
    }

    {
        CHILD_TIMER(timer_special_q, "process_bucket_region outer container");
        TIMER_CATEGORY(timer_special_q, sieving_mixed());
        if (ws.toplevel == 1) {
            /* Process bucket regions in parallel */
            process_many_bucket_regions(ws, wc_p, aux_p, pool, 0, w);
        } else {
            // Prepare plattices at internal levels

            // Visit the downsorting tree depth-first.
            // If toplevel = 1, then this is just processing all bucket
            // regions.
            size_t  const(&BRS)[FB_MAX_PARTS] = BUCKET_REGIONS;
            static_assert(MAX_TOPLEVEL == 3);
            for (int i = 0; i < ws.nb_buckets[ws.toplevel]; i++) {
                /* Dividing by BRS[1] is actually correct if we want to
                 * fill the first_region0_index parameter. Of course we
                 * must make sure that for the recursive downsort, this
                 * doesn't entail an extra multiplication (e.g. by
                 * BRS[2]/BRS[1]. XXX we must check this!
                 */
                switch (ws.toplevel) {
#if MAX_TOPLEVEL >= 2
                    case 2:
                        downsort_tree<1>(ws, wc_p, aux_p, pool,
                                i, i*BRS[2]/BRS[1],
                                precomp_plattices, w);
                        break;
#endif

#if MAX_TOPLEVEL >= 3
                    case 3:
                        downsort_tree<2>(ws, wc_p, aux_p, pool, i,
                                i*BRS[3]/BRS[1],
                                precomp_plattices, w);
                        break;
#endif
                    default:
                        ASSERT_ALWAYS(0);
                }
            }
        }
    }

    BOOKKEEPING_TIMER(timer_special_q);
    /* This ensures proper serialization of stuff that is in queue 0.
     * Maybe we could be looser about this.
     */
    if (sync_at_special_q) {
        pool.drain_all_queues();
    } else {
        pool.drain_queue(0);
    }
}/*}}}*/

/* This returns false if the special-q was discarded */
static bool
do_one_special_q(
        las_info & las,
        nfs_work & ws,
        std::shared_ptr<nfs_aux> const & aux_p,
        thread_pool & pool)/*{{{*/
{
    nfs_aux & aux(*aux_p);

    check_whether_special_q_is_root(las.cpoly, aux.doing);

    timetree_t& timer_special_q(aux.rt.timer);
    las_report& rep(aux.rt.rep);

    per_special_q_banner(aux.doing);

    SIBLING_TIMER(timer_special_q, "skew Gauss");
    TIMER_CATEGORY(timer_special_q, bookkeeping());

    /* By choosing the sieve area, we trigger the setup of most of the
     * per-special-q fields: In particular, this sets ws.conf, ws.Q, and
     * ws.J.
     *
     */
    if (!choose_sieve_area(las, aux_p, aux.doing, ws.conf, ws.Q, ws.J))
        return false;

    check_whether_q_above_large_prime_bound(ws.conf, aux.doing);

    BOOKKEEPING_TIMER(timer_special_q);

    ws.prepare_for_new_q(las, &aux.doing);

    /* the where_am_I structure is store in nfs_aux. We have a few
     * adjustments to make, and we want to make sure that the threads,
     * which also have their where_am_I object in nfs_aux, have it too.
     */
    WHERE_AM_I_UPDATE(aux.w, logI, ws.conf.logI);
    WHERE_AM_I_UPDATE(aux.w, pQ, &ws.Q);
    WHERE_AM_I_UPDATE(aux.w, sides, decltype(aux.w->sides){ws.sides.size()});
    for (size_t i = 0; i < ws.sides.size(); ++i) {
        WHERE_AM_I_UPDATE(aux.w, sides[i].fbs, ws.sides[i].fbs);
    }
    for(auto & t : aux.th) t.w = aux.w;


    /* Currently we assume that we're doing sieving + resieving on
     * both sides, or we're not. In the latter case, we expect to
     * complete the factoring work with batch cofactorization */
    ASSERT_ALWAYS(las.batch || las.batch_print_survivors.filename || ws.conf.needs_resieving());

    std::shared_ptr<nfs_work_cofac> wc_p;

    {
        wc_p = std::make_shared<nfs_work_cofac>(las, ws);

        rep.total_logI += ws.conf.logI;
        rep.total_J += ws.J;

        std::string extra;
        if (ws.task->depth())
            extra = fmt::format(" # within descent, currently at depth ", ws.task->depth());

        /* should stay consistent with DUPECHECK line printed in
         * sq_finds_relation() */
        verbose_fmt_print(0, 2,
                "# Sieving {}; I={}; J={};{}\n",
                ws.Q,
                1U << ws.conf.logI, ws.J, extra);

        if (!las.allow_composite_q() && !ws.Q.doing.is_prime()) {
            verbose_fmt_print(0, 1,
                    "# Warning, q={} is not prime\n",
                    ws.Q.doing.p);
        }
    }

    unsigned int sublat_bound = ws.Q.sublat.m;
    if (sublat_bound == 0)
        sublat_bound = 1;

    for (unsigned int i_cong = 0; i_cong < sublat_bound; ++i_cong) {
        for (unsigned int j_cong = 0; j_cong < sublat_bound; ++j_cong) {
            if (ws.Q.sublat.m) {
                if (i_cong == 0 && j_cong == 0)
                    continue;
                ws.Q.sublat.i0 = i_cong;
                ws.Q.sublat.j0 = j_cong;
                verbose_fmt_print(0, 1,
                        "# Sublattice (i,j) == ({}, {}) mod {}\n",
                        ws.Q.sublat.i0, ws.Q.sublat.j0, ws.Q.sublat.m);
            }
            do_one_special_q_sublat(ws, wc_p, aux_p, pool);
        }
    }

    /* It's better than before, but still. We're going to keep this data
     * around for longer than we think if we get an exception above. 
     */
    for(auto & wss : ws.sides) {
        wss.precomp_plattice_dense_clear();
        wss.ssd->small_sieve_clear();
    }

    return true;
}/*}}}*/

static void prepare_timer_layout_for_multithreaded_tasks(timetree_t & timer,
                                                         int nb_polys)
{
    /* This does nothing. We're just setting up the required
     * empty shell so that all the multithreaded tasks that
     * begin with MARK_TIMER_FOR_SIDE are properly registered
     * under the "sieving on side X" umbrella.
     */
#ifndef DISABLE_TIMINGS
    timetree_t::accounting_child const x(timer, tdict_slot_for_threads);
    for (int side = 0; side < nb_polys; ++side) {
        MARK_TIMER_FOR_SIDE(timer, side);
        TIMER_CATEGORY(timer, sieving(side));
    }
#else
    timer.nop();
#endif
}

static void print_survivors_job(las_info & las)
{
    las_info::batch_print_survivors_t & B(las.batch_print_survivors);
    std::unique_lock<std::mutex> foo(B.mm);
    for (;;) {
        decltype(las.survivors.L) Lloc;

        for( ; las.survivors.get_size() < B.filesize ; )
            B.cv.wait(foo);

        {
            /* we want to pull as much as we can from las.survivors, and
             * store it for local processing. We need the lock on
             * las.survivors to do this. We don't technically need to
             * keep the lock on B.mm, but it won't be of any use while we
             * work, so we might as well keep it.
             */
            size_t locsize = 0;
            const std::lock_guard<std::mutex> dummy(las.survivors.mm);
            for( ; !las.survivors.L.empty() && locsize <= B.filesize ; ) {
                /* See how many survivors I have in the front list, and
                 * pull only as much as we need. Possibly all, but not
                 * necessarily.
                 */
                auto & F = las.survivors.L.front();
                auto it = F.second.begin();
                for( ; locsize <= B.filesize && it != F.second.end() ; ++it)
                    ++locsize;
                decltype(F.second) M;
                M.splice(M.end(), F.second, F.second.begin(), it);
                Lloc.emplace_back(F.first, std::move(M));
                if (F.second.empty())
                    las.survivors.L.pop_front();
            }
            if (locsize == 0)
                return;

            /* if las.survivors.size is SIZE_MAX, we're in drain mode and
             * all printers are either waiting, or will check
             * las.survivors.get_size() before going to wait. In the
             * other (normal) case, it might make sense to wake up
             * another printer thread while we're busy printing here.
             */
            if (las.survivors.size != SIZE_MAX) {
                las.survivors.size -= locsize;
                if (las.survivors.size)
                    B.cv.notify_one();
            }
        }

        auto const f = fmt::format("{}.{}", B.filename, B.counter++);

        /* Now we temporarily unlock foo. */
        foo.unlock();
        auto const f_part = f + ".part";
        auto out = fopen_helper(f_part, "w");
        for(auto const & m : Lloc) {
            special_q const & doing(m.first);
            fmt::print(out.get(), "# q = ({}, {}, {})\n",
                    doing.p, doing.r, doing.side);
            for(auto const & s : m.second) {
                fmt::print(out.get(),
                        "{} {} {}\n", s.a, s.b, join(s.cofactor, " "));
            }
        }
        int const rc = rename(f_part.c_str(), f.c_str());
        WARN_ERRNO_DIAG(rc != 0, "rename(%s, %s)", f_part.c_str(), f.c_str());
        foo.lock();
    }
}


/* This is only used for batch, and also batch-print-survivors. We just
 * collected some survivors in ws.cofac_candidates. We're going to move
 * them to las.L, with the final outcome being either:
 *  - batch cofactoring at the end of las
 *  - printing, via one of the printing threads.
 *
 * So we're going to handle a pair (special q, list of (a,b)'s).
 *
 * ws.cofac_candidates is empty after this function.
 */
static void transfer_local_cofac_candidates_to_global(las_info & las, nfs_work & ws)
{
    if (ws.cofac_candidates.empty())
        return;

    las.survivors.append(ws.Q.doing, std::move(ws.cofac_candidates));

    if (!las.batch_print_survivors.filename)
        return;

    if (las.survivors.size >= las.batch_print_survivors.filesize) {
        /* We need to print. Wake a thread to do it */
        las.batch_print_survivors.cv.notify_one();
    }
}


static void las_subjob(las_info & las, int subjob, report_and_timer & global_rt)/*{{{*/
{
    where_am_I w MAYBE_UNUSED;
    WHERE_AM_I_UPDATE(w, plas, &las);

    las.set_subjob_binding(subjob);

    report_and_timer botched;

    int nwaste = 0;
    int nq = 0;
    double cumulated_wait_time = 0;     /* for this subjob only */
    {
        /* add scoping to control dtor call */
        /* queue 0: main
         * queue 1: ECM
         * queue 2: things that we join almost immediately, but are
         * multithreaded nevertheless: alloc buckets, ...
         */
        thread_pool pool(las.number_of_threads_per_subjob(), cumulated_wait_time, 3, sync_thread_pool);
        nfs_work ws(las);

        /* {{{ Doc on todo list handling
         * The function las_todo_feed behaves in different
         * ways depending on whether we're in q-range mode or in q-list mode.
         *
         * q-range mode: the special-q's to be handled are specified as a
         * range. Then, whenever the todo list almost runs out, it is
         * refilled if possible, up to the limit q1 (-q0 & -rho just gives a
         * special case of this).
         *
         * q-list mode: the special-q's to be handled are always read from a
         * file. Therefore each new special-q to be handled incurs a
         * _blocking_ read on the file, until EOF. This mode is also used for
         * the descent, which has the implication that the read occurs if and
         * only if the todo list is empty. }}} */

        for(;;) {
            if (global_exit_semaphore)
                break;

            main_output->fflush();

            auto * task = las.tree->pull();

            if (!task) break;

            /* maybe examine what we have here in the todo list, and
             * decide on the relevance of creating a new output object */

            /* (non-blocking) join results from detached cofac */
            for(task_result * r ; (r = pool.get_result(1, false)) ; delete r);

            nq++;

            /* We'll convert that to a shared_ptr later on, because this is
             * to be kept by the cofactoring tasks that will linger on quite
             * late.
             * Note that we must construct this object *outside* the try
             * block, because we want this list to be common to all attempts
             * for this q.
             */
            auto rel_hash_p = std::make_shared<nfs_aux::rel_hash_t>();

            for(;;) {
                /*
                 * The lifetime of the nfs_aux object *is* thread-safe,
                 * because the shared_ptr construct provides this
                 * guarantee. Helgrind, however, is not able to detect it
                 * properly, and sees stuff happening in the nfs_aux dtor
                 * as conflicting with what is done in the try{} block.
                 * This is a false positive. cado-nfs.supp has an
                 * explicit suppression for this.
                 */

                /* ready to start over if we encounter an exception */
                try {
                    /* The nfs_aux ctor below starts the special-q timer.
                     * However we must not give it a category right now,
                     * since it is an essential property ot the timer trees
                     * that the root of the trees must not have a nontrivial
                     * category */
                    auto aux_p = std::make_shared<nfs_aux>(las, *task, rel_hash_p, las.number_of_threads_per_subjob());
                    nfs_aux & aux(*aux_p);
                    las_report & rep(aux.rt.rep);
                    timetree_t & timer_special_q(aux.rt.timer);
                    /* in case we get an exception */
                    aux.dest_rt = &botched;

                    ACTIVATE_TIMER(timer_special_q);

                    prepare_timer_layout_for_multithreaded_tasks(timer_special_q, las.cpoly->nb_polys);

                    bool const done = do_one_special_q(las, ws, aux_p, pool);

                    if (!done) {
                        /* Then we don't even keep track of the time, it's
                         * totally insignificant.
                         */
                        rep.nr_sq_discarded++;
                        break;
                    }

                    /* At this point we no longer risk an exception,
                     * therefore it is safe to tinker with the todo list
                     */
                    las.tree->postprocess(task, las.config_pool.max_descent_attempts_allowed(), timer_special_q);

                    transfer_local_cofac_candidates_to_global(las, ws);

                    aux.complete = true;
                    aux.dest_rt = &global_rt;

                    if (exit_after_rel_found > 1 && rep.reports > 0)
                        break;

                    break;
                } catch (buckets_are_full const & e) {
                    nwaste++;
                    verbose_fmt_print (2, 1,
                            "# redoing {} because {} buckets are full\n"
                            "# {}\n",
                            task->sq(),
                            bkmult_specifier::printkey(e.key),
                            e.what());

                    /* reason on the bk_multiplier that we used when we
                     * did the allocation ! It is set by
                     * prepare_for_new_q, called from do_one_special_q.
                     */
                    double const old_value = ws.bk_multiplier.get(e.key);
                    auto ratio = double_ratio(e.reached_size, e.theoretical_max_size) * 1.05;
                    double new_value = old_value * ratio;
                    double las_value;
                    if (!las.grow_bk_multiplier(e.key, ratio, new_value, las_value)) {

                        verbose_fmt_print(0, 1, "# Global {} bucket multiplier has already grown to {:.3f}. Not updating, since this will cover {:.3f}*{}/{}*1.05={:.3f}\n",
                                bkmult_specifier::printkey(e.key),
                                las_value,
                                old_value,
                                e.reached_size,
                                e.theoretical_max_size,
                                new_value
                                );
                    } else {
                        verbose_fmt_print(0, 1, "# Updating {} bucket multiplier to {:.3f}*{}/{}*1.05, ={:.3f}\n",
                                bkmult_specifier::printkey(e.key),
                                old_value,
                                e.reached_size,
                                e.theoretical_max_size,
                                new_value
                                );
                    }
                    if (las.config_pool.default_config_ptr) {
                        expected_memory_usage(las.config_pool.base, las, true, base_memory);
                    }
                    continue;
                }
                break;
            }

        } // end of loop over special q ideals.

        /* we delete the "pool" and "ws" variables at this point. */
        /* The dtor for "pool" is a synchronization point */
    }

    {
        std::lock_guard<std::mutex> const lock(global_rt.mm);

        global_rt.rep.nr_sq_processed += nq;
        global_rt.rep.nwaste += nwaste;
        global_rt.rep.cumulated_wait_time += cumulated_wait_time;
        global_rt.rep.waste += botched.timer.total_counted_time();
    }

    verbose_fmt_print(0, 1, "# subjob {} done ({} special-q's), now waiting for other jobs\n", subjob, nq);
}/*}}}*/

static std::string relation_cache_subdir_name(std::vector<unsigned long> const & splits, std::vector<unsigned long> const & split_q)/*{{{*/
{
    std::string d;
    /* find the file */
    for(unsigned int i = 0 ; i + 1 < split_q.size() ; i++) {
        int l = 0;
        for(unsigned long s = 1 ; splits[i] > s ; s*=10, l++);
        d += fmt::format("/{:0{}}", split_q[i], l);
    }
    return d;
}/*}}}*/

static std::string relation_cache_find_filepath_inner(std::string const & d, unsigned long qq)/*{{{*/
{
    std::string filepath;
    DIR * dir = opendir(d.c_str());
    DIE_ERRNO_DIAG(dir == nullptr, "opendir(%s)", d.c_str());
    for(struct dirent * ent ; (ent = readdir(dir)) != nullptr ; ) {
        unsigned long q0, q1;
        if (sscanf(ent->d_name, "%lu-%lu", &q0, &q1) != 2) continue;
        if (qq < q0 || qq >= q1) continue;
        filepath = d + "/" + ent->d_name;
        break;
    }
    closedir(dir);

    return filepath;
}/*}}}*/

static std::string relation_cache_find_filepath(std::string const & cache_path, std::vector<unsigned long> const & splits, cxx_mpz q)/*{{{*/
{
    std::vector<std::string> searched;

    /* write q in the variable basis given by the splits */
    cxx_mpz oq = q;
    std::vector<unsigned long> split_q = splits;
    for(unsigned int i = splits.size() ; i-- ; ) {
        split_q[i] = mpz_fdiv_ui(q, splits[i]);
        mpz_fdiv_q_ui(q, q, splits[i]);
    }
    if (mpz_cmp_ui(q, 0) != 0) {
        fmt::print(stderr, "# q is too large for relation cache\n", oq);
        exit(EXIT_FAILURE);
    }

    std::string d = cache_path + relation_cache_subdir_name(splits, split_q);

    std::string filepath = relation_cache_find_filepath_inner(d, split_q.back());

    if (filepath.empty() && split_q.size() > 1) {
        searched.push_back(d);

        /* Try the previous directory, if qranges cross the
         * boundaries at powers of ten */
        split_q[split_q.size() - 2] -= 1;
        split_q[split_q.size() - 1] += splits[splits.size() - 1];
        d = cache_path + relation_cache_subdir_name(splits, split_q);
        filepath = relation_cache_find_filepath_inner(d, split_q.back());
    }

    if (filepath.empty()) {
        searched.push_back(d);
        fmt::print(stderr, "# no file found in relation cache for q={} (searched directories: {})\n", q, join(searched, " "));
        exit(EXIT_FAILURE);
    }

    return filepath;
}
/*}}}*/

static void quick_subjob_loop_using_cache(las_info & las)/*{{{*/
{
    std::vector<unsigned long> splits;

    try {
        /* recover the list of splits from the config file */
        json dirinfo;
        if (!(std::ifstream(las.relation_cache + "/dirinfo.json") >> dirinfo))
            throw std::exception();
        for(size_t i = 0 ; i < dirinfo["splits"].size() ; i++) {
            splits.push_back((long) dirinfo["splits"][i]);
        }
    } catch (std::exception const & e) {
        fmt::print(stderr, "# Cannot read relation cache, or dirinfo.json in relation cache\n");
        exit(EXIT_FAILURE);
    }

    /* the inner mechanism of the descent loop entails breaking early on
     * when a relation is found. This doesn't interact well with what
     * we're doing here.
     */
    ASSERT_ALWAYS(!dlp_descent);

    double ct0 = seconds();
    double wt0 = wct_seconds();
    unsigned long nreports = 0;
    int nq = 0;

    for(;; nq++) {
        main_output->fflush();
        auto * task = las.tree->pull();
        if (!task) break;

        auto const & doing(task->sq());

        nq++;

        siever_config conf;
        qlattice_basis Q;
        uint32_t J;

        check_whether_special_q_is_root(las.cpoly, doing);
        per_special_q_banner(doing);
        if (!choose_sieve_area(las, *task, conf, Q, J)) continue;
        check_whether_q_above_large_prime_bound(conf, doing);

        verbose_fmt_print(0, 2, "# Sieving {}; I={}; J={};\n",
                Q, 1U << conf.logI, J);

        std::string const filepath = relation_cache_find_filepath(las.relation_cache, splits, doing.p);

        std::ifstream rf(filepath);
        DIE_ERRNO_DIAG(!rf, "open(%s)", filepath.c_str());
        for(std::string line ; getline(rf, line) ; ) {
            if (line.empty()) continue;
            if (line[0] == '#') continue;
            std::istringstream is(line);
            relation rel;
            if (!(is >> rel)) {
                fmt::print(stderr, "# parse error in relation\n");
                exit(EXIT_FAILURE);
            }
            if (!sq_finds_relation(las, doing, conf, Q, J, rel))
                continue;

            std::string prefix;
            nreports++;

            if (las.suppress_duplicates) {
                if (relation_is_duplicate(rel, doing, las)) {
                    prefix = "# DUPE ";
                    nreports--;
                }
            }
            verbose_fmt_print(0, 1, "{}{}\n", prefix, rel);
        }
        
        verbose_fmt_print (0, 1, "# Time for {}: [not reported in relation-cache mode]\n", Q.doing);
    }

    ct0 = seconds() - ct0;
    wt0 = wct_seconds() - wt0;
    verbose_fmt_print (2, 1, "# Total {} reports [{:1.3g}s/r, {:1.1f}r/sq] in {:1.3g} elapsed s [{:.1f}% CPU]\n",
            nreports,
            nreports ? double_ratio(ct0, nreports) : -1,
            nq ? double_ratio(nreports, nq): -1,
            wt0,
            100.0 * ct0/wt0);

}/*}}}*/

// coverity[root_function]
int main (int argc0, char const * argv0[])/*{{{*/
{
    double t0, wct;
    int argc = argc0;
    char const **argv = argv0;

    setvbuf(stdout, nullptr, _IONBF, 0);
    setvbuf(stderr, nullptr, _IONBF, 0);

    cxx_param_list pl;
    cado_sighandlers_install();

    declare_usage(pl);
    configure_switches(pl);
    configure_aliases(pl);

    argv++, argc--;
    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        /* Could also be a file */
        FILE * f = fopen(argv[0], "r");
        if (f) {
            param_list_read_stream(pl, f, 0);
            fclose(f);
            argv++,argc--;
            continue;
        }
        fmt::print(stderr, "Unhandled parameter {}\n", argv[0]);
        param_list_print_usage(pl, argv0[0], stderr);
        exit(EXIT_FAILURE);
    }

    param_list_parse_int(pl, "trialdiv-first-side", &trialdiv_first_side);
    param_list_parse_int(pl, "exit-early", &exit_after_rel_found);
    if (dlp_descent)
        param_list_parse_double(pl, "grace-time-ratio", &general_grace_time_ratio);
    param_list_parse_int(pl, "log-bucket-region", &LOG_BUCKET_REGION);
    param_list_parse_int(pl, "log-bucket-region-step", &LOG_BUCKET_REGION_step);
    set_LOG_BUCKET_REGION();

    main_output = std::unique_ptr<las_output>(new las_output(pl));

    if (las_production_mode) {
        tdict::global_enable = 0;
    }

    las_info las(pl);    /* side effects: prints cmdline and flags */
#ifdef SAFE_BUCKET_ARRAYS
      verbose_fmt_print(0, 0, "# WARNING: SAFE_BUCKET_ARRAYS is on !\n");
#endif
#ifdef SAFE_BUCKETS_SINGLE
      verbose_fmt_print(0, 0, "# WARNING: SAFE_BUCKETS_SINGLE is on !\n");
#endif

    if (las.cpoly->nb_polys > 2) {
        fmt::print(stderr, "las is only working with poly files with 1 or 2 sides\n");
        return EXIT_FAILURE;
    }

    where_am_I::interpret_parameters(pl);

    base_memory = Memusage() << 10U;

    if (las.tree->todo.print_todo_list_flag) {
        /* printing the todo list takes only a very small amount of ram.
         * In all likelihood, nsubjobs will be total number of cores (or
         * the number of threads that were requested on command line)
         */
        las.set_parallel(pl, double_ratio(base_memory, 1U << 30U));
        las.tree->todo.print_todo_list(pl, las.number_of_threads_total());
        return EXIT_SUCCESS;

    }

    /* First have a guess at our memory usage in single-threaded mode,
     * and see how many threads must be together. This means more memory
     * per job than our estimate, and then, maybe more threads together
     * as a consequence. This will converge eventually (of course, if
     * there is an explicit number of threads that is supplied, the
     * result is immediate).
     *
     * The significant subtlety is that the las_parallel setting matters
     * more than marginally: all subjobs within a given binding will
     * supposedly share the factor base, so we must take that into
     * account, while the las_parallel ctor will only grok a per-subjob
     * estimate.
     */
    try {
        las.set_parallel(pl);
    } catch (las_parallel_desc::needs_job_ram & e) {
        if (las.config_pool.default_config_ptr) {
            verbose_fmt_print(0, 2, "# No --job-memory option given, relying on automatic memory estimate\n");
            siever_config const & sc0(las.config_pool.base);
            const size_t ram0 = expected_memory_usage_per_binding_zone(sc0, las, false);
            for(int z = 1, s = 1 << 30, n = 1, spin=0 ; ; spin++) {
                if (spin > 10) {
                    fmt::print(stderr, "Warning: computation of expected memory does not stabilize after {} attempts, picking the situation as it is\n", spin);
                    break;
                }
                const size_t ram1 = expected_memory_usage_per_subjob_worst_logI(sc0, las, n, false);
                const size_t jobram = (base_memory / z + ram0) / s + ram1;
                /*
                std::ostringstream os;
                os << z << " " << s << " " << n
                    << " " << (double) ram0 / (1 << 30)
                    << " " << (double) ram1 / (1 << 30)
                    << " " << (double) jobram / (1 << 30);
                fmt::print(stderr, "{}\n", os.str());
                */
                las.set_parallel(pl, double_ratio(jobram, 1U << 30U));
                const int nz = las.number_of_memory_binding_zones();
                const int ns = las.number_of_subjobs_per_memory_binding_zone();
                const int nn = las.number_of_threads_per_subjob();
                if (s == ns && n == nn && z == nz)
                    break;
                z = nz; s = ns; n = nn;
            }
        } else {
            verbose_fmt_print(0, 0, "# No --job-memory option given. Job placement needs either an explicit placement with -t, a complete siever config with a factor base to allow automatic ram estimates, or a --job-memory option\n");
            exit(EXIT_FAILURE);
        }
    } catch (las_parallel_desc::bad_specification & e) {
        verbose_fmt_print(0, 0, "# Error reported by the cpu binding layer: {}\n", e.what());
        verbose_fmt_print(0, 0, "# The parallelism specification for this job and/or the specifics of the hardware make it difficult for us to decide on what to do on an automatic basis with respect to CPU binding. Please stick to simple \"-t <number of threads>\". More advanced specifications like \"-t auto\" cannot be supported for this hardware.\n");
        exit(EXIT_FAILURE);
    }


    /* These are sometimes looked up a bit late in the process */
    sieve_shared_data::lookup_parameters(pl, las.cpoly->nb_polys);
    batch_side_config::lookup_parameters(pl, las.cpoly->nb_polys);

    param_list_warn_unused(pl);

    /* In the random-sample + relation cache case, we're going to proceed
     * through a special case, as this will spare us the need to load the
     * factor base.
     * The computation of the needed memory (which in itself is cheap)
     * also seems overkill, but we have a slight difficulty: if doing
     * duplicate checking, we need cofac strategies. And these are placed
     * in a cache, which is accessed depending on the memory binding.
     * Note that building strategies is done only when needed.
     * Nevertheless, this strategies stuff can easily become the dominant
     * factor in cached relations processing.
     */
    if (!las.relation_cache.empty()) {
        quick_subjob_loop_using_cache(las);
        return EXIT_SUCCESS;
    }

    ASSERT_ALWAYS(las.number_of_threads_total() > 0);

    las.display_binding_info();
    if (las.config_pool.default_config_ptr)
        expected_memory_usage(las.config_pool.base, las, true, base_memory);

    // already done below set_parallel
    // las.prepare_sieve_shared_data(pl);
    las.load_factor_base(pl);

    /* The global timer and global report structures are not active for
     * now.  Most of the accounting is done per special_q, and this will
     * be summarized once we're done with the multithreaded run.
     */
    report_and_timer global_rt;

    t0 = seconds ();
    wct = wct_seconds();

    try {
        if (las.batch_print_survivors.filename) {
            for(int i = 0 ; i < las.batch_print_survivors.number_of_printers ; i++) {
                las.batch_print_survivors.printer_threads.emplace_back(
                        print_survivors_job, std::ref(las));
            }
        }

        std::vector<std::thread> subjobs;
        /* In theory we would be able to to multiple descents in parallel, of
         * course, but how we should proceed with the todo list, our brace
         * mechanism, and the descent tree thing is altogether not obvious
         */
        const int nsubjobs = /* dlp_descent ? 1 : */ las.number_of_subjobs_total();
        subjobs.reserve(nsubjobs);
        for(int subjob = 0 ; subjob < nsubjobs ; ++subjob) {
            /* when references are passed through variadic template arguments
             * as for the std::thread ctor, we have automatic decaying unless
             * we use std::ref.
             */
            subjobs.emplace_back(
                    las_subjob,
                        std::ref(las),
                        subjob,
                        std::ref(global_rt)
                    );
        }
        for(auto & t : subjobs) t.join();

        las.tree->display_summary(0, 0);

        las.set_loose_binding();

        if (las.batch_print_survivors.filename) {
            las.survivors.mark_drain();
            las.batch_print_survivors.cv.notify_all();
            for(auto & x : las.batch_print_survivors.printer_threads)
                x.join();
        }

        if (las.batch)
          {
              int const nsides = las.cpoly->nb_polys;

              timetree_t batch_timer;
              auto z = call_dtor([&]() {
                      std::lock_guard<std::mutex> const lock(global_rt.mm);
                      global_rt.timer += batch_timer;
                      });
              ACTIVATE_TIMER(batch_timer);

            /* We need to access lim[01] and lpb[01] */
            siever_config const & sc0(las.config_pool.base);
            CHILD_TIMER(batch_timer, "batch cofactorization (time is wrong because of openmp)");
            TIMER_CATEGORY(batch_timer, batch_mixed());
            double extra_time = 0;

            std::vector<cxx_mpz> batchP(nsides);
            auto lpb = siever_side_config::collect_lpb(sc0.sides);
            auto batchlpb = batch_side_config::collect_batchlpb(las.bsides);
            auto batchmfb = batch_side_config::collect_batchmfb(las.bsides);
            auto batchfilename = batch_side_config::collect_batchfilename(las.bsides);
            for(int side = 0 ; side < nsides ; side++) {
                create_batch_file (batchfilename[side],
                        batchP[side],
                        sc0.sides[side].lim,
                        1UL << batchlpb[side],
                        las.cpoly->pols[side],
                        main_output->output,
                        las.number_of_threads_loose(),
                        extra_time);
            }

            double tcof_batch = seconds ();

            /* This one uses openmp, and forks from the current thread (well,
             * I believe so -- it's not entirely clear how openmp deals with
             * cpu binding. At least I presume that it does nothing before
             * the first pragma omp statement.)
             */
            find_smooth (las.survivors.L,
                    batchP, batchlpb, lpb, batchmfb,
                    main_output->output,
                    las.number_of_threads_loose(),
                    extra_time);

            /* We may go back to our general thread placement at this point.
             * Currently the code below still uses openmp */

            int ncurves = 0;
            for(int side = 0 ; side < nsides ; side++) {
                ncurves = MAX(ncurves, sc0.sides[side].ncurves);
                // Possible issue: if lpb=batchlp, ECM is still used for finding
                // the sieved primes in order to print the smooth relations.
                // In that case, we need enough curves to find them.
                if (sc0.sides[side].lpb == las.bsides[side].batchlpb)
                    ncurves = MAX(ncurves, 30);
            }


            if (ncurves <= 0)
                ncurves = 50; // use the same default as finishbatch

            std::list<std::pair<special_q, std::list<relation>>> rels;

            for(auto const & x : las.survivors.L) {
                rels.emplace_back(x.first, factor (x.second,
                    las.cpoly,
                    x.first,
                    batchlpb,
                    lpb,
                    ncurves,
                    main_output->output,
                    las.number_of_threads_loose(),
                    extra_time,
                    1));
            }
            verbose_fmt_print (0, 1, "# batch reported time for additional threads: {:.2f}\n", extra_time);
            batch_timer.add_foreign_time(extra_time);

            verbose_output_start_batch();
            nfs_aux::rel_hash_t rel_hash;
            size_t nondup = 0;
            for(auto const & rq : rels) {
                verbose_fmt_print(0, 1, "# {} relations for {}\n",
                        rq.second.size(),
                        rq.first);
                for(auto const & rel : rq.second) {
                    nfs_aux::abpair_t const ab(rel.a, rel.b);
                    bool const is_new_rel = rel_hash.insert(ab).second;
                    /* if !is_new_rel, it means that we had this (a,b)
                     * pair twice, probably because of a failed attempt,
                     * that was aborted because of an exception. (occurs
                     * only with 2-level sieving)
                     */
                    nondup += is_new_rel;
                    verbose_fmt_print(0, 1, "{}{}\n",
                            is_new_rel ? "" : "# DUP ", rel);
                }
            }
            verbose_output_end_batch();
            global_rt.rep.reports = nondup;

            tcof_batch = seconds () - tcof_batch;
          }
    } catch (std::exception const & e) {
        verbose_fmt_print(0, 0, "\n\n# Program aborted on fatal error\n{}\n", e.what());
        verbose_fmt_print(1, 0, "\n\n# Program aborted on fatal error\n{}\n", e.what());
        return EXIT_FAILURE;
    }

    t0 = seconds () - t0;
    wct = wct_seconds() - wct;

    if (las.config_pool.base.adjust_strategy < 2) {
        verbose_fmt_print (2, 1,
                "# Average J={:1.0f} for {} special-q's, max bucket fill -bkmult {}\n",
                double_ratio(global_rt.rep.total_J, global_rt.rep.nr_sq_processed),
                global_rt.rep.nr_sq_processed,
                las.get_bk_multiplier().print_all());
    } else {
        verbose_fmt_print (2, 1,
                "# Average logI={:1.1f} for {} special-q's, max bucket fill -bkmult {}\n",
                double_ratio(global_rt.rep.total_logI, global_rt.rep.nr_sq_processed),
                global_rt.rep.nr_sq_processed,
                las.get_bk_multiplier().print_all());
    }
    verbose_fmt_print (2, 1, "# Discarded {} special-q's out of {} pushed\n",
            global_rt.rep.nr_sq_discarded, las.tree->todo.created);

    auto D = global_rt.timer.filter_by_category();
    timetree_t::timer_data_type const tcpu = global_rt.timer.total_counted_time();

    if (tdict::global_enable >= 2) {
        verbose_fmt_print (0, 1, "#\n# Hierarchical timings:\n{}",
                global_rt.timer.display());

        verbose_fmt_print (0, 1, "#\n# Categorized timings (total counted time {:.2f}):\n", tcpu);
        for(auto const &c : D) {
            verbose_fmt_print (0, 1, "# {}: {:.2f}\n", 
                    coarse_las_timers::explain(c.first), c.second);
        }
        verbose_fmt_print (0, 1, "# total counted time: {:.2f}\n#\n", tcpu);
    }
    global_rt.rep.display_survivor_counters();


    if (main_output->verbose)
        facul_print_stats (main_output->output);

    /*{{{ Display tally */
    display_bucket_prime_stats();

    if (las_production_mode) {
        verbose_fmt_print (2, 1, "# Total cpu time {:1.2f}s [remove -production flag for timings]\n", t0);
    } else {
        verbose_fmt_print (2, 1, "# Wasted cpu time due to {} bkmult adjustments: {:1.2f}\n", global_rt.rep.nwaste, global_rt.rep.waste);
        verbose_fmt_print(0, 1, "# Cumulated wait time over all threads {:.2f}\n", global_rt.rep.cumulated_wait_time);
        verbose_fmt_print (2, 1, "# Total cpu time {:1.2f}s, useful {:1.2f}s [norm {:1.2f}+{:1.1f}, sieving {:1.1f}"
            " ({:1.1f}+{:1.1f} + {:1.1f}),"
            " factor {:1.1f} ({:1.1f}+{:1.1f} + {:1.1f}),"
            " rest {:1.1f}], wasted+waited {:1.2f}s, rest {:1.2f}s\n",
            t0,
            tcpu,
            D[coarse_las_timers::norms(0)],
            D[coarse_las_timers::norms(1)],

            D[coarse_las_timers::sieving(0)]+
            D[coarse_las_timers::sieving(1)]+
            D[coarse_las_timers::search_survivors()]+
            D[coarse_las_timers::sieving_mixed()],

            D[coarse_las_timers::sieving(0)],
            D[coarse_las_timers::sieving(1)],
            D[coarse_las_timers::search_survivors()]+
            D[coarse_las_timers::sieving_mixed()],

            D[coarse_las_timers::cofactoring(0)]+
            D[coarse_las_timers::cofactoring(1)]+
            D[coarse_las_timers::cofactoring_mixed()],

            D[coarse_las_timers::cofactoring(0)],
            D[coarse_las_timers::cofactoring(1)],
            D[coarse_las_timers::cofactoring_mixed()],

            D[coarse_las_timers::bookkeeping()],
            global_rt.rep.waste+global_rt.rep.cumulated_wait_time,
            t0-tcpu-global_rt.rep.waste-global_rt.rep.cumulated_wait_time
            );
    }

    verbose_fmt_print (2, 1, "# Total elapsed time {:1.2f}s, per special-q {:g}, per relation {:g}\n",
                 wct,
                 double_ratio(wct, global_rt.rep.nr_sq_processed),
                 double_ratio(wct, global_rt.rep.reports));

    /* memory usage */
    if (main_output->verbose >= 1 && las.config_pool.default_config_ptr) {
        expected_memory_usage(las.config_pool.base, las, true, base_memory);
    }
    const size_t peakmem = PeakMemusage();
    if (peakmem > 0)
        verbose_fmt_print (2, 1, "# PeakMemusage (MB) = {}\n",
                peakmem >> 10);
    if (las.suppress_duplicates) {
        verbose_fmt_print(2, 1, "# Total number of eliminated duplicates: {}\n", global_rt.rep.duplicates);
    }
    verbose_fmt_print (2, 1, "# Total {} reports [{:1.3g}s/r, {:1.1f}r/sq] in {:1.3g} elapsed s [{:.1f}% CPU]\n",
            global_rt.rep.reports,
            double_ratio(t0, global_rt.rep.reports),
            double_ratio(global_rt.rep.reports, global_rt.rep.nr_sq_processed),
            wct, 100*t0/wct);


    /*}}}*/

    las.cofac_stats.print();

    return EXIT_SUCCESS;
}/*}}}*/

