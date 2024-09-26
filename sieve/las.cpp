/* the macro above is for #include <cmath> -- however it must happen
 * first, because it may well be that one of the intermediary headers
 * pull stuff that is dependent on this flag.
 */
#include "cado.h" // IWYU pragma: keep

/* This compilation units reacts to TRACK_CODE_PATH and uses macros
 * such as WHERE_AM_I_UPDATE.
 * This compilation unit _must_ produce different object files depending
 * on the value of TRACK_CODE_PATH.
 * The WHERE_AM_I_UPDATE macro itself is defined in las-where-am-i.hpp
 */

#include <cstdint>     /* AIX wants it first (it's a bug) */
#include <cinttypes>                      // for PRId64, PRIu64
#include <climits>                        // for ULONG_MAX
#include <cmath>                          // for log, pow, sqrt
#include <algorithm>                      // for min, max, max_element
#include <array>                          // for array, array<>::value_type
#include <condition_variable>             // for condition_variable
#include <cstdio>                         // for size_t, fprintf, stderr
#include <cstdlib>                        // for exit, EXIT_FAILURE, EXIT_SU...
#include <cstdarg>             // IWYU pragma: keep
#include <functional>                     // for ref
#include <iomanip>                        // for operator<<, setprecision
#include <list>                           // for list, _List_iterator
#include <map>                            // for map
#include <memory>                         // for allocator, shared_ptr, make...
#include <mutex>                          // for mutex, lock_guard, unique_lock
#include <ostream>                        // for operator<<, ostringstream
#include <istream>                        // for operator>>
#include <fstream>                        // for ifstream
#include <string>                         // for string, basic_string, opera...
#include <thread>                         // for thread
#include <type_traits>                    // for remove_reference<>::type
#include <utility>                        // for move, pair
#include <vector>                         // for vector<>::iterator, vector
#include <gmp.h>                          // for mpz_srcptr, gmp_vfprintf
#include "bucket.hpp"                     // for bucket_slice_alloc_defaults
#include "cado_poly.h"
#include "cxx_mpz.hpp"
#include "ecm/batch.hpp"                      // for cofac_list, cofac_candidate
#include "ecm/facul.hpp"                      // for facul_print_stats
#include "fb-types.h"                     // for fbprime_t, sublat_t, slice_...
#include "fb.hpp"                         // for fb_factorbase::key_type
#include "las-auxiliary-data.hpp"         // for report_and_timer, nfs_aux
#include "las-bkmult.hpp"                 // for buckets_are_full, bkmult_sp...
#include "las-choose-sieve-area.hpp"      // for choose_sieve_area, never_di...
#include "las-cofactor.hpp"               // for cofactorization_statistics
#include "las-config.h"                   // for LOG_BUCKET_REGIONS, BUCKET_...
#include "las-descent-trees.hpp"          // descent_tree
#include "las-descent.hpp"                // for postprocess_specialq_descent
#include "las-divide-primes.hpp"          // for display_bucket_prime_stats
#include "las-dlog-base.hpp"              // IWYU pragma: keep
#include "las-fill-in-buckets.hpp"        // for downsort_tree, fill_in_buck...
#include "las-globals.hpp"                // for main_output, base_memory
#include "las-info.hpp"                   // for las_info, las_info::batch_p...
#include "las-multiobj-globals.hpp"       // for dlp_descent
#include "las-norms.hpp"                  // for lognorm_smart, ADJUST_STRAT...
#include "las-output.hpp"                 // for las_output
#include "las-parallel.hpp"               // for las_parallel_desc, las_para...
#include "las-plattice.hpp"               // for plattice_enumerator
#include "las-process-bucket-region.hpp"  // for process_bucket_region_spawn
#include "las-qlattice.hpp"               // for qlattice_basis, operator<<
#include "las-report-stats.hpp"           // for las_report, coarse_las_timers
#include "las-siever-config.hpp"          // for siever_config::side_config
#include "cado-sighandlers.h"
#include "las-smallsieve.hpp"             // for small_sieve_activate_many_s...
#include "las-threads-work-data.hpp"      // for nfs_work, nfs_work::side_data
#include "las-todo-entry.hpp"             // for las_todo_entry
#include "las-todo-list.hpp"              // for las_todo_list
#include "las-where-am-i-proxy.hpp"            // for where_am_I
#include "las-where-am-i.hpp"             // for where_am_I, WHERE_AM_I_UPDATE
#include "lock_guarded_container.hpp"     // for lock_guarded_container
#include "macros.h"                       // for ASSERT_ALWAYS, MAX, iceildiv
#include "memusage.h"   // PeakMemusage
#include "misc.h"          // size_disp
#include "mpz_poly.h"
#include "multityped_array.hpp"           // for multityped_array
#include "params.h"
#include "relation.hpp"                   // for relation, operator<<
#include "tdict.hpp"                      // for timetree_t, slot, SIBLING_T...
#include "threadpool.hpp"                 // for thread_pool, thread_pool::s...
#include "timing.h"             // for seconds
#include "utils_cxx.hpp"
#include "verbose.h"
#include "json.hpp"
#include "ecm/facul_strategies_stats.hpp"
#include "fmt/format.h"
#include <sys/types.h>
#include <dirent.h>
#include "las-duplicate.hpp"




/*************************** main program ************************************/

static void configure_aliases(cxx_param_list & pl)
{
    las_info::configure_aliases(pl);
    param_list_configure_alias(pl, "log-bucket-region", "B");
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

    siever_config::declare_usage(pl);

    param_list_decl_usage(pl, "exit-early", "once a relation has been found, go to next special-q (value==1), or exit (value==2)");
    param_list_decl_usage(pl, "file-cofact", "provide file with strategies for the cofactorization step");
    param_list_decl_usage(pl, "prepend-relation-time", "prefix all relation produced with time offset since beginning of special-q processing");
    param_list_decl_usage(pl, "sync", "synchronize all threads at each special-q");
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
    int hush = print ? 0 : 3;
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

    for(int side = 0 ; side < 2 ; side++) {
        if (!sc.sides[side].lim) continue;
        double p1 = sc.sides[side].lim;
        double p0 = 2;
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
        int d = cpoly->pols[side]->deg;
        double ideals_per_prime = 1;
        double fac=1;
        for(int k = 1 ; k <= d ; k++) {
            fac *= -k;
            ideals_per_prime += 1/fac;
        }
        ideals_per_prime = 1/(1-ideals_per_prime);
        size_t nideals = nprimes_interval(p0, p1);
        size_t nprimes = nideals / ideals_per_prime;
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

        verbose_output_print(0, 3 + hush,
                "# side %d, lim=%lu, %zu fb primes"
                " (d=%d, %f roots per p if G=S_d): %zuMB\n",
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
    char buf[20];
    int hush = print ? 0 : 3;
    bkmult_specifier bkmult = las.get_bk_multiplier();

    /* FIXME: I think that this code misses the case of sublat. */

    /* sc.instantiate_thresholds() depends on sc.logI */
    fb_factorbase::key_type K[2] {
        sc.instantiate_thresholds(0),
        sc.instantiate_thresholds(1)
    };

    bool do_resieve = sc.sides[0].lim && sc.sides[1].lim;

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
        verbose_output_print(0, 3 + hush, "# %d threads: %zuMB\n",
                nthreads,
                (more = nthreads * 0x4801000) >> 20);
        memory += more;
    }
#endif

    // toplevel is computed by fb_factorbase::slicing::slicing, based on
    // thresholds in fbK
    int toplevel = -1;
    for(int side = 0 ; side < 2 ; side++) {
        int m;
        for(m = 0 ; m < FB_MAX_PARTS && K[side].thresholds[m] < sc.sides[side].lim; ++m);
        if (m > toplevel)
            toplevel = m;
    }

    if (toplevel == 0) toplevel++;

    ASSERT_ALWAYS(toplevel == 1 || toplevel == 2);

    /* the code path is radically different depending on toplevel. */

    int nba = NUMBER_OF_BAS_FOR_THREADS(nthreads);

    double m1s, m1l, m2s;
    size_t s1s, s1l, s2s;

    struct round_me {
        slice_index_t initial;
        slice_index_t increase;
        slice_index_t operator()(slice_index_t y) const {
            return std::max(initial, increase * iceildiv(y, increase));
        }
    };

    round_me round1s, round2s, round1l;

    if (do_resieve) {
        typedef bucket_update_t<1, shorthint_t> T1s;
        typedef bucket_update_t<2, shorthint_t> T2s;
        typedef bucket_update_t<1, longhint_t> T1l;
        typedef bucket_slice_alloc_defaults<1, shorthint_t> W1s;
        typedef bucket_slice_alloc_defaults<2, shorthint_t> W2s;
        typedef bucket_slice_alloc_defaults<1, longhint_t> W1l;
        s2s=sizeof(T2s); m2s=bkmult.get<T2s>();
        s1s=sizeof(T1s); m1s=bkmult.get<T1s>();
        s1l=sizeof(T1l); m1l=bkmult.get<T1l>();
        round1s = round_me { W1s::initial, W1s::increase };
        round2s = round_me { W2s::initial, W2s::increase };
        round1l = round_me { W1l::initial, W1l::increase };
    } else {
        typedef bucket_update_t<1, emptyhint_t> T1s;
        typedef bucket_update_t<2, emptyhint_t> T2s;
        typedef bucket_update_t<1, logphint_t> T1l;
        typedef bucket_slice_alloc_defaults<1, emptyhint_t> W1s;
        typedef bucket_slice_alloc_defaults<2, emptyhint_t> W2s;
        typedef bucket_slice_alloc_defaults<1, logphint_t> W1l;
        s2s=sizeof(T2s); m2s=bkmult.get<T2s>();
        s1s=sizeof(T1s); m1s=bkmult.get<T1s>();
        s1l=sizeof(T1l); m1l=bkmult.get<T1l>();
        round1s = round_me { W1s::initial, W1s::increase };
        round2s = round_me { W2s::initial, W2s::increase };
        round1l = round_me { W1l::initial, W1l::increase };
    }

    if (toplevel == 2) {
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

        for(int side = 0 ; side < 2 ; side++) {
            if (!sc.sides[side].lim) continue;
            /* In truth, I sort of know it isn't valid. We've built most
             * of the stuff on the idea that there's a global "toplevel"
             * notion, but that barely applies when one of the factor
             * bases happens to be much smaller than the other one */
            ASSERT_ALWAYS(K[side].thresholds[2] == sc.sides[side].lim);
            double p1 = K[side].thresholds[2];
            double p0 = K[side].thresholds[1];
            p0 = std::min(p1, p0);
            size_t nprimes = nprimes_interval(p0, p1);
            double w = (std::log(std::log(p1)) - std::log(std::log(p0)));

            /* we duplicate code that is found in allocate_memory. TODO:
             * refactor that */
            size_t nreg = 1UL << (sc.logA - LOG_BUCKET_REGIONS[2]);
            size_t nup_per_reg = 0.25 * w * BUCKET_REGIONS[2] / nba;
            /* assume LOG_BUCKET_REGIONS[2] > logI */
            nup_per_reg *= 3;
            nup_per_reg += NB_DEVIATIONS_BUCKET_REGIONS * sqrt(nup_per_reg);
            size_t nupdates = nup_per_reg * nreg * nba;
            {
                verbose_output_print(0, 3 + hush,
                        "# level 2, side %d:"
                        " %zu primes,"
                        " room for %zu 2-updates [2s] in %d arrays:"
                        " %s\n",
                        side, nprimes, nupdates,
                        nba,
                        size_disp(more = m2s * nupdates * s2s, buf));
                memory += more;
            }
            {
                /* Count the slice_start pointers as well. We need to know
                 * how many slices will be processed in each bucket
                 * array. A rough rule of thumb probably works.
                 */
                size_t nslices_estim = iceildiv(nprimes >> 16, nba);
                size_t nslices_alloc = round2s(nslices_estim);
                std::ostringstream os;
                size_t waste = (nslices_alloc - nslices_estim) * nba * nreg * sizeof(void*);
                if (waste > (100<<20))
                    os << " [note: using coarse-grain value of " << nslices_alloc << " slices instead; " << 100.0*(nslices_alloc-nslices_estim)/nslices_alloc << "% waste ("<<(waste>>20)<<" MB) !]";
                verbose_output_print(0, 3 + hush,
                        "# level 2, side %d:"
                        " expect %zu slices per array,"
                        " %zu pointers each, in %d arrays:"
                        " %s%s\n",
                        side,
                        nslices_estim,
                        nreg,
                        nba,
                        size_disp(more = nba * nreg * MAX(nslices_estim, nslices_alloc) * sizeof(void*), buf),
                        os.str().c_str());
                memory += more;
            }
            {
                // how many downsorted updates are alive at a given point in
                // time ?
                size_t nupdates_D = nupdates >> 8;
                verbose_output_print(0, 3 + hush,
                        "# level 1, side %d:"
                        " %zu downsorted 1-updates [1l]: %s\n",
                        side, nupdates_D,
                        size_disp(more = m1l * nupdates_D * s1l, buf));
                memory += more;
            }
            {
                size_t nslices_estim = 1;
                size_t nslices_alloc = round1l(nslices_estim);
                size_t nreg = 1 << (LOG_BUCKET_REGIONS[2] - LOG_BUCKET_REGIONS[1]);
                std::ostringstream os;
                size_t waste = (nslices_alloc - nslices_estim) * nba * nreg * sizeof(void*);
                if (waste > (100<<20))
                    os << " [note: using coarse-grain value of " << nslices_alloc << " slices instead; " << 100.0*(nslices_alloc-nslices_estim)/nslices_alloc << "% waste ("<<(waste>>20)<<" MB) !]";
                verbose_output_print(0, 3 + hush,
                        "# level 1, side %d:"
                        " expect %zu slices per array,"
                        " %zu pointers each, in %d arrays: %s%s\n",
                        side,
                        nslices_estim,
                        nreg,
                        nba,
                        size_disp(more = nba * nreg * MAX(nslices_estim, nslices_alloc) * sizeof(void*), buf),
                        os.str().c_str());
                memory += more;
            }
        }

        for(int side = 0 ; side < 2 ; side++) {
            if (!sc.sides[side].lim) continue;
            double p1 = K[side].thresholds[1];
            double p0 = K[side].thresholds[0];
            size_t nprimes = nprimes_interval(p0, p1);
            double w = (std::log(std::log(p1)) - std::log(std::log(p0)));

            size_t nreg = 1 << (LOG_BUCKET_REGIONS[2] - LOG_BUCKET_REGIONS[1]);
            size_t nup_per_reg = 0.25 * w * BUCKET_REGIONS[1] / nba;
            /* assume LOG_BUCKET_REGIONS[1] > logI -- if it's not the
             * case, the count will not be too wrong anyway. */
            nup_per_reg *= 3;
            nup_per_reg += NB_DEVIATIONS_BUCKET_REGIONS * sqrt(nup_per_reg);
            size_t nupdates = nup_per_reg * nreg * nba;
            verbose_output_print(0, 3 + hush,
                    "# level 1, side %d:"
                    " %zu primes,"
                    " %zu 1-updates [1s] in %d arrays:"
                    " %s\n",
                    side, nprimes, nupdates, nba,
                    size_disp(more = m1s * nupdates * s1s, buf));
            memory += more;
            verbose_output_print(0, 3 + hush,
                    "# level 1, side %d:"
                    " %zu primes => precomp_plattices: %s\n",
                    side, nprimes,
                    size_disp(more = nprimes * sizeof(plattice_enumerator), buf));
            memory += more;

            {
                /* Count the slice_start pointers as well. We need to know
                 * how many slices will be processed in each bucket
                 * array. A rough rule of thumb probably works.
                 */
                size_t nslices_estim = iceildiv(nprimes >> 16, nba);
                size_t nslices_alloc = round1s(nslices_estim);
                std::ostringstream os;
                size_t waste = (nslices_alloc - nslices_estim) * nba * nreg * sizeof(void*);
                if (waste > (100<<20))
                    os << " [note: using coarse-grain value of " << nslices_alloc << " slices instead; " << 100.0*(nslices_alloc-nslices_estim)/nslices_alloc << "% waste ("<<(waste>>20)<<" MB) !]";
                verbose_output_print(0, 3 + hush,
                        "# level 1, side %d:"
                        " expect %zu slices per array,"
                        " %zu pointers each, in %d arrays:"
                        " %s%s\n",
                        side,
                        nslices_estim,
                        nreg,
                        nba,
                        size_disp(more = nba * nreg * MAX(nslices_estim, nslices_alloc) * sizeof(void*), buf),
                        os.str().c_str());
                memory += more;
            }
        }
    } else if (toplevel == 1) {
        // *ALL* bucket updates are computed in one go as
        // bucket_update_t<1, shorthint_t>
        for(int side = 0 ; side < 2 ; side++) {
            if (!sc.sides[side].lim) continue;
            ASSERT_ALWAYS(K[side].thresholds[1] == sc.sides[side].lim);
            double p1 = K[side].thresholds[1];
            double p0 = K[side].thresholds[0];
            size_t nprimes = nprimes_interval(p0, p1);
            double w = (std::log(std::log(p1)) - std::log(std::log(p0)));

            /* we duplicate code that is found in allocate_memory. TODO:
             * refactor that */
            size_t nreg = 1UL << (sc.logA - LOG_BUCKET_REGIONS[1]);
            size_t nup_per_reg = 0.25 * w * BUCKET_REGIONS[1] / nba;
            /* assume LOG_BUCKET_REGIONS[1] > logI -- if it's not the
             * case, the count will not be too wrong anyway. */
            nup_per_reg *= 3;
            nup_per_reg += NB_DEVIATIONS_BUCKET_REGIONS * sqrt(nup_per_reg);
            size_t nupdates = nup_per_reg * nreg * nba;
            nupdates += NB_DEVIATIONS_BUCKET_REGIONS * sqrt(nupdates);
            verbose_output_print(0, 3 + hush,
                    "# level 1, side %d:"
                    " %zu primes,"
                    " %zu 1-updates [1s] in %d arrays:"
                    " %s\n",
                    side, nprimes, nupdates, nba,
                    size_disp(more = m1s * nupdates * s1s, buf));
            memory += more;
            {
                /* Count the slice_start pointers as well. We need to know
                 * how many slices will be processed in each bucket
                 * array. A rough rule of thumb probably works.
                 */
                size_t nslices_estim = iceildiv(nprimes >> 16, nba);
                size_t nslices_alloc = round1s(nslices_estim);
                std::ostringstream os;
                size_t waste = (nslices_alloc - nslices_estim) * nba * nreg * sizeof(void*);
                if (waste > (100<<20))
                    os << " [note: using coarse-grain value of " << nslices_alloc << " slices instead; " << 100.0*(nslices_alloc-nslices_estim)/nslices_alloc << "% waste ("<<(waste>>20)<<" MB) !]";
                verbose_output_print(0, 3 + hush,
                        "# level 1, side %d:"
                        " expect %zu slices per array,"
                        " %zu pointers each, in %d arrays:"
                        " %s%s\n",
                        side,
                        nslices_estim,
                        nreg,
                        nba,
                        size_disp(more = nba * nreg * MAX(nslices_estim, nslices_alloc) * sizeof(void*), buf),
                        os.str().c_str());
                memory += more;
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
    int hush = print ? 0 : 3;
    char buf[20];
    /* do the estimate based on the typical config stuff provided.
     * This is most often going to give a reasonable rough idea anyway.
     */
    siever_config sc = sc0;
    /* See for which I we'll have the most expensive setting */
    int logImin, logImax;
    if (las.adjust_strategy == 2) {
        logImin = (1+sc.logA)/2 - ADJUST_STRATEGY2_MAX_SQUEEZE;
        logImax = (1+sc.logA)/2 - ADJUST_STRATEGY2_MIN_SQUEEZE;
    } else {
        logImin = logImax = (1+sc.logA)/2;
    }
    size_t max_memory = 0;
    int logI_max_memory = 0;
    for(int logI = logImin ; logI <= logImax ; logI++) {
        sc.logI = logI;
        verbose_output_print(0, 3 + hush,
                "# Expected memory usage per subjob for logI=%d [%d threads]:\n",
                sc.logI, nthreads);

        size_t memory = expected_memory_usage_per_subjob(sc, las, nthreads, print);

        verbose_output_print(0, 2 + hush,
                "# Expected memory usage per subjob for logI=%d: %s\n",
                sc.logI, size_disp(memory, buf));

        if (memory > max_memory) {
            logI_max_memory = sc.logI;
            max_memory = memory;
        }
    }
    if (logImin != logImax || main_output.verbose < 2 + hush)
        verbose_output_print(0, 0 + hush,
                "# Expected memory use per subjob (max reached for logI=%d):"
                " %s\n",
                logI_max_memory, size_disp(max_memory, buf));
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
    int hush = print ? 0 : 3;
    char buf[20];

    size_t fb_memory = expected_memory_usage_per_binding_zone(sc, las, print);
    verbose_output_print(0, 0 + hush,
            "# Expected memory usage per binding zone for the factor base: %s\n",
            size_disp(fb_memory, buf));

    size_t subjob_memory = expected_memory_usage_per_subjob_worst_logI(sc, las, las.number_of_threads_per_subjob(), print);

    size_t memory;
    memory = subjob_memory;
    memory *= las.number_of_subjobs_per_memory_binding_zone();
    memory += fb_memory;
    memory *= las.number_of_memory_binding_zones();
    memory += base_memory;

    verbose_output_print(0, 0 + hush,
            "# Expected memory use for %d binding zone(s) and %d %d-threaded jobs per zone, counting %zu MB of base footprint: %s\n",
            las.number_of_memory_binding_zones(),
            las.number_of_subjobs_per_memory_binding_zone(),
            las.number_of_threads_per_subjob(),     /* per subjob, always */
            base_memory >> 20,
            size_disp(memory, buf));

    return memory;

}/*}}}*/

void check_whether_q_above_lare_prime_bound(siever_config const & conf, las_todo_entry const & doing)/*{{{*/
{
    /* Check whether q is larger than the large prime bound.
     * This can create some problems, for instance in characters.
     * By default, this is not allowed, but the parameter
     * -allow-largesq is a by-pass to this test.
     */
    if (allow_largesq) return;

    if (mpz_sizeinbase(doing.p, 2) > conf.sides[doing.side].lpb) {
        fprintf(stderr, "ERROR: The special q (%d bits) is larger than the "
                "large prime bound on side %d (%d bits).\n",
                (int) mpz_sizeinbase(doing.p, 2),
                doing.side,
                conf.sides[doing.side].lpb);
        fprintf(stderr, "       You can disable this check with "
                "the -allow-largesq argument,\n");
        fprintf(stderr, "       It is for instance useful for the "
                "descent.\n");
        fprintf(stderr, "       Use tasks.sieve.allow_largesq=true.\n");
        exit(EXIT_FAILURE);
    }
}
/*}}}*/

void check_whether_special_q_is_root(cado_poly_srcptr cpoly, las_todo_entry const & doing)/*{{{*/
{
    cxx_mpz const & p(doing.p);
    cxx_mpz const & r(doing.r);
    ASSERT_ALWAYS(mpz_poly_is_root(cpoly->pols[doing.side], r, p));
}
/*}}}*/
void per_special_q_banner(las_todo_entry const & doing)
{
    // arrange so that we don't have the same header line as the one
    // which prints the q-lattice basis
    verbose_output_print(0, 2, "#\n");
    std::ostringstream os;
    os << doing;
    verbose_output_print(0, 1, "# Now sieving %s\n", os.str().c_str());
}

/* This is the core of the sieving routine. We do fill-in-buckets,
 * downsort, apply-buckets, lognorm computation, small sieve computation,
 * and survivor search and detection, all from here.
 */
static void do_one_special_q_sublat(nfs_work & ws, std::shared_ptr<nfs_work_cofac> wc_p, std::shared_ptr<nfs_aux> aux_p, thread_pool & pool)/*{{{*/
{
    int nsides = ws.sides.size();

    nfs_aux & aux(*aux_p);
    timetree_t & timer_special_q(aux.rt.timer);
    where_am_I & w(aux.w);

    /* essentially update the fij polynomials and the max log bounds */
    if (main_output.verbose >= 2) {
        verbose_output_start_batch();
        verbose_output_print (0, 1, "# f_0'(x) = ");
        mpz_poly_fprintf(main_output.output, ws.sides[0].lognorms.fij);
        verbose_output_print (0, 1, "# f_1'(x) = ");
        mpz_poly_fprintf(main_output.output, ws.sides[1].lognorms.fij);
        verbose_output_end_batch();
    }

    where_am_I::begin_special_q(ws);

    /* TODO: is there a way to share this in sublat mode ? */

    multityped_array<precomp_plattice_t, 1, FB_MAX_PARTS> precomp_plattice(nsides);

    {
        {
            /* allocate_bucket_regions is probably ridiculously cheap in
             * comparison to allocate_buckets */
            // ws.allocate_bucket_regions();
            ws.allocate_buckets(*aux_p, pool);
        }

        for(int side = 0 ; side < 2 ; side++) {
            nfs_work::side_data & wss(ws.sides[side]);
            if (wss.no_fb()) continue;

            fill_in_buckets_toplevel(ws, aux, pool, side, w);

            fill_in_buckets_prepare_plattices(ws, pool, side, precomp_plattice);

        }

        /*
         * Mixing-and-matching threads here with the fill-in-buckets threads
         * might lead to unbalance.
         */
        BOOKKEEPING_TIMER(timer_special_q);

        for(int side = 0 ; side < 2 ; side++) {
            nfs_work::side_data & wss(ws.sides[side]);
            if (wss.no_fb()) continue;
            pool.add_task_lambda([&ws,aux_p,side](worker_thread * worker,int){
                    timetree_t & timer(aux_p->get_timer(worker));
                    ENTER_THREAD_TIMER(timer);
                    MARK_TIMER_FOR_SIDE(timer, side);

                    SIBLING_TIMER(timer, "prepare small sieve");

                    nfs_work::side_data & wss(ws.sides[side]);
                    // if (wss.no_fb()) return;

                    small_sieve_init(wss.ssd,
                            wss.fbs->small_sieve_entries.resieved,
                            wss.fbs->small_sieve_entries.rest,
                            ws.conf.logI,
                            side,
                            wss.fbK,
                            ws.Q,
                            wss.lognorms.scale);

                    small_sieve_info("small sieve", side, wss.ssd);

                    if (ws.toplevel == 1) {
                        /* when ws.toplevel > 1, this start_many call is done
                         * several times.
                         */
                        SIBLING_TIMER(timer, "small sieve start positions ");
                        small_sieve_prepare_many_start_positions(wss.ssd,
                                0,
                                std::min(SMALL_SIEVE_START_POSITIONS_MAX_ADVANCE, ws.nb_buckets[1]),
                                ws.conf.logI, ws.Q.sublat);
                        small_sieve_activate_many_start_positions(wss.ssd);
                    }
            },0);
        }

        /* Note: we haven't done any downsorting yet ! */

        pool.drain_queue(0);

        ws.check_buckets_max_full<shorthint_t>(ws.toplevel);
        ws.check_buckets_max_full<emptyhint_t>(ws.toplevel);
        auto exc = pool.get_exceptions<buckets_are_full>(0);
        if (!exc.empty())
            throw *std::max_element(exc.begin(), exc.end());
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
            size_t (&BRS)[FB_MAX_PARTS] = BUCKET_REGIONS;
            for (int i = 0; i < ws.nb_buckets[ws.toplevel]; i++) {
                switch (ws.toplevel) {
                    case 2:
                        downsort_tree<1>(ws, wc_p, aux_p, pool,
                                i, i*BRS[2]/BRS[1],
                                precomp_plattice, w);
                        break;
                    case 3:
                        downsort_tree<2>(ws, wc_p, aux_p, pool, i,
                                i*BRS[3]/BRS[1],
                                precomp_plattice, w);
                        break;
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
static bool do_one_special_q(las_info & las, nfs_work & ws, std::shared_ptr<nfs_aux> aux_p, thread_pool & pool)/*{{{*/
{
    nfs_aux & aux(*aux_p);
    ws.Q.doing = aux.doing;     /* will be set by choose_sieve_area anyway */

    check_whether_special_q_is_root(las.cpoly, aux.doing);

    timetree_t& timer_special_q(aux.rt.timer);
    las_report& rep(aux.rt.rep);

    per_special_q_banner(aux.doing);

    SIBLING_TIMER(timer_special_q, "skew Gauss");
    TIMER_CATEGORY(timer_special_q, bookkeeping());

    /* By choosing the sieve area, we trigger the setup of most of the
     * per-special-q fields.
     */
    if (!choose_sieve_area(las, aux_p, aux.doing, ws.conf, ws.Q, ws.J))
        return false;

    check_whether_q_above_lare_prime_bound(ws.conf, aux.doing);

    BOOKKEEPING_TIMER(timer_special_q);

    ws.prepare_for_new_q(las);

    /* the where_am_I structure is store in nfs_aux. We have a few
     * adjustments to make, and we want to make sure that the threads,
     * which also have their where_am_I object in nfs_aux, have it too.
     */
    WHERE_AM_I_UPDATE(aux.w, logI, ws.conf.logI);
    WHERE_AM_I_UPDATE(aux.w, pQ, &ws.Q);
    WHERE_AM_I_UPDATE(aux.w, sides[0].fbs, ws.sides[0].fbs);
    WHERE_AM_I_UPDATE(aux.w, sides[1].fbs, ws.sides[1].fbs);
    for(auto & t : aux.th) t.w = aux.w;


    /* Currently we assume that we're doing sieving + resieving on
     * both sides, or we're not. In the latter case, we expect to
     * complete the factoring work with batch cofactorization */
    ASSERT_ALWAYS(las.batch || las.batch_print_survivors.filename || (ws.conf.sides[0].lim && ws.conf.sides[1].lim));

    std::shared_ptr<nfs_work_cofac> wc_p;

    {
        wc_p = std::make_shared<nfs_work_cofac>(las, ws);

        rep.total_logI += ws.conf.logI;
        rep.total_J += ws.J;

        std::ostringstream extra;
        if (ws.Q.doing.depth)
            extra << " # within descent, currently at depth " << ws.Q.doing.depth;

        /* should stay consistent with DUPECHECK line printed in
         * sq_finds_relation() */
        std::ostringstream os;
        os << ws.Q;
        verbose_output_vfprint(0, 2, gmp_vfprintf,
                "# "
                "Sieving %s; I=%u; J=%u;%s\n",
                os.str().c_str(),
                1u << ws.conf.logI, ws.J, extra.str().c_str());

        if (!las.allow_composite_q && !ws.Q.doing.is_prime()) {
            verbose_output_vfprint(0, 1, gmp_vfprintf,
                    "# Warning, q=%Zd is not prime\n",
                    (mpz_srcptr) ws.Q.doing.p);
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
                verbose_output_print(0, 1, "# Sublattice (i,j) == (%u, %u) mod %u\n",
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
        small_sieve_clear(wss.ssd);
    }

    return true;
}/*}}}*/

static void prepare_timer_layout_for_multithreaded_tasks(timetree_t & timer)
{
    /* This does nothing. We're just setting up the required
     * empty shell so that all the multithreaded tasks that
     * begin with MARK_TIMER_FOR_SIDE are properly registered
     * under the "sieving on side X" umbrella.
     */
#ifndef DISABLE_TIMINGS
    timetree_t::accounting_child x(timer, tdict_slot_for_threads);
    for (int side = 0; side < 2; ++side) {
        MARK_TIMER_FOR_SIDE(timer, side);
        TIMER_CATEGORY(timer, sieving(side));
    }
#else
    timer.nop();
#endif
}

struct ps_params {
    std::shared_ptr<cofac_list> M;
    las_info & las;
    ps_params(std::shared_ptr<cofac_list> M, las_info & las) : M(M), las(las) {}
};

void las_info::batch_print_survivors_t::doit()
{
    for(std::unique_lock<std::mutex> foo(mm);!todo.empty() || !done;) {
        cv.wait(foo);

        /* This is both for spurious wakeups and for the finish condition */
        for( ; !todo.empty() ; ) {

            /* We have the lock held at this point */

            std::string f = std::string(filename) + "." + std::to_string(counter++);
            std::string f_part = f + ".part";

            cofac_list M = std::move(todo.front());

            todo.pop_front();

            /* Now we temporarily unlock foo. */
            foo.unlock();

            FILE * out = fopen(f_part.c_str(), "w");
            las_todo_entry const * curr_sq = NULL;
            for (auto const &s : M) {
                if (s.doing_p != curr_sq) {
                    curr_sq = s.doing_p;
                    gmp_fprintf(out,
                            "# q = (%Zd, %Zd, %d)\n",
                            (mpz_srcptr) s.doing_p->p,
                            (mpz_srcptr) s.doing_p->r,
                            s.doing_p->side);
                }
                gmp_fprintf(out,
                        "%" PRId64 " %" PRIu64 " %Zd %Zd\n", s.a, s.b,
                        (mpz_srcptr) s.cofactor[0],
                        (mpz_srcptr) s.cofactor[1]);
            }
            fclose(out);
            int rc = rename(f_part.c_str(), f.c_str());
            WARN_ERRNO_DIAG(rc != 0, "rename(%s, %s)", f_part.c_str(), f.c_str());

            foo.lock();
        }
    }
}

static void print_survivors_job(las_info & las)
{
    las.batch_print_survivors.doit();
}



static void las_subjob(las_info & las, int subjob, las_todo_list & todo, report_and_timer & global_rt)/*{{{*/
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
        thread_pool pool(las.number_of_threads_per_subjob(), cumulated_wait_time, 3);
        nfs_work workspaces(las);

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
            {
                std::lock_guard<std::mutex> foo(protect_global_exit_semaphore);
                if (global_exit_semaphore)
                    break;
            }
            main_output.fflush();
            las_todo_entry * doing_p = todo.feed_and_pop(las.rstate);
            if (!doing_p) break;
            las_todo_entry& doing(*doing_p);

            if (dlp_descent) {
                /* If the next special-q to try is a special marker, it means
                 * that we're done with a special-q we started before, including
                 * all its spawned sub-special-q's. Indeed, each time we start a
                 * special-q from the todo list, we replace it by a special
                 * marker. But newer special-q's may enver the todo list in turn
                 * (pushed with las_todo_push_withdepth).
                 */
                if (todo.is_closing_brace(doing)) {
                    las.tree.done_node();
                    if (las.tree.depth() == 0) {
                        if (recursive_descent) {
                            /* BEGIN TREE / END TREE are for the python script */
                            fprintf(main_output.output, "# BEGIN TREE\n");
                            las.tree.display_last_tree(main_output.output);
                            fprintf(main_output.output, "# END TREE\n");
                        }
                        las.tree.visited.clear();
                    }
                    continue;
                }

                /* pick a new entry from the stack, and do a few sanity checks */
                todo.push_closing_brace(doing.depth);

                las.tree.new_node(doing);
            }

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
                    auto aux_p = std::make_shared<nfs_aux>(las, doing, rel_hash_p, las.number_of_threads_per_subjob());
                    nfs_aux & aux(*aux_p);
                    las_report & rep(aux.rt.rep);
                    timetree_t & timer_special_q(aux.rt.timer);
                    /* in case we get an exception */
                    aux.dest_rt = &botched;

                    ACTIVATE_TIMER(timer_special_q);

                    prepare_timer_layout_for_multithreaded_tasks(timer_special_q);

                    bool done = do_one_special_q(las, workspaces, aux_p, pool);

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
                    if (dlp_descent)
                        postprocess_specialq_descent(las, todo, doing, timer_special_q);

                    if (!workspaces.cofac_candidates.empty()) {
                        {
                            std::lock_guard<std::mutex> foo(las.L.mutex());
#if 0
                            las.L.reserve(las.L.size() +
                                    workspaces.cofac_candidates.size());
                            /* do we fear that this move could be
                             * expensive ?  We're single-threaded at this
                             * point, so many threads could be waiting
                             * idly while we're moving pointers around.
                             * Worse, we have a mutex on las.L, so other
                             * subjobs might be waiting too. */
                            for(auto & x : workspaces.cofac_candidates)
                                las.L.emplace_back(std::move(x));
                            /* we can release the mutex now */
#else
                            las.L.splice(las.L.end(), workspaces.cofac_candidates);
                            if (las.batch_print_survivors.filename) {
                                while (las.L.size() >= las.batch_print_survivors.filesize) {
                                    // M is another list containing the
                                    // elements that are going to be printed
                                    // by another thread.
                                    cofac_list M;
                                    auto it = las.L.begin();
                                    for (uint64_t i = 0; i < las.batch_print_survivors.filesize; ++i) {
                                        ++it;
                                    }
                                    M.splice(M.end(), las.L, las.L.begin(), it);
                                    {
                                    std::lock_guard<std::mutex> dummy(las.batch_print_survivors.mm);
                                    las.batch_print_survivors.todo.push_back(std::move(M));
                                    }
                                    las.batch_print_survivors.cv.notify_one();

#if 0
                                    /* temporarily release the lock while
                                     * we're doing print_survivors, and
                                     * especially while we're waiting for
                                     * the previous print_survivors to
                                     * complete. */
                                    struct bar {
                                        std::mutex& m;
                                        bar(std::mutex & m) : m(m) { m.unlock(); }
                                        ~bar() { m.lock(); }
                                    };
                                    bar dummy2(las.L.mutex());
                                    print_survivors(M, las);
#endif
                                }
                            }
#endif
                        }
                        workspaces.cofac_candidates.clear();
                    }

                    aux.complete = true;
                    aux.dest_rt = &global_rt;

                    if (exit_after_rel_found > 1 && rep.reports > 0)
                        break;

                    break;
                } catch (buckets_are_full const & e) {
                    nwaste++;
                    verbose_output_vfprint (2, 1, gmp_vfprintf,
                            "# redoing q=%Zd, rho=%Zd because %s buckets are full\n"
                            "# %s\n",
                            (mpz_srcptr) doing.p, (mpz_srcptr) doing.r,
                            bkmult_specifier::printkey(e.key).c_str(),
                            e.what());

                    /* reason on the bk_multiplier that we used when we
                     * did the allocation ! It is set by
                     * prepare_for_new_q, called from do_one_special_q.
                     */
                    double old_value = workspaces.bk_multiplier.get(e.key);
                    double ratio = (double) e.reached_size / e.theoretical_max_size * 1.05;
                    double new_value = old_value * ratio;
                    double las_value;
                    if (!las.grow_bk_multiplier(e.key, ratio, new_value, las_value)) {

                        verbose_output_print(0, 1, "# Global %s bucket multiplier has already grown to %.3f. Not updating, since this will cover %.3f*%d/%d*1.05=%.3f\n",
                                bkmult_specifier::printkey(e.key).c_str(),
                                las_value,
                                old_value,
                                e.reached_size,
                                e.theoretical_max_size,
                                new_value
                                );
                    } else {
                        verbose_output_print(0, 1, "# Updating %s bucket multiplier to %.3f*%d/%d*1.05, =%.3f\n",
                                bkmult_specifier::printkey(e.key).c_str(),
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

        /* we delete the "pool" and "workspaces" variables at this point. */
        /* The dtor for "pool" is a synchronization point */
    }

    {
        std::lock_guard<std::mutex> lock(global_rt.mm);

        global_rt.rep.nr_sq_processed += nq;
        global_rt.rep.nwaste += nwaste;
        global_rt.rep.cumulated_wait_time += cumulated_wait_time;
        global_rt.rep.waste += botched.timer.total_counted_time();
    }

    verbose_output_print(0, 1, "# subjob %d done (%d special-q's), now waiting for other jobs\n", subjob, nq);
}/*}}}*/

static std::string relation_cache_subdir_name(std::vector<unsigned long> const & splits, std::vector<unsigned long> const & split_q)/*{{{*/
{
    std::string d;
    /* find the file */
    for(unsigned int i = 0 ; i + 1 < split_q.size() ; i++) {
        int l = 0;
        for(unsigned long s = 1 ; splits[i] > s ; s*=10, l++);
        d += fmt::format(FMT_STRING("/{:0{}}"), split_q[i], l);
    }
    return d;
}/*}}}*/

static std::string relation_cache_find_filepath_inner(std::string const & d, unsigned long qq)/*{{{*/
{
    std::string filepath;
    DIR * dir = opendir(d.c_str());
    DIE_ERRNO_DIAG(dir == NULL, "opendir(%s)", d.c_str());
    for(struct dirent * ent ; (ent = readdir(dir)) != NULL ; ) {
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
        gmp_fprintf(stderr, "# q is too large for relation cache\n",
                (mpz_srcptr) oq);
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
        std::ostringstream os;
        for(auto const & s : searched) os << " " << s;
        gmp_fprintf(stderr, "# no file found in relation cache for q=%Zd (searched directories:%s)\n", (mpz_srcptr) q, os.str().c_str());
        exit(EXIT_FAILURE);
    }

    return filepath;
}
/*}}}*/

static void quick_subjob_loop_using_cache(las_info & las, las_todo_list & todo)/*{{{*/
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
        fprintf(stderr, "# Cannot read relation cache, or dirinfo.json in relation cache\n");
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
        main_output.fflush();
        las_todo_entry * doing_p = todo.feed_and_pop(las.rstate);
        if (!doing_p) break;

        nq++;

        // auto rel_hash_p = std::make_shared<nfs_aux::rel_hash_t>();
        // nfs_aux aux(las, *doing_p, rel_hash_p, 1);
        struct { las_todo_entry doing; } aux;
        aux.doing = *doing_p;

        siever_config conf;
        qlattice_basis Q;
        uint32_t J;

        check_whether_special_q_is_root(las.cpoly, aux.doing);
        per_special_q_banner(aux.doing);
        if (!choose_sieve_area(las, aux.doing, conf, Q, J)) continue;
        check_whether_q_above_lare_prime_bound(conf, aux.doing);

        {
            std::ostringstream os;
            os << Q;
            verbose_output_vfprint(0, 2, gmp_vfprintf,
                    "# "
                    "Sieving %s; I=%u; J=%u;\n",
                    os.str().c_str(),
                    1u << conf.logI, J);
        }

        std::string filepath = relation_cache_find_filepath(las.relation_cache, splits, aux.doing.p);

        std::ifstream rf(filepath);
        DIE_ERRNO_DIAG(!rf, "open(%s)", filepath.c_str());
        for(std::string line ; getline(rf, line) ; ) {
            if (line.empty()) continue;
            if (line[0] == '#') continue;
            std::istringstream is(line);
            relation rel;
            if (!(is >> rel)) {
                gmp_fprintf(stderr, "# parse error in relation\n");
                exit(EXIT_FAILURE);
            }
            if (!sq_finds_relation(las, aux.doing, conf, Q, J, rel))
                continue;
            std::ostringstream os;

            nreports++;

            if (las.suppress_duplicates) {
                if (relation_is_duplicate(rel, aux.doing, las)) {
                    os << "# DUPE ";
                    nreports--;
                }
            }
            os << rel << "\n";

            verbose_output_start_batch();     /* unlock I/O */
            verbose_output_print(0, 1, "%s", os.str().c_str());
            verbose_output_end_batch();     /* unlock I/O */
        }
        
        {
            std::ostringstream os;
            os << Q.doing;
            verbose_output_print (0, 1, "# Time for %s: [not reported in relation-cache mode]\n", os.str().c_str());
        }
    }

    ct0 = seconds() - ct0;
    wt0 = wct_seconds() - wt0;
    verbose_output_print (2, 1, "# Total %lu reports [%1.3gs/r, %1.1fr/sq] in %1.3g elapsed s [%.1f%% CPU]\n",
            nreports,
            nreports ? ct0 / nreports : -1,
            (double) (nq ? nreports / nq: -1),
            wt0,
            100.0 * ct0/wt0);

}/*}}}*/

// coverity[root_function]
int main (int argc0, char *argv0[])/*{{{*/
{
    double t0, wct;
    int argc = argc0;
    char **argv = argv0;

    setbuf(stdout, NULL);
    setbuf(stderr, NULL);

    cxx_param_list pl;
    cado_sighandlers_install();

    declare_usage(pl);
    configure_switches(pl);
    configure_aliases(pl);

    argv++, argc--;
    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        /* Could also be a file */
        FILE * f;
        if ((f = fopen(argv[0], "r")) != NULL) {
            param_list_read_stream(pl, f, 0);
            fclose(f);
            argv++,argc--;
            continue;
        }
        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        param_list_print_usage(pl, argv0[0], stderr);
        exit(EXIT_FAILURE);
    }

    param_list_parse_int(pl, "trialdiv-first-side", &trialdiv_first_side);
    param_list_parse_int(pl, "exit-early", &exit_after_rel_found);
    if (dlp_descent)
        param_list_parse_double(pl, "grace-time-ratio", &general_grace_time_ratio);
    param_list_parse_int(pl, "log-bucket-region", &LOG_BUCKET_REGION);
    set_LOG_BUCKET_REGION();

    main_output.set(pl);

    if (las_production_mode) {
        tdict::global_enable = 0;
    }

    las_info las(pl);    /* side effects: prints cmdline and flags */
#ifdef SAFE_BUCKET_ARRAYS
      verbose_output_print(0, 0, "# WARNING: SAFE_BUCKET_ARRAYS is on !\n");
#endif
#ifdef SAFE_BUCKETS_SINGLE
      verbose_output_print(0, 0, "# WARNING: SAFE_BUCKETS_SINGLE is on !\n");
#endif

    las_todo_list todo(las.cpoly, pl);

    /* If qmin is not given, use lim on the special-q side by default.
     * This makes sense only if the relevant fields have been filled from
     * the command line.
     *
     * This is a kludge, really.
     */
    if (todo.sqside >= 0 && las.dupqmin[todo.sqside] == ULONG_MAX)
        las.dupqmin[todo.sqside] = las.config_pool.base.sides[todo.sqside].lim;

    where_am_I::interpret_parameters(pl);

    base_memory = Memusage() << 10;

    if (todo.print_todo_list_flag) {
        /* printing the todo list takes only a very small amount of ram.
         * In all likelihood, nsubjobs will be total number of cores (or
         * the number of threads that were requested on command line)
         */
        las.set_parallel(pl, base_memory / (double) (1 << 30));
        todo.print_todo_list(pl, las.rstate, las.number_of_threads_total());
        main_output.release();
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
            verbose_output_print(0, 2, "# No --job-memory option given, relying on automatic memory estimate\n");
            siever_config const & sc0(las.config_pool.base);
            size_t ram0 = expected_memory_usage_per_binding_zone(sc0, las, false);
            for(int z = 1, s = 1 << 30, n = 1, spin=0 ; ; spin++) {
                if (spin > 10) {
                    fprintf(stderr, "Warning: computation of expected memory does not stabilize after %d attempts, picking the situation as it is\n", spin);
                    break;
                }
                size_t ram1 = expected_memory_usage_per_subjob_worst_logI(sc0, las, n, false);
                size_t jobram = (base_memory / z + ram0) / s + ram1;
                /*
                std::ostringstream os;
                os << z << " " << s << " " << n
                    << " " << (double) ram0 / (1 << 30)
                    << " " << (double) ram1 / (1 << 30)
                    << " " << (double) jobram / (1 << 30);
                fprintf(stderr, "%s\n", os.str().c_str());
                */
                las.set_parallel(pl, (double) jobram / (1 << 30));
                int nz = las.number_of_memory_binding_zones();
                int ns = las.number_of_subjobs_per_memory_binding_zone();
                int nn = las.number_of_threads_per_subjob();
                if (s == ns && n == nn && z == nz)
                    break;
                z = nz; s = ns; n = nn;
            }
        } else {
            verbose_output_print(0, 0, "# No --job-memory option given. Job placement needs either an explicit placement with -t, a complete siever config with a factor base to allow automatic ram estimates, or a --job-memory option\n");
            exit(EXIT_FAILURE);
        }
    } catch (las_parallel_desc::bad_specification & e) {
        verbose_output_print(0, 0, "# Error reported by the cpu binding layer: %s\n", e.what());
        verbose_output_print(0, 0, "# The parallelism specification for this job and/or the specifics of the hardware make it difficult for us to decide on what to do on an automatic basis with respect to CPU binding. Please stick to simple \"-t <number of threads>\". More advanced specifications like \"-t auto\" cannot be supported for this hardware.\n");
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
        quick_subjob_loop_using_cache(las, todo);
        main_output.release();
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

    if (las.batch_print_survivors.filename) {
        for(int i = 0 ; i < las.batch_print_survivors.number_of_printers ; i++) {
            las.batch_print_survivors.printers.push_back(
                    std::thread(print_survivors_job, std::ref(las)));
        }
    }

    std::vector<std::thread> subjobs;
    /* In theory we would be able to to multiple descents in parallel, of
     * course, but how we should proceed with the todo list, our brace
     * mechanism, and the descent tree thing is altogether not obvious
     */
    int nsubjobs = dlp_descent ? 1 : las.number_of_subjobs_total();
    for(int subjob = 0 ; subjob < nsubjobs ; ++subjob) {
        /* when references are passed through variadic template arguments
         * as for the std::thread ctor, we have automatic decaying unless
         * we use std::ref.
         */
        subjobs.push_back(
                std::thread(las_subjob,
                    std::ref(las),
                    subjob,
                    std::ref(todo),
                    std::ref(global_rt)
                ));
    }
    for(auto & t : subjobs) t.join();

    if (dlp_descent && recursive_descent) {
        verbose_output_print(0, 1, "# Now displaying again the results of all descents\n");
        las.tree.display_all_trees(main_output.output);
    }

    las.set_loose_binding();

    if (las.batch_print_survivors.filename) {
        las.batch_print_survivors.mm.lock();
        las.batch_print_survivors.done = true;
        las.batch_print_survivors.todo.push_back(std::move(las.L));
        las.batch_print_survivors.mm.unlock();
        las.batch_print_survivors.cv.notify_all();
        for(auto & x : las.batch_print_survivors.printers)
            x.join();
    }

    if (las.batch)
      {
          int nsides = las.cpoly->nb_polys;

          timetree_t batch_timer;
          auto z = call_dtor([&]() {
                  std::lock_guard<std::mutex> lock(global_rt.mm);
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
                    main_output.output,
                    las.number_of_threads_loose(),
                    extra_time);
        }

	double tcof_batch = seconds ();

        /* This one uses openmp, and forks from the current thread (well,
         * I believe so -- it's not entirely clear how openmp deals with
         * cpu binding. At least I presume that it does nothing before
         * the first pragma omp statement.)
         */
	find_smooth (las.L,
                batchP, batchlpb, lpb, batchmfb,
                main_output.output,
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

        std::list<relation> rels = factor (las.L,
                las.cpoly,
                batchlpb,
                lpb,
                ncurves,
		main_output.output,
                las.number_of_threads_loose(),
                extra_time,
                1);
        verbose_output_print (0, 1, "# batch reported time for additional threads: %.2f\n", extra_time);
        batch_timer.add_foreign_time(extra_time);

        verbose_output_start_batch();
        nfs_aux::rel_hash_t rel_hash;
        size_t nondup = 0;
        for(auto const & rel : rels) {
            std::ostringstream os;
            nfs_aux::abpair_t ab(rel.a, rel.b);
            bool is_new_rel = rel_hash.insert(ab).second;
            if (!is_new_rel) {
                /* we had this (a,b) pair twice, probably because of a
                 * failed attempt, that was aborted because of an
                 * exception. (occurs only with 2-level sieving) */
                os << "# DUP ";
            } else {
                nondup++;
            }
            os << rel;
            verbose_output_print(0, 1, "%s\n", os.str().c_str());
        }
        verbose_output_end_batch();
        global_rt.rep.reports = nondup;

	tcof_batch = seconds () - tcof_batch;
      }

    t0 = seconds () - t0;
    wct = wct_seconds() - wct;

    if (las.adjust_strategy < 2) {
        verbose_output_print (2, 1, "# Average J=%1.0f for %lu special-q's, max bucket fill -bkmult %s\n",
                global_rt.rep.total_J / (double) global_rt.rep.nr_sq_processed, global_rt.rep.nr_sq_processed, las.get_bk_multiplier().print_all().c_str());
    } else {
        verbose_output_print (2, 1, "# Average logI=%1.1f for %lu special-q's, max bucket fill -bkmult %s\n",
                global_rt.rep.total_logI / (double) global_rt.rep.nr_sq_processed, global_rt.rep.nr_sq_processed, las.get_bk_multiplier().print_all().c_str());
    }
    verbose_output_print (2, 1, "# Discarded %lu special-q's out of %u pushed\n",
            global_rt.rep.nr_sq_discarded, todo.nq_pushed);

    auto D = global_rt.timer.filter_by_category();
    timetree_t::timer_data_type tcpu = global_rt.timer.total_counted_time();

    if (tdict::global_enable >= 2) {
        verbose_output_print (0, 1, "#\n# Hierarchical timings:\n%s", global_rt.timer.display().c_str());

        std::ostringstream os;
        os << std::fixed << std::setprecision(2) << tcpu;
        verbose_output_print (0, 1, "#\n# Categorized timings (total counted time %s):\n", os.str().c_str());
        for(auto const &c : D) {
            std::ostringstream xos;
            xos << std::fixed << std::setprecision(2) << c.second;
            verbose_output_print (0, 1, "# %s: %s\n", 
                    coarse_las_timers::explain(c.first).c_str(),
                    xos.str().c_str());
        }
        verbose_output_print (0, 1, "# total counted time: %s\n#\n", os.str().c_str());
    }
    global_rt.rep.display_survivor_counters();


    if (main_output.verbose)
        facul_print_stats (main_output.output);

    /*{{{ Display tally */
    display_bucket_prime_stats();

    if (las_production_mode) {
        verbose_output_print (2, 1, "# Total cpu time %1.2fs [remove -production flag for timings]\n", t0);
    } else {
        verbose_output_print (2, 1, "# Wasted cpu time due to %d bkmult adjustments: %1.2f\n", global_rt.rep.nwaste, global_rt.rep.waste);
        verbose_output_print(0, 1, "# Cumulated wait time over all threads %.2f\n", global_rt.rep.cumulated_wait_time);
        verbose_output_print (2, 1, "# Total cpu time %1.2fs, useful %1.2fs [norm %1.2f+%1.1f, sieving %1.1f"
            " (%1.1f+%1.1f + %1.1f),"
            " factor %1.1f (%1.1f+%1.1f + %1.1f),"
            " rest %1.1f], wasted+waited %1.2fs, rest %1.2fs\n",
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

    verbose_output_print (2, 1, "# Total elapsed time %1.2fs, per special-q %gs, per relation %gs\n",
                 wct, wct / (double) global_rt.rep.nr_sq_processed, wct / (double) global_rt.rep.reports);

    /* memory usage */
    if (main_output.verbose >= 1 && las.config_pool.default_config_ptr) {
        expected_memory_usage(las.config_pool.base, las, true, base_memory);
    }
    const long peakmem = PeakMemusage();
    if (peakmem > 0)
        verbose_output_print (2, 1, "# PeakMemusage (MB) = %ld \n",
                peakmem >> 10);
    if (las.suppress_duplicates) {
        verbose_output_print(2, 1, "# Total number of eliminated duplicates: %lu\n", global_rt.rep.duplicates);
    }
    verbose_output_print (2, 1, "# Total %lu reports [%1.3gs/r, %1.1fr/sq] in %1.3g elapsed s [%.1f%% CPU]\n",
            global_rt.rep.reports, t0 / (double) global_rt.rep.reports,
            (double) global_rt.rep.reports / (double) global_rt.rep.nr_sq_processed,
            wct,
            100*t0/wct);


    /*}}}*/

    las.cofac_stats.print();

    /* In essence, almost a dtor, but we want it to be before the pl dtor */
    main_output.release();

    return EXIT_SUCCESS;
}/*}}}*/

