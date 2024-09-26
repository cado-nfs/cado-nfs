#include "cado.h" // IWYU pragma: keep

#include <cinttypes>                 // for PRId64, PRIu64
#include <cstdint>                    // for uint8_t
#include <cstdio>                     // for NULL
#include <mutex>                      // for lock_guard, mutex
#include <ostream>                    // for operator<<, ostringstream, basi...
#include <string>                     // for char_traits, basic_string
#include <utility>                    // for pair
#include <vector>                     // for vector
#include <cstdarg>             // IWYU pragma: keep

#include <gmp.h>                      // for mpz_srcptr, gmp_vfprintf

#include "cxx_mpz.hpp"
#include "las-detached-cofac.hpp"
#include "las-auxiliary-data.hpp"     // for nfs_aux, nfs_aux::rel_hash_t
#include "las-cofactor.hpp"           // for cofactorization_statistics
#include "las-duplicate.hpp"          // for relation_is_duplicate
#include "las-galois.hpp"             // for add_relations_with_galois
#include "las-globals.hpp"            // for prepend_relation_time, tt_qstart
#include "las-info.hpp"               // for las_info
#include "las-multiobj-globals.hpp"     // for dlp_descent
#include "las-output.hpp"             // for TRACE_CHANNEL
#include "las-report-stats.hpp"       // for las_report, TIMER_CATEGORY, las...
#include "las-threads-work-data.hpp"  // for nfs_work_cofac
#include "relation.hpp"               // for relation, operator<<
#include "tdict.hpp"                  // for timetree_t, SIBLING_TIMER
#include "timing.h"                 // for seconds
#include "utils_cxx.hpp"        // call_dtor
#include "verbose.h"

/* asynchronous cofactorization */

detached_cofac_result * detached_cofac_inner(worker_thread * worker, detached_cofac_parameters * param)
{
    /* Import some contextual stuff. */
    int id = worker->rank();
    nfs_work_cofac & wc(*param->wc_p);
    nfs_aux & aux(*param->aux_p);
    nfs_aux::thread_data & taux(aux.th[id]);
    las_info const & las(wc.las);
    las_report & rep(taux.rep);
    timetree_t & timer(aux.get_timer(worker));

    nfs_aux::rel_hash_t& rel_hash(aux.get_rel_hash());

    cofac_standalone & cur(*param);

    int nsides = las.cpoly->nb_polys;

    std::vector<int> cof_bitsize(nsides, 0);
    las.cofac_stats.call(cur.norm, cof_bitsize);

    SIBLING_TIMER(timer, "cofactoring"); // aka factor_both_leftover_norms
    TIMER_CATEGORY(timer, cofactoring_mixed());

    int pass = cur.factor_both_leftover_norms(wc);
    rep.survivors.cofactored += (pass != 0);

    auto res = new detached_cofac_result;

    if (cur.trace_on_spot() && pass == 0) {
        verbose_output_print(TRACE_CHANNEL, 0,
                "# factor_both_leftover_norm failed for (%" PRId64 ",%" PRIu64 "), ", cur.a, cur.b);
        verbose_output_vfprint(TRACE_CHANNEL, 0, gmp_vfprintf,
                "remains %Zd, %Zd unfactored\n",
                (mpz_srcptr) cur.norm[0],
                (mpz_srcptr) cur.norm[1]);
    }
    if (pass <= 0) {
        /* a factor was > 2^lpb, or some
           factorization was incomplete */
        return res;
    }

    rep.survivors.smooth++;

    /* yippee: we found a relation! */
    SIBLING_TIMER(timer, "print relations");
    TIMER_CATEGORY(timer, bookkeeping());

    las.cofac_stats.success(cof_bitsize);

    relation rel = cur.get_relation(aux.doing);

    if (cur.trace_on_spot()) {
        verbose_output_print(TRACE_CHANNEL, 0, "# Relation for (%"
                PRId64 ",%" PRIu64 ") printed\n", cur.a, cur.b);
    }

    {
        int do_check = las.suppress_duplicates;

        /* note that if we have large primes which don't fit in
         * an unsigned long, then the duplicate check will
         * quickly return "no".
         */

        nfs_aux::abpair_t ab(cur.a, cur.b);
        bool is_new_rel;
        {
            std::lock_guard<std::mutex> foo(rel_hash.mutex());
            is_new_rel = rel_hash.insert(ab).second;
        }

        const char * dup_comment = NULL;

        if (do_check && relation_is_duplicate(rel, wc.doing, wc.las)) {
            dup_comment = "# DUPE ";
            rep.duplicates ++;
        } else {
            if (!is_new_rel) {
                /* the relation was already printed in a prior attempt,
                 * that was aborted because of an exception. */
                dup_comment = "# DUP ";
            }
            /* Even if the relation was already printed, the las_report
             * object is (now) specific to the current attempt, and has
             * no memory of the number of reports for the failed
             * attempt. Therefore we need to count the relation as a
             * report no matter what.
             */
            rep.reports ++;
            /* Not clear what gives when we have Galois relations.  */
        }

        if (!dup_comment) dup_comment = "";

        std::ostringstream os;

        if (prepend_relation_time)
            os << "(" << seconds() - tt_qstart << ") ";

        // verbose_output_print(0, 3, "# i=%d, j=%u, lognorms = %hhu, %hhu\n", i, j, cur.S[0], cur.S[1]);

        os << dup_comment << rel << "\n";

        if(las.galois != NULL) {
            // adding relations on the fly in Galois cases
            // once filtering is ok for all Galois cases, 
            // this entire block would have to disappear
            add_relations_with_galois(las.galois, os, dup_comment,
                    &rep.reports, rel);
        }

        /* print all in one go */
        verbose_output_start_batch();     /* unlock I/O */
        verbose_output_print(0, 1, "%s", os.str().c_str());
        verbose_output_end_batch();     /* unlock I/O */
        if (dlp_descent)
            res->rel_p = std::make_shared<relation>(std::move(rel));
    }
    return res;
}

task_result * detached_cofac(worker_thread * worker, task_parameters * _param, int) /* {{{ */
{
    auto clean_param = call_dtor([_param]() { delete _param; });

    /* We must exit by cleaning the param structure we've been given. But
     * everything we do with the objects whose life is dependent on our
     * param structure must of course be completed at this point. This
     * holds as well for the timer. Yet ACTIVATE_TIMER below registers
     * some stuff to be done at dtor time, so it's important that we
     * clean up the parameters *after* the timer cleans up.
     */
    detached_cofac_parameters *param = static_cast<detached_cofac_parameters *>(_param);

    /* Import some contextual stuff. */
    int id = worker->rank();
    nfs_aux & aux(*param->aux_p);
    nfs_aux::thread_data & taux(aux.th[id]);
    las_report & rep(taux.rep);
    timetree_t & timer(aux.get_timer(worker));
    /* The timer is normally not running, as we're in a thread task.
     * However, in descent mode, this is called synchronously, and then
     * the situation is different since the timer has already been
     * activated above.
     */
    cofac_standalone & cur(*param);
    detached_cofac_result * res;
    if (dlp_descent) {
        CHILD_TIMER(timer, __func__);
        res = detached_cofac_inner(worker, param);
    } else {
        ENTER_THREAD_TIMER(timer);
        res = detached_cofac_inner(worker, param);
    }

    /* Build histogram of lucky S[x] values. Not sure it still works... */
    rep.mark_report(cur.S[0], cur.S[1]);

    return (task_result*) res;
}

/* }}} */


