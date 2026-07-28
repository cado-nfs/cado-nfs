#include "cado.h" // IWYU pragma: keep

#include <ostream>
#include <string>
#include <utility>
#include <vector>
#include <sstream>

#include "las-detached-cofac.hpp"
#include "las-auxiliary-data.hpp"
#include "las-cofactor.hpp"
#include "las-cofac-standalone.hpp"
#include "las-duplicate.hpp"
#include "las-galois.hpp"
#include "las-globals.hpp"
#include "las-info.hpp"
#include "las-multiobj-globals.hpp"
#include "las-output.hpp"
#include "las-report-stats.hpp"
#include "las-threads-work-data.hpp"
#include "relation.hpp"
#include "tdict.hpp"
#include "timing.h"
#include "utils_cxx.hpp"
#include "threadpool.hpp"
#include "verbose.hpp"

/* asynchronous cofactorization */

static relation detached_cofac_inner(worker_thread * worker,
        nfs_work_cofac & wc,
        nfs_aux & aux,
        cofac_standalone && cur)
{
    /* Import some contextual stuff. */
    int const id = worker->rank();
    nfs_aux::thread_data & taux(aux.th[id]);
    las_info const & las(wc.las);
    las_report & rep(taux.rep);
    timetree_t & timer(aux.get_timer(worker));

    nfs_aux::rel_hash_t& rel_hash(aux.get_rel_hash());

    auto tt = timer.trace(id, chronograms::ECM());

    int const nsides = las.cpoly.nsides();

    std::vector<int> cof_bitsize(nsides, 0);
    las.cofac_stats.call(cur.norm, cof_bitsize);

    SIBLING_TIMER(timer, "cofactoring"); // aka factor_leftover_norms
    TIMER_CATEGORY(timer, cofactoring_mixed());

    int const pass = cur.factor_leftover_norms(wc);
    rep.survivors.cofactored += (pass != 0);

    if (cur.trace_on_spot() && pass == 0) {
        verbose_fmt_print(TRACE_CHANNEL, 0,
                "# factor_leftover_norm failed for ({},{}),"
                " remains {} unfactored\n",
                cur.a, cur.b,
                join(cur.norm, " "));
    }
    if (pass <= 0) {
        /* a factor was > 2^lpb, or some
           factorization was incomplete */
        return {};
    }

    rep.survivors.smooth++;

    /* yippee: we found a relation! */
    SIBLING_TIMER(timer, "print relations");
    TIMER_CATEGORY(timer, bookkeeping());

    las.cofac_stats.success(cof_bitsize);

    relation rel = cur.get_relation(aux.doing);

    if (cur.trace_on_spot()) {
        verbose_fmt_print(TRACE_CHANNEL, 0, "# Relation for ({}) printed\n",
                rel.ab());
    }

    {
        int const do_check = las.suppress_duplicates;

        /* note that if we have large primes which don't fit in
         * an unsigned long, then the duplicate check will
         * quickly return "no".
         */

        nfs_aux::abpair_t const ab(cur.a, cur.b);
        const bool is_new_rel = rel_hash.locked()->insert(ab).second;

        const char * dup_comment = nullptr;

        auto dupcheck = [&]() {
          auto tt = timer.trace(id, chronograms::DUPCHECK());
          return relation_is_duplicate(rel, wc.doing, wc.las);
        };

        if (do_check && dupcheck()) {
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

        // verbose_fmt_print(0, 3, "# i={}, j={}, lognorms = {}, {}\n", i, j, cur.S[0], cur.S[1]);

        os << dup_comment << rel << "\n";

        if (las.galois != nullptr) {
            // adding relations on the fly in Galois cases
            // once filtering is ok for all Galois cases, 
            // this entire block would have to disappear
            add_relations_with_galois(las.galois, os, dup_comment,
                    &rep.reports, rel);
        }

        /* print all in one go */
        verbose_fmt_print(0, 1, "{}", os.str());

        if (dlp_descent)
            return rel;
    }

    /* Build histogram of lucky S[x] values. Not sure it still works... */
    // rep.mark_report(cur.S);

    return {};
}

relation detached_cofac(worker_thread * worker,
        nfs_work_cofac & wc,
        nfs_aux & aux,
        cofac_standalone && cur) /* {{{ */
{
    /* We must exit by cleaning the param structure we've been given. But
     * everything we do with the objects whose life is dependent on our
     * param structure must of course be completed at this point. This
     * holds as well for the timer. Yet ACTIVATE_TIMER below registers
     * some stuff to be done at dtor time, so it's important that we
     * clean up the parameters *after* the timer cleans up.
     */

    /* Import some contextual stuff. */
    timetree_t & timer(aux.get_timer(worker));
    /* The timer is normally not running, as we're in a thread task.
     * However, in descent mode, this is called synchronously, and then
     * the situation is different since the timer has already been
     * activated above.
     */

    relation res;
    if (dlp_descent) {
        CHILD_TIMER(timer, __func__);
        res = detached_cofac_inner(worker, wc, aux, std::move(cur));
    } else {
        ENTER_THREAD_TIMER(timer);
        res = detached_cofac_inner(worker, wc, aux, std::move(cur));
    }

    return res;
}

/* }}} */


