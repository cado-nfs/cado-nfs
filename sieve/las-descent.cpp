#include "cado.h" // IWYU pragma: keep
#include <algorithm>              // for max
#include <cmath>                  // for isfinite, INFINITY
#include <memory>                 // for allocator_traits<>::value_type
#include <utility>                // for pair, make_pair
#include <vector>                 // for vector
#include <cstdio>              // IWYU pragma: keep
#include <cstdarg>             // IWYU pragma: keep
#include <gmp.h>                  // for mpz_srcptr, gmp_vfprintf, mpz_sizei...
#include "cxx_mpz.hpp"
#include "las-descent.hpp"
#include "las-descent-trees.hpp"  // for descent_tree::candidate_relation
#include "las-dlog-base.hpp"      // for las_dlog_base
#include "las-globals.hpp"        // for general_grace_time_ratio, recursive...
#include "las-report-stats.hpp"   // for TIMER_CATEGORY
#include "las-siever-config.hpp"  // for siever_config_pool, siever_config_p...
#include "las-todo-entry.hpp"     // for las_todo_entry
#include "las-todo-list.hpp"      // for las_todo_list
#include "relation.hpp"           // for relation::pr, relation
#include "typedefs.h"             // p_r_values_t
#include "verbose.h"

/* This returns true only if this descent node is now done, either based
 * on the new relation we have registered, or because the previous
 * relation is better anyway */
bool register_contending_relation(las_info const & las, las_todo_entry const & doing, relation & rel)/*{{{*/
{
    if (las.tree.must_avoid(rel)) {
        verbose_output_vfprint(0, 1, gmp_vfprintf, "# [descent] Warning: we have already used this relation, avoiding\n");
        return true;
    }

    /* compute rho for all primes, even on the rational side */
    rel.fixup_r(true);

    descent_tree::candidate_relation contender;
    contender.rel = rel;
    double time_left = 0;

    for(int side = 0 ; side < 2 ; side++) {
        for(unsigned int i = 0 ; i < rel.sides[side].size() ; i++) {
            relation::pr const & v(rel.sides[side][i]);
            if (mpz_cmp(doing.p, v.p) == 0)
                continue;
            p_r_values_t p = mpz_get_ui(v.p);
            if (v.p.fits<p_r_values_t>()) {
                p_r_values_t r = mpz_get_ui(v.r);
                if (las.dlog_base.is_known(side, p, r))
                    continue;
            }

            unsigned int n = mpz_sizeinbase(v.p, 2);
            siever_config_pool::key_type K(side, n);
            double e = las.config_pool.hint_expected_time(K);
            if (e < 0) {
                /* This is not worrysome per se. We just do
                 * not have the info in the descent hint table,
                 * period.
                 */
                verbose_output_vfprint(0, 1, gmp_vfprintf, "# [descent] Warning: cannot estimate refactoring time for relation involving %d@%d (%Zd,%Zd)\n", n, side, (mpz_srcptr) v.p, (mpz_srcptr) v.r);
                time_left = INFINITY;
            } else {
                if (std::isfinite(time_left))
                    time_left += e;
            }
            contender.outstanding.push_back(std::make_pair(side, v));
        }
    }
    verbose_output_print(0, 1, "# [descent] This relation entails an additional time of %.2f for the smoothing process (%zu children)\n",
            time_left, contender.outstanding.size());

    /* when we're re-examining this special-q because of a previous
     * failure, there's absolutely no reason to hurry up on a relation */
    contender.set_time_left(time_left, doing.iteration ? INFINITY : general_grace_time_ratio);

    return las.tree.new_candidate_relation(contender);
}/*}}}*/

void postprocess_specialq_descent(las_info & las, las_todo_list & todo, las_todo_entry const & doing, timetree_t & timer_special_q)/*{{{*/
{
    SIBLING_TIMER(timer_special_q, "descent");
    TIMER_CATEGORY(timer_special_q, bookkeeping());

    descent_tree::candidate_relation const & winner(las.tree.current_best_candidate());
    if (winner) {
        /* Even if not going for recursion, store this as being a
         * winning relation. This is useful for preparing the hint
         * file, and also for the initialization of the descent.
         */
        las.tree.take_decision();
        verbose_output_start_batch();
        FILE * output;
        for (size_t i = 0;
                (output = verbose_output_get(0, 0, i)) != NULL;
                i++) {
            winner.rel.print(output, "Taken: ");
        }
        verbose_output_end_batch();
        {
            unsigned int n = mpz_sizeinbase(doing.p, 2);
            verbose_output_start_batch();
            verbose_output_print (0, 1, "# taking path: ");
            for(int i = 0 ; i < doing.depth ; i++) {
                verbose_output_print (0, 1, " ");
            }
            verbose_output_print (0, 1, "%d@%d ->", n, doing.side);
            for(unsigned int i = 0 ; i < winner.outstanding.size() ; i++) {
                int side = winner.outstanding[i].first;
                relation::pr const & v(winner.outstanding[i].second);
                unsigned int n = mpz_sizeinbase(v.p, 2);
                verbose_output_print (0, 1, " %d@%d", n, side);
            }
            if (winner.outstanding.empty()) {
                verbose_output_print (0, 1, " done");
            }
            verbose_output_vfprint (0, 1, gmp_vfprintf, " \t%d %Zd %Zd\n", doing.side,
                    (mpz_srcptr) doing.p,
                    (mpz_srcptr) doing.r);
            verbose_output_end_batch();
        }
        if (recursive_descent) {
            /* reschedule the possibly still missing large primes in the
             * todo list */
            for(unsigned int i = 0 ; i < winner.outstanding.size() ; i++) {
                int side = winner.outstanding[i].first;
                relation::pr const & v(winner.outstanding[i].second);
                unsigned int n = mpz_sizeinbase(v.p, 2);
                verbose_output_vfprint(0, 1, gmp_vfprintf, "# [descent] " HILIGHT_START "pushing side-%d (%Zd,%Zd) [%d@%d]" HILIGHT_END " to todo list (now size %zu)\n", side, (mpz_srcptr) v.p, (mpz_srcptr) v.r, n, side, todo.size() + 1);
                todo.push_withdepth(v.p, v.r, side, doing.depth + 1);
            }
        }
    } else {
        las.tree.mark_try_again(doing.iteration + 1);
        unsigned int n = mpz_sizeinbase(doing.p, 2);
        verbose_output_print (0, 1, "# taking path: %d@%d -> loop (#%d)", n, doing.side, doing.iteration + 1);
        verbose_output_vfprint (0, 1, gmp_vfprintf, " \t%d %Zd %Zd\n", doing.side,
                (mpz_srcptr) doing.p,
                (mpz_srcptr) doing.r);
        verbose_output_vfprint(0, 1, gmp_vfprintf, "# [descent] Failed to find a relation for " HILIGHT_START "side-%d (%Zd,%Zd) [%d@%d]" HILIGHT_END " (iteration %d). Putting back to todo list.\n", doing.side,
                (mpz_srcptr) doing.p,
                (mpz_srcptr) doing.r, n, doing.side, doing.iteration);
        todo.push_withdepth(doing.p, doing.r, doing.side, doing.depth + 1, doing.iteration + 1);
    }
}/*}}}*/
