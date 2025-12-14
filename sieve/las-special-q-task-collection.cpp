#include "cado.h" // IWYU pragma: keep

#include <list>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>

#include "fmt/base.h"
#include "fmt/format.h"
#include "fmt/std.h"    // IWYU pragma: keep    // fmt std::thread_id

#include "cado_poly.h"
#include "special-q.hpp"
#include "params.h"
#include "las-globals.hpp"
#include "las-report-stats.hpp"
#include "las-special-q-task.hpp"
#include "las-special-q-task-collection.hpp"
#include "las-special-q-task-tree.hpp"
#include "las-special-q-task-simple.hpp"
#include "portability.h"
#include "macros.h"
#include "timing.h"
#include "tdict.hpp"
#include "verbose.h"
#include "utils_cxx.hpp"
#include "las-multiobj-globals.hpp"

std::unique_ptr<special_q_task_collection_base> special_q_task_collection_base::create(cxx_cado_poly const & cpoly, cxx_param_list & pl)
{
    special_q_task_collection_base * t;
    if (dlp_descent)
        t = new special_q_task_collection_tree(cpoly, pl);
    else
        t = new special_q_task_collection_simple(cpoly, pl);
    return std::unique_ptr<special_q_task_collection_base>(t);
}

    /* TODO: display stats... */
#if 0
    /* (since we've dropped the dependence on siever_config, we no
     * longer do this)
     */
    /* We also use this to provide an updated stats file (whose
     * results are based on actual descents, and thus less
     * speculative) */
    
    list<string> new_hints;
    for(sit_t it = stats.begin() ; it != stats.end() ; it++) {
        /* Collect and print stats for this sq size */
        fprintf(o, "# Stats for %s\n", it->first().c_str());
        int n = it->second.size();
        int nok = 0;
        double t1 = 0;
        double t2 = 0;
        double tmin = 999999;
        double tmax = 0;
        int d1 = 0;
        int dmin = 0;
        int dmax = 0;
        int w1 = 0;
        int wmin = 0;
        int wmax = 0;
        using lit_t = list<collected_stats>::iterator;
        for(lit_t i = it->second.begin() ; i != it->second.end() ; i++) {
            nok += i->ok;
            t1 += i->t;
            t2 += i->t * i->t;
            d1 += i->d;
            w1 += i->w;
            tmin = min(tmin, i->t);
            tmax = max(tmax, i->t);
            dmin = min(dmin, i->d);
            dmax = max(dmax, i->d);
            wmin = min(wmin, i->w);
            wmax = max(wmax, i->w);
        }
        fprintf(o, "#   success %.2f (%d/%d)\n", (double) nok / n, nok, n);
        fprintf(o, "#   depth avg %.1f, min %d, max %d\n",
                (double) d1/n, dmin, dmax);
        fprintf(o, "#   weight avg %.1f, min %d, max %d\n",
                (double) w1/n, wmin, wmax);
        double tt1 = (double) t1/n;
        double tt2 = (double) t2/n;
        fprintf(o, "#   time avg %.3f, sdev %.3f, min %.3f, max %.3f\n",
                tt1, sqrt(tt2-tt1*tt1), tmin, tmax);

        ostringstream os;
        os << it->first()
            << " " << tt1
            << " " << (double) nok / n 
            << " I=" << it->first.sc->logI;
        for(int i = 0 ; i < 2 ; i++) {
            os << " " << it->first.sc->sides[i]->lim
                << "," << it->first.sc->sides[i]->lpb
                << "," << it->first.sc->sides[i]->mfb
                << "," << it->first.sc->sides[i]->lambda;
        }
        new_hints.push_back(os.str());
    }
    fprintf(o, "# The following data is an _example_ which can be used to provide a hint file\n");
    for(list<string>::iterator i = new_hints.begin() ; i != new_hints.end() ; i++) {
        fprintf(o, "# %s\n", i->c_str());
    }
#endif

void special_q_task_tree::update_child_status(special_q_task_tree * item, status_code before, status_code after)
{
    auto previous = item->status;
    ASSERT_ALWAYS(previous == before);
    item->status = after;
    children_by_status[after].insert(item);
    auto removed = children_by_status[before].erase(item);
    if (!removed)
        throw std::runtime_error(fmt::format(
                    "Cannot find \"{}\" tree node for {}"
                    " among the children of parent {}",
                    previous,
                    item->sq(),
                    item->parent->sq()));
}

/* This creates a tree for this special-q, optionally below a given
 * parent node. It's created as a "pending" node, meaning that:
 *  - it's among the children_pending of its parent
 *  - the timer "spent" is not started yet.
 *  - a pointer to this tree is added to the all_pending list.
 *
 * parent=nullptr means that we create a top-level entry (which, in
 * effect, means that it's a child of the forest element, but that's
 * an implementation detail).
 */
void special_q_task_collection_tree::new_node_unlocked(special_q const & doing, special_q_task_tree * parent)
{
    if (!parent)
        parent = &forest;

    std::unique_ptr<special_q_task_tree> kid(new special_q_task_tree(doing, parent));
    all_pending.push_back(kid.get());
    parent->children_by_status[special_q_task_tree::PENDING].insert(kid.get());
    parent->children.emplace_back(std::move(kid));
    created++;

    /* when this moves from children_pending to children_inprogress,
     * then we can do kid->spent = -seconds()
     */
}

void special_q_task_collection_tree::display_summary(int channel, int verbose)
{
    verbose_fmt_print(channel, verbose, "# BEGIN SUMMARY TREE"
            " [npending={}, created={} pulled={} done={} abandoned={}]\n",
            all_pending.size(), created, pulled, done, abandoned);
    verbose_fmt_print(channel, verbose, "{}",
            special_q_task_tree::prefixed { .t=&forest, .prefix="# " });
    verbose_fmt_print(channel, verbose, "# END SUMMARY TREE\n");
    verbose_fmt_print(channel, verbose, "# The status of the root node is {}\n", forest.status);
}


void special_q_task_collection_tree::done_node_unlocked(special_q_task_tree * item)
{
    constexpr auto PENDING = special_q_task::status_code::PENDING;
    constexpr auto IN_PROGRESS = special_q_task::status_code::IN_PROGRESS;
    constexpr auto IN_RECURSION = special_q_task::status_code::IN_RECURSION;
    constexpr auto DONE = special_q_task::status_code::DONE;

    ASSERT_ALWAYS(item->children_by_status[PENDING].empty());
    ASSERT_ALWAYS(item->children_by_status[IN_PROGRESS].empty());
    ASSERT_ALWAYS(item->children_by_status[IN_RECURSION].empty());

    if (item != &forest) {
        ASSERT_ALWAYS(item->status == IN_PROGRESS || item->status == IN_RECURSION);
        item->update_status(item->status, DONE);
        done++;
    }

    item->spent += seconds();

    if (item->parent == &forest) {
        work_to_do.notify_all();
        verbose_fmt_print(0, 0, "# BEGIN TREE\n");
        verbose_fmt_print(0, 0, "{}",
                special_q_task_tree::prefixed { .t=item, .prefix="# " });
        verbose_fmt_print(0, 0, "# END TREE\n");
        // XXX we used to do las.tree.visited.clear() at this point,
        // because the top level special-qs were sequence points anyway
        // (we had only one subjob). It wouldn't work as is, but still,
        // we might want to do something. See also the question in
        // abandon_node.
    } else if (item == &forest) {
        display_summary(0, 2);
    }

    /* We need to mark the parent as DONE if it has reached this status */
    if (item->parent
            && item->parent->children_by_status[PENDING].empty()
            && item->parent->children_by_status[IN_PROGRESS].empty()
            && item->parent->children_by_status[IN_RECURSION].empty()
            )
        done_node_unlocked(item->parent);
}

void special_q_task_collection_tree::recurse_node_unlocked(special_q_task_tree * item)
{
    constexpr auto IN_PROGRESS = special_q_task::status_code::IN_PROGRESS;
    constexpr auto IN_RECURSION = special_q_task::status_code::IN_RECURSION;

    item->update_status(IN_PROGRESS, IN_RECURSION);
}

special_q_task_tree * special_q_task_collection_tree::pull_internal()
{
    for( ; !all_pending.empty() ; ) {
        special_q_task_tree * item = all_pending.front();
        all_pending.pop_front();
        pulled++;
        if (item->status != special_q_task_tree::ABANDONED)
            return item;
        verbose_fmt_print(0, 0,
                "# **NOT** sieving {} (parent already abandoned)\n",
                item->sq());
        /* this->abandoned has already been increased */
    }
    auto q = todo.feed_and_pop();
    if (!q) {
        verbose_fmt_print(0, 3,
                "# {} gets NULL [npending={}, created={} pulled={} done={} abandoned={}]\n",
                std::this_thread::get_id(),
                all_pending.size(),
                created, pulled, done, abandoned);
        return nullptr;
    } else {
        pulled++;
        new_node_unlocked(q, nullptr);
        verbose_fmt_print(0, 3,
                "# {} gets fresh {} [created={} pulled={} done={} abandoned={}]\n",
                std::this_thread::get_id(), q,
                created, pulled, done, abandoned);
        special_q_task_tree * item = all_pending.front();
        all_pending.pop_front();
        ASSERT_ALWAYS(item->sq() == q);
        ASSERT_ALWAYS(item->status == special_q_task_tree::PENDING);
        return item;
    }
}

special_q_task * special_q_task_collection_tree::pull()
{
    std::unique_lock<std::mutex> lock(tree_lock);

    special_q_task_tree * item;
    /* In the case of a tree structure, we must be prepared for the case
     * where more special q's enter the todo list later.
     */
    for( ; (item = pull_internal()) == nullptr && (done + abandoned) < created ; ) {
        work_to_do.wait(lock);
    }

    if (!item)
        return nullptr;

    constexpr auto PENDING = special_q_task::status_code::PENDING;
    constexpr auto IN_PROGRESS = special_q_task::status_code::IN_PROGRESS;

    item->update_status(PENDING, IN_PROGRESS);
    item->spent = -seconds();

    return item;
}

void special_q_task_collection_tree::take_decision(special_q_task_tree * item)
{
    const std::lock_guard<std::mutex> lock(tree_lock);

    verbose_fmt_print(0, 0, "# Taking decision on {} {}\n",
            item->shortname(), item->sq());

    constexpr auto PENDING = special_q_task::status_code::PENDING;
    constexpr auto IN_PROGRESS = special_q_task::status_code::IN_PROGRESS;
    constexpr auto ABANDONED = special_q_task::status_code::ABANDONED;

    if (item->status == ABANDONED) {
        verbose_fmt_print(0, 0, "# {}: relation search successful,"
                " but node aborted because of earlier failure!\n",
                item->sq());
        return;
    }

    /* If this node was tried earlier and reached a recursion that failed
     * because a child node failed many times, it is absolutely possible
     * that we already have children in the ABANDONED or maybe in DONE
     * state. It does not matter.
     */
    ASSERT_ALWAYS(item->children_by_status[PENDING].empty());
    ASSERT_ALWAYS(item->children_by_status[IN_PROGRESS].empty());

    visited.insert(item->contender.rel);

    verbose_fmt_print(0, 0, "Taken: {}\n", item->contender.rel);

    verbose_fmt_print(0, 1, "# taking path: {:<{}} -> {}\t{}\n",
            item->shortname(), item->depth(),
            item->contender.outstanding.empty() ? "done" :
            join(item->contender.outstanding,
                " ",
                [&](special_q const & x) { return x.shortname(); }
                ),
            item->sq());

    if (recursive_descent) {
        recurse_node_unlocked(item);

        for(auto const & qq : item->contender.outstanding) {
            special_q child(qq);
            verbose_fmt_print(0, 1, 
                    "# [descent] pushing {} [{}]"
                    " to todo list (now size {})\n",
                    child, child.shortname(), all_pending.size());
            new_node_unlocked(child, item);
        }

        /* In the happy case, we can jump to the "done" status right away! */
        if (item->contender.outstanding.empty()) {
            done_node_unlocked(item);
        } else {
            /* it makes sense to wake up waiters, if there are any */
            work_to_do.notify_all();
        }
    } else {
        /* do not register the new nodes as pending. This mode is mostly
         * useless.
         */
        done_node_unlocked(item);
    }

    verbose_fmt_print(0, 2,
            "# CURRENT STATUS (after decision on {} {})"
            " [npending={}, created={} pulled={} done={} abandoned={}]\n"
            "{}"
            "# END CURRENT STATUS\n",
            item->shortname(), item->sq(),
            all_pending.size(),
            created, pulled, done, abandoned,
            special_q_task_tree::prefixed { .t=&forest, .prefix="# "});
}

void special_q_task_collection_tree::abandon_node(special_q_task_tree * item, int max_descent_attempts_allowed)
{
    constexpr auto IN_RECURSION = special_q_task::status_code::IN_RECURSION;
    constexpr auto IN_PROGRESS = special_q_task::status_code::IN_PROGRESS;
    constexpr auto DONE = special_q_task::status_code::DONE;
    constexpr auto PENDING = special_q_task::status_code::PENDING;

    const std::lock_guard<std::mutex> lock(tree_lock);

    verbose_fmt_print (0, 1, 
            "# taking path: {} -> loop (#{})\t{}\n",
            item->shortname(), item->try_again, item->sq());
    verbose_fmt_print (0, 1,
            "# [descent] Failed to find a relation"
            " for {} [{}] (iteration {}). Putting back to todo list.\n",
            item->sq(), item->shortname(), item->try_again);

    verbose_fmt_print (0, 1,
        "# Aborting {} {} ;"
        " children: {} pending, {} in progress, {} in recursion, {} done\n",
        item->shortname(),
        item->sq(),
        item->children_by_status[PENDING].size(),
        item->children_by_status[IN_PROGRESS].size(),
        item->children_by_status[IN_RECURSION].size(),
        item->children_by_status[DONE].size());

    abandon_node_unlocked(item, max_descent_attempts_allowed);

    verbose_fmt_print(0, 2,
            "# CURRENT STATUS (after abort of {} {})"
            " [npending={}, created={} pulled={} done={} abandoned={}]\n"
            "{}"
            "# END CURRENT STATUS\n",
            item->shortname(), item->sq(),
            all_pending.size(),
            created, pulled, done, abandoned,
            special_q_task_tree::prefixed { .t=&forest, .prefix="# "});
}

/* abandoning a node happens when all relations that we tried for this
 * node have led to at least one failure. The course of action is:
 *
 *  - cancel all children that are still in progress or pending
 *
 *  - if the t->try_again count is small, redo sieving with increased
 *    parameters (typically with increased A). The relations that we
 *    found so far will be skipped on the further pass because of the
 *    special_q_task_collection_tree::must_avoid mechanism.
 *
 *  - if no further option is left, abandon the parent.
 *
 * abandoning is a difficult operation in the multithreaded context,
 * because some children or grandchildren might still be in progress.
 */
void special_q_task_collection_tree::abandon_node_unlocked(special_q_task_tree * item, int max_descent_attempts_allowed, bool recursive_failure)
{
    constexpr auto IN_RECURSION = special_q_task::status_code::IN_RECURSION;
    constexpr auto IN_PROGRESS = special_q_task::status_code::IN_PROGRESS;
    constexpr auto DONE = special_q_task::status_code::DONE;
    constexpr auto ABANDONED = special_q_task::status_code::ABANDONED;
    constexpr auto PENDING = special_q_task::status_code::PENDING;

    if (item->status == ABANDONED)
        return;

    /* We only have IN_PROGRESS or IN_RECURSION children if the job is
     * IN_RECURSION.
     * An IN_PROGRESS node is still within relation search.
     * An IN_RECURSION node that we're abandoning just had a child fail.
     */

    /* The question of the DONE nodes is a bit of a heart breaker. If
     * they're done, then we may have a preference towards making these
     * available for further descents, in case we happen to encounter
     * these sqs again. Or at least not forbid the relations that led to
     * these splittings. Our behavior here, on the other hand, is the
     * most extreme one: even though these splittings have proved useful,
     * we just burn them, and forbid their further use.
     */
    for(const status_code s : { PENDING, IN_PROGRESS, IN_RECURSION }) {
        /* We can't loop on children_by_status since abandoning the
         * children nodes will modify the list! */
        auto & L = item->children_by_status[s];
        for( ; !L.empty() ; ) 
            abandon_node_unlocked(*L.begin(), max_descent_attempts_allowed, true);
    }

    /* special case for DONE nodes. We won't abandon them in the
     * recursive_failure branch below. By definition, all their children
     * are also DONE, so we won't update any status there. Bottom line,
     * the branch below implies that we don't want to touch any DONE node
     * at all. (the loop above would be an infinite loop, see #30127)
     */

    if (recursive_failure) {
        /* Here we don't update the timings.
         *  - if the task is IN_PROGRESS, the spent time will be tallied
         *    later.
         *  - if it's IN_RECURSION or DONE, it's already correctly
         *    counted.
         *  - if it's PENDING, there's no spent time to update.
         *  - it probably cannot be ABANDONED in significant cases (the
         *    only case I see is that when this failure ripples to upper
         *    layers, which then see abandoned children, but that is an
         *    implementation detail.
         */
        if (item->status != DONE) {
            /* Do not change DONE nodes to ABANDONED status
             */
            abandoned += item != &forest;
            item->update_status(item->status, ABANDONED);
        }
    } else {
        ASSERT_ALWAYS(item->status == IN_PROGRESS || item->status == IN_RECURSION);
        if (!item->sq()) {
            abandoned += item != &forest;
            item->update_status(item->status, ABANDONED);
            verbose_fmt_print (0, 1, "# The root node fails. This is sad.\n");
        } else if (item->try_again < max_descent_attempts_allowed) {
            item->spent += seconds();
            item->try_again++;
            all_pending.push_back(item);
            item->update_status(item->status, PENDING);
        } else {
            ASSERT_ALWAYS(item != &forest);
            item->spent += seconds();
            abandoned++;
            item->update_status(item->status, ABANDONED);
            verbose_fmt_print (0, 1, 
                    "# {} [{}] -> failed {} times,"
                    " now failing the parent node {}\n",
                    item->shortname(), 
                    item->sq(),
                    item->try_again,
                    item->parent->sq());
            abandon_node_unlocked(item->parent, max_descent_attempts_allowed);
        }
    }
}

void special_q_task_collection_tree::postprocess(special_q_task * task, int max_descent_attempts_allowed, timetree_t & timer_special_q)/*{{{*/
{
    auto * tt = dynamic_cast<special_q_task_tree *>(task);

    if (!tt) return;

    SIBLING_TIMER(timer_special_q, "descent");
    TIMER_CATEGORY(timer_special_q, bookkeeping());

    if (tt->contender.rel) {
        /* Even if not going for recursion, store this as being a
         * winning relation. This is useful for preparing the hint
         * file, and also for the initialization of the descent.
         */
        take_decision(tt);
    } else {
        abandon_node(tt, max_descent_attempts_allowed);
    }
}/*}}}*/

special_q_task * special_q_task_collection_simple::pull()
{
    auto q = todo.feed_and_pop();
    if (!q)
        return nullptr;
    created++;
    pulled++;

    auto * kid = new special_q_task_simple(q);

    const std::lock_guard<std::mutex> lock(history_lock);
    history.emplace_back(std::unique_ptr<special_q_task_simple>(kid));

    return kid;
}

