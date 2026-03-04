#ifndef SIEVE_LAS_SPECIAL_Q_TASK_COLLECTION_HPP_
#define SIEVE_LAS_SPECIAL_Q_TASK_COLLECTION_HPP_

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <condition_variable>
#include <list>
#include <memory>
#include <mutex>
#include <set>
#include <deque>

#include "cado_poly.h"
#include "params.h"
#include "special-q.hpp"
#include "relation.hpp"
#include "las-todo-list.hpp"
#include "tdict.hpp"
#include "timing.h"

#include "las-special-q-task.hpp"
#include "las-special-q-task-tree.hpp"
#include "las-special-q-task-simple.hpp"

#include "sieve-methods.hpp"

/* A "task collection" designates the set of special-q's that a single
 * las process handles over its lifetime.
 *
 * If we're doing a special-q descent, we use the complicated
 * special_q_task_collection_tree object. Otherwise (normal case) we use
 * the much simpler special_q_task_collection_simple.
 *
 * In most cases, special q tasks are created as unique_ptr's, and used
 * via the base interface special_q_task * (which can point to an
 * instance of either special_q_task_simple or special_q_task_tree).
 *
 * Within the special_q_task_collection, a special q can exist in several
 * different states, which are encoded by special_q_task::status_code.
 *
 *  - PENDING: this task has been created, but no las thread currently
 *    works on this special-q. Either it was very recently allocated, or
 *    it is sitting in a storage area, waiting to be picked up (pulled)
 *    by a thread willing to do work.
 *  - IN_PROGRESS: this task was pulled by a thread, which is currently
 *    actively collecting relation for this special-q
 *  - IN_RECURSION (for the descent only): relation collection for this
 *    task is done, but some children tasks were created out of it. These
 *    tasks are not done yet, so the final state of this task depends on
 *    them.
 *  - ABANDONED (for the descent only): we could not find any useful
 *    relation for this task, or all the relations that we found have
 *    failed.  Some extra attempts are allowed, but once this is over,
 *    the task is marked as failed (and the parent task has to try
 *    another relation). 
 *  - DONE: everything is good.
 *
 * The state machine is as follows. Most of what is here only applies to
 * the recursive descent case. The "normal", non-descent case is simple
 * enough and does not need much further detail.
 *
 * Creation of a task is done either from the las_todo_list object, or
 * antonomously from the recursion (in the descent case). When a task is
 * created, the counter `created` is increased by one.
 *
 * Pulling a task (via the `pull()` method) either picks a task from the
 * `all_pending` structure/
 */

struct special_q_task_collection_base {

    /* this is the source from which we pull the work to be done */
    std::unique_ptr<todo_list_base> todo;

    /* done+abandoned is always <= created
     * pulled - (done + abandoned) is the number of retries
     */

    size_t created = 0;
    size_t pulled = 0;
    size_t done = 0;
    size_t abandoned = 0;

    template<sieve_method Algo>
    special_q_task_collection_base(
            cxx_cado_poly const & cpoly,
            cxx_param_list & pl,
            Algo)
        : todo(todo_list_base::create<typename Algo::todo_list>(cpoly, pl))
    { }
    special_q_task_collection_base(special_q_task_collection_base const& t) = delete;
    special_q_task_collection_base(special_q_task_collection_base && t) = delete;
    special_q_task_collection_base& operator=(special_q_task_collection_base const& t) = delete;
    special_q_task_collection_base& operator=(special_q_task_collection_base && t) = delete;
    virtual ~special_q_task_collection_base() = default;

    public:
    virtual bool must_avoid(relation_ab const& ab) const = 0;

    virtual void new_candidate_relation(las_info const &, special_q_task *, relation &) = 0;

    public:
    virtual special_q_task * pull() = 0;
    virtual void postprocess(special_q_task * task, int, timetree_t & timer_special_q) = 0;

    virtual void display_summary(int, int) {}

    template<sieve_method Algo>
    static std::unique_ptr<special_q_task_collection_base> create(
            cxx_cado_poly const & cpoly,
            cxx_param_list & pl);
};

struct special_q_task_collection_simple : public special_q_task_collection_base {
    typedef special_q_task::status_code status_code;

    std::mutex history_lock;
    std::list<std::unique_ptr<special_q_task_simple>> history;

    special_q_task_collection_simple(
            cxx_cado_poly const & cpoly,
            cxx_param_list & pl,
            sieve_method auto tag)
        : special_q_task_collection_base(cpoly, pl, tag)
    { }

    public:
    bool must_avoid(relation_ab const&) const override { return false; }
    void new_candidate_relation(las_info const &, special_q_task *, relation &) override {}

    public:
    special_q_task * pull() override;
    void postprocess(special_q_task *, int, timetree_t &) override {};
};

struct special_q_task_collection_tree : public special_q_task_collection_base {
    private:
        /* we have const members which need to lock the mutex */
        mutable std::mutex tree_lock;
        std::condition_variable work_to_do;
    public:
    /* A "descent tree" is really a set of descent trees, one for each
     * top-level special-q we've considered in the main loop in las. The
     * top-level roots are found in [forest].
     */

    typedef special_q_task::status_code status_code;

    /* We abuse the tree structure to make one of them a parent of all
     * trees that we'll consider. By doing so, we'll be able to
     * distinguish between pending, in-progress, and done special-qs.
     */
    special_q_task_tree forest;

    std::set<relation_ab> visited;

    std::deque<special_q_task_tree *> all_pending;

    special_q_task_collection_tree(
            cxx_cado_poly const & cpoly,
            cxx_param_list & pl)
        : special_q_task_collection_base(cpoly, pl, NFS{})
    {
        forest.status = status_code::IN_RECURSION;
        forest.spent = -seconds();
        // created++;
        // pulled++;
    }

    private:
    void new_node_unlocked(special_q const & doing, special_q_task_tree * parent);

    void done_node_unlocked(special_q_task_tree * item);
    void recurse_node_unlocked(special_q_task_tree * item);

    void abandon_node_unlocked(special_q_task_tree *, int, bool only_down = false);

    /***/

    void new_node(special_q const & doing, special_q_task_tree * parent) {
        const std::lock_guard<std::mutex> lock(tree_lock);
        new_node_unlocked(doing, parent);
    }

    void done_node(special_q_task_tree * item) {
        const std::lock_guard<std::mutex> lock(tree_lock);
        done_node_unlocked(item);
    }

    void recurse_node(special_q_task_tree * item) {
        const std::lock_guard<std::mutex> lock(tree_lock);
        recurse_node_unlocked(item);
    }

    void take_decision(special_q_task_tree * t);

    void abandon_node(special_q_task_tree * item, int);

    special_q_task_tree * pull_internal();

    public:
    bool must_avoid(relation_ab const& ab) const override {
        const std::lock_guard<std::mutex> lock(tree_lock);
        return visited.find(ab) != visited.end();
    }

    void new_candidate_relation(las_info const & las, special_q_task * task, relation & rel) override {
        dynamic_cast<special_q_task_tree *>(task)->new_candidate_relation(las, rel, tree_lock);
    }

    public:
    special_q_task * pull() override;
    void postprocess(special_q_task *, int, timetree_t &) override;
    void display_summary(int, int) override;
};




#endif	/* SIEVE_LAS_SPECIAL_Q_TASK_COLLECTION_HPP_ */
