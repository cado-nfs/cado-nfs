#ifndef SIEVE_LAS_WORK_LOGBOOK_HPP_
#define SIEVE_LAS_WORK_LOGBOOK_HPP_

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
#include "las-special-q.hpp"
#include "relation.hpp"
#include "las-todo-list.hpp"
#include "tdict.hpp"
#include "timing.h"

#include "las-special-q-task.hpp"
#include "las-special-q-task-tree.hpp"
#include "las-special-q-task-simple.hpp"

struct special_q_task_collection_base {

    /* this is the source from which we pull the work to be done */
    las_todo_list todo;

    /* done+abandoned is always <= created
     * pulled - (done + abandoned) is the number of retries
     */

    size_t created = 0;
    size_t pulled = 0;
    size_t done = 0;
    size_t abandoned = 0;

    special_q_task_collection_base(cxx_cado_poly const & cpoly, cxx_param_list & pl)
        : todo(cpoly, pl)
    { }
    special_q_task_collection_base(special_q_task_collection_base const& t) = delete;
    special_q_task_collection_base(special_q_task_collection_base && t) = delete;
    special_q_task_collection_base& operator=(special_q_task_collection_base const& t) = delete;
    special_q_task_collection_base& operator=(special_q_task_collection_base && t) = delete;
    virtual ~special_q_task_collection_base() = default;

    virtual special_q_task * pull() = 0;
    virtual void postprocess(special_q_task * task, timetree_t & timer_special_q) = 0;
    virtual bool must_avoid(relation_ab const& ab) const = 0;

    static std::unique_ptr<special_q_task_collection_base> create(cxx_cado_poly const & cpoly, cxx_param_list & pl);
};

struct special_q_task_collection_simple : public special_q_task_collection_base {
    typedef special_q_task::status_code status_code;

    std::mutex history_lock;
    std::list<std::unique_ptr<special_q_task_simple>> history;

    special_q_task_collection_simple(cxx_cado_poly const & cpoly, cxx_param_list & pl)
        : special_q_task_collection_base(cpoly, pl)
    { }

    special_q_task * pull() override;

    bool must_avoid(relation_ab const&) const override { return false; }

    void postprocess(special_q_task *, timetree_t &) override {};
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

    special_q_task_collection_tree(cxx_cado_poly const & cpoly, cxx_param_list & pl)
        : special_q_task_collection_base(cpoly, pl)
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

    void abandon_node_unlocked(special_q_task_tree *, bool only_down = false);

    void abandon_node(special_q_task_tree * item);

    special_q_task_tree * pull_internal();

    public:

    special_q_task * pull() override;

    bool must_avoid(relation_ab const& ab) const override {
        const std::lock_guard<std::mutex> lock(tree_lock);
        return visited.find(ab) != visited.end();
    }

    void postprocess(special_q_task *, timetree_t &) override;
};




#endif	/* SIEVE_LAS_WORK_LOGBOOK_HPP_ */
