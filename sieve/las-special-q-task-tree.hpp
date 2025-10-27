#ifndef SIEVE_LAS_SPECIAL_Q_TASK_TREE_HPP_
#define SIEVE_LAS_SPECIAL_Q_TASK_TREE_HPP_

#include <list>
#include <map>
#include <memory>
#include <mutex>
#include <ostream>
#include <set>
#include <string>
#include <utility>

#include "fmt/base.h"
#include "fmt/ostream.h"

#include "las-descent-candidate-relation.hpp"
#include "las-special-q-task.hpp"
#include "special-q.hpp"
#include "macros.h"
#include "relation.hpp"

struct las_info;

/* a special_q_task_tree is actually one node (and pointers to its
 * children) in a tree that is considered in the descent process. The
 * general access interface is special_q_task_collection_tree, which is
 * an instance of special_q_task_collection
 */

/* the state machine of a tree node is as follows:
 *
 * pending -> in progress -> in recursion -> done
 *            in progress <- in recursion
 *            {in progress, in recursion} -> abandoned
 *
 * going from IN_PROGRESS to DONE must be in two steps, via
 * IN_RECURSION first (the two can be done in lockstep)
 */
struct special_q_task_tree : public special_q_task {
    int root_depth = 0;
    double spent = 0;
    descent_candidate_relation contender;

    std::list<relation> alternatives;
    std::list<std::unique_ptr<special_q_task_tree>> children;

    struct sort_subtrees_by_size {
        bool operator()(special_q_task_tree const * a, special_q_task_tree const * b) const {
            special_q const & aq(*a);
            special_q const & bq(*b);
            return aq < bq;
        }
    };

    std::map<status_code, std::set<special_q_task_tree *, sort_subtrees_by_size>> children_by_status;

    special_q_task_tree * parent = nullptr;

    int try_again = 0;

    private:
    friend struct special_q_task_collection_tree;

    void vivify_lists() {
        /* vivify all lists */
        children_by_status[PENDING];
        children_by_status[IN_PROGRESS];
        children_by_status[IN_RECURSION];
        children_by_status[DONE];
        children_by_status[ABANDONED];
    }

    special_q_task_tree() {
        vivify_lists();
    }
    public:


    explicit special_q_task_tree(special_q q, special_q_task_tree * parent = nullptr)
        : special_q_task(std::move(q))
          , root_depth(parent ? (1 + parent->root_depth) : 0)
          , parent(parent)
          {
              vivify_lists();
          }

    int depth() const {
        int d = 0;
        for(auto const * p = parent ; p ; p = p->parent, d++);
        return d;
    }
    void update_child_status(special_q_task_tree *, status_code, status_code);
    void update_status(status_code before, status_code after) override {
        if (parent) {
            parent->update_child_status(this, before, after);
        } else {
            /* then it degrades to something very straightforward, of
             * course.
             */
            ASSERT_ALWAYS(status == before);
            status = after;
        }
    }

    struct prefixed {
        special_q_task_tree const * t;
        std::string prefix;
    };
    bool must_take_decision() const override {
        return contender.wins_the_game();
    }
    bool new_candidate_relation(las_info const & las, relation & rel, std::mutex & mm);
};

std::ostream& operator<<(std::ostream&, special_q_task_tree::prefixed const &);
std::ostream& operator<<(std::ostream&, special_q_task_tree::status_code const &);

namespace fmt {
    template <> struct formatter<special_q_task_tree::prefixed>: ostream_formatter {};
    template <> struct formatter<special_q_task_tree::status_code>: ostream_formatter {};
}


#endif	/* SIEVE_LAS_SPECIAL_Q_TASK_TREE_HPP_ */
