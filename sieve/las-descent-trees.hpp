#ifndef CADO_LAS_DESCENT_TREES_HPP
#define CADO_LAS_DESCENT_TREES_HPP

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <algorithm>
#include <list>
#include <mutex>
#include <set>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <gmp.h>
#include "fmt/format.h"

#include "las-todo-entry.hpp"
#include "macros.h"
#include "relation.hpp"
#include "timing.h"
#include "verbose.h"

#ifdef isfinite
/* isfinite is c99 and std::isfinite is c++11 ; it's not totally clear
 * that #include <cmath> + accessing std::isfinite works.
 *
 * Under some conditions, we can get #define'd C functions, which
 * obviously invalidate the C++ prototype (icc version 16.0.3 based on
 * gcc-6.1.0 on CentOS 7.2.1511)
 */
#undef isfinite
#endif

struct descent_tree {
    private:
        /* we have const members which need to lock the mutex */
        mutable std::mutex tree_lock;
    public:
    struct tree_label {
        int side = 0;
        relation::pr pr;
        tree_label() = default;
        tree_label(int side, relation::pr pr )
            : side(side)
            , pr(std::move(pr)) {}
        tree_label(int side, mpz_srcptr p, mpz_srcptr r) :side(side), pr(p, r) {}
        std::string shortname() const {
            return fmt::format("{}@{}", mpz_sizeinbase(pr.p, 2), side);
        }
        std::string fullname() const {
            return fmt::format("{} {} {}", side, pr.p, pr.r);
        }
        bool operator<(const tree_label& o) const {
            if (pr < o.pr) return true;
            if (o.pr < pr) return false;
            return side < o.side;
        }
    };
    /* For descent mode: we compute the expected time to finish given the
     * factor sizes, and deduce a deadline.  Assuming that not all
     * encountered factors are below the factor base bound, if we expect
     * an additional time T to finish the decomposition, we keep looking
     * for a better decomposition for a grace time which is computed as
     * x*T, for some configurable ratio x (one might think of x=0.2 for
     * instance. x is the grace_time_ratio member), which defines a
     * ``deadline'' for next step.  [If all factors happen to be smooth,
     * the deadline is immediate, of course.] If within the grace period,
     * a new relation is found, with an earlier implied deadline, the
     * deadline is updated. We say that the "decision is taken" when the
     * deadline passes, and the las machinery is told to decide that it
     * should proceed with the descent, and stop processing the current
     * special-q.
     */
    struct candidate_relation {
        relation rel;
        std::vector<std::pair<int, relation::pr> > outstanding;
        double time_left = INFINITY;
        double deadline = INFINITY;
        // bool marked_taken;      /* false until we take the decision */
        candidate_relation() = default;
        candidate_relation& operator=(candidate_relation const& o) {
            if (this == &o) return *this;
            /* nothing very fancy, except that we keep the old deadline.
             * */
            rel = o.rel;
            outstanding = o.outstanding;
            time_left = o.time_left;
            deadline = std::min(deadline, o.deadline);
            return *this;
        }
        candidate_relation(candidate_relation const & o) = default;
        candidate_relation& operator=(candidate_relation && o) noexcept {
            rel = o.rel;
            outstanding = o.outstanding;
            time_left = o.time_left;
            deadline = std::min(deadline, o.deadline);
            return *this;
        }
        candidate_relation(candidate_relation && o) noexcept = default;
        ~candidate_relation() = default;

        bool operator<(candidate_relation const& b) const
        {
            if (!rel) return false;
            if (!b.rel) return true;
            if (std::isfinite(time_left)) { return time_left < b.time_left; }
            return outstanding.size() < b.outstanding.size();
        }
        explicit operator bool() const { return (bool) rel; }
        bool wins_the_game() const {
            return (bool) rel && (outstanding.empty() || seconds() >= deadline);
        }
        void set_time_left(double t, double grace_time_ratio) {
            time_left = t;
            if (outstanding.empty()) {
                deadline = seconds();
            } else {
                deadline = seconds() + grace_time_ratio * t;
            }
        }
    };
    struct tree {
        tree_label label;
        int root_depth = 0;
        double spent = 0;
        candidate_relation contender;
        std::list<std::unique_ptr<tree>> children;
        bool try_again = false;
        explicit tree(tree_label label, int root_depth = 0)
            : label(std::move(label))
            , root_depth(root_depth)
        {}
        bool is_successful() const {
            if (!contender.rel && !try_again)
                return false;
            for(auto const & c : children) {
                if (!c->is_successful())
                    return false;
            }
            return true;
        }
        /* this returns the tree depth _below this level_ only. Add
         * this->root_depth to get the overall depth (that this subtree
         * contributes to) */
        int depth() const {
            int d = 0;
            for(auto const & c : children)
                d = std::max(d, 1 + c->depth());
            return d;
        }
        /* this returns the tree weight _below this level_ only. */
        int weight() const {
            int w = 1;
            for(auto const & c :children)
                w += c->weight();
            return w;
        }
    };

    /* A "descent tree" is really a set of descent trees, one for each
     * top-level special-q we've considered in the main loop in las. The
     * top-level roots are found in [forest].
     *
     * At a given point in the process, the list [current] contains
     * pointers to the different levels of the tree that form the lineage
     * of the last special-q being considered. More precisely,
     * current.back()->label is the special-q that is being considered,
     * current.back()->children() is empty, and
     * (*----current.end())->label is the one whose consideration led us
     * to consider it, and so on.
     *
     * forest (and its children) owns the pointers. current does not.
     *
     * pointers in current are all within forest.back()
     */
    std::list<std::unique_ptr<tree>> forest;
    std::list<tree *> current;       /* stack of trees */

    /* This is an ugly temporary hack */
    typedef std::set<relation_ab> visited_t;
    visited_t visited;


    ~descent_tree() = default;

    /* designing a copy ctor for this structure would be tedious */
    descent_tree() = default;
    descent_tree(descent_tree const& t) = delete;
    descent_tree(descent_tree && t) = delete;
    descent_tree& operator=(descent_tree const& t) = delete;
    descent_tree& operator=(descent_tree && t) = delete;

    void new_node(las_todo_entry const & doing) {
        const std::lock_guard<std::mutex> lock(tree_lock);
        const int level = doing.depth;
        ASSERT_ALWAYS(level == (int) current.size());
        std::unique_ptr<tree> kid(new tree(tree_label(doing.side, doing.p, doing.r)));
        kid->spent = -seconds();
        if (current.empty()) {
            current.push_back(kid.get());
            forest.push_back(std::move(kid));
        } else {
            auto & tail = *current.back();
            current.push_back(kid.get());
            tail.children.push_back(std::move(kid));
        }
    }

    void ditch_node() {
        tree * kid = current.back();
        current.pop_back();
        if (!current.empty()) {
            ASSERT_ALWAYS(current.back()->children.back().get() == kid);
            current.back()->children.pop_back();
        } else {
            ASSERT_ALWAYS(forest.back().get() == kid);
            forest.pop_back();
        }
        /* no deletion is necessary. It's taken care of by the two
         * pop_back() operations above */
    }
    void done_node() {
        current.back()->spent += seconds();
        current.pop_back();
    }

    candidate_relation const& current_best_candidate() const {
        /* outside multithreaded context */
        return current.back()->contender;
    }

    void mark_try_again(int i) {
        current.back()->try_again = i;
    }

    /* return true if the decision to go to the next step of the descent
     * should be taken now, and register this relation has having been
     * taken */
    bool must_take_decision() {
        const std::lock_guard<std::mutex> lock(tree_lock);
        return current.back()->contender.wins_the_game();
    }

    void take_decision() {
        /* must be called outside multithreaded context */
        // current.back()->contender.marked_taken = true;
        visited.insert(current.back()->contender.rel);
    }

    /* this returns true if the decision should be taken now */
    bool new_candidate_relation(candidate_relation& newcomer)
    {
        const std::lock_guard<std::mutex> lock(tree_lock);
        candidate_relation & defender(current.back()->contender);
        if (newcomer < defender) {
            if (newcomer.outstanding.empty()) {
                verbose_output_print(0, 1, "# [descent] Yiippee, splitting done\n");
            } else if (std::isfinite(defender.deadline)) {
                /* This implies that newcomer.deadline is also finite */
                const double delta = defender.time_left-newcomer.time_left;
                verbose_output_print(0, 1, "# [descent] Improved ETA by %.2f\n", delta);
            } else if (defender) {
                /* This implies that we have fewer outstanding
                 * special-q's */
                verbose_output_print(0, 1, "# [descent] Improved number of children to split from %u to %u\n",
                        (unsigned int) defender.outstanding.size(),
                        (unsigned int) newcomer.outstanding.size());
            }
            defender = newcomer;
            if (!defender.outstanding.empty()) {
                verbose_output_print(0, 1, "# [descent] still searching for %.2f\n", defender.deadline - seconds());
            }
        }
        return defender.wins_the_game();
    }

    bool must_avoid(relation_ab const& ab) const {
        const std::lock_guard<std::mutex> lock(tree_lock);
        return visited.find(ab) != visited.end();
    }
    bool must_avoid(relation const& rel) const {
        return must_avoid((relation_ab const &) rel);
    }
    int depth() const {
        return (int) current.size();
    }
    int display_tree(FILE* o, tree const * t, std::string const& prefix);
    void display_last_tree(FILE * o);

    void display_all_trees(FILE * o);
};
#endif	/* CADO_LAS_DESCENT_TREES_HPP */
