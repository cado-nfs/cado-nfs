#ifndef TREE_STATS_HPP_
#define TREE_STATS_HPP_

#include <string>
#include <vector>
#include <map>
#include <climits>
#include <ostream>
#include "macros.h"
#include "timing.h"

/* This structure is meant to help the timing of a recursive tree-like
 * algorithms like the linear generator algorithm we use in block
 * Wiedemann. The features are the following.
 *
 *  - the algorithms is made of (recursive) levels. Level 0 is the top of
 *    the recursion, called with an argument that we refer to as the
 *    total breadth.
 *  - each level is allowed to call possibly different functions (in
 *    principle, at most two, byt we don't make any difference).
 *  - each function is allowed to present its timings as divided into
 *    several (hierarchical) small steps, which are either functions or
 *    fragments of functions.
 *
 * While the algorithm is running, the tree_stats structure maintains
 * several timers. Some are "hot", some are "cold". "Hot" timers are
 * timers that are still counting.
 *
 * We have:
 *  - tree_stats::levels: a vector, indexed by the recursive level, of
 *    all functions we have encountered so far. More precisely, each item
 *    in this vector is a map indexed by the function name.
 *  - the value_type of the preceding map is a generic step_time object,
 *    wrapped as a function_stats object that carries some additional data.
 *  - likewise, we have a stack of all hot timers (from depth 0 down to
 *    the current depth). Those are step_time objects again, but wrapped
 *    into running_stats object that carry more data.
 *  - step_time timer objects an have recursive structure, and contain
 *    sub-timers.
 */

/* 
 *  - TODO: maybe replace map by simple vectors -- this will allow a
 *    better ordering of the reported timings. Maybe this will also make
 *    it possible to remove or simplify curstack
 */

class tree_stats {
    
    struct step_time {
        double real = 0;
        /* We have an ncalled field for each small step, because we _can_
         * have a substep called several times per enclosing level.
         * Therefore, while at first sight we would think that all steps
         * at a given level are called the same number of times, it need
         * not be the case.
         *
         * Note that there is an important difference between the ncalled
         * field for functions (that may still be hot at the time of
         * printing) and for small steps (that are _never_ hot when we
         * print). small steps register their ncalled field right when
         * they begin execution, while functions do that only on
         * completion.
         */
        unsigned int ncalled = 0;
        bool is_hot() const { return real < 0; }
        void heat_up(double wct = wct_seconds()) { real -= wct; }
        void cool_down(double wct = wct_seconds()) { real += wct; }

        // expected time per call. This is set by tree_stats::plan_smallstep()
        // (actually it is _increased_ by plan_smallstep)
        double theoretical=0;

        typedef std::map<std::string, step_time> steps_t;
        steps_t steps;
        step_time & operator+=(step_time const & x) {
            ASSERT_ALWAYS(!is_hot());
            if (!x.is_hot()) {
                real += x.real;
                ncalled += x.ncalled;
            }
            theoretical += x.theoretical;
            for(auto const & s : x.steps)
                steps[s.first] += s.second;
            return *this;
        }
    };

    struct running_stats : public step_time {
        /* inputsize is tracked only for very top level step_time
         * objects. It's (currently) unused for small steps */
        unsigned int inputsize;
        bool leaf;
        double time_children;
        running_stats(unsigned int inputsize, bool leaf)
            : inputsize(inputsize)
            , leaf(leaf)
            , time_children(0)
        {}
        std::vector<steps_t::iterator> nested_substeps;
        step_time& current_substep() { return nested_substeps.back()->second; }
    };

    struct function_stats : public step_time {
        unsigned int min_inputsize = UINT_MAX;
        unsigned int max_inputsize = 0;
        unsigned int sum_inputsize = 0;
        /* 8d85180cc : how much of the total input size for this function
         * is actually _not_ processed recursively.
         * This value is either sum_inputsize or 0. (but two functions at
         * the same level may conceivably have different behaviours).
         */
        unsigned int trimmed = 0;
        /* projected_calls and projected_time are cached value, computed
         * by level_stats::projected_time()
         */
        double projected_time = 0;
        unsigned int projected_calls = 0;
        function_stats & operator+=(running_stats const&);
    };

    struct level_stats : public std::map<std::string, function_stats> {
        /* compute the projected time at this level, knowing the total
         * called size at this level as well as the total size that was
         * trimmed in the calls below.
         */
        double projected_time(unsigned int total_breadth, unsigned int trimmed_breadth);
        double last_printed_projected_time;
        /* cached value, computed by projected_time. It is typically used
         * also by loops that call projected_time(), since the
         * trimmed_breadth argument of the projected_time() call is
         * typically the sum of the trimmed_here value for the levels of
         * smaller depth. */
        unsigned int trimmed_here;       /* trimmed at this level ! */
    };

    std::vector<level_stats> levels;
    typedef std::vector<std::pair<std::string, running_stats>> curstack_t;
    curstack_t curstack;

    unsigned int tree_total_breadth;
    double last_print_time = 0;
    std::pair<unsigned int, unsigned int> last_print_position { 0,0 };

    static std::ostream& recursively_print_substeps_at_depth(
        std::ostream & os,
        step_time::steps_t const & FS,
        unsigned int ncalled,
        unsigned int projected_calls,
        int nesting);

    void print(unsigned int level);

public:
    unsigned int depth = 0;
    
    void enter(std::string const & func, unsigned int inputsize, bool leaf = false); 
    void leave();
    struct sentinel {
        tree_stats & stats;
        sentinel(sentinel const&) = delete;
        sentinel(tree_stats & stats,
                std::string const & func, unsigned int inputsize, bool leaf = false)
            : stats(stats) { stats.enter(func, inputsize, leaf); }
        ~sentinel() { stats.leave(); }
    };

    void begin_smallstep(std::string const & func, unsigned int ncalls=1);
    void end_smallstep();
    struct smallstep_sentinel {
        tree_stats & stats;
        smallstep_sentinel(smallstep_sentinel const&) = delete;
        smallstep_sentinel(tree_stats & stats,
                std::string const & func, unsigned int ncalls=1)
            : stats(stats) { stats.begin_smallstep(func, ncalls); }
        ~smallstep_sentinel() { stats.end_smallstep(); }
    };


    void begin_plan_smallstep(std::string const & func, double theory=0);
    void end_plan_smallstep();
    struct plan_smallstep_sentinel {
        tree_stats & stats;
        plan_smallstep_sentinel(plan_smallstep_sentinel const&) = delete;
        plan_smallstep_sentinel(tree_stats & stats,
                std::string const & func, double theory=0)
            : stats(stats) { stats.begin_plan_smallstep(func, theory); }
        ~plan_smallstep_sentinel() { stats.end_plan_smallstep(); }
    };
    inline void plan_smallstep(std::string const & func, double theory=0)
    {
        plan_smallstep_sentinel dummy(*this, func, theory);
    }

    void final_print();
};

#endif	/* TREE_STATS_HPP_ */
