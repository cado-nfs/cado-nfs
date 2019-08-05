#ifndef TREE_STATS_HPP_
#define TREE_STATS_HPP_

#include <string>
#include <vector>
#include <map>
#include <climits>
#include <cmath>
#include <ostream>
#include "macros.h"
#include "timing.h"
#include "params.h"

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

extern size_t lingen_round_operand_size(size_t x, int bits = 2);

/* 
 *  - TODO: maybe replace map by simple vectors -- this will allow a
 *    better ordering of the reported timings. Maybe this will also make
 *    it possible to remove or simplify curstack
 */

class tree_stats {
    static int max_nesting;
    public:
    static void interpret_parameters(cxx_param_list & pl);
    static void declare_usage(cxx_param_list & pl);
    private:
    struct step_time {
        std::string name;
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
        protected:
        /* If negative, then this is a transition level only */
        int _total_ncalls = 0;
        public:
        void set_total_ncalls(int t) { _total_ncalls = t; }
        unsigned int total_ncalls() const { return std::abs(_total_ncalls); }
        bool is_transition_level() const { return _total_ncalls < 0; }

        inline double projected_time() const {
            return ncalled ? real * total_ncalls() / (double) ncalled : NAN;
        }
        bool is_hot() const { return real < 0; }
        void heat_up(double wct = wct_seconds()) { real -= wct; }
        void cool_down(double wct = wct_seconds()) { real += wct; }


        // expected time per call. This is set by tree_stats::plan_smallstep()
        // (actually it is _increased_ by plan_smallstep)
        double planned_time = 0;
        unsigned int planned_calls = 0;
        unsigned int items_per_call = 1;
        // items_pending is only non-zero when we're doing the small
        // steps. It doesn't make sense to take it into account for +=.
        unsigned int items_pending = 0;

        step_time(std::string const & name) : name(name) {}
        step_time(step_time const &) = default;
        step_time(step_time &&) = default;

        typedef std::map<std::string, step_time> steps_t;
        steps_t steps;

        step_time const * find_pending_smallstep() const {
            if (items_pending) return this;
            for(auto const & x : steps) {
                auto y = x.second.find_pending_smallstep();
                if (y) return y;
            }
            return nullptr;
        }
        bool has_pending_smallsteps() const {
            return find_pending_smallstep() != nullptr;
        }

        step_time & operator+=(step_time const & x) {
            ASSERT_ALWAYS(!is_hot());
            ASSERT_ALWAYS(_total_ncalls == 0 || _total_ncalls == x._total_ncalls);
            _total_ncalls = x._total_ncalls;
            ASSERT_ALWAYS(items_per_call == 1 || items_per_call == x.items_per_call);
            ASSERT_ALWAYS(ncalled <= planned_calls);
            ASSERT_ALWAYS(ncalled + 1 >= planned_calls);
            ASSERT_ALWAYS(x.ncalled <= x.planned_calls);
            ASSERT_ALWAYS(x.ncalled + 1 >= x.planned_calls);
            items_per_call = x.items_per_call;
            if (!x.is_hot()) {
                real += x.real;
                ncalled += x.ncalled;
            }
            planned_time += x.planned_time;
            planned_calls += x.planned_calls;
            ASSERT_ALWAYS(ncalled <= planned_calls);
            ASSERT_ALWAYS(ncalled + 1 >= planned_calls);
            for(auto const & s : x.steps) {
                auto itb = steps.emplace(s);
                step_time & N(itb.first->second);
                ASSERT_ALWAYS(N.name == s.second.name);
                if (!itb.second) N += s.second;
                // steps[s.first] += s.second;
            }
            return *this;
        }
    };

    struct function_with_input_size {
        std::string func;
        unsigned int inputsize;
        function_with_input_size(std::string const & func, unsigned int inputsize)
            : func(func)
            , inputsize(inputsize)
        {
        }
        bool operator<(function_with_input_size const& a) const {
            if (func < a.func) return true;
            if (func > a.func) return false;
            return lingen_round_operand_size(inputsize) < lingen_round_operand_size(a.inputsize);
        }
        bool operator==(function_with_input_size const& a) const {
            return func == a.func && lingen_round_operand_size(inputsize) == lingen_round_operand_size(a.inputsize);
        }
    };

    struct running_stats : public step_time {
        /* inputsize is tracked only for very top level step_time
         * objects. It's (currently) unused for small steps */
        unsigned int inputsize;
        bool leaf;
        double time_children;
        running_stats(std::string const & name, unsigned int inputsize, bool leaf)
            : step_time(name)
            , inputsize(inputsize)
            , leaf(leaf)
            , time_children(0)
        { }
        std::vector<steps_t::iterator> nested_substeps;
        step_time& current_substep() { return nested_substeps.back()->second; }
        step_time const & current_substep() const { return nested_substeps.back()->second; }
    };

    struct function_stats : public step_time {
        unsigned int min_inputsize = UINT_MAX;
        unsigned int max_inputsize = 0;
        function_stats & operator+=(running_stats const&);
        function_stats (running_stats const& r) : step_time(r) {
            ASSERT_ALWAYS(!r.is_hot());
            min_inputsize = max_inputsize = r.inputsize;
        }
    };

    struct level_stats : public std::map<function_with_input_size, function_stats> {
        /* compute the projected time at this level */
        double projected_time();
        double last_printed_projected_time;
    };

    std::vector<level_stats> levels;
    typedef std::vector<std::pair<function_with_input_size, running_stats>> curstack_t;
    curstack_t curstack;

    unsigned int tree_total_breadth;
    double last_print_time = 0;
    std::pair<unsigned int, unsigned int> last_print_position { 0,0 };

    static std::ostream& recursively_print_substeps_at_depth(
        std::ostream & os,
        step_time::steps_t const & FS,
        unsigned int ncalled,
        unsigned int total_ncalls,
        int nesting);

    void print(unsigned int level);

private:
    void enter(std::string const & func, unsigned int inputsize, int total_ncalls, bool leaf = false); 
    void leave();
    unsigned int depth = 0;
    unsigned int transition_levels_in_depth = 0;
public:
    unsigned int non_transition_depth() const {
        return depth - transition_levels_in_depth;
    }
    struct sentinel {
        tree_stats & stats;
        sentinel(sentinel const&) = delete;
        sentinel(tree_stats & stats,
                std::string const & func, unsigned int inputsize, unsigned int total_ncalls, bool leaf = false)
            : stats(stats) { stats.enter(func, inputsize, total_ncalls, leaf); }
        ~sentinel() { stats.leave(); }
    };
    struct transition_sentinel {
        tree_stats & stats;
        transition_sentinel(transition_sentinel const&) = delete;
        transition_sentinel(tree_stats & stats,
                std::string const & func, unsigned int inputsize, unsigned int total_ncalls, bool leaf = false)
            : stats(stats) { stats.transition_levels_in_depth++; stats.enter(func, inputsize, -total_ncalls, leaf); }
        ~transition_sentinel() { stats.leave(); stats.transition_levels_in_depth--; }
    };

    void begin_smallstep(std::string const & func, unsigned int ncalls=1);
    void end_smallstep();
    /* skip = begin+end, but do nothing inbetween */
    void skip_smallstep(std::string const & func, unsigned int ncalls=1);
    struct smallstep_sentinel {
        tree_stats & stats;
        smallstep_sentinel(smallstep_sentinel const&) = delete;
        smallstep_sentinel(tree_stats & stats,
                std::string const & func, unsigned int ncalls=1)
            : stats(stats) { stats.begin_smallstep(func, ncalls); }
        ~smallstep_sentinel() { stats.end_smallstep(); }
    };
    bool local_smallsteps_done() const;


    void begin_plan_smallstep(std::string const & func, weighted_double const &);
    void end_plan_smallstep();
    struct plan_smallstep_sentinel {
        tree_stats & stats;
        plan_smallstep_sentinel(plan_smallstep_sentinel const&) = delete;
        plan_smallstep_sentinel(tree_stats & stats,
                std::string const & func, weighted_double const & theory)
            : stats(stats) { stats.begin_plan_smallstep(func, theory); }
        ~plan_smallstep_sentinel() { stats.end_plan_smallstep(); }
    };
    inline void plan_smallstep(std::string const & func, weighted_double const & theory)
    {
        plan_smallstep_sentinel dummy(*this, func, theory);
    }
    inline void plan_smallstep(std::string const & func, double theory)
    {
        plan_smallstep(func, { 1, theory });
    }

    void final_print();
};

#endif	/* TREE_STATS_HPP_ */
