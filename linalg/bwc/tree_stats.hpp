#ifndef TREE_STATS_HPP_
#define TREE_STATS_HPP_
// IWYU pragma: no_include <ext/alloc_traits.h>
#include <cstdlib>                       // for abs
#include <utility>   // for pair
#include <string>
#include <vector>
#include <map>
#include <climits>
#include <cmath>
#include <istream> // IWYU pragma: keep
#include <ostream> // IWYU pragma: keep
#include "macros.h"
#include "timing.h"
#include "lingen_round_operand_size.hpp"
struct cxx_param_list; // IWYU pragma: keep


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
    static int max_nesting;
    public:
    static void interpret_parameters(cxx_param_list & pl);
    static void declare_usage(cxx_param_list & pl);
    private:

    // for clarity, fwd-decl. We really _must_ keep the fwd-decl of
    // level_stats, despite what iwyu says (and as a matter of fact, it
    // does not seem to abide by its very own pragma)
    // struct step_time; // IWYU pragma: keep
    struct function_with_input_size; // IWYU pragma: keep
    struct running_stats; // IWYU pragma: keep
    // struct function_stats; // IWYU pragma: keep
    struct level_stats; // IWYU pragma: keep

    typedef std::vector<std::pair<function_with_input_size, running_stats>> curstack_t;

    struct step_time {/*{{{*/
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
        // expected time per call. This is set by tree_stats::plan_smallstep()
        // (actually it is _increased_ by plan_smallstep)
        double planned_time = 0;
        unsigned int planned_calls = 0;
        unsigned int items_per_call = 1;
        // items_pending is only non-zero when we're doing the small
        // steps. It doesn't make sense to take it into account for +=.
        unsigned int items_pending = 0;

        typedef std::map<std::string, step_time> steps_t;
        steps_t steps;


        void set_total_ncalls(int t) { _total_ncalls = t; }
        unsigned int total_ncalls() const { return std::abs(_total_ncalls); }
        bool is_transition_level() const { return _total_ncalls < 0; }

        inline double projected_time() const {
            return ncalled ? real * total_ncalls() / (double) ncalled : NAN;
        }
        bool is_hot() const { return real < 0; }
        void heat_up(double wct = wct_seconds()) { real -= wct; }
        void cool_down(double wct = wct_seconds()) { real += wct; }

        std::ostream& serialize(std::ostream& os) const {/*{{{*/
            os << " " << name;
            if (real < 0) {
                double delta = -(real + wct_seconds());
                if (!(delta < 0)) delta = -0.001;
                os << " " << delta;
            } else {
                os << " " << real;
            }
            os  << " " << ncalled
                << " " << _total_ncalls
                << " " << planned_time
                << " " << planned_calls
                << " " << items_per_call
                << " " << items_pending;
            os << "\n" << steps.size() << "\n";
            for(auto const & x : steps) {
                os << " " << x.first;
                x.second.serialize(os);
            }
            os << "\n";
            return os;
        }/*}}}*/
        std::istream& unserialize(std::istream& is) {
            step_time res(is);
            if (!is) return is;
            *this = std::move(res);
            return is;
        }
        protected:
        step_time(std::istream& is) {/*{{{*/
            is  >> name;
            double delta;
            is >>  delta;
            if (delta < 0) {
                real = -(wct_seconds() + delta);
            } else {
                real = delta;
            }
            is  >> ncalled
                >> _total_ncalls
                >> planned_time
                >> planned_calls
                >> items_per_call
                >> items_pending;
            size_t n;
            for(is >> n ; is && n-- ; ) {
                std::string k;
                is >> k;
                auto it = steps.find(k);
                if (it != steps.end()) {
                    it->second = step_time(is);
                } else {
                    steps.emplace(k, step_time(is));
                }
            }
        }/*}}}*/
        public:

        step_time(std::string const & name) : name(name) {}
        step_time(step_time const &) = default;
        step_time& operator=(step_time const &) = default;
        step_time& operator=(step_time &&) = default;

        step_time const * find_pending_smallstep() const {/*{{{*/
            if (items_pending) return this;
            for(auto const & x : steps) {
                auto y = x.second.find_pending_smallstep();
                if (y) return y;
            }
            return nullptr;
        }/*}}}*/
        bool has_pending_smallsteps() const {
            return find_pending_smallstep() != nullptr;
        }

        std::ostream& debug_print(std::ostream& s, std::string indent) const;
        step_time & operator+=(step_time const & x);
    };/*}}}*/

    struct function_with_input_size {/*{{{*/
        std::string func;
        unsigned int inputsize;
        function_with_input_size(std::string const & func, unsigned int inputsize)
            : func(func)
            , inputsize(inputsize)
        {
        }
        std::ostream& serialize(std::ostream& os) const {
            return os << " " << func
                      << " " << inputsize;
        }
        std::istream& unserialize(std::istream& is) {
            function_with_input_size res(is);
            if (!is) return is;
            *this = std::move(res);
            return is;
        }
        protected:
        function_with_input_size(std::istream& is) {
            is >> func >> inputsize;
        }
        friend struct level_stats;
        friend class tree_stats;
        public:
        bool operator<(function_with_input_size const& a) const {
            if (func < a.func) return true;
            if (func > a.func) return false;
            return lingen_round_operand_size(inputsize) < lingen_round_operand_size(a.inputsize);
        }
        bool operator==(function_with_input_size const& a) const {
            return func == a.func && lingen_round_operand_size(inputsize) == lingen_round_operand_size(a.inputsize);
        }
    };/*}}}*/

    struct running_stats : public step_time {/*{{{*/
        /* inputsize is tracked only for very top level step_time
         * objects. It's (currently) unused for small steps */
        unsigned int inputsize;
        bool leaf;
        double time_children;
        std::vector<steps_t::iterator> nested_substeps;
        running_stats(std::string const & name, unsigned int inputsize, bool leaf)
            : step_time(name)
            , inputsize(inputsize)
            , leaf(leaf)
            , time_children(0)
        { }
        std::ostream& serialize(std::ostream& os /* , curstack_t::const_iterator place */) const {
            step_time::serialize(os)
                    << " " << inputsize
                    << " " << leaf
                    << " " << time_children;
            /* not sure we ever need to serialize that, in fact */
            ASSERT_ALWAYS(nested_substeps.empty());
#if 0
            os << nested_substeps.size();
            for(auto const & x : nested_substeps) {
                /* serializing an iterator is tricky, as we need to know
                 * the parent structure...
                 */
                auto const & S = place->second.steps;
                ASSERT_ALWAYS(S.find(x->first) != S.end());
                os << x->first;
            }
#endif
            return os;
        }
        std::istream& unserialize(std::istream& is /* , curstack_t::iterator place */) {
            running_stats res(is);
            if (!is) return is;
            *this = std::move(res);
            return is;
        }
        step_time& current_substep() { return nested_substeps.back()->second; }
        step_time const & current_substep() const { return nested_substeps.back()->second; }
        protected:
        running_stats(std::istream& is /*, curstack_t::iterator place */) : step_time(is) {
            is >> inputsize
               >> leaf
               >> time_children;
#if 0
            size_t n;
            is >> n;
            for( ; n-- ; ) {
                /* serializing an iterator is tricky, as we need to know
                 * the parent structure...
                 */
                auto & S = place->second.steps;
                steps_t::key_type x;
                is >> x;
                auto it = S.find(x);
                ASSERT_ALWAYS(it != S.end());
                nested_substeps.push_back(it);
            }
#endif
        }
        friend struct level_stats;
        friend class tree_stats;
    };/*}}}*/

    struct function_stats : public step_time {/*{{{*/
        unsigned int min_inputsize = UINT_MAX;
        unsigned int max_inputsize = 0;
        function_stats & operator+=(running_stats const&);
        function_stats (running_stats const& r) : step_time(r) {
            ASSERT_ALWAYS(!r.is_hot());
            min_inputsize = max_inputsize = r.inputsize;
        }
        std::ostream& serialize(std::ostream& os) const {
            return step_time::serialize(os)
                << " " << min_inputsize
                << " " << max_inputsize
                << "\n";
        }
        std::istream& unserialize(std::istream& is) {
            function_stats res(is);
            if (!is) return is;
            *this = std::move(res);
            return is;
        }
        protected:
        function_stats(std::istream& is) : step_time(is) {
            is >> min_inputsize >> max_inputsize;
        }
        friend struct level_stats;
    };/*}}}*/

    struct level_stats : public std::map<function_with_input_size, function_stats> {/*{{{*/
        /* compute the projected time at this level */
        double projected_time();
        double last_printed_projected_time = 0;
        std::ostream& serialize(std::ostream& os) const {
            os << "\n" << size() << "\n";
            for(auto const & x : *this) {
                x.first.serialize(os);
                x.second.serialize(os);
            }
            return os;
        }
        std::istream& unserialize(std::istream& is) {
            level_stats res(is);
            if (!is) return is;
            *this = std::move(res);
            return is;
        }
        level_stats() = default;
        protected:
        level_stats(std::istream& is) {
            size_t n;
            for(is >> n ; is && n-- ; ) {
                key_type K(is);
                mapped_type M(is);
                emplace(std::move(K), std::move(M));
            }
        }
        friend class tree_stats;
    };/*}}}*/

    std::vector<level_stats> levels;
    curstack_t curstack;

    unsigned int tree_total_breadth = 0;
    double last_print_time = 0;
    std::pair<unsigned int, unsigned int> last_print_position { 0,0 };
    private:
    unsigned int depth = 0;
    unsigned int transition_levels_in_depth = 0;
    public:
    std::ostream& serialize(std::ostream& os) const {
        os << "\n" << levels.size() << "\n";
        for(auto const & x : levels) x.serialize(os);
        os << "\n" << curstack.size() << "\n";
        for(auto x = curstack.cbegin() ; x != curstack.cend() ; ++x) {
            x->first.serialize(os);
            os << " ";
            x->second.serialize(os /* , x */);
        }
        os << " " << tree_total_breadth
            << " " << depth
            << " " << transition_levels_in_depth;
        os << "\n";
        return os;
    }
    tree_stats() = default;
    protected:
    tree_stats(std::istream& is) {
        size_t n;
        for(is >> n ; is && n-- ; ) levels.push_back(level_stats(is));
        for(is >> n ; is && n-- ; ) {
            function_with_input_size K(is);
            running_stats M(is);
            curstack.emplace_back(std::move(K), std::move(M));
        }
        is  >> tree_total_breadth
            >> depth
            >> transition_levels_in_depth;
    }

    public:
    std::istream& unserialize(std::istream& is) {
        tree_stats res(is);
        if (!is) return is;
        *this = std::move(res);
        return is;
    }

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
        // we're in a dtor, hence nothrow, yet we have
        // ASSERT_ALWAYS...
        // coverity[exn_spec_violation]
        ~sentinel() {
            stats.leave();
        }
    };
    struct transition_sentinel {
        tree_stats & stats;
        transition_sentinel(transition_sentinel const&) = delete;
        transition_sentinel(tree_stats & stats,
                std::string const & func, unsigned int inputsize, unsigned int total_ncalls, bool leaf = false)
            : stats(stats) { stats.transition_levels_in_depth++; stats.enter(func, inputsize, -total_ncalls, leaf); }
        // we're in a dtor, hence nothrow, yet we have
        // ASSERT_ALWAYS...
        // coverity[exn_spec_violation]
        ~transition_sentinel() {
            stats.leave();
            stats.transition_levels_in_depth--;
        }
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
        // coverity[exn_spec_violation]
        ~smallstep_sentinel() { stats.end_smallstep(); }
    };
    bool local_smallsteps_done(bool compulsory = false) const;


    void begin_plan_smallstep(std::string const & func, weighted_double const &);
    void begin_plan_smallstep_microsteps(std::string const & func);
    void end_plan_smallstep();
    struct plan_smallstep_sentinel {
        tree_stats & stats;
        plan_smallstep_sentinel(plan_smallstep_sentinel const&) = delete;
        plan_smallstep_sentinel(tree_stats & stats,
                std::string const & func, weighted_double const & theory)
            : stats(stats) { stats.begin_plan_smallstep(func, theory); }
        // coverity[exn_spec_violation]
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
    std::ostream& debug_print(std::ostream&) const;
};

inline std::ostream& operator<<(std::ostream& os, tree_stats const & a)
{
    return a.serialize(os);
}
inline std::istream& operator>>(std::istream& is, tree_stats & a)
{
    return a.unserialize(is);
}

#endif	/* TREE_STATS_HPP_ */
