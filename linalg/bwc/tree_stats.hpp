#ifndef TREE_STATS_HPP_
#define TREE_STATS_HPP_

#include <string>
#include <vector>
#include <map>
#include <climits>

class tree_stats {
    struct small_step_time {
        double real = 0;
        double artificial = 0;  // only for deprecated draft mode.
        double theoretical=0;   // expected time to the current point
        unsigned int ncalled = 0;
        small_step_time operator+=(small_step_time const & x) {
            real += x.real;
            artificial += x.artificial;
            theoretical += x.theoretical;
            ncalled += x.ncalled;
            return *this;
        }
    };
    struct function_stats {
        unsigned int ncalled = 0;
        unsigned int min_inputsize = UINT_MAX;
        unsigned int max_inputsize = 0;
        unsigned int sum_inputsize = 0;
        unsigned int trimmed = 0;       /* either sum_inputsize or 0 */
        /* The "spent" time also includes artificial time. So does the
         * projected time. */
        double spent = 0;
        double projected_time = 0;
        unsigned int projected_calls = 0;
        std::map<std::string, small_step_time> small_steps;
    };

    struct level_stats : public std::map<std::string, function_stats> {
        double projected_time(unsigned int, unsigned int);
        double last_printed_projected_time;
        unsigned int trimmed_here;       /* trimmed at this level ! */
    };


    struct running_stats {
        std::string func;
        unsigned int inputsize = 0;
        unsigned int trimmed = 0;
        double time_self = 0;
        /* this time must not be subtracted from the parent time (because
         * it was not spent for real), but we must take it into account
         * for projected timings.
         */
        double time_artificial = 0;
        double time_children = 0;
        std::map<std::string, small_step_time> small_steps;
        small_step_time * substep = NULL;
        inline void add_artificial_time(double t) {
            if (substep) substep->artificial += t;
            /* it does get counted inside our timing anyway */
            time_artificial += t;
        }
    };

    std::vector<level_stats> stack;
    std::vector<running_stats> curstack;

    unsigned int tree_total_breadth;
    double last_print_time;
    std::pair<unsigned int, unsigned int> last_print_position { 0,0 };

    double begin;       /* stored as a wct_seconds() return value */


    void print(unsigned int level);

    /* max number of seconds to use as a bound before we consider that
     * the timing is settled for good.
     */
    int draft = 0;
public:
    unsigned int depth = 0;
    inline int is_draft_mode() const { return draft; }
    inline void set_draft_mode(int d) { draft = d; }
    
    void enter(const char * func, unsigned int inputsize, bool recurse = true); 
    void leave();

    void final_print();

    void plan_smallstep(const char * func, double theory=0);
    void begin_smallstep(const char * func);
    void end_smallstep();

    struct sentinel {
        tree_stats & stats;
        sentinel(sentinel const&) = delete;
        sentinel(tree_stats & stats,
                const char * func, unsigned int inputsize, bool recurse = true)
            : stats(stats) { stats.enter(func, inputsize, recurse); }
        ~sentinel() { stats.leave(); }
    };

    void add_artificial_time(double t);

    /* This returns the time spent on this function for all calls at this
     * level, including the artificial time that has been reported so far by
     * these calls.
     */
    double spent_so_far();
};

#endif	/* TREE_STATS_HPP_ */
