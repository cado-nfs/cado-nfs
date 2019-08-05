#include "cado.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <cctype>
#include <algorithm>
#include <sstream>

#include "fmt/format.h"
#include "fmt/printf.h"
#include "portability.h"
#include "select_mpi.h"
#include "utils.h"
#include "tree_stats.hpp"

using namespace std;

size_t lingen_round_operand_size(size_t x, int bits) {/*{{{*/
    /* round x up to the next size that has all but its six most significant
     * bits set to 0.
     */
    if (x == 0) return x;
    x -= 1;
    size_t y = x >> bits;
    for(int i = 1 ; y ; y >>= i) x |= y;
    x += 1;
    return x;
}/*}}}*/

int tree_stats::max_nesting = 0;

void tree_stats::interpret_parameters(cxx_param_list & pl)
{
    param_list_parse_int(pl, "tree_stats_max_nesting", &max_nesting);
}


void tree_stats::declare_usage(cxx_param_list & pl)
{
    param_list_decl_usage(pl, "tree_stats_max_nesting", "max nesting level of small steps to display within plingen");
}

double tree_stats::level_stats::projected_time()
{
    /* Now count the contribution of each sub-function */
    double contrib = 0;
    for(auto & x : *this)
        contrib += x.second.projected_time();
    return contrib;
}

std::ostream& tree_stats::recursively_print_substeps_at_depth(
        std::ostream & os,
        step_time::steps_t const & FS,
        unsigned int ncalled,
        unsigned int total_ncalls,
        int nesting)
{
    std::string prefix(2 + nesting, ' ');

    if (nesting >= max_nesting)
        return os;

    for(auto const & y : FS) {
        std::string const & func(y.first);
        unsigned int n = y.second.ncalled;
        double t = y.second.real;
        unsigned int items_per_call = y.second.items_per_call;
        double th = y.second.planned_time;
        double th_n = y.second.planned_calls;
        os << prefix
           << "(" << func;
        /* For a function that has been entered k times (each time going
         * through a tree of small steps), we have the following values
         * for the different step_time objects (here, t0 is some
         * theoretical value, t1 is a measured value)
         *
         * top-level (parent of the function_stats) object, **NOT**
         * printed by this function:
         *      ncalled = k
         *      real = k * t1 (without children)
         *      total_ncalls = should be a power of two.
         *      projected_time = total_ncalls/ncalled * real
         *
         * small steps:
         *      items_per_call = X (X can be the contribution of multiple
         *          calls to begin_smallstep, with ncalls arguments summing
         *          up to X).
         *      ncalled = k * X
         *      planned_time = k * X * t0
         *      real = k * X * t1
         */
        size_t projected_n = (size_t) total_ncalls * items_per_call;
        os << " "
            << fmt::sprintf("%u/%zu", n, projected_n);
        if (n) {
            if (th)
                os << fmt::sprintf(" [%.1f%%]", 100.0 * (t/n) / (th/th_n));
            os << " " << fmt::sprintf("%.2g -> %.1f",
                    t / n, t / n * (double) total_ncalls);
        } else if (th) {
            /* Since n == 0, begin_smallstep has never been called.
             * We only had plan_smallstep, and we had it only once,
             * so that th corresponds t 
             * that th_n corresponds to one call of the function, not
             * more */
            ASSERT_ALWAYS(th_n == 1);
            os << fmt::sprintf(" %.2g -> %.1f",
                    (th/th_n), (th/th_n) * (double) total_ncalls);
        }
        os << ")\n";
        recursively_print_substeps_at_depth(os, y.second.steps,
                ncalled, total_ncalls, nesting+1);
    }
    return os;
}



void tree_stats::print(unsigned int)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank) return;

    double sum = 0;
    double time_to_go = 0;

    std::ostringstream os;

    int total_transition_levels = 0;
    for(unsigned int k = 0 ; k < levels.size() ; k++) {
        /* At depth k, we have:
         *  - cold timers for several functions, in levels[k]
         *  - if k < curstack.size(), one hot timer, which may actually
         *  contain info that is really interesting to report, so we'll
         *  strive to do this. They will maybe 
         */
        level_stats & u(levels[k]);

        curstack_t::value_type const * wip = nullptr;
        if (k < curstack.size())
            wip = &curstack[k];

        /* does this have any use ? XXX */
        // if (k > level && u.empty()) break;

        /* Compute projected time for this level based on the total
         * breadth, and the calls which we had so far.
         */
        double pt = u.projected_time();
        u.last_printed_projected_time = pt;

        char mixedlevel_code = (u.size() <= 1) ? '\0' : 'a';

        bool all_functions_are_transitions = true;

        for(auto x = u.cbegin(); x != u.cend(); ++x, mixedlevel_code++) {
            function_with_input_size const& fi(x->first);
            function_stats const& F(x->second);
            sum += F.projected_time();
            time_to_go += F.projected_time() - F.real;
            /* when printing at this level, we must take into account the
             * running_stats pointed to by "extra" as well.
             */
            if (!F.is_transition_level()) {
                os << k-total_transition_levels;
                if (mixedlevel_code) os << mixedlevel_code;
                all_functions_are_transitions = false;
            } else {
                os << "--";
            }
            os << " " << fmt::sprintf("[%u-%u, %s]", F.min_inputsize, F.max_inputsize, fi.func)
               << " " << fmt::sprintf("%u/%u", F.ncalled, F.total_ncalls())
               << " " << fmt::sprintf("%.2g -> %.1f",
                       F.real / F.ncalled, F.projected_time())
               << " " << fmt::sprintf(" (total: %.1f)", sum)
               << "\n";

            step_time FS = F;
            if (wip && wip->first == fi) {
                FS += wip->second;
                wip = nullptr;
            }

            recursively_print_substeps_at_depth(os, FS.steps, F.ncalled, F.total_ncalls(), 0);
        }

        /* Since the running stats are a different type, with different
         * semantics, we must do things a bit differently.
         */
        if (wip) {
            function_with_input_size const & fi(wip->first);
            running_stats const & r(wip->second);
            double level_th = 0;
            for(auto const & y : r.steps) {
                double t = y.second.real;
                double th = y.second.planned_time;
                unsigned int n = y.second.ncalled;
                level_th += n ? t : th;
                time_to_go += th * (r.total_ncalls() - n);
            }
            if (!r.is_transition_level()) {
                os << k-total_transition_levels;
                all_functions_are_transitions = false;
            } else {
                os << "--";
            }
            sum += level_th * r.total_ncalls();
            os << fmt::sprintf("* [%u, %s] 0/%u %.2g -> %.1f (total: %.1f)\n",
                    r.inputsize,
                    fi.func,
                    r.total_ncalls(),
                    level_th, level_th * r.total_ncalls(),
                    sum);
            recursively_print_substeps_at_depth(os, r.steps, 1, r.total_ncalls(), 0);
        }
        if (all_functions_are_transitions) total_transition_levels++;
    }

    puts(os.str().c_str());

    /* Note that time_to_go is only relative to the levels for which we
     * have got at least one data point */

    {
        /* print ETA */
        time_t eta[1];
        char eta_string[32] = "not available yet\n";
        *eta = wct_seconds() + time_to_go;
#ifdef HAVE_CTIME_R
        ctime_r(eta, eta_string);
#else
        strncpy(eta_string, ctime(eta), sizeof(eta_string));
#endif
        unsigned int s = strlen(eta_string);
        for( ; s && isspace((int)(unsigned char)eta_string[s-1]) ; eta_string[--s]='\0') ;

        printf("lingen ETA: %s\n", eta_string);
    }
}

void tree_stats::enter(std::string const & func, unsigned int inputsize, int total_ncalls, bool leaf)
{
    int rank;
    if (depth == 0)
        tree_total_breadth = inputsize;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank) { ++depth; return; }
    running_stats s(func, inputsize, leaf);
    s.set_total_ncalls(total_ncalls);
    s.heat_up();
    // ASSERT_ALWAYS(curstack.empty() || !curstack.back().second.in_substep());
    curstack.push_back({{func, inputsize}, s});
    ++depth;
    // we used to have the following, but now we allow ourselves to
    // insert empty shells between levels (so that intermediate
    // operations such as the mpi gather and scatter steps can be
    // identified as small steps of something...)
    // ASSERT_ALWAYS(depth == curstack.size());
    //
    // Bottom line, we can only assert the weaker:
    ASSERT_ALWAYS(depth <= curstack.size());
    if (curstack.size() > levels.size())
        levels.insert(levels.end(), curstack.size() - levels.size(), level_stats());
}

tree_stats::function_stats& tree_stats::function_stats::operator+=(tree_stats::running_stats const & s)
{
    (tree_stats::step_time&)*this += (tree_stats::step_time const&) s;
    if (s.inputsize < min_inputsize) min_inputsize = s.inputsize;
    if (s.inputsize > max_inputsize) max_inputsize = s.inputsize;
    return *this;
}

void tree_stats::leave()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank) { --depth; return; }
    double now = wct_seconds();
    auto sback = std::move(curstack.back());
    function_with_input_size const & sfunc = sback.first;
    running_stats & s = sback.second;
    curstack.pop_back();

    s.cool_down(now);
    if (!curstack.empty())
        curstack.back().second.time_children += s.real;
    s.real -= s.time_children;
    unsigned int level = --depth;
    ASSERT_ALWAYS(depth == curstack.size());
    ASSERT_ALWAYS(level < levels.size());
    ASSERT_ALWAYS(!s.has_pending_smallsteps());

    /* merge our running stats into the level_stats */
    function_with_input_size fi = sfunc;
    //function_stats & F(levels[level][fi]);
    // F += s;
    auto itb = levels[level].emplace(fi, s);
    function_stats & F(itb.first->second);
    if (!itb.second) F += s;
    F.planned_calls++;  /* just for consistency */
    F.ncalled++;

    /* Is it any useful to print something new ? */
    if (now < last_print_time + 2) return;

    int needprint = 0;
    for(unsigned int k = 0 ; !needprint && k < levels.size() ; k++) {
        level_stats & u(levels[k]);
        double t = u.projected_time();
        double t0 = u.last_printed_projected_time;
        needprint = (t < 0.98 * t0) || (t > 1.02 * t0);
    }

    if (!needprint) return;

    last_print_time = now;
    last_print_position = make_pair(level, F.ncalled);

    print(level);
}

void tree_stats::final_print()
{
    ASSERT_ALWAYS(depth == 0);
    if (last_print_position != make_pair(0u, 1u))
        print(0);
    {
        /* print ETA */
        time_t eta[1];
        char eta_string[32];
        *eta = wct_seconds();
#ifdef HAVE_CTIME_R
        ctime_r(eta, eta_string);
#else
        strncpy(eta_string, ctime(eta), sizeof(eta_string));
#endif
        unsigned int s = strlen(eta_string);
        for( ; s && isspace((int)(unsigned char)eta_string[s-1]) ; eta_string[--s]='\0') ;

        printf("lingen done at: %s\n", eta_string);
    }
}

void tree_stats::begin_plan_smallstep(std::string const & func, weighted_double const & theory)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank) return;
    ASSERT_ALWAYS(!curstack.empty());
    running_stats& s(curstack.back().second);

    step_time::steps_t * where = &s.steps;
    if (!s.nested_substeps.empty())
        where = & s.current_substep().steps;
    running_stats::steps_t::iterator it =where->insert({func, step_time(func)}).first;
    s.nested_substeps.push_back(it);

    step_time & S(it->second);
    S.set_total_ncalls(s.total_ncalls());
    S.items_per_call = theory.n;
    S.planned_time += theory.t * theory.n;
    ASSERT_ALWAYS(S.ncalled == S.planned_calls);
    S.planned_calls++;
}
void tree_stats::end_plan_smallstep()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank) return;
    ASSERT_ALWAYS(!curstack.empty());
    running_stats& s(curstack.back().second);
    s.nested_substeps.pop_back();
}

void tree_stats::begin_smallstep(std::string const & func, unsigned int ncalls)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank) return;
    ASSERT_ALWAYS(!curstack.empty());
    running_stats& s(curstack.back().second);
    // At first thought, we never have two substeps of the same name at a given
    // level. Alas, this is not always right. One example is at the mpi
    // threshold in lingen. We have 2 gather and 2 scatter steps.
    // In effect, the way we're proceeding sets the substep pointer to the
    // current smallstep, so that it's counted as "in progress", although
    // that does not properly acknowledge the fact that one instance of
    // that sub-step has already run. In that situation:
    //  - we're failing to list something for which we do have info.
    //  - at the end of the day, the "number of calls" will consider both
    //    instances as one single call.
    //
    // This is considered harmless, given that:
    //  - the cut-off point (MPI threshold, in our current use) is early
    //    enough that timings are displayed.
    //  - it is possible to embed in the substep name some info hinting
    //    at the fact that we have 2 substeps (e.g. "foo(1+2)" or
    //    "foo(L+R)").
    //
    // ASSERT_ALWAYS(ssi.second);
    step_time::steps_t * where = &s.steps;
    if (!s.nested_substeps.empty())
        where = & s.current_substep().steps;
    running_stats::steps_t::iterator it =where->insert({func, step_time(func)}).first;
    it->second.set_total_ncalls(s.total_ncalls());
    s.nested_substeps.push_back(it);
    step_time & S(s.current_substep());
    S.items_pending += ncalls;
    ASSERT_ALWAYS(S.items_pending <= S.items_per_call);
    S.heat_up();
}

void tree_stats::end_smallstep()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank) return;
    ASSERT_ALWAYS(!curstack.empty());
    running_stats& s(curstack.back().second);
    step_time & S(s.current_substep());
    S.cool_down();
    ASSERT_ALWAYS(S.items_pending <= S.items_per_call);
    if (S.items_pending == S.items_per_call) {
        S.ncalled++;
        S.items_pending = 0;
    }
    s.nested_substeps.pop_back();
}

void tree_stats::skip_smallstep(std::string const & func, unsigned int ncalls)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank) return;
    ASSERT_ALWAYS(!curstack.empty());
    running_stats& s(curstack.back().second);
    step_time::steps_t * where = &s.steps;
    if (!s.nested_substeps.empty())
        where = & s.current_substep().steps;
    running_stats::steps_t::iterator it = where->insert({func, step_time(func)}).first;
    it->second.set_total_ncalls(s.total_ncalls());
    s.nested_substeps.push_back(it);
    step_time & S(s.current_substep());
    S.items_pending += ncalls;
    ASSERT_ALWAYS(S.items_pending <= S.items_per_call);
    if (S.items_pending == S.items_per_call) {
        S.ncalled++;
        S.items_pending = 0;
    }
    s.nested_substeps.pop_back();
}

bool tree_stats::local_smallsteps_done() const
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank) return true;
    ASSERT_ALWAYS(!curstack.empty());
    running_stats const & s(curstack.back().second);
    step_time::steps_t const * where = &s.steps;
    if (!s.nested_substeps.empty())
        where = & s.current_substep().steps;
    for(auto const & s : *where) {
        if (s.second.items_pending) return false;
    }
    return true;
}
