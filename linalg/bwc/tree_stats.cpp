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

double tree_stats::level_stats::projected_time(unsigned int total_breadth, unsigned int trimmed_breadth)
{
    // return total_breadth * spent / sum_inputsize;
    unsigned int sum_inputsize = 0;
    // unsigned int ncalled = 0;
    trimmed_here = 0;
    for(auto const & x : *this) {
        function_stats const& F(x.second);
        sum_inputsize += F.sum_inputsize;
        trimmed_here += F.trimmed;
    }
    /* Now count the contribution of each sub-function */
    double contrib = 0;
    for(auto & x : *this) {
        function_stats & F(x.second);
        double r = (double) (total_breadth - trimmed_breadth) / sum_inputsize;
        ASSERT_ALWAYS(sum_inputsize <= (total_breadth - trimmed_breadth));
        F.projected_calls = round(r * F.ncalled);
        F.projected_time = r * F.real;
        contrib += F.projected_time;
    }
    return contrib;
}

std::ostream& tree_stats::recursively_print_substeps_at_depth(
        std::ostream & os,
        step_time::steps_t const & FS,
        unsigned int ncalled,
        unsigned int projected_calls,
        int nesting)
{
    std::string prefix(2 + nesting, ' ');

    for(auto const & y : FS) {
        std::string const & func(y.first);
        double t = y.second.real;
        double th = y.second.theoretical;
        unsigned int n = y.second.ncalled;
        os << prefix
           << "(" << func;
        if (n) {
            unsigned int projected_n = n * projected_calls / ncalled;
            os << " "
               << fmt::sprintf("%u/%u", n, projected_n);
        }
        if (n) {
            if (th)
                os << fmt::sprintf(" [%.1f%%]", 100.0*t / th);
            os << " " << fmt::sprintf("%.2g -> %.1f",
                    t / n, t * (double) projected_calls / ncalled);
        } else if (th) {
            /* At this point we don't know in advance how many calls of
             * this small step we will have per level. This should be
             * added to plan_smallstep XXX
             */
            int call_ratio = 1;
            os << " " << fmt::sprintf("%.2g -> %.1f",
                    th, th * call_ratio * (double) projected_calls);
        }
        os << ")\n";
        recursively_print_substeps_at_depth(os, y.second.steps,
                ncalled, projected_calls, nesting+1);
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
    unsigned int tree_trimmed_breadth = 0;

    std::ostringstream os;

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
        double pt = u.projected_time(tree_total_breadth, tree_trimmed_breadth);
        u.last_printed_projected_time = pt;

        char mixedlevel_code = (u.size() <= 1) ? '\0' : 'a';

        for(auto x = u.cbegin(); x != u.cend(); ++x, mixedlevel_code++) {
            string const& func(x->first);
            function_stats const& F(x->second);
            sum += F.projected_time;
            time_to_go += F.projected_time - F.real;
            /* when printing at this level, we must take into account the
             * running_stats pointed to by "extra" as well.
             */
            os << k;
            if (mixedlevel_code) os << mixedlevel_code;
            os << " " << fmt::sprintf("[%u-%u, %s]", F.min_inputsize, F.max_inputsize, func)
               << " " << fmt::sprintf("%u/%u", F.ncalled, F.projected_calls)
               << " " << fmt::sprintf("%.2g -> %.1f",
                       F.real / F.ncalled, F.projected_time)
               << " " << fmt::sprintf(" (total: %.1f)", sum)
               << "\n";

            step_time FS = F;
            if (wip && wip->first == func) {
                FS += wip->second;
                wip = nullptr;
            }

            recursively_print_substeps_at_depth(os, FS.steps, F.ncalled, F.projected_calls, 0);
        }

        /* Since the running stats are a differetn type, with different
         * semantics, we must do things a bit differently.
         */
        if (wip) {
            std::string const & func(wip->first);
            running_stats const & r(wip->second);
            unsigned int exp_ncalls = round((double) tree_total_breadth / r.inputsize);
            double level_th = 0;
            for(auto const & y : r.steps) {
                double t = y.second.real;
                double th = y.second.theoretical;
                unsigned int n = y.second.ncalled;
                level_th += n ? t : th;
            }
            sum += level_th * exp_ncalls;
            os << fmt::sprintf("%u * [%u, %s] 0/%u %.2g -> %.1f (total: %.1f)\n",
                    k,
                    r.inputsize,
                    func.c_str(),
                    exp_ncalls,
                    level_th, level_th * exp_ncalls,
                    sum);
            recursively_print_substeps_at_depth(os, r.steps, 1, exp_ncalls, 0);
        }

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

void tree_stats::enter(const char * func, unsigned int inputsize, bool recurse)
{
    int rank;
    if (depth == 0)
        tree_total_breadth = inputsize;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank) { ++depth; return; }
    running_stats s(inputsize, !recurse);
    s.heat_up();
    // ASSERT_ALWAYS(curstack.empty() || !curstack.back().second.in_substep());
    curstack.push_back({func, s});
    ++depth;
    ASSERT_ALWAYS(depth == curstack.size());
    if (curstack.size() > levels.size())
        levels.insert(levels.end(), curstack.size() - levels.size(), level_stats());
}

tree_stats::function_stats& tree_stats::function_stats::operator+=(tree_stats::running_stats const & s)
{
    (tree_stats::step_time&)*this += (tree_stats::step_time const&) s;
    if (s.inputsize < min_inputsize) min_inputsize = s.inputsize;
    if (s.inputsize > max_inputsize) max_inputsize = s.inputsize;
    sum_inputsize += s.inputsize;
    if (s.leaf) trimmed += s.inputsize;
    return *this;
}

void tree_stats::leave()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank) { --depth; return; }
    double now = wct_seconds();
    auto sback = std::move(curstack.back());
    std::string const & sfunc = sback.first;
    running_stats & s = sback.second;
    curstack.pop_back();

    s.cool_down(now);
    if (!curstack.empty())
        curstack.back().second.time_children += s.real;
    s.real -= s.time_children;
    unsigned int level = --depth;
    ASSERT_ALWAYS(depth == curstack.size());
    ASSERT_ALWAYS(level < levels.size());

    /* merge our running stats into the level_stats */
    function_stats & F(levels[level][sfunc]);
    F += s;
    F.ncalled++;

    /* Is it any useful to print something new ? */
    if (now < last_print_time + 2) return;

    int needprint = 0;
    unsigned int trimmed_breadth = 0;
    for(unsigned int k = 0 ; !needprint && k < levels.size() ; k++) {
        level_stats & u(levels[k]);
        double t = u.projected_time(tree_total_breadth, trimmed_breadth);
        double t0 = u.last_printed_projected_time;
        needprint = (t < 0.98 * t0) || (t > 1.02 * t0);
        trimmed_breadth += u.trimmed_here;
    }

    if (!needprint)
        return;

    last_print_time = now;
    last_print_position = make_pair(level, F.sum_inputsize);

    print(level);
}

void tree_stats::final_print()
{
    ASSERT_ALWAYS(depth == 0);
    if (last_print_position != make_pair(0u, tree_total_breadth))
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

void tree_stats::plan_smallstep(std::string const & func, double theory)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank) return;
    ASSERT_ALWAYS(!curstack.empty());
    running_stats& s(curstack.back().second);
    s.steps[func].theoretical += theory;
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
    running_stats::steps_t::iterator it =where->insert({func, step_time()}).first;
    s.nested_substeps.push_back(it);
    s.current_substep().ncalled += ncalls;
    s.current_substep().heat_up();
}

void tree_stats::end_smallstep()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank) return;
    ASSERT_ALWAYS(!curstack.empty());
    running_stats& s(curstack.back().second);
    s.current_substep().cool_down();
    // s.current_substep().ncalled++;
    s.nested_substeps.pop_back();
}
