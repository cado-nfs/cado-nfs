#include "cado.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <cctype>

#include "select_mpi.h"
#include "tree_stats.hpp"
#include "timing.h"     // wct_seconds

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
        // ncalled += ncalled;
    }
    // unsigned int expected_total_calls = total_breadth / (sum_inputsize / (double) ncalled);
    /* Now count the contribution of each sub-function */
    double contrib = 0;
    for(auto & x : *this) {
        function_stats & F(x.second);
        // expected_calls = expected_total_calls * (double) ncalled / ncalled;
        double r = (double) (total_breadth - trimmed_breadth) / sum_inputsize;
        ASSERT_ALWAYS(sum_inputsize <= (total_breadth - trimmed_breadth));
        // double contrib = n * spent / ncalled;
        // total_breadth / sum_inputsize * (double) ncalled * (double) ncalled / ncalled * spent / ncalled;
        F.projected_calls = round(r * F.ncalled);
        F.projected_time = r * F.spent;

        contrib += F.projected_time;
    }
    return contrib;
}


void tree_stats::print(unsigned int level)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank) return;

    double sum = 0;
    double time_to_go = 0;
    int nstars=0;
    int nok=0;
    double firstok = 0;
    double complement = 0;
    unsigned int tree_trimmed_breadth = 0;
    const char * prefix = draft ? "##DRAFT## ":"";

    for(unsigned int k = 0 ; k < stack.size() ; k++) {
        level_stats & u(stack[k]);

        if (k > level && u.empty())
            break;

        if (!u.empty()) {
            /* Compute projected time for this level based on the total
             * breadth, and the calls which we had so far.
             */
            double pt = u.projected_time(tree_total_breadth, tree_trimmed_breadth);

            if (nstars && nok == 0) {
                firstok = pt;
            } else if (nstars && nok == 1) {
                double ratio = firstok / pt;
                complement = ratio * firstok * (pow(ratio, nstars) - 1) / (ratio - 1);
            }
            nok++;

            char code[2]={'\0', '\0'};
            if (u.size() > 1) code[0] = 'a';
            for(auto const & x : u) {
                string const& key(x.first);
                function_stats const& F(x.second);
                sum += F.projected_time;
                time_to_go += F.projected_time - F.spent;
                printf("%s%u%s [%u-%u, %s] %u/%u %.2g -> %.1f (total: %.1f wct)\n",
                        prefix,
                        k, (const char*) code,
                        F.min_inputsize, F.max_inputsize,
                        key.c_str(),
                        F.ncalled, F.projected_calls,
                        F.spent / F.ncalled, F.projected_time, sum);
                if (code[0]) code[0]++;
                for(auto const & y : F.small_steps) {
                    double t = y.second.real + y.second.artificial;
                    unsigned int n = F.ncalled;
                    if (k < curstack.size()) {
                        running_stats const& r(curstack[k]);
                        if (r.func == key) {
                            /* shall we count one extra call which has alreay
                             * been done ?? */
                            auto z = r.small_steps.find(y.first);
                            if (z != r.small_steps.end() && &(z->second) != r.substep) {
                                /* also count the artificial time */
                                t += z->second.real;
                                t += z->second.artificial;
                                n++;
                            }
                        }
                    }
                    printf("%s   (%s %u/%u %.2g -> %.1f)\n",
                            prefix,
                            y.first.c_str(),
                            n, F.projected_calls,
                            t / n,
                            t * (double) F.projected_calls / n);
                }
            }
            tree_trimmed_breadth += u.trimmed_here;
            u.last_printed_projected_time = pt;
        } else {
            /* We're not in a situation where we can control exactly
             * what's going to happen, because none of the running_stats
             * has issued a leave() command. The only thing we have is
             * the currently running stats, which nevertheless has some
             * data which can be printed.
             */
            ASSERT_ALWAYS(k < curstack.size());
            running_stats const& r(curstack[k]);
            unsigned int exp_ncalls = round((double) tree_total_breadth / r.inputsize);
            printf("%s%u * [%u, %s] 0/%u\n",
                    prefix,
                    k,
                    r.inputsize,
                    r.func.c_str(),
                    exp_ncalls);
            for(auto const & y : r.small_steps) {
                unsigned int n = 0;
                ASSERT_ALWAYS(&(y.second) != r.substep);
                /* also count the artificial time */
                double t = y.second.real + y.second.artificial;
                n++;
                printf("%s   (%s %u/%u %.2g -> %.1f)\n",
                        prefix,
                        y.first.c_str(),
                        n, exp_ncalls,
                        t / n,
                        t * (double) exp_ncalls / n);
            }
            nstars++;
            continue;
        }

    }

    if (nstars && nok >= 2) {
        printf("%sexpected time for levels 0-%u: %.1f (total: %.1f)\n",
                prefix,
                nstars-1, complement, sum + complement);
    }

    /* Note that time_to_go is only relative to the levels for which we
     * have got at least one data point */

    if (draft) {
        /* draft mode means that we're going to start over anyway ! */
        time_to_go = sum;
    }

    {
        /* print ETA */
        time_t eta[1];
        char eta_string[32] = "not available yet\n";
        *eta = wct_seconds() + time_to_go + complement;
#ifdef HAVE_CTIME_R
        ctime_r(eta, eta_string);
#else
        strncpy(eta_string, ctime(eta), sizeof(eta_string));
#endif
        unsigned int s = strlen(eta_string);
        for( ; s && isspace((int)(unsigned char)eta_string[s-1]) ; eta_string[--s]='\0') ;

        if (draft) {
            printf("%slingen expected duration: %f s (ETA from now: %s)\n", prefix, sum + complement, eta_string);
        } else {
            printf("lingen ETA: %s\n", eta_string);
        }
    }
}

void tree_stats::add_artificial_time(double t) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank) return;
    curstack.back().add_artificial_time(t);
}

void tree_stats::enter(const char * func, unsigned int inputsize, bool recurse)
{
    unsigned int trimmed = recurse ? 0 : inputsize;

    int rank;
    if (depth == 0)
        tree_total_breadth = inputsize;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank) { ++depth; return; }
    if (depth == 0) {
        begin = wct_seconds();
    }
    running_stats s;
    s.time_self = -wct_seconds();
    s.func = func;
    s.inputsize = inputsize;
    s.trimmed = trimmed;
    if (!curstack.empty()) {
        ASSERT_ALWAYS(!curstack.back().substep);
    }
    curstack.push_back(s);
    ++depth;
    ASSERT_ALWAYS(depth == curstack.size());
    if (curstack.size() > stack.size())
        stack.insert(stack.end(), curstack.size() - stack.size(), level_stats());
}

/* This returns the time spent on this function for all calls at this
 * level, including the artificial time that has been reported so far by
 * these calls.
 */
double tree_stats::spent_so_far()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank) return 0;
    double now = wct_seconds();
    running_stats s = curstack.back();
    unsigned int level = depth-1;
    ASSERT_ALWAYS(level < stack.size());
    double t = s.time_self + now + s.time_artificial;
    auto fi = stack[level].find(s.func);
    if (fi != stack[level].end()) t += fi->second.spent;
    return t;
}

void tree_stats::leave()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank) { --depth; return; }
    double now = wct_seconds();
    running_stats s = curstack.back();
    curstack.pop_back();

    s.time_self += now;
    if (!curstack.empty())
        curstack.back().time_children += s.time_self;
    // comparing timings is asking for (non-reproducible) trouble.
    // ASSERT_ALWAYS(s.time_children <= s.time_self);
    s.time_self -= s.time_children;
    unsigned int level = --depth;
    ASSERT_ALWAYS(depth == curstack.size());
    ASSERT_ALWAYS(level < stack.size());

    /* merge our running stats into the level_stats */

    function_stats & F(stack[level][s.func]);
    F.ncalled++;
    if (s.inputsize < F.min_inputsize)
        F.min_inputsize = s.inputsize;
    if (s.inputsize > F.max_inputsize)
        F.max_inputsize = s.inputsize;
    F.sum_inputsize += s.inputsize;
    F.trimmed += s.trimmed;
    F.spent += s.time_self + s.time_artificial;
    for(auto const & x : s.small_steps) {
        F.small_steps[x.first] += x.second;
    }

    /* Is it any useful to print something new ? */

    if (now < last_print_time + 2) return;

    int needprint = 0;
    unsigned int trimmed_breadth = 0;
    for(unsigned int k = 0 ; !needprint && k < stack.size() ; k++) {
        level_stats & u(stack[k]);
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

        if (!draft) printf("lingen done at: %s\n", eta_string);
    }
}

void tree_stats::begin_smallstep(const char * func)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank) return;
    ASSERT_ALWAYS(!curstack.empty());
    running_stats& s(curstack.back());
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
    s.substep = &(s.small_steps[func]);
    s.substep->real -= wct_seconds();
}

void tree_stats::end_smallstep()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank) return;
    ASSERT_ALWAYS(!curstack.empty());
    running_stats& s(curstack.back());
    s.substep->real += wct_seconds();
    s.substep = NULL;
}
