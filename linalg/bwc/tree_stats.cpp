#include "cado.h" // IWYU pragma: keep
// IWYU pragma: no_include <ext/alloc_traits.h>
// IWYU pragma: no_include <memory>

#include <cstdio>
#include <cstring>
#include <ctime>
#include <cctype>

#include <map>             // for map<>::iterator, _Rb_tree_iterator, _Rb_tr...
#include <string>          // for operator<<, string, basic_string, operator+
#include <utility>         // for pair, make_pair, move, operator!=
#include <vector>          // for vector, vector<>::value_type
#include <sstream> // IWYU pragma: keep
#include <ostream> // ostream operator<<
#include <stdexcept>     // for runtime_error

#include "fmt/base.h"
#include "fmt/format.h"

#include "params.h"      // for cxx_param_list, param_list_decl_usage, param...
#include "select_mpi.h"
#include "tree_stats.hpp"
#include "timing.h"     // wct_seconds
#include "macros.h"

int tree_stats::max_nesting = 0;

void tree_stats::interpret_parameters(cxx_param_list & pl)
{
    param_list_parse_int(pl, "tree_stats_max_nesting", &max_nesting);
}


void tree_stats::declare_usage(cxx_param_list & pl)
{
    param_list_decl_usage(pl, "tree_stats_max_nesting", "max nesting level of small steps to display within lingen");
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
    std::string const prefix(2 + nesting, ' ');

    if (nesting >= max_nesting)
        return os;

    for(auto const & y : FS) {
        std::string const & func(y.first);
        unsigned int const n = y.second.ncalled;
        double const t = y.second.real;
        unsigned int const items_per_call = y.second.items_per_call;
        double const th = y.second.planned_time;
        double const th_n = y.second.planned_calls;
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
        os << " "
            << fmt::format("{}/{}",
                    items_per_call * (size_t) n,
                    items_per_call * (size_t) total_ncalls);
        if (n) {
            if (th > 0)
                os << fmt::format(" [{:.1f}%%]", 100.0 * (t/n) / (th/th_n));
            os << " " << fmt::format("{:.2g} -> {:.1f}",
                    t / n / items_per_call, t / n * (double) total_ncalls);
        } else if (th > 0) {
            /* Since n == 0, begin_smallstep has never been called.
             * We only had plan_smallstep, and we had it only once,
             * so that th corresponds t 
             * that th_n corresponds to one call of the function, not
             * more */
            ASSERT_ALWAYS(th_n == 1);
            os << fmt::format(" {:.2g} -> {:.1f}",
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
    bool has_untimed = false;

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
        double const pt = u.projected_time();
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
            os << " " << fmt::format("[{}-{}, {}]", F.min_inputsize, F.max_inputsize, fi.func)
               << " " << fmt::format("{}/{}", F.ncalled, F.total_ncalls())
               << " " << fmt::format("{:.2g} -> {:.1f}",
                       F.real / F.ncalled, F.projected_time())
               << " " << fmt::format(" (total: {:.1f} wct)", sum)
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
                double const t = y.second.real;
                double const th = y.second.planned_time;
                unsigned int const n = y.second.ncalled;
                level_th += n ? t : th;
                if (th > 0)
                    time_to_go += th * (r.total_ncalls() - n);
                else
                    has_untimed = true;
            }
            if (!r.is_transition_level()) {
                os << k-total_transition_levels;
                all_functions_are_transitions = false;
            } else {
                os << "--";
            }
            os << fmt::format("* [{}, {}] 0/{}",
                    r.inputsize,
                    fi.func,
                    r.total_ncalls()
                    );
            if (level_th > 0) {
                sum += level_th * r.total_ncalls();
                os << fmt::format(" {:.2g} -> {:.1f} (total: {:.1f} wct)",
                        level_th, level_th * r.total_ncalls(),
                        sum);
            }
            os << "\n";
            recursively_print_substeps_at_depth(os, r.steps, 1, r.total_ncalls(), 0);
        }
        if (all_functions_are_transitions) total_transition_levels++;
    }

    puts(os.str().c_str());

    /* Note that time_to_go is only relative to the levels for which we
     * have got at least one data point */

    {
        /* print ETA */
        char eta_string[32] = "not available yet\n";
        if (!has_untimed) {
            time_t eta[1];
            *eta = wct_seconds() + time_to_go;
#ifdef HAVE_CTIME_R
            ctime_r(eta, eta_string);
#else
            strncpy(eta_string, ctime(eta), sizeof(eta_string));
#endif
        }
        unsigned int s = strlen(eta_string);
        for( ; s && isspace((int)(unsigned char)eta_string[s-1]) ; eta_string[--s]='\0') ;

        fmt::print("lingen ETA: {}\n", eta_string);
    }
}

void tree_stats::enter(std::string const & func, unsigned int inputsize, int total_ncalls, bool leaf)
{
    if (func.find(" ") != std::string::npos)
        abort();
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
    min_inputsize = std::min(min_inputsize, s.inputsize);
    max_inputsize = std::max(max_inputsize, s.inputsize);
    return *this;
}

void tree_stats::leave()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank) { --depth; return; }
    double const now = wct_seconds();
    auto sback = std::move(curstack.back());
    function_with_input_size const & sfunc = sback.first;
    running_stats & s = sback.second;
    curstack.pop_back();

    s.cool_down(now);
    if (!curstack.empty())
        curstack.back().second.time_children += s.real;
    s.real -= s.time_children;
    unsigned int const level = --depth;
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

#define xxxFORCE_PRINT_ALWAYS

#ifndef FORCE_PRINT_ALWAYS
    /* Is it any useful to print something new ? */
    if (now < last_print_time + 5) return;

    int needprint = 0;
    for(unsigned int k = 0 ; !needprint && k < levels.size() ; k++) {
        level_stats & u(levels[k]);
        double const t = u.projected_time();
        double const t0 = u.last_printed_projected_time;
        needprint = (t < 0.9 * t0) || (t > 1.1 * t0);
    }

    if (!needprint) return;
#endif

    last_print_time = now;
    last_print_position = std::make_pair(level, F.ncalled);

    print(level);
}

void tree_stats::final_print()
{
    ASSERT_ALWAYS(depth == 0);
    if (last_print_position != std::make_pair(0u, 1u))
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

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (!rank)
            fmt::print("lingen done at: {}\n", eta_string);
    }
}

void tree_stats::begin_plan_smallstep(std::string const & func, weighted_double const & theory)
{
    try {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank) return;
        ASSERT_ALWAYS(!curstack.empty());
        running_stats& s(curstack.back().second);

        step_time::steps_t * where = &s.steps;
        if (!s.nested_substeps.empty())
            where = & s.current_substep().steps;
        running_stats::steps_t::iterator const it =where->insert({func, step_time(func)}).first;
        s.nested_substeps.push_back(it);

        step_time & S(it->second);
        S.set_total_ncalls(s.total_ncalls());
        S.items_per_call = theory.n;
        S.planned_time += theory.t * theory.n;
        ASSERT_ALWAYS(S.ncalled == S.planned_calls);
        S.planned_calls++;
    } catch (std::runtime_error const & e) {
        std::stringstream os;
        os << fmt::format("Exception at {}({}, {:.3g}, {})\n",
                __func__, func, theory.t, theory.n);
        os << "State of *this\n";
        debug_print(os);
        os << "Error message: " << e.what() << "\n";
        throw std::runtime_error(os.str());
    }
}
void tree_stats::end_plan_smallstep()
{
    try {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank) return;
        ASSERT_ALWAYS(!curstack.empty());
        running_stats& s(curstack.back().second);
        s.nested_substeps.pop_back();
    } catch (std::runtime_error const & e) {
        std::stringstream os;
        // now that we have compile-time checking of format strings, at
        // least with c++17, we can safely silence this false positive.
        // coverity[fun_call_w_exception]
        os << fmt::format("Exception at {}()\n", __func__);
        os << "State of *this\n";
        debug_print(os);
        os << "Error message: " << e.what() << "\n";
        throw std::runtime_error(os.str());
    }
}

/* The only subtlety here is that we expect that the planning for the
 * smallstep is already done, and that we're now interested in planning
 * for the micro steps only.
 */
void tree_stats::begin_plan_smallstep_microsteps(std::string const & func)
{
    try {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank) return;
        ASSERT_ALWAYS(!curstack.empty());
        running_stats& s(curstack.back().second);

        step_time::steps_t * where = &s.steps;
        if (!s.nested_substeps.empty())
            where = & s.current_substep().steps;
        running_stats::steps_t::iterator const it =where->insert({func, step_time(func)}).first;
        s.nested_substeps.push_back(it);

    } catch (std::runtime_error const & e) {
        std::stringstream os;
        // see above
        // coverity[fun_call_w_exception]
        os << fmt::format("Exception at {}\n", __func__);
        os << "State of *this\n";
        debug_print(os);
        os << "Error message: " << e.what() << "\n";
        throw std::runtime_error(os.str());
    }
}

void tree_stats::begin_smallstep(std::string const & func, unsigned int ncalls)
{
    try {
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
        running_stats::steps_t::iterator const it =where->insert({func, step_time(func)}).first;
        it->second.set_total_ncalls(s.total_ncalls());
        s.nested_substeps.push_back(it);
        step_time & S(s.current_substep());
        S.items_pending += ncalls;
        ASSERT_ALWAYS(S.items_pending <= S.items_per_call);
        S.heat_up();
    } catch (std::runtime_error const & e) {
        std::stringstream os;
        // see above
        // coverity[fun_call_w_exception]
        os << fmt::format("Exception at {}({},{})\n", __func__, func, ncalls);
        os << "State of *this\n";
        debug_print(os);
        os << "Error message: " << e.what() << "\n";
        throw std::runtime_error(os.str());
    }
}

void tree_stats::end_smallstep()
{
    try {
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
    } catch (std::runtime_error const & e) {
        std::stringstream os;
        // see above
        // coverity[fun_call_w_exception]
        os << fmt::format("Exception at {}()\n", __func__);
        os << "State of *this\n";
        debug_print(os);
        os << "Error message: " << e.what() << "\n";
        throw std::runtime_error(os.str());
    }
}

void tree_stats::skip_smallstep(std::string const & func, unsigned int ncalls)
{
    try {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank) return;
        ASSERT_ALWAYS(!curstack.empty());
        running_stats& s(curstack.back().second);
        step_time::steps_t * where = &s.steps;
        if (!s.nested_substeps.empty())
            where = & s.current_substep().steps;
        running_stats::steps_t::iterator const it = where->insert({func, step_time(func)}).first;
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
    } catch (std::runtime_error const & e) {
        std::stringstream os;
        // see above
        // coverity[fun_call_w_exception]
        os << fmt::format("Exception at {}({}, {})\n", __func__, func, ncalls);
        os << "State of *this\n";
        debug_print(os);
        os << "Error message: " << e.what() << "\n";
        throw std::runtime_error(os.str());
    }
}

bool tree_stats::local_smallsteps_done(bool compulsory) const
{
    try {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank) return true;
        ASSERT_ALWAYS(!curstack.empty());
        running_stats const & s(curstack.back().second);
        step_time::steps_t const * where = &s.steps;
        if (!s.nested_substeps.empty())
            where = & s.current_substep().steps;
        for(auto const & s : *where) {
            if (s.second.items_pending) {
                if (compulsory) {
                    throw std::runtime_error(s.second.name + " has pending items");
                }
                return false;
            }
        }
        return true;
    } catch (std::runtime_error const & e) {
        std::stringstream os;
        // see above
        // coverity[fun_call_w_exception]
        os << fmt::format("Exception at {}()\n", __func__);
        os << "State of *this\n";
        debug_print(os);
        os << "Error message: " << e.what() << "\n";
        throw std::runtime_error(os.str());
    }
}

std::ostream& tree_stats::debug_print(std::ostream& os) const
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank) return os;
    ASSERT_ALWAYS(!curstack.empty());
    running_stats const & s(curstack.back().second);
    s.debug_print(os, "");
    return os;
}

std::ostream& tree_stats::step_time::debug_print(std::ostream& os, std::string const & indent) const {
    os << indent << fmt::format("name={}\n", name);
    os << indent << fmt::format("items_pending={}\n", items_pending);
    os << indent << fmt::format("items_per_call={}\n", items_per_call);
    os << indent << fmt::format("ncalled={}\n", ncalled);
    os << indent << fmt::format("planned_calls={}\n", planned_calls);
    os << indent << fmt::format("planned_time={}\n", planned_time);
    os << indent << fmt::format("total_ncalls={}\n", total_ncalls());
    os << indent << fmt::format("is_transition={}\n", is_transition_level());
    os << indent << fmt::format("real={}\n", real);
    for(auto const & s : steps) {
        s.second.debug_print(os, indent + "  ");
    }
    return os;
}


tree_stats::step_time & tree_stats::step_time::operator+=(step_time const & x)
{
    try {
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
    } catch (std::runtime_error const & e) {
        std::stringstream os;
        os << "Exception while adding counters with id " << name << "\n";
        os << "State of *this\n";
        debug_print(os, "(*this)  ");
        os << "State of x\n";
        x.debug_print(os, "(x)      ");
        os << "Error message: " << e.what() << "\n";
        throw std::runtime_error(os.str());
    }
} 
