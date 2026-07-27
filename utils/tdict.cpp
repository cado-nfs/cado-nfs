#include "cado.h" // IWYU pragma: keep

#include <cstdint>

#include <ostream>
#include <algorithm>
#include <limits>
#include <vector>
#include <string>
#include <set>

#include "nlohmann/json.hpp"
#include "fmt/format.h"
#include "fmt/base.h"

#include "tdict.hpp"
#include "utils_cxx.hpp"
#include "params.hpp"
#include "verbose.hpp"
#include "fstream_maybe_compressed.hpp"

namespace chronograms {
    static int enable = 0;
    static std::string chronogram_file;

    void configure_switches(cxx_param_list &) {
        // pl.configure_switch("-gantt");
    }
    void declare_usage(cxx_param_list & pl) {
        pl.declare_usage("chronogram", "Output data to generate a time chart and save it to the given file");
    }
    void interpret_parameters(cxx_param_list & pl) {
        enable = pl.parse("-chronogram", chronogram_file);
    }
    bool is_enabled() { return enable; }
} /* namespace chronograms */

namespace tdict {

    static int production_mode = 0;
    static int global_enable = 0;

    int is_production_mode() { return production_mode; }
    int is_enabled() { return global_enable; }
#ifndef DISABLE_TIMINGS

void declare_usage(cxx_param_list & pl)
{
    pl.declare_usage("T",   "enable fine-grain timings (use twice to get them for each q)");
    pl.declare_usage("production", "Sort of an opposite to -v. Disable all diagnostics except the cheap or critical ones. See #21688 and #21825.");
}

void configure_aliases(cxx_param_list &)
{
}

void configure_switches(cxx_param_list & pl)
{
    pl.configure_switch("-production");
    pl.configure_switch("-T");
}

void interpret_parameters(cxx_param_list & pl)
{
    pl.parse("-production", production_mode);
    pl.parse("-T", global_enable);
    if (is_production_mode() && is_enabled())
        pl.fail("-T and -production are incompatible");
}

std::ostream& operator<<(std::ostream & o, timer_seconds_thread_and_wct::type const & a) {
    o << a.t << " (" << a.w << " wct";
    if (a.w) o << ", " << a.t/a.w*100.0 << "% cpu";
    o << ")";
    return o;
}

#else

void declare_usage(cxx_param_list &) {}
void configure_aliases(cxx_param_list &) {}
void configure_switches(cxx_param_list &) {}
void interpret_parameters(cxx_param_list & pl) {}

#endif

} /* namespace tdict */

std::string chronograms::format_as(const bubble_info& info)
{
    using kind_t = chronograms::bubble_info::kind_t;
    auto const & D = info.data;
    // NOLINTBEGIN(cppcoreguidelines-pro-type-union-access)
    switch (info.kind) {
        case kind_t::SSS:
            return fmt::format(
                    "SSS side {} level {}",
                    D.sss.side, D.sss.level);
        case kind_t::FIB:
            return fmt::format(
                    "FIB side {} level {} B {} slice {}",
                    D.fib.side, D.fib.level, D.fib.B, D.fib.slice);
        case kind_t::DS:
            return fmt::format(
                    "DS side {} level {} B {}",
                    D.ds.side, D.ds.level, D.ds.B);
        case kind_t::PCLAT:
            return fmt::format(
                    "PCLAT side {} level {} slice {}",
                    D.fib.side, D.fib.level, D.fib.slice);
        case kind_t::PBR:
            return fmt::format(
                    "PBR M {} B {}",
                    D.pbr.M, D.pbr.B);
        case kind_t::INIT: return "INIT";
        case kind_t::QLATTICE: return "QLATTICE";
        case kind_t::SLICING: return "SLICING";
        case kind_t::ALLOC: return "ALLOC";
        case kind_t::AB: return "AB";
        case kind_t::ECM: return "ECM";
        case kind_t::DUPCHECK: return "DUPCHECK";
        case kind_t::BOTCHED: return "BOTCHED";
    }
    // NOLINTEND(cppcoreguidelines-pro-type-union-access)
    return {};
}

void timetree_t::display_chart() const
{
    if (!chronograms::is_enabled())
        return;

    verbose_fmt_print (0, 0,
            "# Chronogram info ({} entries) will go to {}\n",
            chart.size(),
            chronograms::chronogram_file);
    
    ofstream_maybe_compressed out(chronograms::chronogram_file);;

    std::set<int> thread_set;

    auto time_min = std::numeric_limits<uint64_t>::max();
    auto time_max = std::numeric_limits<uint64_t>::min();
    auto thr_min = std::numeric_limits<int>::max();
    auto thr_max = std::numeric_limits<int>::min();

    for (const auto& b : chart) {
        time_min = std::min(time_min, b.t0);
        time_max = std::max(time_max, b.t1);
        thr_min = std::min(thr_min, b.thread);
        thr_max = std::max(thr_max, b.thread);

        thread_set.insert(b.thread);
    }

    verbose_fmt_print (0, 0,
            "# Chronogram info has data for {} threads over {:.2f} seconds\n",
            thr_max - thr_min + 1,
            double_ratio(time_max - time_min, 1.0e9));

    using json = nlohmann::json;

    json J;
    J["win_start"] = time_min;
    J["win_end"] = time_max;
    J["categories"] = chronograms::bubble_info::kind_names;
    J["events"] = json();

    for (const auto& b : chart) {
        json E;
        E.push_back(b.t0 - time_min);
        E.push_back(b.t1 - b.t0);
        E.push_back(static_cast<int>(b.info.kind));
        E.push_back(b.on_cpu);
        E.push_back(fmt::format("{}", b.info));
        J["events"][std::to_string(b.thread)].push_back(E);
    }
    out << J;

    verbose_fmt_print (0, 0,
            "# Chronogram info can be viewed"
            " using https://cado-nfs.inria.fr/chronogram.html"
            " or ./scripts/chronograms/chronograms.html\n");
}
