#include "cado.h" // IWYU pragma: keep

#include <ostream>
#include <algorithm>
#include <limits>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <algorithm>

#include "fmt/format.h"
#include "fmt/base.h"
#include "fmt/ostream.h"

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
        case kind_t::SKEWGAUSS: return "SKEWGAUSS";
        case kind_t::ADJUST: return "ADJUST";
        case kind_t::SLICING: return "SLICING";
        case kind_t::ALLOC: return "ALLOC";
        case kind_t::AB: return "AB";
        case kind_t::ECM: return "ECM";
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

    double time_min = std::numeric_limits<double>::max();
    double time_max = std::numeric_limits<double>::min();
    int thr_min = std::numeric_limits<int>::max();
    int thr_max = std::numeric_limits<int>::min();

    for (const auto& b : chart) {
        double t0 = b.tv_get.tv_sec + b.tv_get.tv_usec * 1e-6;
        double t1 = b.tv_put.tv_sec + b.tv_put.tv_usec * 1e-6;
        time_min = std::min(time_min, t0);
        time_max = std::max(time_max, t1);
        thr_min = std::min(thr_min, b.thread);
        thr_max = std::max(thr_max, b.thread);

        thread_set.insert(b.thread);
    }

    verbose_fmt_print (0, 0,
            "# Chronogram info has data for {} threads over {:.2f} seconds\n",
            thr_max - thr_min + 1, time_max - time_min);

    out << fmt::format("{{\n"
            "\"win_start\":{:.9f},\n"
            "\"win_end\":{:.9f},\n",
            time_min, time_max);

    out << "\"categories\":[";
    out << join(chronograms::bubble_info::kind_names, ",",
            [](auto c) { return fmt::format("\"{}\"", c); });
    out << "],\n\"events\":{\n";

    // Write Events Grouped by Thread
    std::map<int, std::vector<const chronograms::bubble*>> thread_bubbles;
    for (const auto& b : chart) thread_bubbles[b.thread].push_back(&b);

    size_t t_count = 0;
    for (const auto& [thread_id, bubbles] : thread_bubbles) {
        out << fmt::format("\"{}\":[", thread_id);
        for (size_t i = 0; i < bubbles.size(); ++i) {
            const auto& b = *bubbles[i];
            double t0 = b.tv_get.tv_sec + b.tv_get.tv_usec * 1e-6;
            double t1 = b.tv_put.tv_sec + b.tv_put.tv_usec * 1e-6;
            t0 -= time_min;
            t1 -= time_min;

            std::string cat = fmt::format("{}", b.info);
            int cat_idx = static_cast<int>(b.info.kind);

            // Output tuple: [rel_start, duration, cat_idx, metric, "desc"]
            out << fmt::format("[{:.9f},{:.9f},{},{},\"{}\"]{}",
                               t0, t1 - t0,
                               cat_idx,
                               0, cat,
                               (i + 1 < bubbles.size() ? ",\n" : ""));
        }
        out << "]" << (++t_count < thread_bubbles.size() ? ",\n" : "\n");
    }
    out << "}\n}\n";

    verbose_fmt_print (0, 0,
            "# Chronogram info can be viewed"
            " using https://cado-nfs.inria.fr/chronogram.html"
            " or ./scripts/chronograms/chronograms.html\n");
}

#if 0   /* why do I even need this??? */

template class std::map<tdict::key, tdict::slot_base const *>;

template struct tdict::tree<tdict::timer_seconds_thread>;
template class std::map<tdict::key, tdict::tree<tdict::timer_seconds_thread> >;
// template struct std::pair<tdict::key const, tdict::slot_base const *>;
template struct tdict::tree<tdict::timer_seconds_thread>::accounting_child_meta<tdict::tree<tdict::timer_seconds_thread>::accounting_base>;

#ifdef  HAVE_GCC_STYLE_AMD64_INLINE_ASM
template struct tdict::tree<tdict::timer_ticks>;
template class std::map<tdict::key, tdict::tree<tdict::timer_ticks> >;
template struct tdict::tree<tdict::timer_ticks>::accounting_child_meta<tdict::tree<tdict::timer_ticks>::accounting_base>;
#else
template struct tdict::tree<tdict::timer_none>;
template class std::map<tdict::key, tdict::tree<tdict::timer_none> >;
template struct tdict::tree<tdict::timer_none>::accounting_child_meta<tdict::tree<tdict::timer_none>::accounting_base>;
#endif

#endif

