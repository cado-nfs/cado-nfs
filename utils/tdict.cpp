#include "cado.h" // IWYU pragma: keep

#include <string>
#include <ostream>
#include <algorithm>
#include <limits>

#include "fmt/format.h"
#include "fmt/base.h"
#include "fmt/ostream.h"

#include "tdict.hpp"
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
    
    ofstream_maybe_compressed chronogram_stream(chronograms::chronogram_file);;

    double time_min = std::numeric_limits<double>::max();
    double time_max = std::numeric_limits<double>::min();
    int thr_min = std::numeric_limits<int>::max();
    int thr_max = std::numeric_limits<int>::min();

    for(auto const & T : chart) {
        double t0 = T.tv_get.tv_sec + 1.0e-6 * T.tv_get.tv_usec;
        double t1 = T.tv_put.tv_sec + 1.0e-6 * T.tv_put.tv_usec;

        fmt::print(chronogram_stream, "t {} {:.9f} {:.9f} {} {}\n",
                T.thread, t0, t1 - t0, T.on_cpu, T.info);
        time_min = std::min(time_min, t0);
        time_max = std::max(time_max, t1);
        thr_min = std::min(thr_min, T.thread);
        thr_max = std::max(thr_max, T.thread);
    }
    verbose_fmt_print (0, 0, "# Chronogram info has data for {} threads over {:.2f} seconds\n",
            thr_max - thr_min + 1, time_max - time_min);
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

