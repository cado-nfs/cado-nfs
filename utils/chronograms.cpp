#include "cado.h" // IWYU pragma: keep

#include <algorithm>
#include <limits>
#include <string>
#include <vector>

#include "fmt/base.h"
#include "fmt/format.h"

#include "params.hpp"
#include "chronograms.hpp"
#include "verbose.hpp"
#include "fstream_maybe_compressed.hpp"
#include "utils_cxx.hpp"
#include "nlohmann/json.hpp"

namespace chronograms {
    static int enable = 0;
    static std::string chronogram_file;

    void configure_switches(cxx_param_list &) {
    }
    void declare_usage(cxx_param_list & pl) {
        pl.declare_usage("chronogram", "Output data to generate a time chart and save it to the given file");
    }
    void interpret_parameters(cxx_param_list & pl) {
        enable = pl.parse("-chronogram", chronogram_file);
    }
    bool is_enabled() { return enable; }

    std::string format_as(const bubble_info& info)
    {
        using kind_t = bubble_info::kind_t;
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

            /* we should never see NONE */
            case kind_t::NONE: return "NONE";
        }
        // NOLINTEND(cppcoreguidelines-pro-type-union-access)
        return {};
    }

    void display(std::map<size_t, std::vector<chronograms::bubble>> const & M)
    {
        if (!chronograms::is_enabled())
            return;

        size_t total_entries = 0;
        for(auto const & [ k, v ] : M)
            total_entries += v.size();

        verbose_fmt_print (0, 0,
                "# Chronogram info ({} threads, {} entries) will go to {}\n",
                M.size(), total_entries,
                chronograms::chronogram_file);

        ofstream_maybe_compressed out(chronograms::chronogram_file);;

        auto time_min = std::numeric_limits<uint64_t>::max();
        auto time_max = std::numeric_limits<uint64_t>::min();
        auto thr_min = std::numeric_limits<size_t>::max();
        auto thr_max = std::numeric_limits<size_t>::min();

        for(auto const & [ k, v ] : M) {
            for (const auto& b : v) {
                time_min = std::min(time_min, b.t0);
                time_max = std::max(time_max, b.t1);
            }
            thr_min = std::min(thr_min, k);
            thr_max = std::max(thr_max, k);
        }

        verbose_fmt_print (0, 0,
                "# Chronogram info has data for {:.2f} seconds\n",
                double_ratio(time_max - time_min, 1.0e9));

        using json = nlohmann::json;

        json J;
        J["format"] = 20260727;
        J["win_start"] = time_min;
        J["win_end"] = time_max;
        J["categories"] = chronograms::bubble_info::kind_names;
        J["events"] = json();

        for(auto const & [ k, v ] : M) {
            for (const auto& b : v) {
                json E;
                E.push_back(b.t0 - time_min);
                E.push_back(b.t1 - b.t0);
                E.push_back(static_cast<int>(b.info.kind));
                E.push_back(b.on_cpu);
                E.push_back(fmt::format("{}", b.info));
                J["events"][std::to_string(k)].push_back(E);
            }
        }
        out << J;

        verbose_fmt_print (0, 0,
                "# Chronogram info can be viewed"
                " using https://cado-nfs.inria.fr/chronogram.html"
                " or ./scripts/chronograms/chronograms.html\n");
    }
} /* namespace chronograms */

