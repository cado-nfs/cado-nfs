#include "cado.h" // IWYU pragma: keep

#include <ostream>
#include <string>

#include "fmt/base.h"

#include "tdict.hpp"
#include "params.hpp"

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

