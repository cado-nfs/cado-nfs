#include "cado.h" // IWYU pragma: keep
#include <ostream>    // std::ostream // IWYU pragma: keep
#include <pthread.h>
#include "tdict.hpp"
#include "params.h"

namespace tdict {

/* default value is set in configure_switches */ 
int global_enable;

#ifndef DISABLE_TIMINGS

pthread_mutex_t slot_base::m = PTHREAD_MUTEX_INITIALIZER;

void declare_usage(cxx_param_list & pl)
{
    param_list_decl_usage(pl, "T",   "enable fine-grain timings (use twice to get them for each q)");
}

void configure_aliases(cxx_param_list &)
{
}

void configure_switches(cxx_param_list & pl)
{
    param_list_configure_switch(pl, "-T", &global_enable);
    /* We now rely on fine-grain switches to provide _all_ timings */
    global_enable = 1;
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

#endif

};

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

