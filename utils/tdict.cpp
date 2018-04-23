#include "cado.h"
#include <iomanip> // for std::fixed and std::setprecision
#include <pthread.h>
#include <sys/time.h>
#include "verbose.h"
#include "tdict.hpp"
#include "params.h"

int tdict::global_enable = 0;
int time_bubble_chaser::enable = 0;

pthread_mutex_t tdict::slot_base::m = PTHREAD_MUTEX_INITIALIZER;

void tdict_decl_usage(param_list pl)
{
    param_list_decl_usage(pl, "T",   "enable fine-grain timings (use twice to get them for each q)");
}

void tdict_configure_switch(param_list pl)
{
    param_list_configure_switch(pl, "-T", &tdict::global_enable);
}

time_bubble_chaser::time_bubble_chaser(int thread, kind_t kind, id_t const & id): thread(thread), kind(kind), id(id) {
    if (!enable) return;
    gettimeofday(&tv_get, NULL);
    on_cpu = - (double) microseconds_thread();
}
time_bubble_chaser& time_bubble_chaser::put() {
    if (!enable) return *this;
    gettimeofday(&tv_put, NULL);
    on_cpu += microseconds_thread();
    return *this;
}
void timetree_t::display_chart() const
{
    verbose_output_print (0, 2, "# time chart has %zu entries\n", chart.size());
    for(auto const & T : chart) {
        std::ostringstream os;

        os << "thread " << T.thread;

        /* This switch should really be passed as a function pointer */
        switch(T.kind) {
            case time_bubble_chaser::FIB:
                os << " FIB"
                    << " side " << T.id[0]
                    << " level " << T.id[1]
                    << " B " << T.id[2]
                    << " slice " << T.id[3];
                break;
            case time_bubble_chaser::DS:
                os << " DS"
                    << " side " << T.id[0]
                    << " level " << T.id[1]
                    << " B " << T.id[2];
                break;
            case time_bubble_chaser::AB:
                os << " AB";
                for(auto const & x : T.id)
                    os << " x" << (&x-begin(T.id)) << " " << x;
                break;
            case time_bubble_chaser::PBR:
                os << " PBR";
                os  << " M " << T.id[1]
                    << " B " << T.id[2];
                //for(auto const & x : T.id)
                    //os << " x" << (&x-begin(T.id)) << " " << x;
                break;
            case time_bubble_chaser::PCLAT:
                os << " PCLAT"
                    << " side " << T.id[0]
                    << " level " << T.id[1]
                    << " slice " << T.id[3];
                break;
        };

        double t0 = T.tv_get.tv_sec + 1.0e-6 * T.tv_get.tv_usec;
        double t1 = T.tv_put.tv_sec + 1.0e-6 * T.tv_put.tv_usec;
        double t = T.on_cpu;
        os  << std::fixed << std::setprecision(9)
            << " t0 " << t0
            << " t1 " << t1
            << " time " << t;

        verbose_output_print (0, 2, "#  %s\n", os.str().c_str());
    }
}
