#include "cado.h"
#include <pthread.h>
#include "tdict.hpp"
#include "params.h"

int tdict::global_enable = 0;

pthread_mutex_t tdict::slot_base::m = PTHREAD_MUTEX_INITIALIZER;

void tdict_decl_usage(param_list pl)
{
    param_list_decl_usage(pl, "T",   "enable fine-grain timings (use twice to get them for each q)");
}

void tdict_configure_switch(param_list pl)
{
    param_list_configure_switch(pl, "-T", &tdict::global_enable);
}
