#include "cado.h" // IWYU pragma: keep

#include "las-globals.hpp"

#include "las-config.h"
#include "las-output.hpp"


size_t base_memory = 0;
int recursive_descent = 0;
int prepend_relation_time = 0;
int exit_after_rel_found = 0;
std::mutex protect_global_exit_semaphore;
int global_exit_semaphore = 0;

int allow_largesq = 0;
int sync_at_special_q = 0;
int trialdiv_first_side = 0;

double general_grace_time_ratio = DESCENT_DEFAULT_GRACE_TIME_RATIO;

double tt_qstart;

las_output main_output;


