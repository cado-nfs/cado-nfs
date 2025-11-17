#include "cado.h" // IWYU pragma: keep

#include <cstddef>

#include <memory>
#include <atomic>

#include "las-globals.hpp"

#include "las-config.hpp"
#include "las-output.hpp"

// NOLINTBEGIN(cppcoreguidelines-avoid-non-const-global-variables)
size_t base_memory = 0;
int recursive_descent = 0;
int prepend_relation_time = 0;
int exit_after_rel_found = 0;

std::atomic<bool> global_exit_semaphore;

int allow_largesq = 0;
int sync_at_special_q = 0;
int sync_thread_pool = 0;
int trialdiv_first_side = 0;

double general_grace_time_ratio = DESCENT_DEFAULT_GRACE_TIME_RATIO;

double tt_qstart;

std::unique_ptr<las_output> main_output;

// NOLINTEND(cppcoreguidelines-avoid-non-const-global-variables)

