#ifndef CADO_LAS_GLOBALS_HPP
#define CADO_LAS_GLOBALS_HPP

#include <cstddef>
#include <mutex>
#include <memory>

#include "las-output.hpp"

/* these are influenced by las's command line parameters, for example.
 * Most of the timers here are contenders for removal or better
 * integration in some other mechanism.
 */

// NOLINTBEGIN(cppcoreguidelines-avoid-non-const-global-variables)
extern size_t base_memory;
extern int recursive_descent;
extern int prepend_relation_time;
extern int exit_after_rel_found;
extern std::mutex protect_global_exit_semaphore;
extern int global_exit_semaphore;
extern int allow_largesq;
extern int sync_at_special_q;
extern int trialdiv_first_side;
extern double general_grace_time_ratio;
extern double tt_qstart;
extern std::unique_ptr<las_output> main_output;
// NOLINTEND(cppcoreguidelines-avoid-non-const-global-variables)

#endif	/* CADO_LAS_GLOBALS_HPP */
