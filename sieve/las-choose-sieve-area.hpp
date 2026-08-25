#ifndef CADO_LAS_CHOOSE_SIEVE_AREA_HPP
#define CADO_LAS_CHOOSE_SIEVE_AREA_HPP

#include <cstdint>       // for uint32_t
#include <memory>        // for shared_ptr
#include "las-info.hpp"  // for las_info
#include "sieve-methods.hpp"
#include "threadpool.hpp"

class nfs_aux;
struct special_q;
struct siever_config;

extern int never_discard;

bool choose_sieve_area(
        las_info const & las,
        thread_pool & pool,
        special_q_task const &,
        siever_config & conf,
        special_q_data_class auto & Q,
        uint32_t & J);

/* This second version is single-threaded, and does not define internal
 * timing measurements
 */
bool choose_sieve_area(
        las_info const & las,
        special_q_task const &,
        siever_config & conf,
        special_q_data_class auto & Q,
        uint32_t & J);

#endif	/* CADO_LAS_CHOOSE_SIEVE_AREA_HPP */
