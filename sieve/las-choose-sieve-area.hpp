#ifndef LAS_CHOOSE_SIEVE_AREA_HPP_
#define LAS_CHOOSE_SIEVE_AREA_HPP_

#include "las-types.hpp"
#include "verbose.h"

extern int never_discard;

bool choose_sieve_area(las_info const & las,
        std::shared_ptr<nfs_aux> aux_p,
        int adjust_strategy, las_todo_entry const & doing, siever_config & conf, qlattice_basis & Q, uint32_t & J);

/* This second version is single-threaded, and does not define internal
 * timing measurements
 */
bool choose_sieve_area(las_info const & las,
        int adjust_strategy, las_todo_entry const & doing, siever_config & conf, qlattice_basis & Q, uint32_t & J);

#endif	/* LAS_CHOOSE_SIEVE_AREA_HPP_ */
