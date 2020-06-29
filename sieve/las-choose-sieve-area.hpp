#ifndef LAS_CHOOSE_SIEVE_AREA_HPP_
#define LAS_CHOOSE_SIEVE_AREA_HPP_

#include <cstdint>       // for uint32_t
#include <memory>        // for shared_ptr
#include "las-info.hpp"  // for las_info
class nfs_aux;
struct las_todo_entry;
struct qlattice_basis;
struct siever_config;

extern int never_discard;

bool choose_sieve_area(las_info const & las,
        std::shared_ptr<nfs_aux> aux_p,
        las_todo_entry const & doing,
        siever_config & conf, qlattice_basis & Q, uint32_t & J);

/* This second version is single-threaded, and does not define internal
 * timing measurements
 */
bool choose_sieve_area(las_info const & las,
        las_todo_entry const & doing,
        siever_config & conf, qlattice_basis & Q, uint32_t & J);

#endif	/* LAS_CHOOSE_SIEVE_AREA_HPP_ */
