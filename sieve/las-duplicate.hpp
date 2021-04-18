#ifndef LAS_DUPLICATE_HPP_
#define LAS_DUPLICATE_HPP_

#include "las-info.hpp"  // for las_info
struct las_todo_entry;
struct relation;
struct siever_config;
struct qlattice_basis;

/* Return true if the relation is found when sieving [doing]. conf, Q, J
 * are as returned by choose_sieve_area.
 */
bool
sq_finds_relation(las_info const & las,
        las_todo_entry const & doing,
        siever_config const & conf,
        qlattice_basis const & Q,
        uint32_t J,
        relation const& rel);

int
relation_is_duplicate(relation const& rel,
        las_todo_entry const & doing,
        las_info const& las);

#endif
