#ifndef CADO_LAS_DUPLICATE_HPP
#define CADO_LAS_DUPLICATE_HPP

#include "las-info.hpp"  // for las_info
#include "sieve-methods.hpp"

struct special_q;
struct relation;
struct siever_config;

/* Return true if the relation is found when sieving [doing]. conf, Q, J
 * are as returned by choose_sieve_area.
 */
bool
sq_finds_relation(las_info const & las,
        special_q const & doing,
        siever_config const & conf,
        special_q_data_class auto const & Q,
        uint32_t J,
        relation const& rel);

int
relation_is_duplicate(relation const& rel,
        special_q const & doing,
        las_info const& las);

#endif
