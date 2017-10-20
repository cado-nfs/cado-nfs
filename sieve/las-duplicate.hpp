#ifndef LAS_DUPLICATE_HPP_
#define LAS_DUPLICATE_HPP_

#include "las-types.hpp"
#include "relation.h"

sieve_info *
fill_in_sieve_info(las_todo_entry const& doing,
                   uint32_t I, uint32_t J,
                   cado_poly_ptr cpoly, siever_config const & conf, int nb_threas);

/* We take a non-const reference because we're (temporarily) sharing the
 * pointers used for strategies and such.
 */
int relation_is_duplicate(relation const&, int, sieve_info &);

#endif
