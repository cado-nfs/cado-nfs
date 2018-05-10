#ifndef LAS_DUPLICATE_HPP_
#define LAS_DUPLICATE_HPP_

#include "las-types.hpp"
#include "relation.hpp"
#include <memory>

sieve_info *
fill_in_sieve_info(las_todo_entry const& doing,
                   uint32_t I, uint32_t J,
                   cxx_cado_poly const &, siever_config const &);

int
relation_is_duplicate(relation const& rel,
        las_todo_entry const & doing,
        las_info const& las,
        siever_config const & sc,
        std::shared_ptr<facul_strategies_t> old_strategies,
        int adjust_strategy);

#endif
