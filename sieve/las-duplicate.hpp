#ifndef LAS_DUPLICATE_HPP_
#define LAS_DUPLICATE_HPP_

#include "las-info.hpp"  // for las_info
struct las_todo_entry;
struct relation;

int
relation_is_duplicate(relation const& rel,
        las_todo_entry const & doing,
        las_info const& las);

#endif
