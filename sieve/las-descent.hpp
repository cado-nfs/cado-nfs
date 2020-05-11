#ifndef LAS_DESCENT_HPP_
#define LAS_DESCENT_HPP_

#include "las-info.hpp"
#include "las-todo-list.hpp"
#include "las-todo-entry.hpp"
#include "relation.hpp"

bool register_contending_relation(las_info const & las, las_todo_entry const & doing, relation & rel);
void postprocess_specialq_descent(las_info & las, las_todo_list & todo, las_todo_entry const & doing, timetree_t & timer_special_q);

#endif	/* LAS_DESCENT_HPP_ */
