#ifndef LAS_DESCENT_HPP_
#define LAS_DESCENT_HPP_

#include "las-info.hpp"
#include "tdict.hpp"     // for timetree_t
class las_todo_list;
struct las_todo_entry;
struct relation;

bool register_contending_relation(las_info const & las, las_todo_entry const & doing, relation & rel);
void postprocess_specialq_descent(las_info & las, las_todo_list & todo, las_todo_entry const & doing, timetree_t & timer_special_q);

#endif	/* LAS_DESCENT_HPP_ */
