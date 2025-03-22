#ifndef CONVEX_HULL_HPP
#define CONVEX_HULL_HPP

#include "point.hpp"
#include "tab_point.hpp"

// return index of the new point of convex hull
int select_next_point(tabular_point const & t, point const & pt);

int search_init_point(tabular_point const & t);

tabular_point convex_hull(tabular_point const & t);
;

#endif /* CONVEX_HULL_HPP */
