#ifndef CONVEX_HULL_H
#define CONVEX_HULL_H

#include "point.h"      // point_t
#include "tab_point.h"

#ifdef __cplusplus
extern "C" {
#endif

//return index of the new point of convex hull
int select_next_point(tabular_point_t * t, point_t * pt);

int search_init_point(tabular_point_t * t);

tabular_point_t *convex_hull(tabular_point_t * t);

#ifdef __cplusplus
}
#endif

#endif				/* CONVEX_HULL_H */
