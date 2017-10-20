#ifndef TAB_POINT
#define TAB_POINT

#include "point.h"

typedef struct tabular_point {
    point_t **tab;
    int index;
    int size;
} tabular_point_t;

#ifdef __cplusplus
extern "C" {
#endif

tabular_point_t *tabular_point_create(void);

void tabular_point_free(tabular_point_t * t);

void tabular_point_realloc(tabular_point_t * t);

int tabular_point_get_index(tabular_point_t * t);

point_t *tabular_point_get_point(tabular_point_t * t, int index);

void tabular_point_add_point(tabular_point_t * t, point_t * pt);

void tabular_point_add(tabular_point_t * t, int numero, double x, double y);

#ifdef __cplusplus
}
#endif

#endif				/* tab_fm */
