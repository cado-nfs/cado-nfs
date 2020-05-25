#include "cado.h" // IWYU pragma: keep
#include "tab_point.h"
#include "macros.h"

#include <stdlib.h>

tabular_point_t *tabular_point_create(void)
{
    tabular_point_t *t = malloc(sizeof(*t));
    ASSERT(t != NULL);

    t->index = 0;
    t->size = 2;

    t->tab = malloc(t->size * sizeof(point_t *));
    ASSERT(t->tab != NULL);

    return t;
}

void tabular_point_free(tabular_point_t * t)
{
    for (int i = 0; i < t->index; i++)	//size
	point_free(t->tab[i]);
    free(t->tab);
    free(t);
}

void tabular_point_realloc(tabular_point_t * t)
{
    t->tab = realloc(t->tab, t->size * 2 * (sizeof(point_t *)));
    ASSERT(t->tab != NULL);
    t->size *= 2;
}

int tabular_point_get_index(tabular_point_t * t)
{
    return t->index;
}

point_t *tabular_point_get_point(tabular_point_t * t, int index)
{
    ASSERT(index < t->index);
    return t->tab[index];
}

void tabular_point_add_point(tabular_point_t * t, point_t * pt)
{
    tabular_point_add(t, pt->number, pt->x, pt->y);
}

void tabular_point_add(tabular_point_t * t, int numero, double x, double y)
{
    if (t->index >= t->size)
	tabular_point_realloc(t);
    t->tab[t->index] = point_create(numero, x, y);
    ASSERT(t->tab[t->index] != NULL);
    t->index++;
}
