#include "cado.h"
#include <stdlib.h>
#include "convex_hull.h"
#include <math.h>

/*
  The following model begin from the most left point and finish by the
  most right point.
*/


static int cmp_double (double a, double b)
{
    double diff = a - b;
    double precision = 0.0000000001;

    if ( diff < precision && diff > -1*precision)
	return 0;
    else if (diff < precision)
	return -1;
    else
	return 1;
}

static double convertEnDegre(const double angle)
{
    return angle * 180 / M_PI;
}

double CST_MUL = 1;

//we suppose that the coordinates of pt2 are highter than those of pt1
static double compute_angle(point_t * pt1, point_t * pt2)
{
    double op = pt2->y - pt1->y;
    double adj = (pt2->x - pt1->x)* CST_MUL;

    if (adj == 0) {
	if (op > 0)
	    return 90;
	else
	    return 360;
    } else if (op == 0) {
	if (adj > 0)
	    return 0;
	else
	    return 180;
    }
    double angle = convertEnDegre(atan(op / adj));
    return angle;
}

//return index of the new point of convex hull
int select_next_point(tabular_point_t * t, point_t * pt)
{
    int res = -1;
    int len = tabular_point_get_index(t);
    point_t *elem;
    double angle_min = 90;
    double angle;

    for (int i = 0; i < len; i++) {
	elem = tabular_point_get_point(t, i);
	if (cmp_double(pt->x, elem->x) == -1) {
	    angle = compute_angle(pt, elem);
	    if (angle < angle_min) {
		angle_min = angle;
		res = i;
	    }
	}
    }
    return res;
}

int search_init_point(tabular_point_t * t)
{
    int res = 0;
    int len = tabular_point_get_index(t);
    point_t *pt_min = tabular_point_get_point(t, res);
    point_t *elem;
    for (int i = 1; i < len; i++) {
	elem = tabular_point_get_point(t, i);
	if (cmp_double(pt_min->x, elem->x) > 0
	    || (cmp_double(pt_min->x, elem->x) == 0
		&& cmp_double(pt_min->y, elem->y) > 0)){
	    res = i;
	    pt_min = elem;
	}
    }
    return res;
}

tabular_point_t *convex_hull(tabular_point_t * t)
{
    //max et min!
    int len = t->index;
    double minx = INFINITY;
    double maxx = 0;
    double miny = INFINITY;
    double maxy = 0;
    for (int i =0; i < len; i++)
	{
	    point_t* p = tabular_point_get_point (t, i);
	    if (p->x < minx)
		minx = p->x;
	    if (p->x > maxx)
		maxx = p->x;
	    if (p->y < miny)
		miny = p->y;
	    if (p->y > maxy)
		maxy = p->y;
	}
    double scaley = log(maxy - miny)/log(10);
    double scalex = log(maxx - minx)/log(10);
    //to keep in the same scale between x and y!
    CST_MUL = pow(10, scaley - scalex);


    tabular_point_t *convex_hull = tabular_point_create();
    int index = search_init_point(t);
    while (index != -1) {
	tabular_point_add_point(convex_hull, t->tab[index]);
	index = select_next_point(t, t->tab[index]);
    }

    return convex_hull;
}
