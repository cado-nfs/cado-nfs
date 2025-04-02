#include "cado.h" // IWYU pragma: keep

#include <cmath>

#include <algorithm>

#include "convex_hull.hpp"
#include "point.hpp"
#include "tab_point.hpp"

/*
  The following model begin from the most left point and finish by the
  most right point.
*/

static int cmp_double(double a, double b)
{
    double diff = a - b;
    double precision = 1e-10;

    if (diff < precision && diff > -1 * precision)
        return 0;
    else if (diff < precision)
        return -1;
    else
        return 1;
}

static double CST_MUL = 1;

// we suppose that the coordinates of p2 are highter than those of p1
static double compute_angle(point const & p1, point const & p2)
{
    double op = p2.y - p1.y;
    double adj = (p2.x - p1.x) * CST_MUL;

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
    return atan(op / adj);
}

// return index of the new point of convex hull
int select_next_point(tabular_point const & t, point const & pt)
{
    int res = -1;
    double angle_min = 90;
    double angle;

    for (int i = 0; i < (int)t.size(); i++) {
        auto const & elem = t[i];
        if (cmp_double(pt.x, elem.x) == -1) {
            angle = compute_angle(pt, elem);
            if (angle < angle_min) {
                angle_min = angle;
                res = i;
            }
        }
    }
    return res;
}

int search_init_point(tabular_point const & t)
{
    int res = 0;
    point pt_min = t[0];
    for (int i = 0; i < (int)t.size(); i++) {
        auto const & elem = t[i];
        if (cmp_double(pt_min.x, elem.x) > 0 ||
            (cmp_double(pt_min.x, elem.x) == 0 &&
             cmp_double(pt_min.y, elem.y) > 0)) {
            res = i;
            pt_min = elem;
        }
    }
    return res;
}

tabular_point convex_hull(tabular_point const & t)
{
    double minx = INFINITY;
    double maxx = 0;
    double miny = INFINITY;
    double maxy = 0;
    for (auto const & p: t) {
        minx = std::min(minx, p.x);
        miny = std::min(miny, p.y);
        maxx = std::max(maxx, p.x);
        maxy = std::max(maxy, p.y);
    }
    double scaley = log(maxy - miny) / log(10);
    double scalex = log(maxx - minx) / log(10);
    // to keep in the same scale between x and y!
    CST_MUL = pow(10, scaley - scalex);

    tabular_point convex_hull;
    int index = search_init_point(t);
    while (index != -1) {
        convex_hull.push_back(t[index]);
        index = select_next_point(t, t[index]);
    }

    return convex_hull;
}
