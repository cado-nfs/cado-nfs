#include "cado.h" // IWYU pragma: keep
                  //
#include "double_poly.h"
#include "double_poly_complex_roots.h"
#include "polyroots.h"

/* This code cannot be in double_poly.cpp because that would trigger a
 * warning about _Complex being deprecated post-c++17. At least g++-15
 * moans.
 *
 * And anyway, this is all slated for removal eventually. Once we get the
 * generic polynomial<double> type right, we should be able to reap the
 * benefits of it, and get a replacement for much of the code here.
 */

int double_poly_complex_roots(double _Complex *roots, double_poly_srcptr f)
{
    return poly_roots_double(f->coeff, f->deg, roots);
}

int double_poly_complex_roots_long(long double _Complex *roots, double_poly_srcptr f)
{
    return poly_roots_longdouble(f->coeff, f->deg, roots);
}

