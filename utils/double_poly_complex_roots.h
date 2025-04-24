#ifndef CADO_UTILS_DOUBLE_POLY_COMPLEX_ROOTS_H
#define CADO_UTILS_DOUBLE_POLY_COMPLEX_ROOTS_H

#include <complex.h>    // IWYU pragma: keep
#include "double_poly.h"

/* This interface uses C99 _Complex types, and these don't go well with
 * c++ std::complex. So most probably new code should keep away from
 * these.
 */
#ifdef __cplusplus
extern "C" {
#endif

/* we use _Complex and not complex below, because it's mandatory in order
 * to be nice to C++ code */
int double_poly_complex_roots(double _Complex *roots, double_poly_srcptr);
int double_poly_complex_roots_long(long double _Complex *roots, double_poly_srcptr);

#ifdef __cplusplus
}
#endif

#endif	/* CADO_UTILS_DOUBLE_POLY_COMPLEX_ROOTS_H */
