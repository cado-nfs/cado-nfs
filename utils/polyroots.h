#ifndef UTILS_POLYROOTS_H_
#define UTILS_POLYROOTS_H_

/* Do not include this file directly. The way to access these functions
 * is via double_poly.h
 */
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define MAX_POLYROOTS_ROOTFINDER_DEGREE   10

uint32_t poly_roots_double(const double *poly, uint32_t degree, double _Complex *roots);
uint32_t poly_roots_longdouble(const double *poly, uint32_t degree, long double _Complex *roots);


#ifdef __cplusplus
}
#endif

#endif	/* UTILS_POLYROOTS_H_ */
