#ifndef GALOIS_UTILS_H_
#define GALOIS_UTILS_H_

#include "mod_ul.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Galois automorphisms
   autom2.1: 1/x
   autom2.2: -x
   autom3.1: 1-1/x
   autom3.2: -1-1/x
   autom4.1: -(x+1)/(x-1)
   autom6.1: -(2x+1)/(x-1)
*/

void automorphism_init(int *ord, int mat[4], const char *galois_autom);
unsigned long automorphism_apply(residueul_t mat[4], unsigned long r, const modulusul_t mm, const unsigned long qq);

#ifdef __cplusplus
}
#endif


#endif /* GALOIS_UTILS_H_ */
