#ifndef BBLAS_LEVEL3B_HPP_
#define BBLAS_LEVEL3B_HPP_

#include "bblas.hpp"

/**********************************************************************/
/* level 3b: matrix multiplications
 *      mul_6464_6464   (several variants)
 *      addmul_6464_6464_fragment       TODO: missing test
 */

/* implemented here:
 *    - mul_6464_6464 ; self-explanatory
 *    - addmul_6464_6464_fragment: takes only a subset of rows of the
 *      first input, and a subset of columns of the second input.
 */

/* implementation details, variants */
#if defined(HAVE_SSE2) && ULONG_BITS == 64
void mul_6464_6464_sse(mat64_ptr C, mat64_srcptr A, mat64_srcptr B);
#endif

void mul_6464_6464_v2(mat64_ptr C, mat64_srcptr A, mat64_srcptr B);

void addmul_6464_6464_fragment_lookup4(uint64_t *C,
                   const uint64_t *A,
                   const uint64_t *B,
                   unsigned int i0,
                   unsigned int i1,
                   unsigned int yi0,
                   unsigned int yi1);

/* final exported choices */
void mul_6464_6464(mat64 C, mat64 A, mat64 B);
void addmul_6464_6464_fragment(uint64_t *C,
                   const uint64_t *A,
                   const uint64_t *B,
                   unsigned int i0,
                   unsigned int i1,
                   unsigned int yi0,
                   unsigned int yi1);


#endif	/* BBLAS_LEVEL3B_HPP_ */
