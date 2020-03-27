#ifndef BBLAS_LEVEL2A_HPP_
#define BBLAS_LEVEL2A_HPP_

#include "bblas.hpp"

/**********************************************************************/
/* level 2a: rank-1 updates.
 *      addmul_To64_o64 (several variants)
 */

/* implemented here:
 *    - addmul_To64_o64 ; add to a 64*64 matrix the rank-1 product
 *      obtained my multiplying two vectors together.
 *
 */

/* implementation details, variants */
void addmul_To64_o64_lsb(uint64_t * r, uint64_t a, uint64_t w);
void addmul_To64_o64_msb(uint64_t * r, uint64_t a, uint64_t w);
void addmul_To64_o64_lsb_packof2(uint64_t * r, uint64_t a, uint64_t w);

#if defined(HAVE_SSE2) && ULONG_BITS == 64
void addmul_To64_o64_lsb_sse_v1(uint64_t * r, uint64_t a, uint64_t w);
#endif

/* final exported choices */
void addmul_To64_o64(uint64_t * r, uint64_t a, uint64_t w);

#endif	/* BBLAS_LEVEL2A_HPP_ */
