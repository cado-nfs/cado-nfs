#ifndef BBLAS_LEVEL2A_HPP_
#define BBLAS_LEVEL2A_HPP_

// IWYU pragma: private, include "bblas.hpp"
// IWYU pragma: friend ".*/bblas.*"

#include "cado_config.h"    // for HAVE_SSE2, ULONG_BITS
#include <cstdint>         // for uint64_t
#include "bblas_mat64.hpp"

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
void addmul_To64_o64_lsb(mat64 & r, uint64_t a, uint64_t w);
void addmul_To64_o64_msb(mat64 & r, uint64_t a, uint64_t w);
void addmul_To64_o64_lsb_packof2(mat64 & r, uint64_t a, uint64_t w);

#if defined(HAVE_SSE2) && ULONG_BITS == 64
void addmul_To64_o64_lsb_sse_v1(mat64 & r, uint64_t a, uint64_t w);
#endif

/* final exported choices */
void addmul_To64_o64(mat64 & r, uint64_t a, uint64_t w);

#endif	/* BBLAS_LEVEL2A_HPP_ */
