#ifndef BBLAS_LEVEL2B_HPP_
#define BBLAS_LEVEL2B_HPP_

// IWYU pragma: private, include "bblas.hpp"
// IWYU pragma: friend ".*/bblas.*"

#include <cstdint>         // for uint64_t
#include "bblas_mat64.hpp"

/**********************************************************************/
/* level 2b: vector times matrices.
 *      mul_o64_6464    (several variants)
 *      mul_o64_T6464   (several variants)
 */

/* implemented here:
 *    - mul_o64_6464 ; multiply a vector by a matrix
 *    - mul_o64_T6464 ; multiply a vector by the transpose of a matrix
 */

/* implementation details, variants */
void mul_o64_6464_C_lsb(uint64_t * r, uint64_t a, mat64 const & w);
void mul_o64_6464_C_msb(uint64_t *r, uint64_t a, mat64 const & w);
void mul_o64_T6464_C_parity(uint64_t * w, uint64_t a, mat64 const & b);
void mul_o64_T6464_C_parity3(uint64_t * w, uint64_t a, mat64 const & b);

/* final exported choices. */
void mul_o64_6464(uint64_t * r, uint64_t a, mat64 const & w);
void mul_o64_T6464(uint64_t * w, uint64_t a, mat64 const & b);

#endif	/* BBLAS_LEVEL2B_HPP_ */
