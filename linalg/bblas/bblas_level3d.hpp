#ifndef BBLAS_LEVEL3D_HPP_
#define BBLAS_LEVEL3D_HPP_

// IWYU pragma: private, include "bblas.hpp"
// IWYU pragma: friend ".*/bblas.*"

#include "bblas_mat64.hpp"

/**********************************************************************/
/* level 3d: solution of triangular linear systems
 *      trsm64
 *      trsm64_general
 */

/* implemented here:
 *    - trsm64: gen two 64*64 bit matrices L and U, with L unit lower
 *      triangular, replace U by L^-1*U.
 *    - trsm64_general: same, but apply only the square submatrix of L
 *      whose diagonal indices are in the given integer interval.
 */

void trsm64_general(mat64 const & L, mat64 & U, unsigned int n0, unsigned int n1);

void trsm64(mat64 const & L, mat64 & U);

#endif	/* BBLAS_LEVEL3D_HPP_ */
