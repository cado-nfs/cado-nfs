#ifndef BBLAS_LEVEL3A_HPP_
#define BBLAS_LEVEL3A_HPP_

#include "bblas.hpp"
#include <gmp.h>

/**********************************************************************/
/* level 3a: basic operations on 64*64 matrices.
 *      mat64_eq
 *      mat64_add
 *      mat64_is_uppertriangular
 *      mat64_is_lowertriangular
 *      mat64_triangular_is_unit
 *      mat64_set_identity
 *      mat64_set_zero
 *      mat64_copy
 *      mat64_transpose     (several variants, one obvious winner)
 */

/* implemented here: see above. names are self-explanatory. */

int mat64_eq(mat64 const & a, mat64 const & b);
void mat64_fill_random(mat64 & w, gmp_randstate_t rstate);
int mat64_is_uppertriangular(mat64 const & u);
int mat64_is_lowertriangular(mat64 const & u);
int mat64_triangular_is_unit(mat64 const & u);
void mat64_set_identity(mat64 & m);
void mat64_copy(mat64 & b, mat64 const & a);
void mat64_set_zero(mat64 & m);

/* implementation details, variants */
void mat64_add_C(mat64 & C, mat64 const & A, mat64 const & B);
void mat64_transpose_recursive_inplace(mat64 & a);
void mat64_transpose_simple_and_stupid(mat64 & dst, mat64 const & src);
void mat64_transpose_recursive(mat64 & dst, mat64 const & src);

/* final exported choices. */
void mat64_add(mat64 & C, mat64 const & A, mat64 const & B);
void mat64_transpose(mat64 & dst, mat64 const & src);

#endif	/* BBLAS_LEVEL3A_HPP_ */
