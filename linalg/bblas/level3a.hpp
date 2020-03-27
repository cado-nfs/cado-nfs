#ifndef BBLAS_LEVEL3A_HPP_
#define BBLAS_LEVEL3A_HPP_

#include "bblas.hpp"

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

int mat64_eq(mat64_srcptr a, mat64_srcptr b);
int mat64_is_uppertriangular(mat64_srcptr u);
int mat64_is_lowertriangular(mat64_srcptr u);
int mat64_triangular_is_unit(mat64_srcptr u);
void mat64_set_identity(mat64_ptr m);
void mat64_copy(mat64_ptr b, mat64_srcptr a);
void mat64_set_zero(mat64_ptr m);

/* implementation details, variants */
void mat64_add_C(mat64_ptr C, mat64_srcptr A, mat64_srcptr B);
void mat64_transpose_recursive_inplace(mat64_ptr a);
void mat64_transpose_simple_and_stupid(mat64_ptr dst, mat64_srcptr src);
void mat64_transpose_recursive(mat64_ptr dst, mat64_srcptr src);

/* final exported choices. */
void mat64_add(mat64_ptr C, mat64_srcptr A, mat64_srcptr B);
void mat64_transpose(mat64_ptr dst, mat64_srcptr src);

#endif	/* BBLAS_LEVEL3A_HPP_ */
