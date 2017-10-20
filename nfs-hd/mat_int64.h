#ifndef MAT_INT64_H
#define MAT_INT64_H

#include <stdint.h>
#include "mat_Z.h"
#include "int64_vector.h"
#include "int64_poly.h"
#include "macros.h"
#include "matrix.h"
#include "list.h"

#ifdef __cplusplus
extern "C" {
#endif


/*
 * WARNING: to correspond to the definition of mat_Z, the first column and line
 *  are unusued.
 */

/*
 * Initialize a matrix.
 *
 * matrix: the matrix.
 * NumRows: number of rows in the matrix.
 * NumCols: number of cols in the matrix.
 */
void mat_int64_init(mat_int64_ptr matrix, unsigned int NumRows,
    unsigned int NumCols);

/*
 * Initialize a matrix and set the coefficients.
 *
 * matrix: the matrix.
 * NumRows: number of rows in the matrix.
 * NumCols: number of cols in the matrix.
 * coeff: an array with the coefficient of the matrix, given in the order to
 *  fill the row.
 */
void mat_int64_init_with_array(mat_int64_ptr matrix, unsigned int NumRows,
    unsigned int NumCols, int64_t * coeff);

/*
 * All the coeffecients of the matrix are equal to 0.
 *
 * matrix: the matrix.
 */
void mat_int64_set_zero(mat_int64_ptr matrix);

/*
 * Set the matrix to the identity matrix.
 *
 * matrix: the matrix.
 */
void mat_int64_set_identity(mat_int64_ptr matrix);

/*
 * Compare to matrices. Return 1 if these are equal, 0 otherwise.
 */
unsigned int mat_int64_equal(mat_int64_srcptr A, mat_int64_srcptr B);

/*
 * Copy A in B. Assume that B and A have the same size.
 *
 * B: the new matrix.
 * A: the root matrix.
 */
void mat_int64_copy(mat_int64_ptr B, mat_int64_srcptr A);

/*
 * Set a coefficient to i, this coefficient is in the row (row) and the column
 *  (col).
 *
 * matrix: the matrix.
 * i: the new coefficient.
 * row: number for row.
 * col: number for column.
 */
static inline void mat_int64_set_coeff(mat_int64_ptr matrix, int64_t i,
    unsigned int row, unsigned int col)
{
  ASSERT(row < matrix->NumRows);
  ASSERT(1 <= row);
  ASSERT(col < matrix->NumCols);
  ASSERT(1 <= col);

  matrix->coeff[row + 1][col + 1] = i;
}

/*
 * Get the coefficient at the row (row) and the column (col).
 *
 * matrix: the matrix.
 * row: number for row.
 * col: number for column.
 */
static inline int64_t mat_int64_get_coeff(mat_int64_srcptr matrix,
      unsigned int row, unsigned int col)
{
  ASSERT(row < matrix->NumRows);
  ASSERT(1 <= row);
  ASSERT(col < matrix->NumCols);
  ASSERT(1 <= col);

  return matrix->coeff[row + 1][col + 1];
}

/*
 * Delete a matrix.
 *
 * matrix: mat_int64_ptr, the matrix.
 */
void mat_int64_clear(mat_int64_ptr matrix);

/*
 * Write a matrix in a file.
 *
 * file: the file.
 * matrix: the matrix.
 */
void mat_int64_fprintf(FILE * file, mat_int64_srcptr matrix);

/*
 * Write a matrix with # before each line.
 *
 * file: the file.
 * matrix: the matrix.
 */
void mat_int64_fprintf_comment(FILE * file, mat_int64_srcptr matrix);

/*
 * Transpose a matrix_src, and set the result in matrix. Assume that matrix have
 *  the good size and is already initialized.
 *
 * matrix: the matrix equals to the transpose of matrix_src.
 * matrix_src: the root matrix.
 */
void mat_int64_transpose(mat_int64_ptr matrix, mat_int64_srcptr matrix_src);

/*
 * Compute naively the product A * B = C. Assume that C have the good size.
 *
 * C: the result.
 * A: a matrix.
 * B: a matrix.
 */
void mat_int64_mul_mat_int64(mat_int64_ptr C, mat_int64_srcptr A,
    mat_int64_srcptr B);

/*
 * Compute the product A * c and set the coefficient in a polynomial a.
 *
 * a: a polynomial.
 * A: a matrix.
 * c: a vector.
 */
void mat_int64_mul_int64_vector_to_int64_poly(int64_poly_ptr a,
    mat_int64_srcptr A, int64_vector_srcptr c);

/*
 * Compute a = A * c.
 *
 * a: result vector.
 * A: matrix.
 * c: vector.
 */
void mat_int64_mul_int64_vector(int64_vector_ptr a, mat_int64_srcptr A,
    int64_vector_srcptr c);

/*
 * Compute the product A * c (the coefficients of c is considered as
 *  coefficients of a vector) and set the coefficient in a polynomial a.
 *
 * a: a polynomial.
 * A: a matrix.
 * c: a polynomial.
 */
void mat_int64_mul_mpz_poly(mpz_poly_ptr a, mat_int64_srcptr A,
    mpz_poly_srcptr c);

/*
 * Convert a mat_Z to a mat_int64.
 *
 * matrix_int: the output matrix.
 * matrix: the root matrix.
 */
void mat_Z_to_mat_int64(mat_int64_ptr matrix_int, mat_Z_srcptr matrix);

/*
 * Extract a submatrix of matrix_in and put it in matrix_out. The point defined
 *  by (ulx, uly) is the left up corner of the submatrix, and (drx - 1, dry -
 *  1) define the right down corner.
 *
 * matrix_out: the output matrix.
 * matrix_in: the input matrix.
 * ulx, uly: coordinate of the left up corner.
 * drx, dry: coordinate of the right down corner.
 */
void mat_int64_extract(mat_int64_ptr matrix_out, mat_int64_srcptr matrix_in,
    unsigned int ulx, unsigned int uly, unsigned int drx, unsigned int dry);

/*
 * Extract the vector of the column col of a matrix.
 *
 * v: the vector.
 * matrix: the matrix.
 * col: index of the column.
 */
void mat_int64_extract_vector(int64_vector_ptr v, mat_int64_srcptr matrix,
    unsigned int col);

/*
 * Matrix is a diagonal matrix, with the coefficients of x.
 */
void mat_int64_set_diag(mat_int64_ptr matrix, int64_vector_srcptr x);

/*
 * Perform LLL on A and store it in C.
 */
void mat_int64_LLL(mat_int64_ptr C, mat_int64_srcptr A);

/*
 * Perform LLL on the transposition of A and store the transposition in C.
 */
void mat_int64_LLL_transpose(mat_int64_ptr C, mat_int64_srcptr A);

/*
 * Perform LLL on A and store the unimodular matrix in C.
 */
void mat_int64_LLL_unimodular(mat_int64_ptr C, mat_int64_srcptr A);

/*
 * Perform LLL on the transposition of A and store the transposition of the
 *  unimodular matrix in C.
 */
void mat_int64_LLL_unimodular_transpose(mat_int64_ptr C, mat_int64_srcptr A);

/*
 * Fill a matrix with the vectors of list in columns.
 */
void mat_int64_vector_from_list_int64(mat_int64_ptr matrix,
    list_int64_vector_srcptr list);

#ifdef __cplusplus
}
#endif
#endif // MAT_INT64_H
