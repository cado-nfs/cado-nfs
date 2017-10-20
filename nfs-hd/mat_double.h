#ifndef MAT_DOUBLE_H
#define MAT_DOUBLE_H 

#include "macros.h"
#include "matrix.h"

/* TODO: STL or something */
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
void mat_double_init(mat_double_ptr matrix, unsigned int NumRows,
    unsigned int NumCols);

/*
 * Initialize a matrix and set the coefficients. Assume there are enough
 *  coefficients.
 *
 * matrix: the matrix.
 * NumRows: number of rows in the matrix.
 * NumCols: number of cols in the matrix.
 * coeff: an array with the coefficient of the matrix, given in the order to
 *  fill the row.
 */
void mat_double_init_with_array(mat_double_ptr matrix, unsigned int NumRows,
    unsigned int NumCols, double * coeff);

/*
 * All the coeffecients of the matrix are equal to 0.
 *
 * matrix: the matrix.
 */
void mat_double_set_zero(mat_double_ptr matrix);

/*
 * Copy A in B. Assume that B and A have the same size.
 *
 * B: the new matrix.
 * A: the root matrix.
 */
void mat_double_copy(mat_double_ptr B, mat_double_srcptr A);

/*
 * Set a coefficient to i, this coefficient is in the row (row) and the column
 *  (col). Row and col are defined without the offset on rows and columns.
 *
 * matrix: the matrix.
 * i: the new coefficient.
 * row: number for row.
 * col: number for column.
 */
static inline void mat_double_set_coeff(mat_double_ptr matrix, double i,
    unsigned int row, unsigned int col)
{
  ASSERT(row < matrix->NumRows);
  ASSERT(1 <= row);
  ASSERT(col < matrix->NumCols);
  ASSERT(1 <= col);

  matrix->coeff[row + 1][col + 1] = i;
}

/*
 * Return the coeffecient at postition [row, col]. Row and col are defined
 *  without the offset on rows and columns.
 */
static inline double mat_double_get_coeff(mat_double_srcptr matrix,
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
 * matrix: mat_double_ptr, the matrix.
 */
void mat_double_clear(mat_double_ptr matrix);

/*
 * Write a matrix in a file.
 *
 * file: the file.
 * matrix: the matrix.
 */
void mat_double_fprintf(FILE * file, mat_double_srcptr matrix);

/*
 * Write a matrix with # before each line.
 *
 * file: the file.
 * matrix: the matrix.
 */
void mat_double_fprintf_comment(FILE * file, mat_double_srcptr matrix);

/*
 * Transpose a matrix_src, and set the result in matrix. Assume that matrix have
 *  the good size and is already initialized.
 *
 * matrix: the matrix equals to the transpose of matrix_src.
 * matrix_src: the root matrix.
 */
void mat_double_transpose(mat_double_ptr matrix, mat_double_srcptr matrix_src);

#ifdef __cplusplus
}
#endif
#endif /* MAT_DOUBLE_H */
