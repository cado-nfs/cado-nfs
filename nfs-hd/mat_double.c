#include "cado.h" // IWYU pragma: keep
#include "macros.h"
#include <stdio.h>
#include <stdlib.h>
#include "mat_double.h"

void mat_double_init(mat_double_ptr matrix, unsigned int NumRows, unsigned int
                NumCols)
{
  ASSERT(1 <= NumRows);
  ASSERT(1 <= NumCols);

  matrix->NumRows = NumRows;
  matrix->NumCols = NumCols;

  matrix->coeff = (double ** ) malloc(sizeof(double *) * (NumRows + 1));
  for (unsigned int row = 0; row < (NumRows + 1); row++) {
    matrix->coeff[row] = (double * ) malloc(sizeof(double) * (NumCols + 1));
  }
}

void mat_double_init_with_array(mat_double_ptr matrix, unsigned int NumRows,
    unsigned int NumCols, double * coeff)
{
  ASSERT(1 <= NumRows);
  ASSERT(1 <= NumCols);

  matrix->NumRows = NumRows;
  matrix->NumCols = NumCols;

  matrix->coeff = (double ** ) malloc(sizeof(double *) * (NumRows + 1));
  for (unsigned int row = 0; row < (NumRows + 1); row++) {
    matrix->coeff[row] = (double * ) malloc(sizeof(double) * (NumCols + 1));
  }

  for (unsigned int row = 1; row < NumRows + 1; row++) {
    for (unsigned int col = 1; col < NumCols + 1; col++) {
      matrix->coeff[row][col] = coeff[(col - 1) + NumCols * (row - 1)];
    }
  }
}

void mat_double_set_zero(mat_double_ptr matrix)
{
  for (unsigned int row = 1; row < matrix->NumRows + 1; row++) {
    for (unsigned int col = 1; col < matrix->NumCols + 1; col++) {
      matrix->coeff[row][col] = 0;
    }
  }
}

void mat_double_copy(mat_double_ptr B, mat_double_srcptr A)
{
  ASSERT(B->NumRows == A->NumRows);
  ASSERT(B->NumCols == A->NumCols);

  for (unsigned int row = 1; row < B->NumRows + 1; row++) {
    for (unsigned int col = 1; col < B->NumCols + 1; col++) {
      B->coeff[row][col] = A->coeff[row][col];
    }
  }
}

void mat_double_clear(mat_double_ptr matrix)
{
  for (unsigned int i = 0; i < (matrix->NumRows + 1); i++) {
    free(matrix->coeff[i]);
  }
  free(matrix->coeff);
}

void mat_double_fprintf(FILE * file, mat_double_srcptr matrix)
{
  ASSERT(matrix->NumRows >= 1);
  ASSERT(matrix->NumCols >= 1);

  fprintf(file, "[");
  for (unsigned int row = 1; row < matrix->NumRows; row++) {
    fprintf(file, "[");
    for (unsigned int col = 1; col < matrix->NumCols; col++) {
      fprintf(file, "%f, ", matrix->coeff[row][col]);
    }
    fprintf(file, "%f],\n", matrix->coeff[row][matrix->NumCols]);
  }
  fprintf(file, "[");
  for (unsigned int col = 1; col < matrix->NumCols; col++) {
    fprintf(file, "%f, ", matrix->coeff[matrix->NumRows][col]);
  }
  fprintf(file, "%f]]\n",
      matrix->coeff[matrix->NumRows][matrix->NumCols]);
}

void mat_double_fprintf_comment(FILE * file, mat_double_srcptr matrix)
{
  ASSERT(matrix->NumRows >= 1);
  ASSERT(matrix->NumCols >= 1);

  fprintf(file, "# [");
  fprintf(file, "[");
  for (unsigned int col = 1; col < matrix->NumCols; col++) {
    fprintf(file, "%f, ", matrix->coeff[1][col]);
  }
  fprintf(file, "%f],\n", matrix->coeff[1][matrix->NumCols]);

  for (unsigned int row = 2; row < matrix->NumRows; row++) {
    fprintf(file, "# [");
    for (unsigned int col = 1; col < matrix->NumCols; col++) {
      fprintf(file, "%f, ", matrix->coeff[row][col]);
    }
    fprintf(file, "%f],\n", matrix->coeff[row][matrix->NumCols]);
  }
  fprintf(file, "# [");
  for (unsigned int col = 1; col < matrix->NumCols; col++) {
    fprintf(file, "%f, ", matrix->coeff[matrix->NumRows][col]);
  }
  fprintf(file, "%f]]\n",
      matrix->coeff[matrix->NumRows][matrix->NumCols]);
}

void mat_double_transpose(mat_double_ptr matrix, mat_double_srcptr matrix_src)
{
  ASSERT(matrix->NumRows == matrix_src->NumCols);
  ASSERT(matrix->NumCols == matrix_src->NumRows);
  
  mat_double_t tmp;
  mat_double_init(tmp, matrix_src->NumRows, matrix_src->NumCols);
  mat_double_copy(tmp, matrix_src);
  for (unsigned int row = 1; row < matrix->NumRows + 1; row++) {
    for (unsigned int col = 1; col < matrix->NumCols + 1; col++) {
      matrix->coeff[col][row] = tmp->coeff[row][col];
    }
  }
  mat_double_clear(tmp);
}
