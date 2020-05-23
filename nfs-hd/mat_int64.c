#include "cado.h" // IWYU pragma: keep
#include "macros.h"
#include "mat_Z.h"
#include <stdio.h>
#include <stdlib.h>
#include "int64_vector.h"
#include "int64_poly.h"
#include "mat_int64.h"
#include <inttypes.h>
#include "utils_int64.h"

void mat_int64_init(mat_int64_ptr matrix, unsigned int NumRows, unsigned int
                NumCols)
{
  ASSERT(1 <= NumRows);
  ASSERT(1 <= NumCols);

  matrix->NumRows = NumRows;
  matrix->NumCols = NumCols;

  matrix->coeff = (int64_t ** ) malloc(sizeof(int64_t *) * (NumRows + 1));
  for (unsigned int row = 0; row < (NumRows + 1); row++) {
    matrix->coeff[row] = (int64_t * ) malloc(sizeof(int64_t) * (NumCols + 1));
  }
}

void mat_int64_init_with_array(mat_int64_ptr matrix, unsigned int NumRows,
    unsigned int NumCols, int64_t * coeff)
{
  ASSERT(1 <= NumRows);
  ASSERT(1 <= NumCols);

  matrix->NumRows = NumRows;
  matrix->NumCols = NumCols;

  matrix->coeff = (int64_t ** ) malloc(sizeof(int64_t *) * (NumRows + 1));
  for (unsigned int row = 0; row < (NumRows + 1); row++) {
    matrix->coeff[row] = (int64_t * ) malloc(sizeof(int64_t) * (NumCols + 1));
  }

  for (unsigned int row = 1; row < NumRows + 1; row++) {
    for (unsigned int col = 1; col < NumCols + 1; col++) {
      matrix->coeff[row][col] = coeff[(col - 1) + NumCols * (row - 1)];
    }
  }
}

void mat_int64_set_zero(mat_int64_ptr matrix)
{
  for (unsigned int row = 1; row < matrix->NumRows + 1; row++) {
    for (unsigned int col = 1; col < matrix->NumCols + 1; col++) {
      matrix->coeff[row][col] = 0;
    }
  }
}

void mat_int64_set_identity(mat_int64_ptr matrix)
{
  for (unsigned int row = 1; row < matrix->NumRows + 1; row++) {
    for (unsigned int col = 1; col < matrix->NumCols + 1; col++) {
      if (row == col) {
        matrix->coeff[row][col] = 1;
      } else {
        matrix->coeff[row][col] = 0;
      }
    }
  }
}

unsigned int mat_int64_equal(mat_int64_srcptr A, mat_int64_srcptr B)
{
  if (A->NumCols != B->NumCols) {
    return 0;
  }
  if (A->NumRows != B->NumRows) {
    return 0;
  }
  for (unsigned int col = 1; col <= A->NumCols; col++) {
    for (unsigned int row = 1; row <= A->NumRows; row++) {
      if (A->coeff[row][col] != B->coeff[row][col]) {
        return 0;
      }
    }
  }
  return 1;
}

void mat_int64_copy(mat_int64_ptr B, mat_int64_srcptr A)
{
  ASSERT(B->NumRows == A->NumRows);
  ASSERT(B->NumCols == A->NumCols);

  for (unsigned int row = 1; row < B->NumRows + 1; row++) {
    for (unsigned int col = 1; col < B->NumCols + 1; col++) {
      B->coeff[row][col] = A->coeff[row][col];
    }
  }
}

void mat_int64_clear(mat_int64_ptr matrix)
{
  for (unsigned int i = 0; i < (matrix->NumRows + 1); i++) {
    free(matrix->coeff[i]);
  }
  free(matrix->coeff);
}

void mat_int64_fprintf(FILE * file, mat_int64_srcptr matrix)
{
  ASSERT(matrix->NumRows >= 1);
  ASSERT(matrix->NumCols >= 1);

  fprintf(file, "[");
  for (unsigned int row = 1; row < matrix->NumRows; row++) {
    fprintf(file, "[");
    for (unsigned int col = 1; col < matrix->NumCols; col++) {
      fprintf(file, "%" PRId64 ", ", matrix->coeff[row][col]);
    }
    fprintf(file, "%" PRId64 "],\n", matrix->coeff[row][matrix->NumCols]);
  }
  fprintf(file, "[");
  for (unsigned int col = 1; col < matrix->NumCols; col++) {
    fprintf(file, "%" PRId64 ", ", matrix->coeff[matrix->NumRows][col]);
  }
  fprintf(file, "%" PRId64 "]]\n",
      matrix->coeff[matrix->NumRows][matrix->NumCols]);
}

void mat_int64_fprintf_comment(FILE * file, mat_int64_srcptr matrix)
{
  ASSERT(matrix->NumRows >= 1);
  ASSERT(matrix->NumCols >= 1);

  fprintf(file, "# [");
  fprintf(file, "[");
  for (unsigned int col = 1; col < matrix->NumCols; col++) {
    fprintf(file, "%" PRId64 ", ", matrix->coeff[1][col]);
  }
  fprintf(file, "%" PRId64 "],\n", matrix->coeff[1][matrix->NumCols]);

  for (unsigned int row = 2; row < matrix->NumRows; row++) {
    fprintf(file, "# [");
    for (unsigned int col = 1; col < matrix->NumCols; col++) {
      fprintf(file, "%" PRId64 ", ", matrix->coeff[row][col]);
    }
    fprintf(file, "%" PRId64 "],\n", matrix->coeff[row][matrix->NumCols]);
  }
  fprintf(file, "# [");
  for (unsigned int col = 1; col < matrix->NumCols; col++) {
    fprintf(file, "%" PRId64 ", ", matrix->coeff[matrix->NumRows][col]);
  }
  fprintf(file, "%" PRId64 "]]\n",
      matrix->coeff[matrix->NumRows][matrix->NumCols]);
}

void mat_int64_transpose(mat_int64_ptr matrix, mat_int64_srcptr matrix_src)
{
  ASSERT(matrix->NumRows == matrix_src->NumCols);
  ASSERT(matrix->NumCols == matrix_src->NumRows);
  
  mat_int64_t tmp;
  mat_int64_init(tmp, matrix_src->NumRows, matrix_src->NumCols);
  mat_int64_copy(tmp, matrix_src);
  for (unsigned int row = 1; row < matrix->NumRows + 1; row++) {
    for (unsigned int col = 1; col < matrix->NumCols + 1; col++) {
      matrix->coeff[col][row] = tmp->coeff[row][col];
    }
  }
  mat_int64_clear(tmp);
}

void mat_int64_mul_mat_int64(mat_int64_ptr C, mat_int64_srcptr A, 
    mat_int64_srcptr B)
{
  ASSERT(A->NumCols == B->NumRows);
  ASSERT(A->NumRows == C->NumRows);
  ASSERT(B->NumCols == C->NumCols);

  mat_int64_set_zero(C);
  mat_int64_t A_tmp;
  mat_int64_init(A_tmp, A->NumRows, A->NumCols);
  mat_int64_copy(A_tmp, A);
  mat_int64_t B_tmp;
  mat_int64_init(B_tmp, B->NumRows, B->NumCols);
  mat_int64_copy(B_tmp, B);
  for (unsigned int row = 1; row < A_tmp->NumRows + 1; row++) {
    for (unsigned int col = 1; col < B_tmp->NumCols + 1; col++) {
      for (unsigned int k = 1; k < A_tmp->NumCols + 1; k++) {
        C->coeff[row][col] = C->coeff[row][col] +
          A_tmp->coeff[row][k] * B_tmp->coeff[k][col];
      }
    }
  }
  mat_int64_clear(A_tmp);
  mat_int64_clear(B_tmp);
}

void mat_int64_mul_int64_vector_to_int64_poly(int64_poly_ptr a,
    mat_int64_srcptr A, int64_vector_srcptr c)
{
  ASSERT(A->NumCols == c->dim);
  
  int64_t tmp = 0;
  for (unsigned int row = 1; row < A->NumRows + 1; row++) {
    tmp = 0;
    for (unsigned int col = 1; col < A->NumCols + 1; col++) {
      tmp = tmp + A->coeff[row][col] * c->c[col - 1];
    }
    int64_poly_setcoeff(a, row - 1, tmp);
  }
}

void mat_int64_mul_int64_vector(int64_vector_ptr a, mat_int64_srcptr A,
    int64_vector_srcptr c)
{
  ASSERT(A->NumCols == c->dim);
  ASSERT(A->NumRows == a->dim);

  int64_vector_t v;
  int64_vector_init(v, A->NumRows);
  int64_t tmp = 0;
  for (unsigned int row = 1; row < A->NumRows + 1; row++) {
    tmp = 0;
    for (unsigned int col = 1; col < A->NumCols + 1; col++) {
      tmp = tmp + A->coeff[row][col] * c->c[col - 1];
    }
    int64_vector_setcoordinate(v, row - 1, tmp);
  }
  int64_vector_set(a, v);
  int64_vector_clear(v);
}

//TODO: problem here if c has a deg not equal to A->NumCols - 1.
void mat_int64_mul_int64_poly(int64_poly_ptr a, mat_int64_srcptr A,
    int64_poly_srcptr c)
{
  ASSERT(A->NumCols >= (unsigned int) (c->deg + 1));
 
  int64_poly_t p_tmp;
  int64_poly_init(p_tmp, c->deg);
  int64_t tmp = 0;
  for (unsigned int row = 1; row < A->NumRows + 1; row++) {
    tmp = 0;
    for (unsigned int col = 1; col < (unsigned int) (c->deg + 2); col++) {
      if (c->coeff[col - 1] != 0) {
        tmp = tmp + A->coeff[row][col] * c->coeff[col - 1];
      }
    }
    int64_poly_setcoeff(p_tmp, row - 1, tmp);
  }
  int64_poly_copy(a, p_tmp);
  int64_poly_clear(p_tmp);
}

void mat_Z_to_mat_int64(mat_int64_ptr matrix_int, mat_Z_srcptr matrix)
{
  ASSERT(matrix_int->NumRows == matrix->NumRows);
  ASSERT(matrix_int->NumCols == matrix->NumCols);

  for (unsigned int row = 1; row < matrix->NumRows + 1; row++) {
    for (unsigned int col = 1; col < matrix->NumCols + 1; col++) {
      matrix_int->coeff[row][col] = mpz_get_si(matrix->coeff[row][col]);
    }
  }
}

void mat_int64_extract(mat_int64_ptr matrix_out, mat_int64_srcptr matrix_in,
    unsigned int ulx, unsigned int uly, unsigned int drx, unsigned int dry)
{
  ASSERT(drx - ulx > 0);
  ASSERT(dry - uly > 0);

  mat_int64_init(matrix_out, drx - ulx, dry - uly);

#ifndef NDEBUG
  mat_int64_set_zero(matrix_out);
#endif // NDEBUG

  for (unsigned int i = 0; i < drx - ulx; i++) {
    for (unsigned int j = 0; j < dry - uly; j++) {
      matrix_out->coeff[i + 1][j + 1] =
        matrix_in->coeff[i + ulx + 1][j + uly + 1];
    }
  }
}

void mat_int64_extract_vector(int64_vector_ptr v, mat_int64_srcptr matrix,
    unsigned int col)
{
  ASSERT(col < matrix->NumCols);
  ASSERT(v->dim == matrix->NumRows);

  for (unsigned int i = 1; i <= matrix->NumRows; i++) {
    v->c[i - 1] = matrix->coeff[i][col + 1];
  }
}

void mat_int64_set_diag(mat_int64_ptr matrix, int64_vector_srcptr x)
{
  ASSERT(matrix->NumRows == matrix->NumCols);
  ASSERT(x->dim == matrix->NumRows);

  mat_int64_set_zero(matrix);
  for (unsigned int i = 0; i < matrix->NumCols; i++) {
    matrix->coeff[i + 1][i + 1] = x->c[i];
  }
}

void mat_int64_LLL(mat_int64_ptr C, mat_int64_srcptr A)
{
  ASSERT(A->NumRows == C->NumRows);
  ASSERT(A->NumCols == C->NumCols);

  mat_Z_t A_Z;
  mat_Z_init(A_Z, A->NumRows, A->NumCols);
  mat_int64_to_mat_Z(A_Z, A);

  mat_Z_t C_Z;
  mat_Z_init(C_Z, C->NumRows, C->NumCols);

  mat_Z_LLL(C_Z, A_Z);

  mat_Z_to_mat_int64(C, C_Z);

  mat_Z_clear(C_Z);
  mat_Z_clear(A_Z);
}

void mat_int64_LLL_transpose(mat_int64_ptr C, mat_int64_srcptr A)
{
  ASSERT(A->NumRows == C->NumRows);
  ASSERT(A->NumCols == C->NumCols);

  mat_int64_t A_tmp;
  mat_int64_init(A_tmp, A->NumRows, A->NumCols);
  mat_int64_transpose(A_tmp, A);

  mat_int64_LLL(C, A_tmp);

  mat_int64_transpose(C, C);
  mat_int64_clear(A_tmp);
}

void mat_int64_LLL_unimodular(mat_int64_ptr C, mat_int64_srcptr A)
{
  ASSERT(A->NumRows == C->NumRows);
  ASSERT(A->NumCols == C->NumCols);

  mat_Z_t A_Z;
  mat_Z_init(A_Z, A->NumRows, A->NumCols);
  mat_int64_to_mat_Z(A_Z, A);

  mat_Z_t C_Z;
  mat_Z_init(C_Z, C->NumRows, C->NumCols);

  mat_Z_LLL_unimodular(C_Z, A_Z);

  mat_Z_to_mat_int64(C, C_Z);

  mat_Z_clear(C_Z);
  mat_Z_clear(A_Z);
}

void mat_int64_LLL_unimodular_transpose(mat_int64_ptr C, mat_int64_srcptr A)
{
  ASSERT(A->NumRows == C->NumRows);
  ASSERT(A->NumCols == C->NumCols);

  mat_int64_t A_tmp;
  mat_int64_init(A_tmp, A->NumRows, A->NumCols);
  mat_int64_transpose(A_tmp, A);

  mat_int64_LLL_unimodular(C, A_tmp);

  mat_int64_transpose(C, C);
  mat_int64_clear(A_tmp);
}


void mat_int64_from_list_int64_vector(mat_int64_ptr matrix,
    list_int64_vector_srcptr list)
{
  ASSERT(list->length == matrix->NumCols);
  ASSERT(list->v[0]->dim == matrix->NumRows);

  for (unsigned int col = 0; col < list->length; col++) {
    for (unsigned int row = 0; row < list->v[0]->dim; row++) {
      matrix->coeff[row + 1][col + 1] = list->v[col]->c[row];
    }
  }
}
