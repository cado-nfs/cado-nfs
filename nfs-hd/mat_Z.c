#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "macros.h"
#include "mat_Z.h"
#include "lll.h"
#include "mpz_vector.h"

void mat_Z_init(mat_Z_ptr matrix, unsigned int NumRows, unsigned int NumCols)
{
  ASSERT(1 <= NumRows);
  ASSERT(1 <= NumCols);

  matrix->NumRows = NumRows;
  matrix->NumCols = NumCols;

  matrix->coeff = (mpz_t ** ) malloc(sizeof(mpz_t *) * (NumRows + 1));
  for (unsigned int row = 0; row < (NumRows + 1); row++) {
    matrix->coeff[row] = (mpz_t * )malloc(sizeof(mpz_t) * (NumCols + 1));
  }

  for (unsigned int row = 0; row < NumRows + 1; row++) {
    for (unsigned int col = 0; col < NumCols + 1; col++) {
      mpz_init(matrix->coeff[row][col]);
    }
  }
}

void mat_Z_init_with_array(mat_Z_ptr matrix, unsigned int NumRows,
    unsigned int NumCols, mpz_t * coeff)
{
  ASSERT(1 <= NumRows);
  ASSERT(1 <= NumCols);

  matrix->NumRows = NumRows;
  matrix->NumCols = NumCols;

  matrix->coeff = (mpz_t ** )malloc(sizeof(mpz_t *) * (NumRows + 1));
  for (unsigned int row = 0; row < (NumRows + 1); row++) {
    matrix->coeff[row] = (mpz_t * ) malloc(sizeof(mpz_t) * (NumCols + 1));
  }

  for (unsigned int row = 0; row < NumRows + 1; row++) {
    for (unsigned int col = 0; col < NumCols + 1; col++) {
      mpz_init(matrix->coeff[row][col]);
    }
  }

  for (unsigned int row = 1; row < NumRows + 1; row++) {
    for (unsigned int col = 1; col < NumCols + 1; col++) {
      mpz_set(matrix->coeff[row][col], coeff[(col - 1) + NumCols * (row - 1)]);
    }
  }
}

void mat_Z_copy(mat_Z_ptr B, mat_Z_srcptr A)
{
  ASSERT(B->NumRows == A->NumRows);
  ASSERT(B->NumCols == A->NumCols);

  for (unsigned int row = 1; row < B->NumRows + 1; row++) {
    for (unsigned int col = 1; col < B->NumCols + 1; col++) {
      mpz_set(B->coeff[row][col], A->coeff[row][col]);
    }
  }
}

void mat_Z_set_coeff(mat_Z_ptr matrix, mpz_t i, unsigned int row,
    unsigned int col)
{
  ASSERT(row <= matrix->NumRows);
  ASSERT(col <= matrix->NumCols);
  ASSERT(1 <= row);
  ASSERT(1 <= col);

  mpz_set(matrix->coeff[row][col], i);
}

void mat_Z_set_coeff_int64(mat_Z_ptr matrix, int64_t i, unsigned int row,
                           unsigned int  col)
{
  ASSERT(row <= matrix->NumRows);
  ASSERT(col <= matrix->NumCols);
  ASSERT(1 <= row);
  ASSERT(1 <= col);

  mpz_set_si(matrix->coeff[row][col], i);
}

void mat_Z_set_coeff_uint64(mat_Z_ptr matrix, uint64_t i, unsigned int row,
                            unsigned int  col)
{
  ASSERT(row <= matrix->NumRows);
  ASSERT(col <= matrix->NumCols);
  ASSERT(1 <= row);
  ASSERT(1 <= col);

  mpz_set_ui(matrix->coeff[row][col], i);
}

void mat_Z_set_zero(mat_Z_ptr matrix)
{
  for (unsigned int i = 0; i < matrix->NumRows + 1; i++) {
    for (unsigned int j = 0; j < matrix->NumCols + 1; j++) {
      mpz_set_ui(matrix->coeff[i][j], 0);
    }
  }
}

void mat_Z_clear(mat_Z_ptr matrix)
{
  for (unsigned int i = 0; i < matrix->NumRows + 1; i++) {
    for (unsigned int j = 0; j < matrix->NumCols + 1; j++) {
      mpz_clear(matrix->coeff[i][j]);
    }
  }

  for (unsigned int i = 0; i < (matrix->NumRows + 1); i++) {
    free(matrix->coeff[i]);
  }
  free(matrix->coeff);
}

void mat_Z_fprintf(FILE * file, mat_Z_srcptr matrix)
{
  fprintf(file, "[");
  for (unsigned int row = 1; row < matrix->NumRows; row++) {
    fprintf(file, "[");
    for (unsigned int col = 1; col < matrix->NumCols; col++) {
      gmp_fprintf(file, "%Zd, ", matrix->coeff[row][col]);
    }
    gmp_fprintf(file, "%Zd],\n", matrix->coeff[row][matrix->NumCols]);
  }
  fprintf(file, "[");
  for (unsigned int col = 1; col < matrix->NumCols; col++) {
    gmp_fprintf(file, "%Zd, ", matrix->coeff[matrix->NumRows][col]);
  }
  gmp_fprintf(file, "%Zd]]\n", matrix->coeff[matrix->NumRows][matrix->NumCols]);
}

void mat_Z_fprintf_comment(FILE * file, mat_Z_srcptr matrix)
{
  fprintf(file, "# [");
  fprintf(file, "[");
  for (unsigned int col = 1; col < matrix->NumCols; col++) {
    gmp_fprintf(file, "%Zd, ", matrix->coeff[1][col]);
  }
  gmp_fprintf(file, "%Zd],\n", matrix->coeff[1][matrix->NumCols]);
  for (unsigned int row = 2; row < matrix->NumRows; row++) {
    fprintf(file, "# [");
    for (unsigned int col = 1; col < matrix->NumCols; col++) {
      gmp_fprintf(file, "%Zd, ", matrix->coeff[row][col]);
    }
    gmp_fprintf(file, "%Zd],\n", matrix->coeff[row][matrix->NumCols]);
  }
  fprintf(file, "# [");
  for (unsigned int col = 1; col < matrix->NumCols; col++) {
    gmp_fprintf(file, "%Zd, ", matrix->coeff[matrix->NumRows][col]);
  }
  gmp_fprintf(file, "%Zd]]\n", matrix->coeff[matrix->NumRows][matrix->NumCols]);
}

void mat_Z_transpose(mat_Z_ptr matrix, mat_Z_srcptr matrix_src)
{
  ASSERT(matrix->NumRows == matrix_src->NumCols);
  ASSERT(matrix->NumCols == matrix_src->NumRows);

  mat_Z_t tmp;
  mat_Z_init(tmp, matrix_src->NumRows, matrix_src->NumCols);
  mat_Z_copy(tmp, matrix_src);
  for (unsigned int row = 1; row < matrix->NumRows + 1; row++) {
    for (unsigned int col = 1; col < matrix->NumCols + 1; col++) {
      mpz_set(matrix->coeff[col][row], tmp->coeff[row][col]);
    }
  }
  mat_Z_clear(tmp);
}

void mat_Z_mul_mat_Z(mat_Z_ptr C, mat_Z_srcptr A, mat_Z_srcptr B)
{
  ASSERT(A->NumCols == B->NumRows);
  ASSERT(A->NumRows == C->NumRows);
  ASSERT(B->NumCols == C->NumCols);

  mat_Z_t tmp;
  mat_Z_init(tmp, C->NumRows, C->NumCols);

  for (unsigned int row = 1; row < A->NumRows + 1; row++) {
    for (unsigned int col = 1; col < B->NumCols + 1; col++) {
      for (unsigned int k = 1; k < A->NumCols + 1; k++) {
        mpz_addmul(tmp->coeff[row][col], A->coeff[row][k], B->coeff[k][col]);
      }
    }
  }
  mat_Z_copy(C, tmp);

  mat_Z_clear(tmp);
}

void mat_Z_mul_mpz_vector_to_mpz_poly(mpz_poly_ptr a, mat_Z_srcptr A,
                                      mpz_vector_srcptr c)
{
  ASSERT(A->NumCols == c->dim);

  for (unsigned int row = 1; row < A->NumRows + 1; row++) {
    mpz_t tmpz;
    mpz_init(tmpz);
    for (unsigned int col = 1; col < A->NumCols + 1; col++) {
      mpz_addmul(tmpz, A->coeff[row][col], c->c[col - 1]);
    }
    mpz_poly_setcoeff(a, row - 1, tmpz);
    mpz_clear(tmpz);
  }
}

void mat_Z_mul_mpz_poly(mpz_poly_ptr a, mat_Z_srcptr A, mpz_poly_srcptr c)
{
  ASSERT(A->NumCols >= (unsigned int) (c->deg + 1));

  for (unsigned int row = 1; row < A->NumRows + 1; row++) {
    mpz_t tmpz;
    mpz_init(tmpz);
    for (unsigned int col = 1; col < (unsigned int) (c->deg + 2); col++) {
      mpz_addmul(tmpz, A->coeff[row][col], c->coeff[col - 1]);
    }
    mpz_poly_setcoeff(a, row - 1, tmpz);
    mpz_clear(tmpz);
  }
}

void mat_Z_LLL(mat_Z_ptr C, mat_Z_srcptr A)
{
  ASSERT(A->NumRows == C->NumRows);
  ASSERT(A->NumCols == C->NumCols);

  mpz_t a;
  mpz_t b;
  mpz_t det;
  mpz_init(a);
  mpz_init(b);
  mpz_init(det);

  mpz_set_ui(a, 3);
  mpz_set_ui(b, 4);

  mat_Z B;
  LLL_init(&B, A->NumRows, A->NumCols);
  for (unsigned int row = 1; row < A->NumRows + 1; row++) {
    for (unsigned int col = 1; col < A->NumCols + 1; col++) {
      mpz_set(B.coeff[row][col], A->coeff[row][col]);
    }
  }

  LLL(det, B, NULL, a, b);

  ASSERT(B.NumRows == (int)C->NumRows);
  ASSERT(B.NumCols == (int)C->NumCols);

  for (unsigned int row = 1; row < A->NumRows + 1; row++) {
    for (unsigned int col = 1; col < A->NumCols + 1; col++) {
      mpz_set(C->coeff[row][col], B.coeff[row][col]);
    }
  }

  mpz_clear(a);
  mpz_clear(b);
  mpz_clear(det);
  LLL_clear(&B);
}

void mat_Z_LLL_transpose(mat_Z_ptr matrix_new, mat_Z_srcptr matrix_old)
{
  mat_Z_transpose(matrix_new, matrix_old);
  mat_Z_LLL(matrix_new, matrix_new);
  mat_Z_transpose(matrix_new, matrix_new);
}

void mat_Z_LLL_unimodular(mat_Z_ptr C, mat_Z_srcptr A)
{
  mpz_t a;
  mpz_t b;
  mpz_t det;
  mpz_init(a);
  mpz_init(b);
  mpz_init(det);

  mpz_set_ui(a, 3);
  mpz_set_ui(b, 4);

  mat_Z B;
  LLL_init(&B, A->NumRows, A->NumCols);
  mat_Z U;
  LLL_init(&U, A->NumRows, A->NumCols);
  for (unsigned int row = 1; row < A->NumRows + 1; row++) {
    for (unsigned int col = 1; col < A->NumCols + 1; col++) {
      mpz_set(B.coeff[row][col], A->coeff[row][col]);
    }
  }

  LLL(det, B, &U, a, b);

  ASSERT(B.NumRows == (int)C->NumRows);
  ASSERT(B.NumCols == (int)C->NumCols);

  for (unsigned int row = 1; row < A->NumRows + 1; row++) {
    for (unsigned int col = 1; col < A->NumCols + 1; col++) {
      mpz_set(C->coeff[row][col], U.coeff[row][col]);
    }
  }

  mpz_clear(a);
  mpz_clear(b);
  mpz_clear(det);
  LLL_clear(&B);
  LLL_clear(&U);
}

void mat_Z_LLL_unimodular_transpose(mat_Z_ptr mat_new, mat_Z_srcptr mat_old)
{
  mat_Z_transpose(mat_new, mat_old);
  mat_Z_LLL_unimodular(mat_new, mat_new);
  mat_Z_transpose(mat_new, mat_new);
}

void mat_int64_to_mat_Z(mat_Z_ptr mat_Z, mat_int64_srcptr mat_int)
{
  ASSERT(mat_Z->NumRows == mat_int->NumRows);
  ASSERT(mat_Z->NumCols == mat_int->NumCols);

  for (unsigned int row = 1; row < mat_int->NumRows + 1; row++) {
    for (unsigned int col = 1; col < mat_int->NumCols + 1; col++) {
      mpz_set_si(mat_Z->coeff[row][col], mat_int->coeff[row][col]);
    }
  }
}

static int compare_last(const void * p0, const void * p1)
{
  const mpz_vector_t * v0 = (const mpz_vector_t * ) p0;
  const mpz_vector_t * v1 = (const mpz_vector_t * ) p1;

  return mpz_cmp((*v0)->c[(*v0)->dim - 1], (*v1)->c[(*v1)->dim - 1]);
}

void mat_Z_sort_last(mat_Z_ptr M_out, mat_Z_srcptr M_in)
{
  mpz_vector_t * v = (mpz_vector_t *) malloc(sizeof(mpz_vector_t) *
      M_in->NumCols);
  for (unsigned int col = 1; col <= M_in->NumCols; col++) {
    mpz_vector_init(v[col - 1], M_in->NumRows);
    for (unsigned int row = 1; row <= M_in->NumRows; row++) {
      mpz_set(v[col - 1]->c[row - 1], M_in->coeff[row][col]);
    }
  }

  qsort(v, M_in->NumCols, sizeof(v[0]), compare_last);

  for (unsigned int col = 1; col <= M_in->NumCols; col++) {
    for (unsigned int row = 1; row <= M_in->NumRows; row++) {
      mpz_set(M_out->coeff[row][col], v[col - 1]->c[row - 1]);
    }
  }

  for (unsigned int col = 0; col < M_in->NumCols; col++) {
    mpz_vector_clear(v[col]);
  }
  free(v);
}

void mat_Z_set_diag(mat_Z_ptr M, mpz_vector_srcptr diag)
{
  ASSERT(diag->dim == M->NumCols);
  ASSERT(diag->dim == M->NumRows);

  mat_Z_set_zero(M);

  for (unsigned int i = 1; i <= diag->dim; i++) {
    mpz_set(M->coeff[i][i], diag->c[i - 1]);
  }
}

void mat_Z_skew_LLL(mat_Z_ptr MSLLL, mat_Z_srcptr M_root,
    mpz_vector_srcptr skewness)
{
  ASSERT(skewness->dim == M_root->NumCols);
  ASSERT(M_root->NumRows == M_root->NumCols);
  ASSERT(MSLLL->NumRows == M_root->NumCols);
  ASSERT(MSLLL->NumRows == M_root->NumCols);

  mat_Z_t I_s;
  mat_Z_init(I_s, M_root->NumRows, M_root->NumCols);
  mat_Z_set_diag(I_s, skewness);

  mat_Z_t M;
  mat_Z_init(M, M_root->NumRows, M_root->NumCols);
  mat_Z_mul_mat_Z(M, I_s, M_root);

  mat_Z_t U;
  mat_Z_init(U, M_root->NumRows, M_root->NumCols);
  mat_Z_LLL_unimodular_transpose(U, M);

  mat_Z_mul_mat_Z(MSLLL, M_root, U);

  mat_Z_clear(U);
  mat_Z_clear(I_s);
  mat_Z_clear(M);
}
