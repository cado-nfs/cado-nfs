#include "cado.h"
#include "mat_int64.h"
#include "timing.h"

int main()
{
  //[[0, 2, 3],[4, 5, 6],[7, 8, 10]]
  mat_Z_t A_Z;
  mat_Z_init(A_Z, 3, 3);
  mat_Z_set_coeff_int64(A_Z, 0, 1, 1);
  mat_Z_set_coeff_int64(A_Z, 2, 1, 2);
  mat_Z_set_coeff_int64(A_Z, 3, 1, 3);
  mat_Z_set_coeff_int64(A_Z, 4, 2, 1);
  mat_Z_set_coeff_int64(A_Z, 5, 2, 2);
  mat_Z_set_coeff_int64(A_Z, 6, 2, 3);
  mat_Z_set_coeff_int64(A_Z, 7, 3, 1);
  mat_Z_set_coeff_int64(A_Z, 8, 3, 2);
  mat_Z_set_coeff_int64(A_Z, 10, 3, 3);
  mat_Z_t C_Z;
  mat_Z_init(C_Z, 3, 3);
  mat_Z_LLL(C_Z, A_Z);
  mat_Z_clear(A_Z);

  mat_int64_t A;
  mat_int64_init(A, 3, 3);
  mat_int64_set_identity(A);
  A->coeff[1][1] = 0;
  A->coeff[1][2] = 2;
  A->coeff[1][3] = 3;
  A->coeff[2][1] = 4;
  A->coeff[2][2] = 5;
  A->coeff[2][3] = 6;
  A->coeff[3][1] = 7;
  A->coeff[3][2] = 8;
  A->coeff[3][3] = 10;

  mat_int64_t C;
  mat_int64_init(C, 3, 3);
  mat_int64_set_zero(C);
  mat_int64_LLL(C, A);

  mat_Z_to_mat_int64(A, C_Z);

  ASSERT_ALWAYS(mat_int64_equal(A, C));

  mat_int64_clear(A);
  mat_int64_clear(C);
  mat_Z_clear(C_Z);

  //[[1, 0, 0],[0, 1, 0],[0, 0, 1]]
  mat_Z_init(A_Z, 3, 3);
  mat_Z_set_coeff_int64(A_Z, 1, 1, 1);
  mat_Z_set_coeff_int64(A_Z, 0, 1, 2);
  mat_Z_set_coeff_int64(A_Z, 0, 1, 3);
  mat_Z_set_coeff_int64(A_Z, 0, 2, 1);
  mat_Z_set_coeff_int64(A_Z, 1, 2, 2);
  mat_Z_set_coeff_int64(A_Z, 0, 2, 3);
  mat_Z_set_coeff_int64(A_Z, 0, 3, 1);
  mat_Z_set_coeff_int64(A_Z, 0, 3, 2);
  mat_Z_set_coeff_int64(A_Z, 1, 3, 3);
  mat_Z_init(C_Z, 3, 3);
  mat_Z_LLL(C_Z, A_Z);
  mat_Z_clear(A_Z);

  mat_int64_init(A, 3, 3);
  mat_int64_set_identity(A);
  A->coeff[1][1] = 1;
  A->coeff[1][2] = 0;
  A->coeff[1][3] = 0;
  A->coeff[2][1] = 0;
  A->coeff[2][2] = 1;
  A->coeff[2][3] = 0;
  A->coeff[3][1] = 0;
  A->coeff[3][2] = 0;
  A->coeff[3][3] = 1;

  mat_int64_init(C, 3, 3);
  mat_int64_set_zero(C);
  mat_int64_LLL(C, A);

  mat_Z_to_mat_int64(A, C_Z);

  ASSERT_ALWAYS(mat_int64_equal(A, C));

  mat_int64_clear(A);
  mat_int64_clear(C);
  mat_Z_clear(C_Z);
}
