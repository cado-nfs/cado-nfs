#ifndef MATRIX_H
#define MATRIX_H 

#include <gmp.h>
#include <stdint.h>

/* TODO: compare with utils/mpz_mat.cpp, see if there's stuff that should
 * be merged */

#ifdef __cplusplus
extern "C" {
#endif


typedef struct {
  mpz_t **coeff;
  unsigned int NumRows;
  unsigned int NumCols;
} s_mat_Z_t;

typedef s_mat_Z_t mat_Z_t[1];
typedef s_mat_Z_t * mat_Z_ptr;
typedef const s_mat_Z_t * mat_Z_srcptr;

typedef struct {
  int64_t **coeff;
  unsigned int NumRows;
  unsigned int NumCols;
} s_mat_int64_t;

typedef s_mat_int64_t mat_int64_t[1];
typedef s_mat_int64_t * mat_int64_ptr;
typedef const s_mat_int64_t * mat_int64_srcptr;

typedef struct {
  double **coeff;
  unsigned int NumRows;
  unsigned int NumCols;
} s_mat_double_t;

typedef s_mat_double_t mat_double_t[1];
typedef s_mat_double_t * mat_double_ptr;
typedef const s_mat_double_t * mat_double_srcptr;

#ifdef __cplusplus
}
#endif
#endif /* MATRIX_H */
