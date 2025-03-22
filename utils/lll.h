#ifndef CADO_UTILS_LLL_H
#define CADO_UTILS_LLL_H

#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    /* TODO: Refactor. Outstanding uses should use mpz_mat_LLL instead.
     * We have some in size_optimization.c
     * */
    mpz_t ** coeff;
    int NumRows, NumCols;
} mat_Z;

void LLL_init(mat_Z *, int, int);
void LLL_clear(mat_Z *);
int LLL(mpz_t det, mat_Z B, mat_Z * U, mpz_srcptr a, mpz_srcptr b);

#ifdef __cplusplus
}
#endif

#endif /* CADO_UTILS_LLL_H */
