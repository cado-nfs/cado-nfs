#include "cado.h"
#include "test_bblas_level3.hpp"
#include "time_bblas_common.hpp"
#include <cstring>

/* level-3 combine matrices together. Most of our focus is on
 * fixed-size matrices, while blas do not fix the matrix size a
 * priori.
 */
test_bblas_base::tags_t test_bblas_level3::level3a_tags { "basic", "l3a", "l3"};/*{{{ */ 
void test_bblas_level3::level3a()
{
    /* copy */
    mat64_copy(r, a);
    mat64_copy(xr, r);

    mat64_copy(r, a);
    ASSERT_ALWAYS(mat64_eq(xr, r));
    TIME1(1, mat64_copy, r, a);

#ifdef  HAVE_M4RI
    mzd_copy(R, A);
    mzd_check_mem(R, xr, 64);
    TIME1(1, mzd_copy, R, A);
#endif				/* HAVE_M4RI */

    /* add */
    mat64_add(r, a, w);
    mat64_copy(xr, r);

    mat64_add_C(r, a, w);
    ASSERT_ALWAYS(mat64_eq(xr, r));
    TIME1(1, mat64_add_C, r, a, w);

#ifdef  HAVE_M4RI
    mzd_add(R, A, W);
    mzd_check_mem(R, xr, 64);
    TIME1(1, mzd_add, R, A, W);
#endif				/* HAVE_M4RI */
}				/*}}} */

test_bblas_base::tags_t test_bblas_level3::transpose_tags { "transpose", "matpoly_polmat", "poly", "l3a", "l3" };/*{{{*/
void test_bblas_level3::transpose() {
    /* There is no BLAS analogue, but it seems that it should be
     * categorized together with the other 3a things.
     */
    mat64_transpose(r, a);
    mat64_transpose(xr, r);
    ASSERT_ALWAYS (mat64_eq(xr, a));
    TIME1(1, mat64_transpose_simple_and_stupid, r, a);
    TIME1(1, mat64_transpose_recursive, r, a);
    TIME1(1, mat64_transpose, r, a);

#ifdef  HAVE_M4RI
    mzd_transpose(R, A);
    mzd_check_mem(R, xr, 64);
    TIME1(1, mzd_transpose, R, A);
#endif  /* HAVE_M4RI */
}/*}}}*/

test_bblas_base::tags_t test_bblas_level3::matpoly_polmat_tags { "matpoly_polmat", "l3a", "l3" };/*{{{*/
void test_bblas_level3::matpoly_polmat() {
    size_t n = 512;
    unsigned int K = 16;
    unsigned int L = 8;
    mat64 * A = (mat64 *) malloc(n * K * L * sizeof(mat64));
    mat64 * B = (mat64 *) malloc(n * K * L * sizeof(mat64));
    size_t datasize = K * 64 * L * 64 * iceildiv(n, 64);
    uint64_t * data = (uint64_t *) malloc(datasize * sizeof(uint64_t));
    uint64_t * data_t = (uint64_t *) malloc(datasize * sizeof(uint64_t));
    memfill_random(A, n * K * L * sizeof(mat64), rstate);
    printf("-- expand %u*%u*%zu bit matrices into %u*%u %zu-bit polynomials\n", 
            K, L, n, K * 64, L * 64, n);
    TIME1(5, binary_polmat_to_matpoly_simple_and_stupid, data, A, K * 64, L * 64, n);
    TIME1(5, binary_polmat_to_matpoly_nested_transpositions, data_t, A, K * 64, L * 64, n);
    ASSERT_ALWAYS(memcmp(data, data_t, datasize * sizeof(uint64_t)) == 0);
    printf("-- retrieve %u*%u*%zu bit matrices from %u*%u %zu-bit polynomials\n", 
            K, L, n, K * 64, L * 64, n);
    TIME1(5, binary_matpoly_to_polmat_simple_and_stupid, B, data, K * 64, L * 64, n);
    for(unsigned int i = 0 ; i < K ; i++ ) {
        for(unsigned int j = 0 ; j < L ; j++ ) {
            for(size_t k = 0 ; k < n ; k++) {
                ASSERT_ALWAYS(mat64_eq(A[(k*K+i)*L+j], B[(k*K+i)*L+j]));
            }
        }
    }
    TIME1(5, binary_matpoly_to_polmat_nested_transpositions, B, data, K * 64, L * 64, n);
    for(unsigned int i = 0 ; i < K ; i++ ) {
        for(unsigned int j = 0 ; j < L ; j++ ) {
            for(size_t k = 0 ; k < n ; k++) {
                ASSERT_ALWAYS(mat64_eq(A[(k*K+i)*L+j], B[(k*K+i)*L+j]));
            }
        }
    }
    free(data);
    free(data_t);
    free(A);
    free(B);
}
/*}}}*/

test_bblas_base::tags_t test_bblas_level3::matmul_tags { "matmul", "l3b", "l3" };/*{{{*/
void test_bblas_level3::matmul() {
    unsigned int n = 64;

    /* BLAS level 3 analogue: cblas_dgemm */

    /* multiplicate of two 64x64 matrices */
    mul_6464_6464(r, a, w);
    memcpy(xr, r, n * sizeof(uint64_t));

#if defined(HAVE_SSE2) && ULONG_BITS == 64
    mul_6464_6464_sse(r, a, w);
    ASSERT_ALWAYS(memcmp(xr, r, n * sizeof(uint64_t)) == 0);
    TIME1(1, mul_6464_6464_sse, r, a, w);
#endif

    mul_6464_6464_v2(r, a, w);
    ASSERT_ALWAYS(memcmp(xr, r, n * sizeof(uint64_t)) == 0);
    TIME1(1, mul_6464_6464_v2, r, a, w);

    /* Functions which can do any n can also do n=64 */
    test_bblas_level3(64).level3c();
}/*}}}*/

test_bblas_base::tags_t test_bblas_level3::rank_n_update_tags { "rank_n_update", "l3c", "l3", };/*{{{*/
void test_bblas_level3::rank_n_update() {
    unsigned int n = nmax;
    TIME1N(1, mul_TN64_N64_addmul, r, a, b, n);
    TIME1N(5, mul_TN32_N64_C, r, (uint32_t*)a, b, n);
    TIME1N(5, mul_TN64_N64_C, r, a, b, n);
}/* }}} */

test_bblas_base::tags_t test_bblas_level3::level3c_tags { "l3c", "l3" };/*{{{*/
void test_bblas_level3::level3c()
{
    unsigned int n = nmax;
    /* multiplication of a vector by a matrix */
    mul_N64_6464_vec(r, a, w, n);
    memcpy(xr, r, n * sizeof(uint64_t));

    mul_N64_6464_vec(r, a, w, n);
    ASSERT_ALWAYS(memcmp(xr, r, n * sizeof(uint64_t)) == 0);
    TIME1N(2, mul_N64_6464_vec, r, a, w, n);

    mul_N64_6464_transB(r, a, w, n);
    ASSERT_ALWAYS(memcmp(xr, r, n * sizeof(uint64_t)) == 0);
    TIME1N(2, mul_N64_6464_transB, r, a, w, n);

#if defined(HAVE_SSE2) && ULONG_BITS == 64
    mul_N64_6464_sse(r, a, w, n);
    ASSERT_ALWAYS(memcmp(xr, r, n * sizeof(uint64_t)) == 0);
    TIME1N(2, mul_N64_6464_sse, r, a, w, n);
#endif

    mul_N64_6464_lookup4(r, a, w, n);
    ASSERT_ALWAYS(memcmp(xr, r, n * sizeof(uint64_t)) == 0);
    TIME1N(2, mul_N64_6464_lookup4, r, a, w, n);

    mul_N64_6464_lookup8(r, a, w, n);
    ASSERT_ALWAYS(memcmp(xr, r, n * sizeof(uint64_t)) == 0);
    TIME1N(2, mul_N64_6464_lookup8, r, a, w, n);

    mul_N64_T6464_vec(r, a, wt, n);
    ASSERT_ALWAYS(memcmp(xr, r, n * sizeof(uint64_t)) == 0);
    TIME1N(2, mul_N64_T6464_vec, r, a, w, n);

    mul_N64_T6464_transB(r, a, wt, n);
    ASSERT_ALWAYS(memcmp(xr, r, n * sizeof(uint64_t)) == 0);
    TIME1N(2, mul_N64_T6464_transB, r, a, w, n);

#ifdef HAVE_M4RI
    mzd_mul_naive(R, A, W);
    mzd_check_mem(R, xr, n);
    TIME1(1, mzd_mul_naive, R, A, W);

    _mzd_mul_naive(R, A, WT, 1);
    mzd_check_mem(R, xr, n);
    TIME1(1, _mzd_mul_naive, R, A, WT, 1);

    mzd_mul_m4rm(R, A, W, 0);
    mzd_check_mem(R, xr, n);
    TIME1(1, mzd_mul_m4rm, R, A, W, 0);

    _mzd_mul_m4rm(R, A, W, 0, 1);
    mzd_check_mem(R, xr, n);
    TIME1(1, _mzd_mul_m4rm, R, A, W, 0, 1);
#endif /* HAVE_M4RI */
}/*}}}*/

test_bblas_base::tags_t test_bblas_level3::trsm_tags { "trsm", "l3d", "l3" };/*{{{*/
void test_bblas_level3::trsm() {
    /* BLAS level 3 analogue: cblas_dtrsm */
    printf("-- solve unit lower triangular systems\n");
    mat64 L, U0, Li, U1;
    memfill_random(L, (64) * sizeof(uint64_t), rstate);
    for(unsigned int i = 0 ; i < 64 ; i++) {
        L[i] &= (UINT64_C(1) << i) - 1;
        L[i] |=  UINT64_C(1) << i;
    }
    full_echelon_6464_imm(U0, Li, L);
    memfill_random(U0, (64) * sizeof(uint64_t), rstate);
    mat64_copy(U1, U0);
    TIME1(1, trsm64, L, U0);
    mat64_copy(U0, U1);
    trsm64(L, U0);
    trsm64(L, U1);
    ASSERT_ALWAYS(mat64_eq(U0, U1));
}
/*}}}*/

