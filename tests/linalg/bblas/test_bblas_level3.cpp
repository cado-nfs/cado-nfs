#include "cado.h"
#include <cstdint>               // for uint64_t, UINT64_C, UINT8_C, uint32_t
#include <cstring>
#include <gmp.h>                  // for gmp_urandomm_ui
#include "bblas_bitmat.hpp"
#include "bblas_mat64.hpp"
#include "bblas_mat8.hpp"
#include "bblas_level3a.hpp"      // for mat64_fill_random, mat64_transpose
#include "bblas_level3a1.hpp"     // for binary_matpoly_to_polmat_nested_tra...
#include "bblas_level3b.hpp"      // for mul_6464_6464, mul_6464_6464_sse
#include "bblas_level3c.hpp"      // for mul_N64_6464_vec, mul_N64_6464_avx2
#include "bblas_level3d.hpp"      // for trsm64_general, trsm64
#include "gmp_aux.h"              // for memfill_random
#include "macros.h"               // for ASSERT_ALWAYS, iceildiv
#include "test_bblas_base.hpp"
#include "test_bblas_level3.hpp"
#include "time_bblas_common.hpp"

/* level-3 combine matrices together. Most of our focus is on
 * fixed-size matrices, while blas do not fix the matrix size a
 * priori.
 */
test_bblas_base::tags_t test_bblas_level3::level3a_tags { "basic", "l3a", "l3"};/*{{{ */ 
void test_bblas_level3::level3a()
{
    printf(" -- straightforward operations --\n");
#ifdef  HAVE_M4RI
    {
        mat64 & XR = * (mat64 *) xr;
        mat64 & A = * (mat64 *) a;
        XR = A;
    }
    mzd_copy(R64, A64);
    mzd_check_mem(R64, xr, 64);
    TIME1(1, mzd_copy, R64, A64);
#endif				/* HAVE_M4RI */

    /* add */
    {
        mat64 & R = * (mat64 *) r;
        mat64 & XR = * (mat64 *) xr;
        mat64 & A = * (mat64 *) a;
        mat64_add(R, A, w);
        XR = R;

        mat64_add_C(R, A, w);
        ASSERT_ALWAYS(XR == R);
        TIME1(1, mat64_add_C, R, A, w);
    }

#ifdef  HAVE_M4RI
    mzd_add(R64, A64, W);
    mzd_check_mem(R64, xr, 64);
    TIME1(1, mzd_add, R64, A64, W);
#endif				/* HAVE_M4RI */
}				/*}}} */

test_bblas_base::tags_t test_bblas_level3::transpose_tags { "transpose", "matpoly_polmat", "poly", "l3a", "l3" };/*{{{*/
void test_bblas_level3::transpose() {
    printf(" -- transposition --\n");
    /* There is no BLAS analogue, but it seems that it should be
     * categorized together with the other 3a things.
     */
    {
        mat64 & R = * (mat64 *) r;
        mat64 & XR = * (mat64 *) xr;
        mat64 & A = * (mat64 *) a;

        mat64_transpose(R, A);
        mat64_transpose(XR, R);
        ASSERT_ALWAYS(XR == A);
        TIME1(1, mat64_transpose_simple_and_stupid, R, A);
        TIME1(1, mat64_transpose_recursive, R, A);
        TIME1(1, mat64_transpose, R, A);
    }

#ifdef  HAVE_M4RI
    mzd_transpose(R64, A64);
    mzd_transpose(R64, R64);
    mzd_check_mem(R64, xr, 64);
    TIME1(1, mzd_transpose, R64, A64);
#endif  /* HAVE_M4RI */
}/*}}}*/

test_bblas_base::tags_t test_bblas_level3::matpoly_polmat_tags { "matpoly_polmat", "l3a", "l3" };/*{{{*/
void test_bblas_level3::matpoly_polmat() {
    size_t n = 512;
    unsigned int K = 16;
    unsigned int L = 8;
    if (test_accel) {
        n = 128;
        K = 3;
        L = 2;
    }
    mat64 * A = mat64::alloc(n * K * L);
    mat64 * B = mat64::alloc(n * K * L);
    size_t datasize = K * 64 * L * 64 * iceildiv(n, ULONG_BITS) * (64 / ULONG_BITS);
    unsigned long * data = new unsigned long[datasize];
    unsigned long * data_t = new unsigned long[datasize];
    memfill_random(A, n * K * L * sizeof(mat64), rstate);
    printf(" -- conversion of %u*%u*%zu bit matrices to/from %u*%u %zu-bit polynomials\n", 
            K, L, n, K * 64, L * 64, n);
    TIME1(5, binary_polmat_to_matpoly_simple_and_stupid, data, A, K * 64, L * 64, n);
    TIME1(5, binary_polmat_to_matpoly_nested_transpositions, data_t, A, K * 64, L * 64, n);
    ASSERT_ALWAYS(memcmp(data, data_t, datasize * sizeof(unsigned long)) == 0);

    TIME1(5, binary_matpoly_to_polmat_simple_and_stupid, B, data, K * 64, L * 64, n);
    for(unsigned int i = 0 ; i < K ; i++ ) {
        for(unsigned int j = 0 ; j < L ; j++ ) {
            for(size_t k = 0 ; k < n ; k++) {
                ASSERT_ALWAYS(A[(k*K+i)*L+j] == B[(k*K+i)*L+j]);
            }
        }
    }
    TIME1(5, binary_matpoly_to_polmat_nested_transpositions, B, data, K * 64, L * 64, n);
    for(unsigned int i = 0 ; i < K ; i++ ) {
        for(unsigned int j = 0 ; j < L ; j++ ) {
            for(size_t k = 0 ; k < n ; k++) {
                ASSERT_ALWAYS(A[(k*K+i)*L+j] == B[(k*K+i)*L+j]);
            }
        }
    }
    delete[] data;
    delete[] data_t;
    mat64::free(A, n*K*L);
    mat64::free(B, n*K*L);
}
/*}}}*/

test_bblas_base::tags_t test_bblas_level3::matmul_tags { "matmul", "l3b", "l3" };/*{{{*/
void test_bblas_level3::matmul() {
    printf(" -- matrix products --\n");
    /* BLAS level 3 analogue: cblas_dgemm */

    mat64 & R = * (mat64 *) r;
    mat64 & XR = * (mat64 *) xr;
    mat64 & A = * (mat64 *) a;

    /* multiplicate of two 64x64 matrices */
    mul_6464_6464(R, A, w);
    XR = R;

    TIME1(1, mul_6464_6464, R, A, w);

#if defined(HAVE_SSE2) && ULONG_BITS == 64
    mul_6464_6464_sse(R, A, w);
    ASSERT_ALWAYS(XR == R);
    TIME1(1, mul_6464_6464_sse, R, A, w);
#endif

    mul_6464_6464_v2(R, A, w);
    ASSERT_ALWAYS(XR == R);
    TIME1(1, mul_6464_6464_v2, R, A, w);

    /* Functions which can do any n can also do n=64 */
    test_bblas_level3(64).level3c_list();

    mat64 L = A;
    L.make_lowertriangular();
    mul_6464_6464(R, L, w);
    XR = R;
    mul_6464lt_6464(R, L, w);
    ASSERT_ALWAYS(XR == R);
    TIME1(1, mul_6464lt_6464, R, L, w);

}/*}}}*/

test_bblas_base::tags_t test_bblas_level3::rank_n_update_tags { "rank_n_update", "l3c", "l3", };/*{{{*/
void test_bblas_level3::rank_n_update() {
    unsigned int n = nmax;
    printf(" -- rank-n updates --\n");

    mat64 & R = * (mat64 *) r;

    TIME1N(1, mul_TN64_N64_addmul, R, a, b, n);
    TIME1N(5, mul_TN32_N64_C, r, (uint32_t*)a, b, n);
    TIME1N(5, mul_TN64_N64_C, R, a, b, n);
}/* }}} */

test_bblas_base::tags_t test_bblas_level3::level3c_tags { "l3c", "l3" };/*{{{*/
void test_bblas_level3::level3c_list()
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

#if defined(HAVE_AVX2)
    mul_N64_6464_avx2(r, a, w, n);
    ASSERT_ALWAYS(memcmp(xr, r, n * sizeof(uint64_t)) == 0);
    TIME1N(2, mul_N64_6464_avx2, r, a, w, n);
#endif

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
}
void test_bblas_level3::level3c()
{
    printf(" -- matrix column times matrix --\n");
    level3c_list();
}/*}}}*/

test_bblas_base::tags_t test_bblas_level3::trsm_tags { "trsm", "l3d", "l3" };/*{{{*/
void test_bblas_level3::trsm() {
    /* BLAS level 3 analogue: cblas_dtrsm */
    printf(" -- solve unit lower triangular systems\n");
    mat64 L, U0, U1;

    /* test consistency first. Extract L from a matrix full of garbage.
     * We'll feed the garbage matrix to trsm64_general, to verify that it
     * copes with intervals well */

    for(unsigned int k = 0 ; k < 1000 ; k++) {
        mat64 Lr;
        unsigned int n0 = gmp_urandomm_ui(rstate, 64);
        unsigned int n1 = gmp_urandomm_ui(rstate, 65 - n0) + n0;
        mat64_fill_random(Lr, rstate);
        L = 1;
        for(unsigned int i = n0 + 1 ; i < n1 ; i++) {
            /* import bits [n0..i-1] from Lr */
            uint64_t m = (-(UINT64_C(1) << n0)) & ((UINT64_C(1) << i)-1);
            L[i] ^= Lr[i] & m;
        }

        mat64_fill_random(U0, rstate);
        U1 = U0;
        trsm64_general(Lr, U0, n0, n1);     /* U0 is now L^-1 * U1 */
        mul_6464_6464(U0, L, U0);
        ASSERT_ALWAYS(U0 == U1);
    }


    mat64_fill_random(L, rstate);
    mat64_fill_random(U0, rstate);
    TIME1(1, trsm64, L, U0);
}
/*}}}*/

test_bblas_base::tags_t test_bblas_level3::m8_tags { "m8", "l3"};/* */ 
void test_bblas_level3::m8()
{
    printf(" -- straightforward operations --\n");

    /* add */
    mat8 & R = * (mat8 *) r;
    mat8 & XR = * (mat8 *) xr;
    mat8 & A = * (mat8 *) a;
    mat8 w8;

    mat8::fill_random(w8, rstate);
    TIME1(1, mat8::add, R, A, w8);

    mat8::transpose(R, A);
    mat8::transpose(XR, R);
    ASSERT_ALWAYS(XR == A);
    TIME1(1, mat8::transpose, R, A);

    TIME1(1, mat8::mul, R, A, w8);

    mat8 L, U0, U1;

    /* test consistency first. Extract L from a matrix full of garbage.
     * We'll feed the garbage matrix to trsm8_general, to verify that it
     * copes with intervals well */

    for(unsigned int k = 0 ; k < 1000 ; k++) {
        mat8 Lr;
        unsigned int n0 = gmp_urandomm_ui(rstate, mat8::width);
        unsigned int n1 = gmp_urandomm_ui(rstate, mat8::width + 1 - n0) + n0;
        mat8::fill_random(Lr, rstate);
        L = 1;
        for(unsigned int i = n0 + 1 ; i < n1 ; i++) {
            /* import bits [n0..i-1] from Lr */
            uint8_t m = (-(UINT8_C(1) << n0)) & ((UINT8_C(1) << i)-1);
            L[i] ^= Lr[i] & m;
        }

        mat8::fill_random(U0, rstate);
        U1 = U0;
        mat8::trsm(Lr, U0, n0, n1);     /* U0 is now L^-1 * U1 */
        mat8::mul(U0, L, U0);
        ASSERT_ALWAYS(U0 == U1);
    }
    mat8::fill_random(L, rstate);
    mat8::fill_random(U0, rstate);
    /* since mat8::trsm is an overloaded function, the templates in TIME1
     * can't infer the type. It's a bit of a pity, but working around it
     * is easy enough
     */
    auto mat8trsm = [&](mat8 const &L, mat8 & U) { mat8::trsm(L, U); };
    TIME1(1, mat8trsm, L, U0);
}
/*}}}*/

