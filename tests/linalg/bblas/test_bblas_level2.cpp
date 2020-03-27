#include "cado.h"

#include "test_bblas_level2.hpp"
#include "test_bblas_level3.hpp"
#include "time_bblas_common.hpp"
#include <cstring>

/* level-2 blas combines vectors and matrices (i.e. among three
 * operands, two of which are inputs and one is output, we have only
 * one two-dimensional operand. It need not necessarily be an input
 * operand).
 */
test_bblas_base::tags_t test_bblas_level2::level2a_tags { "rank_1_update", "l2a", "l2", };/*{{{*/
void test_bblas_level2::level2a() {
    /* BLAS level 2 analogue: cblas_dger (rank-1 update) */
    /* TODO: check correctness */
    uint64_t * r = (uint64_t *) malloc(64 * sizeof(uint64_t));
    uint64_t a = uint64_random(rstate);
    uint64_t w = uint64_random(rstate);
    TIME1(1, addmul_To64_o64_lsb, r, a, w);
    TIME1(1, addmul_To64_o64_msb, r, a, w);
    TIME1(1, addmul_To64_o64_lsb_packof2, r, a, w);
#if defined(HAVE_SSE2) && ULONG_BITS == 64
    TIME1(1, addmul_To64_o64_lsb_sse_v1, r, a, w);
#endif
    TIME1(1, addmul_To64_o64, r, a, w);
    free(r);
} /*}}}*/

test_bblas_base::tags_t test_bblas_level2::level2_tags { "vecmul", "l2b", "l2" };
/*{{{*/
void test_bblas_level2::level2() {
    unsigned int n = 1;

    /* BLAS level 2 analogue: cblas_dgemv */

    /* reference */
    mul_o64_6464(r, *a, w);
    memcpy(xr, r, n * sizeof(uint64_t));

    /* vectors and matrices */
    mul_o64_6464_C_lsb(r, *a, w);
    ASSERT_ALWAYS(memcmp(xr, r, n * sizeof(uint64_t)) == 0);
    TIME1(1, mul_o64_6464_C_lsb, r, *a, w);

    mul_o64_6464_C_msb(r, *a, w);
    ASSERT_ALWAYS(memcmp(xr, r, n * sizeof(uint64_t)) == 0);
    TIME1(1, mul_o64_6464_C_msb, r, *a, w);


    /* multiply vector by transpose of matrix */
    mul_o64_T6464_C_parity(r, *a, wt);
    ASSERT_ALWAYS(memcmp(xr, r, n * sizeof(uint64_t)) == 0);
    TIME1(1, mul_o64_T6464_C_parity, r, *a, wt);

    mul_o64_T6464_C_parity3(r, *a, wt);
    ASSERT_ALWAYS(memcmp(xr, r, n * sizeof(uint64_t)) == 0);
    TIME1(1, mul_o64_T6464_C_parity3, r, *a, wt);

    /* Functions which can do any n can also do n=1 */
    test_bblas_level3(1).level3c();
}/*}}}*/
