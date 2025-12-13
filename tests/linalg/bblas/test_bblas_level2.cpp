#include "cado.h" // IWYU pragma: keep

#include <cstdint>
#include <cstdio>
#include <cstring>

#include "test_bblas_base.hpp"
#include "test_bblas_level2.hpp"
#include "test_bblas_level3.hpp"
#include "time_bblas_common.hpp"
#include "bblas_mat64.hpp"
#include "bblas_level2a.hpp"              // for addmul_To64_o64_lsb_sse_v1 ...
#include "bblas_level2b.hpp"              // for mul_o64_6464, mul_o64_6464_C_lsb
#include "macros.h"               // for ASSERT_ALWAYS

/* level-2 blas combines vectors and matrices (i.e. among three
 * operands, two of which are inputs and one is output, we have only
 * one two-dimensional operand. It need not necessarily be an input
 * operand).
 */
test_bblas_base::tags_t test_bblas_level2::level2a_tags { "rank_1_update", "l2a", "l2", };/*{{{*/
void test_bblas_level2::level2a() {
    /* BLAS level 2 analogue: cblas_dger (rank-1 update) */
    /* TODO: check correctness */
    printf(" -- rank-1 updates --\n");

    mat64  const& R = * (mat64 *) r;

    TIME1(1, addmul_To64_o64, R, *a, *b);

    TIME1(1, addmul_To64_o64_lsb, R, *a, *b);
    TIME1(1, addmul_To64_o64_msb, R, *a, *b);
    TIME1(1, addmul_To64_o64_lsb_packof2, R, *a, *b);
#if defined(HAVE_SSE2) && ULONG_BITS == 64
    TIME1(1, addmul_To64_o64_lsb_sse_v1, R, *a, *b);
#endif
} /*}}}*/

test_bblas_base::tags_t test_bblas_level2::level2_tags { "vecmul", "l2b", "l2" };
/*{{{*/
void test_bblas_level2::level2() {
    unsigned int const n = 1;

    printf(" -- vector times (transpose of) matrix --\n");
    /* BLAS level 2 analogue: cblas_dgemv */
    TIME1(1, mul_o64_6464, r, *a, wt);
    TIME1(1, mul_o64_T6464, r, *a, wt);

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
    test_bblas_level3(1).level3c_list();
}/*}}}*/
