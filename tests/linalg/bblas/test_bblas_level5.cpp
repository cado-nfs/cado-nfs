#include "cado.h" // IWYU pragma: keep
#include <stdint.h>               // for uint64_t, UINT64_C
#include "bblas_mat64.hpp"
#include "bblas_level5.hpp"
#include "gmp_aux.h"              // for memfill_random
#include "test_bblas_base.hpp"
#include "test_bblas_level5.hpp"
#include "time_bblas_common.hpp"
#ifdef HAVE_M4RIE
/* To test against m4rie routines, include a checkout of
 * https://bitbucket.org/malb/m4rie.git under linalg/m4rie, and run
 * "autoreconf -i" there ; cado-nfs cmake logic then detects it and
 * enables the corresponding code here (and in a few other places in
 * test_bblas). */
#include "m4rie/m4rie.h"
#endif

    /* We also have a completely different set of routines that deal with
     * polynomials of matrices */

test_bblas_base::tags_t test_bblas_level5::polmul_tags { "polmul", "poly", "l5" };/*{{{*/
void test_bblas_level5::polmul() {
    size_t n = 64;
    mat64 * A = mat64::alloc(n);
    mat64 * B = mat64::alloc(n);
    mat64 * C = mat64::alloc(2 * n);
    memfill_random(A, n * sizeof(mat64), rstate);
    memfill_random(B, n * sizeof(mat64), rstate);
    printf("-- polynomials (N=%zu) --\n", n);
    TIME1(5, m64pol_mul, C,A,B,n,n);
    TIME1(5, m64pol_mul_kara, C,A,B,n,n);
    mat64::free(A, n);
    mat64::free(B, n);
    mat64::free(C, n);
}/*}}}*/

test_bblas_base::tags_t test_bblas_level5::polblockmul_tags { "polblockmul", "poly", "l5" };/*{{{*/
void test_bblas_level5::polblockmul() {
    size_t n = 64;
    unsigned int K = 2;
    mat64 * A = mat64::alloc(K * K * n);
    mat64 * B = mat64::alloc(K * K * n);
    mat64 * C = mat64::alloc(K * K * 2 * n);
    memfill_random(A, K * K * n * sizeof(mat64), rstate);
    memfill_random(B, K * K * n * sizeof(mat64), rstate);
    printf("-- polynomials, larger matrices (K=%u, N=%zu) --\n", K, n);
    TIME1(5, m64polblock_mul, C,A,B,n,n,2);
    TIME1(5, m64polblock_mul_kara, C,A,B,n,n,K);
    mat64::free(A, K * K * n);
    mat64::free(B, K * K * n);
    mat64::free(C, K * K * 2 * n);
}/*}}}*/

test_bblas_base::tags_t test_bblas_level5::matpolmul_tags = { "matpolmul", "poly", "l5" };/*{{{*/
void test_bblas_level5::matpolmul() {
    size_t n = 128;
    mat64 * A = mat64::alloc(n);
    mat64 * B = mat64::alloc(n);
    mat64 * C = mat64::alloc(n);
    uint64_t * Al = (uint64_t *) A;
    uint64_t * Bl = (uint64_t *) B;
    uint64_t * Cl = (uint64_t *) C;
    memfill_random(A, n * sizeof(mat64), rstate);
    memfill_random(B, n * sizeof(mat64), rstate);
    printf("-- 64x64 matrices over GF(2^64) --\n");
    TIME1(5, m64pol_mul_gf2_64_bitslice, C,A,B);
    TIME1(5, m64pol_mul_gf2_64_nobitslice, Cl,Al,Bl);
    printf("-- 64x64 matrices over GF(2^128) --\n");
    TIME1(5, m64pol_mul_gf2_128_bitslice, C,A,B);
    TIME1(5, m64pol_mul_gf2_128_nobitslice, Cl,Al,Bl);
    mat64::free(A, n);
    mat64::free(B, n);
    mat64::free(C, n);
    /* On Core i5 (magret), it's almost a tie between the two
     * options... */
#if 0
    -- 64x64 matrices over GF(2^64) --
        m64pol_mul_gf2_64_bitslice       6351 times in 0.7889 ms each
        m64pol_mul_gf2_64_nobitslice    4773 times in 1.0497 ms each
        -- 64x64 matrices over GF(2^128) --
        m64pol_mul_gf2_128_bitslice      2067 times in 2.4238 ms each
        m64pol_mul_gf2_128_nobitslice   1521 times in 3.2939 ms each
#endif
        /* Without pclmul, of course the situation is more clear (truffe,
         * Core2 Duo U9400 */
#if 0
        -- 64x64 matrices over GF(2^64) --
        m64pol_mul_gf2_64_bitslice       2695 times in 1.8590 ms each
        m64pol_mul_gf2_64_nobitslice    435 times in 11.5172 ms each
        -- 64x64 matrices over GF(2^128) --
        m64pol_mul_gf2_128_bitslice      871 times in 5.7520 ms each
        m64pol_mul_gf2_128_nobitslice   158 times in 31.8354 ms each

#endif
#ifdef HAVE_M4RIE
        if (0) {
            printf("-- 64x64 matrices over GF(2^64) using M4RIE --\n");
            /* Now try to see if m4rie can improve these timings */
            /* Unfortunately as of version 20111203, m4rie supports only
             * GF(2^n) up until n==160. Which cleary won't do, for our
             * objectives.
             * Even the interface for the defining polynomial only sets
             * it via a uint64_t, degree included. So it seems that we
             * cannot do degree 64, period...
             * Bottom line, the only thing that the code below does is a
             * nifty segfault...
             */
            gf2e * ff = gf2e_init(/* 2^64 + ... */ UINT64_C(0x0000000247F43CB7));
            mzed_t *Az = mzed_init(ff, 64, 64);
            mzed_t *Bz = mzed_init(ff, 64, 64);
            mzed_t *Cz = mzed_init(ff, 64, 64);
            mzed_randomize(Az);
            mzed_randomize(Bz);
            TIME1(5, mzed_mul, Cz, Az, Bz);
            mzed_free(Az);
            mzed_free(Bz);
            mzed_free(Cz);
            gf2e_free(ff);
        }
#endif /* HAVE_M4RIE */

}/*}}}*/

test_bblas_base::tags_t test_bblas_level5::matpolscale_tags = { "matpolscale", "poly", "l5" };/*{{{*/
void test_bblas_level5::matpolscale() {
    /* Now multiplication by a scalar. We'll do both GF(2^64) and
     * GF(2^128), so let's allocate room for both */
    size_t n = 128;
    /* random values with average hamming weight. */
    uint64_t scalar[2] = { UINT64_C(0x8d5511cbd7f0d885), UINT64_C(0x2073a477a8b5dd8a) };
    mat64 * A = mat64::alloc(n);
    mat64 * B = mat64::alloc(n);
    uint64_t * Al = (uint64_t *) A;
    uint64_t * Bl = (uint64_t *) B;
    memfill_random(A, n * sizeof(mat64), rstate);
    memfill_random(B, n * sizeof(mat64), rstate);
    printf("-- 64x64 matrix over GF(2^64), multiplication by scalar --\n");
    TIME1(5, m64pol_scalmul_gf2_64_bitslice, B,A,scalar);
    TIME1(5, m64pol_scalmul_gf2_64_bitslice2, B,A,scalar);
    TIME1(5, m64pol_scalmul_gf2_64_nobitslice, Bl,Al,scalar);
    printf("-- 64x64 matrix over GF(2^128), multiplication by scalar --\n");
    TIME1(5, m64pol_scalmul_gf2_128_bitslice, B,A,scalar);
    TIME1(5, m64pol_scalmul_gf2_128_nobitslice, Bl,Al,scalar);
    mat64::free(A, n);
    mat64::free(B, n);
    /* The bitsliced version sucks. Really.
     * TODO: See if we can do something. Abandon L1 cache focus, and
     * be content with L2 ? */
} /*}}}*/
