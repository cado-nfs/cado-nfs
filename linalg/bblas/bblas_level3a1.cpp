#include "cado.h" // IWYU pragma: keep
#include <cstdint>                         // for uint64_t
#include <cstring>
#include <algorithm>
#include "bblas_mat64.hpp"  // for mat64
#include "bblas_level3a.hpp"  // for mat64_transpose
#include "macros.h"           // for ASSERT_ALWAYS, iceildiv
#include "bblas_level3a1.hpp"

/**********************************************************************/
/* level 3a (extension, matpoly_polmat): conversions.
 *      binary_polmat_to_matpoly    (several variants, one obvious winner)
 *      binary_matpoly_to_polmat    (several variants, one obvious winner)
 *      TODO: better naming.
 */

/* implemented here:
 *    - binary_matpoly_to_polmat ; takes an m*n matrix of polynomials of
 *      length len (all multiples of 64), and returns a length len
 *      polynomial of (m/64)*(n/64) block matrices, each block being a
 *      64*64 matrix
 *    - binary_polmat_to_matpoly ; converse
 */

/* implements binary_matpoly_to_polmat */
void binary_matpoly_to_polmat_simple_and_stupid(mat64 * dst, unsigned long const * src, unsigned int m, unsigned int n, unsigned int len)/*{{{*/
{
    /* dst must have room for len mat64's. */
    /* src is assumed row-major */
    ASSERT_ALWAYS((n%64)==0);
    ASSERT_ALWAYS((m%64)==0);
    size_t stride = iceildiv(len, ULONG_BITS);
    memset((void *) dst, 0, (m/64) * (n/64) * len * sizeof(mat64));
    mat64 * q0 = dst;
    unsigned long const * p0 = src;
    unsigned long mk = 1;
    for(unsigned int k = 0 ; k < len ; k++) {
        uint64_t mj = 1;
        unsigned long const * p = p0;
        for(unsigned int i = 0 ; i < m ; i++) {
            mat64 * q = q0 + (i / 64) * (n / 64);
            for(unsigned int j = 0 ; j < n ; j++) {
                (*q)[i%64] ^= mj & -((*p&mk) != 0);
                mj <<= 1;
                q  += (mj == 0);
                mj += mj == 0;
                p += stride;
            }
        }
        mk <<= 1;
        p0 += (mk == 0);
        q0 += (m/64) * (n/64);
        mk += (mk == 0);
    }
}/*}}}*/

/* implements binary_polmat_to_matpoly */
void binary_polmat_to_matpoly_simple_and_stupid(unsigned long * dst, mat64 const * src, unsigned int m, unsigned int n, unsigned int len)/*{{{*/
{
    /* dst must have room for len mat64's. */
    /* src is assumed row-major */
    ASSERT_ALWAYS((n%64)==0);
    ASSERT_ALWAYS((m%64)==0);
    size_t stride = iceildiv(len, ULONG_BITS);
    memset(dst, 0, m * n * stride * sizeof(unsigned long));
    mat64 const * q0 = src;
    unsigned long * p0 = dst;
    unsigned long mk = 1;
    for(unsigned int k = 0 ; k < len ; k++) {
        uint64_t mj = 1;
        unsigned long * p = p0;
        for(unsigned int i = 0 ; i < m ; i++) {
            mat64 const * q = q0 + (i / 64) * (n / 64);
            for(unsigned int j = 0 ; j < n ; j++) {
                *p ^= mk & -(((*q)[i%64]&mj) != 0);
                mj <<= 1;
                q  += (mj == 0);
                mj += mj == 0;
                p += stride;
            }
        }
        mk <<= 1;
        p0 += (mk == 0);
        q0 += (m/64) * (n/64);
        mk += (mk == 0);
    }
}/*}}}*/

/* implements binary_matpoly_to_polmat */
void binary_matpoly_to_polmat_nested_transpositions(mat64* dst, unsigned long const* src, unsigned int m, unsigned int n, unsigned int len) /*{{{*/
{
    unsigned int M = m / 64;
    unsigned int N = n / 64;
    unsigned int Lu = iceildiv(len, ULONG_BITS);
    static_assert(64 % ULONG_BITS == 0, "ULONG_BITS must divide 64");
    ASSERT_ALWAYS(Lu % (64 / ULONG_BITS) == 0);
    unsigned int L = Lu / (64 / ULONG_BITS);

    uint64_t* temp = new uint64_t[m * n * L];

    uint64_t* q = (uint64_t*)dst;

    /* We have (M*64*N)*(64)*(L)*(1) U-bit words */
    generic_transpose_words(q, (uint64_t const *) src, m * N, 64, L, 1);
    /* We have (M*64*N)*(L)*(64)*(1) 64-bit words */
    for (unsigned int k = 0; k < m * N * L; k++) {
        mat64_transpose(dst[k], dst[k]);
    }
    /* We have (M*64*N)*(L)*(64)*(1) 64-bit words */
    /* We have 1*(M*64*N)*(L*64)*(1) 64-bit words */
    generic_transpose_words_inplace(q, 1, m * N, L * 64, 1, temp);
    /* We have (L*64)*(M)*(64)*(N) 64-bit words */
    generic_transpose_words_inplace(q, L * 64 * M, 64, N, 1, temp);
    /* We have (L*64)*(M*N)*64 64-bit words */

    delete[] temp;
} /*}}}*/

/* implements binary_polmat_to_matpoly */
void binary_polmat_to_matpoly_nested_transpositions(unsigned long * dst_u, mat64 const * src, unsigned int m, unsigned int n, unsigned int len)/*{{{*/
{
    unsigned int M = m / 64;
    unsigned int N = n / 64;
    unsigned int Lu = iceildiv(len, ULONG_BITS);
    static_assert(64 % ULONG_BITS == 0, "ULONG_BITS must divide 64");
    ASSERT_ALWAYS(Lu % (64 / ULONG_BITS) == 0);
    unsigned int L = Lu / (64 / ULONG_BITS);

    uint64_t * temp = new uint64_t[m*n*L];

    uint64_t * dst = (uint64_t *) dst_u;

    /* We have (L*64)*(M*N)*64 64-bit words */
    /* We have (L*64*M)*(N)*(64)*1 64-bit words */
    generic_transpose_words(dst, (uint64_t const *) src, L * 64 * M, N, 64, 1);

    /* We have (L*64*M)*(64)*(N)*1 64-bit words */
    /* We have 1*(L*64)*(M*64*N)*1 64-bit words */
    generic_transpose_words_inplace(dst, 1, L*64, m * N, 1, temp);

    /* We have 1*(M*64*N)*(L*64)*1 64-bit words */
    if (((uintptr_t) dst) % alignof(mat64) == 0) {
        for(unsigned int k = 0 ; k < m*N*L ; k++) {
            mat64 & t = * (mat64 *) (dst + 64 * k);
            mat64_transpose(t, t);
        }
    } else {
        for(unsigned int k = 0 ; k < m*N*L ; k++) {
            mat64 t;
            memcpy((void*) &t, dst + 64 * k, sizeof(mat64));
            mat64_transpose(t, t);
            memcpy(dst + 64 * k, (void*) &t, sizeof(mat64));
        }
    }
    /* We have 1*(M*64*N)*(L*64)*1 64-bit words */

    /* We have (M*64*N)*(L)*(64)*1 64-bit words */
    generic_transpose_words_inplace(dst, m * N, L, 64, 1, temp);
    /* We have (M*64*N)*(64)*(L)*1 64-bit words */

    delete[] temp;
}/*}}}*/

/* implements binary_matpoly_transpose_to_polmat */
void binary_matpoly_transpose_to_polmat_nested_transpositions(mat64* dst, unsigned long const* src, unsigned int m, unsigned int n, unsigned int len) /*{{{*/
{
    unsigned int M = m / 64;
    unsigned int N = n / 64;
    unsigned int Lu = iceildiv(len, ULONG_BITS);
    static_assert(64 % ULONG_BITS == 0, "ULONG_BITS must divide 64");
    ASSERT_ALWAYS(Lu % (64 / ULONG_BITS) == 0);
    unsigned int L = Lu / (64 / ULONG_BITS);

    uint64_t* temp = new uint64_t[m * n * L];

    uint64_t* q = (uint64_t *) dst;

    /* We have 1*(M*64)*(N*64)*(L) 64-bit words */
    generic_transpose_words(q, (uint64_t const *) src, 1, m, n, L);
    /* We have (N*64*M)*(64)*(L)*(1) 64-bit words */
    generic_transpose_words_inplace(q, n * M, 64, L, 1, temp);
    /* We have (N*64*M)*(L)*(64)*(1) 64-bit words */
    for (unsigned int k = 0; k < n * M * L; k++) {
        mat64_transpose(dst[k], dst[k]);
    }
    /* We have (N*64*M)*(L)*(64)*(1) 64-bit words */
    /* We have 1*(N*64*M)*(L*64)*(1) 64-bit words */
    generic_transpose_words_inplace(q, 1, n * M, L * 64, 1, temp);
    /* We have (L*64)*(N)*(64)*(M) 64-bit words */
    generic_transpose_words_inplace(q, L * 64 * N, 64, M, 1, temp);
    /* We have (L*64)*(N*M)*64 64-bit words */

    delete[] temp;
} /*}}}*/

/* implements binary_polmat_to_matpoly_transpose */
void binary_polmat_to_matpoly_transpose_nested_transpositions(unsigned long * dst_u, mat64 const * src, unsigned int m, unsigned int n, unsigned int len)/*{{{*/
{
    unsigned int M = m / 64;
    unsigned int N = n / 64;
    unsigned int Lu = iceildiv(len, ULONG_BITS);
    static_assert(64 % ULONG_BITS == 0, "ULONG_BITS must divide 64");
    ASSERT_ALWAYS(Lu % (64 / ULONG_BITS) == 0);
    unsigned int L = Lu / (64 / ULONG_BITS);

    uint64_t * temp = new uint64_t[m*n*L];
    uint64_t * dst = (uint64_t *) dst_u;

    /* We have (L*64)*(M*N)*64 64-bit words */
    /* We have (L*64*M)*(N)*(64)*1 64-bit words */
    generic_transpose_words(dst, (uint64_t const *) src, L * 64 * M, N, 64, 1);

    /* We have (L*64*M)*(64)*(N)*1 64-bit words */
    /* We have 1*(L*64)*(M*64*N)*1 64-bit words */
    generic_transpose_words_inplace(dst, 1, L*64, m * N, 1, temp);

    /* We have 1*(M*64*N)*(L*64)*1 64-bit words */

    if (((uintptr_t) dst) % alignof(mat64) == 0) {
        for(unsigned int k = 0 ; k < m*N*L ; k++) {
            mat64 & t = * (mat64 *) (dst + 64 * k);
            mat64_transpose(t, t);
        }
    } else {
        for(unsigned int k = 0 ; k < m*N*L ; k++) {
            mat64 t;
            memcpy((void *) &t, dst + 64 * k, sizeof(mat64));
            mat64_transpose(t, t);
            memcpy(dst + 64 * k, (void *) &t, sizeof(mat64));
        }
    }


    /* We have 1*(M*64*N)*(L*64)*1 64-bit words */

    /* We have (M*64*N)*(L)*(64)*1 64-bit words */
    generic_transpose_words_inplace(dst, m * N, L, 64, 1, temp);
    /* We have (M*64*N)*(64)*(L)*1 64-bit words */

    generic_transpose_words_inplace(dst, 1, m, n, L, temp);

    delete[] temp;
}/*}}}*/

/* {{{ final choices -- these are really clear-cut */
void binary_polmat_to_matpoly(unsigned long * dst, mat64 const * src, unsigned int m, unsigned int n, unsigned int len)/*{{{*/
{
    binary_polmat_to_matpoly_nested_transpositions(dst, src, m, n, len);
}
/*}}}*/
void binary_matpoly_to_polmat(mat64 * dst, unsigned long const * src, unsigned int m, unsigned int n, unsigned int len)/*{{{*/
{
    binary_matpoly_to_polmat_nested_transpositions(dst, src, m, n, len);
}
/*}}}*/
void binary_polmat_to_matpoly_transpose(unsigned long * dst, mat64 const * src, unsigned int m, unsigned int n, unsigned int len)/*{{{*/
{
    binary_polmat_to_matpoly_transpose_nested_transpositions(dst, src, m, n, len);
}
/*}}}*/
void binary_matpoly_transpose_to_polmat(mat64 * dst, unsigned long const * src, unsigned int m, unsigned int n, unsigned int len)/*{{{*/
{
    binary_matpoly_transpose_to_polmat_nested_transpositions(dst, src, m, n, len);
}
/*}}}*/
/*}}}*/
