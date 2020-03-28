#include "cado.h"
#include "level3a1.hpp"

#include <cstring>

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
void binary_matpoly_to_polmat_simple_and_stupid(mat64 * dst, uint64_t const * src, unsigned int m, unsigned int n, unsigned int len)/*{{{*/
{
    /* dst must have room for len mat64's. */
    /* src is assumed row-major */
    ASSERT_ALWAYS((n%64)==0);
    ASSERT_ALWAYS((m%64)==0);
    size_t stride = iceildiv(len, 64);
    memset((void *) dst, 0, (m/64) * (n/64) * len * sizeof(mat64));
    mat64 * q0 = dst;
    uint64_t const * p0 = src;
    uint64_t mk = 1;
    for(unsigned int k = 0 ; k < len ; k++) {
        uint64_t mj = 1;
        uint64_t const * p = p0;
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
void binary_polmat_to_matpoly_simple_and_stupid(uint64_t * dst, mat64 const * src, unsigned int m, unsigned int n, unsigned int len)/*{{{*/
{
    /* dst must have room for len mat64's. */
    /* src is assumed row-major */
    ASSERT_ALWAYS((n%64)==0);
    ASSERT_ALWAYS((m%64)==0);
    size_t stride = iceildiv(len, 64);
    memset(dst, 0, m * n * stride * sizeof(uint64_t));
    mat64 const * q0 = src;
    uint64_t * p0 = dst;
    uint64_t mk = 1;
    for(unsigned int k = 0 ; k < len ; k++) {
        uint64_t mj = 1;
        uint64_t * p = p0;
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

/* {{{ generic transposition utilities -- these could have wider use */
/* transform n0*n1*n2*n3 words into n0*n2*n1*n3 words */
void generic_transpose_words(uint64_t * dst, uint64_t const * src, unsigned int n0, unsigned int n1, unsigned int n2, unsigned int n3)/*{{{*/
{
    ASSERT_ALWAYS(dst != src);
    if (n1 == 1 || n2 == 1) {
        memcpy(dst, src, n0*n1*n2*n3*sizeof(uint64_t));
        return;
    }
    for(unsigned int i0 = 0 ; i0 < n0 ; i0++) {
        uint64_t * q = dst + i0 * n1 * n2 * n3;
        uint64_t const * p = src + i0 * n1 * n2 * n3;
        for(unsigned int i2 = 0 ; i2 < n2 ; i2++) {
            for(unsigned int i1 = 0 ; i1 < n1 ; i1++) {
                memcpy(q, p + (i1 * n2 + i2) * n3, n3 * sizeof(uint64_t));
                q += n3;
            }
        }
    }
}/*}}}*/

/* temp needs n1 * n2 * n3 words.  */
void generic_transpose_words_inplace(uint64_t * x, unsigned int n0, unsigned int n1, unsigned int n2, unsigned int n3, uint64_t * temp) /*{{{*/
{
    if (n1 == 1 || n2 == 1) return;
    for(unsigned int i0 = 0 ; i0 < n0 ; i0++) {
        uint64_t * q = x + i0 * n1 * n2 * n3;
        uint64_t * t = temp;
        for(unsigned int i2 = 0 ; i2 < n2 ; i2++) {
            for(unsigned int i1 = 0 ; i1 < n1 ; i1++) {
                memcpy(t, q + (i1 * n2 + i2) * n3, n3 * sizeof(uint64_t));
                t += n3;
            }
        }
        memcpy(q, temp, n1 * n2 * n3 * sizeof(uint64_t));
    }
}/*}}}*/
/* }}} */

/* implements binary_matpoly_to_polmat */
void binary_matpoly_to_polmat_nested_transpositions(mat64* dst, uint64_t const* src, unsigned int m, unsigned int n, unsigned int len) /*{{{*/
{
    unsigned int M = m / 64;
    unsigned int N = n / 64;
    unsigned int L = iceildiv(len, 64);

    uint64_t* temp = new uint64_t[m * n * L];

    uint64_t* q = (uint64_t*)dst;

    /* We have (M*64*N)*(64)*(L)*(1) 64-bit words */
    generic_transpose_words(q, src, m * N, 64, L, 1);
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
void binary_polmat_to_matpoly_nested_transpositions(uint64_t * dst, mat64 const * src, unsigned int m, unsigned int n, unsigned int len)/*{{{*/
{
    unsigned int M = m / 64;
    unsigned int N = n / 64;
    unsigned int L = iceildiv(len, 64);

    uint64_t * temp = new uint64_t[m*n*L];

    /* We have (L*64)*(M*N)*64 64-bit words */
    /* We have (L*64*M)*(N)*(64)*1 64-bit words */
    generic_transpose_words(dst, (uint64_t const *) src, L * 64 * M, N, 64, 1);

    /* We have (L*64*M)*(64)*(N)*1 64-bit words */
    /* We have 1*(L*64)*(M*64*N)*1 64-bit words */
    generic_transpose_words_inplace(dst, 1, L*64, m * N, 1, temp);

    /* We have 1*(M*64*N)*(L*64)*1 64-bit words */
    for(unsigned int k = 0 ; k < m*N*L ; k++) {
        mat64 & t = * (mat64 *) (dst + 64 * k);
        mat64_transpose(t, t);
    }
    /* We have 1*(M*64*N)*(L*64)*1 64-bit words */

    /* We have (M*64*N)*(L)*(64)*1 64-bit words */
    generic_transpose_words_inplace(dst, m * N, L, 64, 1, temp);
    /* We have (M*64*N)*(64)*(L)*1 64-bit words */

    delete[] temp;
}/*}}}*/

/* {{{ final choices -- these are really clear-cut */
void binary_polmat_to_matpoly(uint64_t * dst, mat64 const * src, unsigned int m, unsigned int n, unsigned int len)/*{{{*/
{
    binary_polmat_to_matpoly_nested_transpositions(dst, src, m, n, len);
}
/*}}}*/
void binary_matpoly_to_polmat(mat64 * dst, uint64_t const * src, unsigned int m, unsigned int n, unsigned int len)/*{{{*/
{
    binary_matpoly_to_polmat_nested_transpositions(dst, src, m, n, len);
}
/*}}}*/
/*}}}*/
