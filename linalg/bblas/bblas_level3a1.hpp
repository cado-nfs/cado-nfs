#ifndef BBLAS_LEVEL3A1_HPP_
#define BBLAS_LEVEL3A1_HPP_

// IWYU pragma: private, include "bblas.hpp"
// IWYU pragma: friend ".*/bblas.*"

#include "bblas_mat64.hpp"
#include <algorithm>

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

/* implementation details, variants */
void binary_matpoly_to_polmat_simple_and_stupid(mat64 * dst, unsigned long const * src, unsigned int m, unsigned int n, unsigned int len);
void binary_polmat_to_matpoly_simple_and_stupid(unsigned long * dst, mat64 const * src, unsigned int m, unsigned int n, unsigned int len);
void binary_matpoly_to_polmat_nested_transpositions(mat64 * dst, unsigned long const * src, unsigned int m, unsigned int n, unsigned int len);
void binary_polmat_to_matpoly_nested_transpositions(unsigned long * dst, mat64 const * src, unsigned int m, unsigned int n, unsigned int len);
void binary_matpoly_transpose_to_polmat_nested_transpositions(mat64 * dst, unsigned long const * src, unsigned int m, unsigned int n, unsigned int len);
void binary_polmat_to_matpoly_transpose_nested_transpositions(unsigned long * dst, mat64 const * src, unsigned int m, unsigned int n, unsigned int len);

/* {{{ generic transposition utilities -- these could have wider use.
 * However the code below has really nothing clever. Maybe it makes
 * just as much sense to copy and adapt it where appropriate...
 */
/* transform n0*n1*n2*n3 words into n0*n2*n1*n3 words */
template<typename T>
void generic_transpose_words(T * dst, T const * src, unsigned int n0, unsigned int n1, unsigned int n2, unsigned int n3)/*{{{*/
{
    ASSERT_ALWAYS(dst != src);
    if (n1 == 1 || n2 == 1) {
        std::copy_n(src, n0*n1*n2*n3, dst);
        return;
    }
    for(unsigned int i0 = 0 ; i0 < n0 ; i0++) {
        T * q = dst + i0 * n1 * n2 * n3;
        T const * p = src + i0 * n1 * n2 * n3;
        for(unsigned int i2 = 0 ; i2 < n2 ; i2++) {
            for(unsigned int i1 = 0 ; i1 < n1 ; i1++) {
                std::copy_n(p + (i1 * n2 + i2) * n3, n3, q);
                q += n3;
            }
        }
    }
}/*}}}*/

/* temp needs n1 * n2 * n3 words.  */
template<typename T>
void generic_transpose_words_inplace(T * x, unsigned int n0, unsigned int n1, unsigned int n2, unsigned int n3, T * temp) /*{{{*/
{
    if (n1 == 1 || n2 == 1) return;
    for(unsigned int i0 = 0 ; i0 < n0 ; i0++) {
        T * q = x + i0 * n1 * n2 * n3;
        T * t = temp;
        for(unsigned int i2 = 0 ; i2 < n2 ; i2++) {
            for(unsigned int i1 = 0 ; i1 < n1 ; i1++) {
                std::copy_n(q + (i1 * n2 + i2) * n3, n3, t);
                t += n3;
            }
        }
        std::copy_n(temp, n1*n2*n3, q);
    }
}/*}}}*/
/* }}} */

/* final exported choices */
void binary_polmat_to_matpoly(unsigned long * dst, mat64 const * src, unsigned int m, unsigned int n, unsigned int len);
void binary_matpoly_to_polmat(mat64 * dst, unsigned long const * src, unsigned int m, unsigned int n, unsigned int len);
void binary_polmat_to_matpoly_transpose(unsigned long * dst, mat64 const * src, unsigned int m, unsigned int n, unsigned int len);
void binary_matpoly_transpose_to_polmat(mat64 * dst, unsigned long const * src, unsigned int m, unsigned int n, unsigned int len);

#endif	/* BBLAS_LEVEL3A1_HPP_ */
