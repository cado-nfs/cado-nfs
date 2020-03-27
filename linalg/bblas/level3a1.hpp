#ifndef BBLAS_LEVEL3A1_HPP_
#define BBLAS_LEVEL3A1_HPP_

#include "bblas.hpp"
#include "m64pol.hpp"

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
void binary_matpoly_to_polmat_simple_and_stupid(m64pol_ptr dst, uint64_t const * src, unsigned int m, unsigned int n, unsigned int len);
void binary_polmat_to_matpoly_simple_and_stupid(uint64_t * dst, m64pol_srcptr src, unsigned int m, unsigned int n, unsigned int len);
void binary_matpoly_to_polmat_nested_transpositions(m64pol_ptr dst, uint64_t const * src, unsigned int m, unsigned int n, unsigned int len);
void binary_polmat_to_matpoly_nested_transpositions(uint64_t * dst, m64pol_srcptr src, unsigned int m, unsigned int n, unsigned int len);

/* final exported choices */
void binary_polmat_to_matpoly(uint64_t * dst, m64pol_srcptr src, unsigned int m, unsigned int n, unsigned int len);
void binary_matpoly_to_polmat(m64pol_ptr dst, uint64_t const * src, unsigned int m, unsigned int n, unsigned int len);

#endif	/* BBLAS_LEVEL3A1_HPP_ */
