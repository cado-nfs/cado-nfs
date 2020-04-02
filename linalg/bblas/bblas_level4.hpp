#ifndef BBLAS_LEVEL4_HPP_
#define BBLAS_LEVEL4_HPP_

#include "bblas.hpp"
#include "bblas_mat64.hpp"
#include "bblas_perm_matrix.hpp"
#include <vector>

/**********************************************************************/
/* level 4: factorizations and reductions of matrices
 *      gauss
 *      pluq
 *      lup
 *      echelon
 *      ple
 */

int gauss_6464_C(mat64 & mm, mat64 & e, mat64 const & m);
int gauss_6464_imm(mat64 & mm, mat64 & e, mat64 const & m);
int gauss_128128_C(mat64 * m);


/* PLUQ -- well we're not computing exactly PLUQ 
 * PLUQ says: Any m*n matrix A with rank r , can be written A = P*L*U*Q
 * where P and Q are two permutation matrices, of dimension respectively
 * m*m and n*n, L is m*r unit lower triangular and U is r*n upper
 * triangular.
 *
 * Here we compute p,l,u,q such that p*l*a*transpose(q) = an upper
 * triangular matrix, whose diagonal has r one and n-r zeros.
 */
/* outer routine */
int PLUQ64(perm_matrix_ptr p, mat64 & l, mat64 & u, perm_matrix_ptr q, mat64 const & m);

int PLUQ64_n(int * phi, mat64 & l, mat64 * u, mat64 const * a, int n);

/* This code is here because someday, I vaguely had the idea of using it
 * as a building block for the binary lingen. In fact, the code fragments
 * here for PLUQ and such have never been put in production, so I'm
 * pretty sure they're quite fragile.
 */
int PLUQ128(perm_matrix_ptr p, mat64 * l, mat64 * u, perm_matrix_ptr q, mat64 const * m);

/* Computes l,u,p, such that:
 *  - l is unit lower triangular
 *  - l*a=u
 *  - up=u*transpose(p) is upper triangular,
 *  -   up[k]==0 iff up[k,k] == 0
 *  -   rank(up)=rank(a)=# of diagonal 1's in up.
 */
int LUP64_imm(mat64 & l, mat64 & u, mat64 & p, mat64 const & a);

/* Computes e,mm such that mm=e*m is in row echelon form */
int full_echelon_6464_imm(mat64 & mm, mat64 & e, mat64 const & m);

/* Goal: obtain a PLE decomposition of the matrix X, together with a list
 * of the pivot rows.
 *
 * we assume that X has size 64*m * 64*n, and that blocks are stored
 * row-major (i.e. we have m lists of n consecutive mat64's).
 *
 * The L part is stored inside the matrix X.
 *
 * The permutations are stored implicitly. We know that permutations are
 * formed as (current index i, other index >= i). Hence it is sufficient
 * to store the list of other indices, up to the rank. This information
 * is sufficient to recover the list of pivot rows.
 */

std::vector<unsigned int> binary_blas_PLE(mat64 * X, unsigned int m, unsigned int n);

#endif	/* BBLAS_LEVEL4_HPP_ */
