#include "cado.h"
#include "bblas_level4.hpp"
#include "bblas_level4_ple_internal.hpp"
#include "bblas_mat64.hpp"
#include "bblas_mat8.hpp"

/**********************************************************************/
/* level 4: factorizations and reductions of matrices
 *      gauss
 *      pluq
 *      ple
 */

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

/* Code is in bblas_level4_ple_internal.hpp and
 * bblas_level4_ple_internal_inl.hpp.
 * It gets exposed only via bpack.cpp
 */
