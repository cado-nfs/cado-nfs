#include "cado.h"
#include "level4.hpp"
#include "level4_ple_internal.hpp"

template<typename T>
constexpr const unsigned int PLE<T>::B;

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

std::vector<unsigned int> binary_blas_PLE(mat64 * X, unsigned int m, unsigned int n)
{
    auto ple = PLE<mat64>(X, m, n);
    return ple();
}
