#include "cado.h"
#include "bblas_level4.hpp"
#include "bblas_level4_ple_internal.hpp"
#include "bblas_level4_ple_internal_inl.hpp"

/* Is it sufficient to meet the ODR requirement ? There are places where
 * we do std::min(B, foo). That requires that a definition of B be
 * available somewhere. It's not clear to me that the following kind of
 * template constexpr definition does the trick.
 */
template<typename matrix> constexpr const unsigned int bpack<matrix>::B;

/* propose implementation for the PLE front-ends */
template<typename matrix>
std::vector<unsigned int> bpack<matrix>::ple()
{
    auto ple = PLE<matrix>(*this);
    return ple();
}

template<typename matrix>
void bpack<matrix>::propagate_row_permutations(std::vector<unsigned int> const & p)
{
    /* X has size 64*m * 64*n, and blocks are stored
     * row-major (i.e. we have m lists of n consecutive mat64's).
     * 
     * exchange row 0 with row p[0]
     * exchange row 1 with row p[1]
     * and so on.
     *
     * It turns out that we have a function in ple that does exactly this !
     */
    PLE<matrix>(*this).propagate_row_permutations(p.size(), UINT_MAX, p.begin(), p.end());
}

/* This replaces, in-place, the blocks of the lower triangular part of
 * the matrix X, by the inverse. Blocks X[i*n+j] with j>i are not
 * touched. Blocks with j>n do not exist in X, but are implicitly assumed
 * to correspond to a fragment of the identity matrix (and need not be
 * touched anyway)
 *
 * Note that the upper part of the diagonal blocks is replaced by zeroes.
 * Fixing that wouldn't be terribly hard.
 */
template<typename matrix>
void bpack<matrix>::invert_lower_triangular()
{
    for(unsigned int j = 0 ; j < n && j < m ; j++) {
        matrix Ljj_inv = 1;
        matrix::trsm(X[j * n + j], Ljj_inv);
        X[j * n + j] = Ljj_inv;
    }
    for(unsigned int j = 0 ; j < n ; j++) {
        for(unsigned int i = j + 1 ; i < m ; i++) {
            // we have i > j
            matrix S = 0;
            // L_{i,j} * R_{j,j} + L_{i,j+1} * R{j+1,j} + ... + L_{i,i} * R_{i,j} = 0
            for(unsigned int k = j ; k < i ; k++) {
                matrix::addmul(S, X[i * n + k], X[k * n + k]);
            }
            /* yes, X[i*n+i] is lower triangular. maybe we could save a
             * few cycles here */
            matrix::mul(X[i * n + j], X[i * n + i], S);
        }
    }
}


template struct bpack<mat64>;
template struct bpack<mat8>;

template struct PLE<mat64>;
template struct PLE<mat8>;
