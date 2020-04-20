#include "cado.h"
#include "bblas_level4.hpp"
#include "bblas_level4_ple_internal.hpp"
#include "bblas_level4_ple_internal_inl.hpp"
#include "gmp_aux.h"

/* Is it sufficient to meet the ODR requirement ? There are places where
 * we do std::min(B, foo). That requires that a definition of B be
 * available somewhere. It's not clear to me that the following kind of
 * template constexpr definition does the trick.
 */
template<typename matrix> constexpr const unsigned int bpack_view_base<matrix>::B;

template<typename matrix>
void bpack_view<matrix>::fill_random(gmp_randstate_t rstate) {
    memfill_random((void *) X, mblocks * nblocks * sizeof(matrix), rstate);
}

/* propose implementation for the PLE front-ends */
template<typename matrix>
std::vector<unsigned int> bpack_view<matrix>::ple()
{
    auto ple = PLE<matrix>(*this);
    return ple();
}

template<typename matrix>
void bpack_view<matrix>::propagate_row_permutations(std::vector<unsigned int> const & p)
{
    /* X has size (B*mblocks) * (B*nblocks), and blocks are stored
     * row-major (i.e. we have mblocks lists of nblocks consecutive mat64's).
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
 * the matrix X, by the inverse. Blocks X[i*nblocks+j] with j>i are not
 * touched. Blocks with j>nblocks do not exist in X, but are implicitly assumed
 * to correspond to a fragment of the identity matrix (and need not be
 * touched anyway)
 *
 * Note that the upper part of the diagonal blocks is replaced by zeroes.
 * Fixing that wouldn't be terribly hard.
 */
template<typename matrix>
void bpack_view<matrix>::invert_lower_triangular()
{
    for(unsigned int j = 0 ; j < nblocks && j < mblocks ; j++) {
        matrix Ljj_inv = 1;
        matrix::trsm(X[j * nblocks + j], Ljj_inv);
        X[j * nblocks + j] = Ljj_inv;
    }
    for(unsigned int j = 0 ; j < nblocks ; j++) {
        for(unsigned int i = j + 1 ; i < mblocks ; i++) {
            // we have i > j
            matrix S = 0;
            // L_{i,j} * R_{j,j} + L_{i,j+1} * R{j+1,j} + ... + L_{i,i} * R_{i,j} = 0
            for(unsigned int k = j ; k < i ; k++) {
                matrix::addmul(S, X[i * nblocks + k], X[k * nblocks + k]);
            }
            /* yes, X[i*nblocks+i] is lower triangular. maybe we could save a
             * few cycles here */
            matrix::mul(X[i * nblocks + j], X[i * nblocks + i], S);
        }
    }
}

template<typename matrix>
void bpack<matrix>::set_zero() {
    std::fill_n(X.begin(), mblocks * nblocks, 0);
}
/*
template<typename matrix>
void bpack_ops<matrix>::fill_random(bpack_view<matrix> A, gmp_randstate_t rstate) {
    A.fill_random(rstate);
}
*/

template<typename matrix>
bool bpack_const_view<matrix>::is_lowertriangular() const {
    for(unsigned int j = 0 ; j < nblocks ; j++) {
        for(unsigned int i = 0 ; i < mblocks ; i++) {
            if (j > i) {
                if (X[i*nblocks+j] != 0) return false;
            } else if (j == i) {
                if (!X[i*nblocks+j].is_lowertriangular()) return false;
            }
        }
    }
    return true;
}

template<typename matrix>
bool bpack_const_view<matrix>::is_uppertriangular() const {
    for(unsigned int j = 0 ; j < nblocks ; j++) {
        for(unsigned int i = 0 ; i < mblocks ; i++) {
            if (j < i) {
                if (X[i*nblocks+j] != 0) return false;
            } else if (j == i) {
                if (!X[i*nblocks+j].is_uppertriangular()) return false;
            }
        }
    }
    return true;
}

template<typename matrix>
bool bpack_const_view<matrix>::triangular_is_unit() const {
    for(unsigned int j = 0 ; j < nblocks && j < mblocks ; j++) {
        if (!X[j*nblocks+j].triangular_is_unit()) return false;
    }
    return true;
}

template<typename matrix>
void bpack_view<matrix>::make_lowertriangular() {
    for(unsigned int j = 0 ; j < nblocks ; j++) {
        for(unsigned int i = 0 ; i < mblocks ; i++) {
            if (j > i)
                X[i*nblocks+j] = 0;
            else if (j == i) {
                X[i*nblocks+j].make_lowertriangular();
            }
        }
    }
}

template<typename matrix>
void bpack_view<matrix>::make_uppertriangular() {
    for(unsigned int j = 0 ; j < nblocks ; j++) {
        for(unsigned int i = 0 ; i < mblocks ; i++) {
            if (j < i)
                X[i*nblocks+j] = 0;
            else if (j == i) {
                X[i*nblocks+j].make_uppertriangular();
            }
        }
    }
}

template<typename matrix>
void bpack_view<matrix>::make_unit_lowertriangular() {
    make_lowertriangular();
    triangular_make_unit();
}

template<typename matrix>
void bpack_view<matrix>::make_unit_uppertriangular() {
    make_uppertriangular();
    triangular_make_unit();
}

template<typename matrix>
void bpack_view<matrix>::triangular_make_unit() {
    for(unsigned int j = 0 ; j < nblocks && j < mblocks ; j++) {
        X[j*nblocks+j].triangular_make_unit();
    }
}

template struct bpack_view<mat64>;
template struct bpack_view<mat8>;
template struct bpack_const_view<mat64>;
template struct bpack_const_view<mat8>;
template struct bpack<mat64>;
template struct bpack<mat8>;

template struct PLE<mat64>;
template struct PLE<mat8>;
