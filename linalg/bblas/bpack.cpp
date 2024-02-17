#include "cado.h" // IWYU pragma: keep
#include <climits>                       // for UINT_MAX
#include <cstddef>                       // for size_t
#include <cstdint>
#include "bpack.hpp"
// #include "bblas_level4.hpp"
#include "bblas_bitmat.hpp"  // for bitmat
#include "bblas_level4_ple_internal.hpp" // PLE
#include "bblas_level4_ple_internal_inl.hpp" // IWYU pragma: keep
#include "gmp_aux.h"    // memfill_random
#include "omp_proxy.h" // IWYU pragma: keep

/* Is it sufficient to meet the ODR requirement ? There are places where
 * we do std::min(B, foo). That requires that a definition of B be
 * available somewhere. It's not clear to me that the following kind of
 * template constexpr definition does the trick.
 */
template<typename T> constexpr const unsigned int bpack_view_base<T>::B;
template<typename T> constexpr const unsigned int bpack<T>::B;

template<typename T>
void bpack_view<T>::fill_random(gmp_randstate_t rstate) {
    memfill_random((void *) X, mblocks * nblocks * sizeof(bitmat<T>), rstate);
}

/* propose implementation for the PLE front-ends */
template<typename T>
std::vector<unsigned int> bpack_view<T>::ple(std::vector<unsigned int> const & d)
{
    auto ple = PLE<T>(*this, d);
    return ple();
}

template<typename T>
std::vector<unsigned int> bpack_view<T>::ple()
{
    auto ple = PLE<T>(*this);
    return ple();
}

template<typename T>
void bpack_view<T>::propagate_row_permutations(std::vector<unsigned int> const & p)
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
    PLE<T>(*this).propagate_row_permutations(p.size(), UINT_MAX, p.begin(), p.end());
}

/* This replaces, in-place, the blocks of the lower triangular part of
 * the bitmat<T> X, by the inverse. Blocks cell(i, j) with j>i are not
 * touched. Blocks with j>nblocks do not exist in X, but are implicitly assumed
 * to correspond to a fragment of the identity bitmat<T> (and need not be
 * touched anyway)
 *
 * Note that the upper part of the diagonal blocks is replaced by zeroes.
 * Fixing that wouldn't be terribly hard.
 */
template<typename T>
void bpack_view<T>::invert_lower_triangular()
{
    for(unsigned int j = 0 ; j < nblocks && j < mblocks ; j++) {
        bitmat<T> Ljj_inv = 1;
        bitmat<T>::trsm(cell(j, j), Ljj_inv);
        cell(j, j) = Ljj_inv;
    }
    for(unsigned int j = 0 ; j < nblocks ; j++) {
        for(unsigned int i = j + 1 ; i < mblocks ; i++) {
            // we have i > j
            bitmat<T> S = 0;
            // L_{i,j} * R_{j,j} + L_{i,j+1} * R{j+1,j} + ... + L_{i,i} * R_{i,j} = 0
            for(unsigned int k = j ; k < i && k < nblocks ; k++) {
                /* if k > nblocks, cell(i, k) is zero */
                bitmat<T>::addmul(S, cell(i, k), cell(k, j));
            }
            if (i < nblocks) {
                /* mul_lt_ge is a very shallow win. */
                bitmat<T>::mul_lt_ge(cell(i, j), cell(i, i), S);
            } else {
                cell(i, j) = S;
            }
        }
    }
}

template<typename T>
bpack_view<T>& bpack_view<T>::set(int a)
{
    std::fill_n(X, mblocks * nblocks, 0);
    if (a&1)
        triangular_make_unit();
    return *this;
}

template<typename T>
bpack_view<T>& bpack_view<T>::set(bpack_const_view<T> v)
{
    ASSERT_ALWAYS(mblocks == v.mblocks);
    ASSERT_ALWAYS(nblocks == v.nblocks);
    for(unsigned int j = 0 ; j < nblocks ; j++) {
        for(unsigned int i = 0 ; i < mblocks ; i++) {
            cell(i, j) = v.cell(i, j);
        }
    }
    return *this;
}


template<typename T>
bool bpack_const_view<T>::operator==(int a) const {
    for(unsigned int j = 0 ; j < nblocks ; j++) {
        for(unsigned int i = 0 ; i < mblocks ; i++) {
            if (j == i) {
                if (cell(j, j) != a) return false;
            } else {
                if (cell(i, j) != 0) return false;
            }
        }
    }
    return true;
}
/*
template<typename T>
void bpack_ops<T>::fill_random(bpack_view<T> A, gmp_randstate_t rstate) {
    A.fill_random(rstate);
}
*/

template<typename T>
bool bpack_const_view<T>::is_lowertriangular() const {
    for(unsigned int j = 0 ; j < nblocks ; j++) {
        for(unsigned int i = 0 ; i < mblocks ; i++) {
            if (j > i) {
                if (cell(i, j) != 0) return false;
            } else if (j == i) {
                if (!cell(i, j).is_lowertriangular()) return false;
            }
        }
    }
    return true;
}

template<typename T>
bool bpack_const_view<T>::is_uppertriangular() const {
    for(unsigned int j = 0 ; j < nblocks ; j++) {
        for(unsigned int i = 0 ; i < mblocks ; i++) {
            if (j < i) {
                if (cell(i, j) != 0) return false;
            } else if (j == i) {
                if (!cell(i, j).is_uppertriangular()) return false;
            }
        }
    }
    return true;
}

template<typename T>
bool bpack_const_view<T>::triangular_is_unit() const {
    for(unsigned int j = 0 ; j < nblocks && j < mblocks ; j++) {
        if (!cell(j, j).triangular_is_unit()) return false;
    }
    return true;
}

template<typename T>
void bpack_view<T>::make_lowertriangular() {
    for(unsigned int j = 0 ; j < nblocks ; j++) {
        for(unsigned int i = 0 ; i < mblocks ; i++) {
            if (j > i)
                cell(i, j) = 0;
            else if (j == i) {
                cell(i, j).make_lowertriangular();
            }
        }
    }
}

template<typename T>
void bpack_view<T>::make_uppertriangular() {
    for(unsigned int j = 0 ; j < nblocks ; j++) {
        for(unsigned int i = 0 ; i < mblocks ; i++) {
            if (j < i)
                cell(i, j) = 0;
            else if (j == i) {
                cell(i, j).make_uppertriangular();
            }
        }
    }
}

template<typename T>
void bpack_view<T>::make_unit_lowertriangular() {
    make_lowertriangular();
    triangular_make_unit();
}

template<typename T>
void bpack_view<T>::make_unit_uppertriangular() {
    make_uppertriangular();
    triangular_make_unit();
}

template<typename T>
void bpack_view<T>::triangular_make_unit() {
    for(unsigned int j = 0 ; j < nblocks && j < mblocks ; j++) {
        cell(j, j).triangular_make_unit();
    }
}

template<typename T>
void bpack_ops<T>::mul(bpack_view<T> C, bpack_const_view<T> A, bpack_const_view<T> B)
{
    if (C.overlaps(A) || C.overlaps(B)) {
        bpack<T> CC(A.nrows(), B.ncols());
        mul(CC.view(), A, B);
        C.set(CC);
        return;
    }
    C.set(0);
    ASSERT_ALWAYS(C.nrowblocks() == A.nrowblocks());
    ASSERT_ALWAYS(C.ncolblocks() == B.ncolblocks());
    ASSERT_ALWAYS(A.ncolblocks() == B.nrowblocks());
    for(unsigned int bi = 0 ; bi < A.nrowblocks() ; bi++) {
        for(unsigned int bj = 0 ; bj < B.ncolblocks() ; bj++) {
            for(unsigned int bk = 0 ; bk < A.ncolblocks() ; bk++) {
                bitmat<T>::addmul(C.cell(bi, bj), A.cell(bi, bk), B.cell(bk, bj));
            }
        }
    }
}

template<typename T>
void bpack_ops<T>::mul_lt_ge(bpack_const_view<T> A, bpack_view<T> X)
{
    /* Do X <- A * X
     *
     * A is considered an an implicitly square bitmat<T>. It is
     * completed to the right to as many blocks as is necessary
     * to match the number of row blocks of X
     */
    ASSERT_ALWAYS(A.nrowblocks() == X.nrowblocks());
    ASSERT_ALWAYS(A.ncolblocks() <= A.nrowblocks());
#if 0
    for(unsigned int bi = X.nrowblocks() ; bi-- ; ) {
// #ifdef HAVE_OPENMP
// #pragma omp parallel for
// #endif
        for(unsigned int cbj = 0 ; cbj < X.ncolblocks() ; cbj++) {
            unsigned int bj = X.ncolblocks() - 1 - cbj;
            /* (AX)_{i,j} is the sum for 0<=k<=i of A_{i,k}X_{k,j}
             * (when X is also lt we even have j<=k).
             *
             * except for bi>=A.ncolblocks, where it is:
             * X_{i,j}+sum for 0<=k<A.ncolblocks A_{i,k}X_{k,j}
             *
             */
            bitmat<T> & S = X.cell(bi, bj);
            if (bi < A.ncolblocks())
                bitmat<T>::mul(S, A.cell(bi, bi), S);
            for(unsigned int bk = 0 ; bk < bi && bk < A.ncolblocks() ; bk++) {
                bitmat<T>::addmul(S, A.cell(bi, bk), X.cell(bk, bj));
            }
        }
    }
#else
    /* This approach is significantly faster when the multiplication code
     * benefits from doing precomputations on its right-hand side.
     */

    /*
    fprintf(stderr, "mul_lt_ge %u %u %u\n",
            A.nrowblocks(),
            A.ncolblocks(),
            X.ncolblocks());
     */

    /* We're going to process X in vertical bands of B columns (or,
     * more precisely, B/bitmat<T>::width blocks. In fact, we may
     * choose the number B to our liking (any multiple of
     * bitmat<T>::width), but it seems that the fewer the better.
     */
    constexpr const unsigned int B = bitmat<T>::width;
#ifdef HAVE_OPENMP
#pragma omp parallel num_threads(X.ncols() / B)
#endif
    {
        /* Use a temp variable */
        bpack<T> Y(X.nrows(), B);
        size_t A_stride = &A.cell(1,0) - &A.cell(0,0);
        size_t T_stride = &Y.cell(1,0) - &Y.cell(0,0);
#ifdef HAVE_OPENMP
#pragma omp for
#endif
    for(unsigned int bj = 0 ; bj < X.ncolblocks() ; bj += Y.ncolblocks()) {
        /* refresh our temp variable */
        Y = 0;
        unsigned int ndbj = std::min(Y.ncolblocks(), X.ncolblocks() - bj);
        for(unsigned int bk = 0 ; bk < A.ncolblocks() ; bk++) {
            for(unsigned int dbj = 0 ; dbj < ndbj ; dbj++) {
                bitmat<T>::addmul_blocks(&Y.cell(0,dbj), &A.cell(0, bk), X.cell(bk, bj + dbj), A.nrowblocks(), T_stride, A_stride);
            }
        }
        /* We could conceivably parallelize the loops below, but openmp
         * won't let us do it, and I don't know how I can work around
         * this limitation (it would require cooperation from the
         * enclosing loop). Y is the trouble maker here.
         */
        for(unsigned int bi = 0 ; bi < A.ncolblocks() ; bi++) {
            for(unsigned int dbj = 0 ; dbj < ndbj ; dbj++) {
                X.cell(bi, bj + dbj) = Y.cell(bi, dbj);
            }
        }
        for(unsigned int bi = A.ncolblocks() ; bi < X.nrowblocks() ; bi++) {
            for(unsigned int dbj = 0 ; dbj < ndbj ; dbj++) {
                bitmat<T>::add(X.cell(bi, bj + dbj), X.cell(bi, bj + dbj), Y.cell(bi, dbj));
            }
        }
    }
    }
#endif
}

template<typename T>
void bpack_ops<T>::extract_uppertriangular(bpack_view<T> a, bpack_const_view<T> const b)
{
    ASSERT_ALWAYS(a.nrowblocks() == b.nrowblocks());
    ASSERT_ALWAYS(a.ncolblocks() == b.ncolblocks());
    for(unsigned int bi = 0 ; bi < b.nrowblocks() ; bi++) {
        for(unsigned int bj = 0 ; bj < b.ncolblocks() ; bj++) {
            if (bi > bj)
                a.cell(bi, bj) = 0;
            else if (bi == bj)
                bitmat<T>::extract_uppertriangular(a.cell(bi, bj), b.cell(bi, bj));
            else
                a.cell(bi, bj) = b.cell(bi, bj);
        }
    }
}

template<typename T>
void bpack_ops<T>::extract_lowertriangular(bpack_view<T> a, bpack_const_view<T> const b)
{
    ASSERT_ALWAYS(a.nrowblocks() == b.nrowblocks());
    ASSERT_ALWAYS(a.ncolblocks() == b.ncolblocks());
    for(unsigned int bi = 0 ; bi < b.nrowblocks() ; bi++) {
        for(unsigned int bj = 0 ; bj < b.ncolblocks() ; bj++) {
            if (bi < bj)
                a.cell(bi, bj) = 0;
            else if (bi == bj)
                bitmat<T>::extract_lowertriangular(a.cell(bi, bj), b.cell(bi, bj));
            else
                a.cell(bi, bj) = b.cell(bi, bj);
        }
    }
}

template<typename T>
void bpack_ops<T>::extract_LU(bpack_view<T> L, bpack_view<T> U)
{
    ASSERT_ALWAYS(L.nrowblocks() == U.nrowblocks());
    ASSERT_ALWAYS(L.ncolblocks() == U.ncolblocks());
    for(unsigned int bi = 0 ; bi < U.nrowblocks() ; bi++) {
        for(unsigned int bj = 0 ; bj < U.ncolblocks() ; bj++) {
            if (bi < bj) {
                L.cell(bi, bj) = 0;
                /* U unchanged */
            } else if (bi == bj) {
                bitmat<T>::extract_LU(L.cell(bi, bj), U.cell(bi, bj));
            } else {
                L.cell(bi, bj) = U.cell(bi, bj);
                U.cell(bi, bj) = 0;
            }
        }
    }
}

template<typename T>
bool bpack_const_view<T>::operator==(bpack_const_view<T> v) const
{
    if (mblocks != v.mblocks) return false;
    if (nblocks != v.nblocks) return false;
    for(unsigned int bi = 0 ; bi < mblocks ; bi++) {
        for(unsigned int bj = 0 ; bj < nblocks ; bj++) {
            if (cell(bi, bj) != v.cell(bi, bj)) return false;
        }
    }
    return true;
}

template struct bpack_ops<uint64_t>;
template struct bpack_ops<uint8_t>;
template struct bpack_view<uint64_t>;
template struct bpack_view<uint8_t>;
template struct bpack_const_view<uint64_t>;
template struct bpack_const_view<uint8_t>;
template struct bpack<uint64_t>;
template struct bpack<uint8_t>;

template struct PLE<uint64_t>;
template struct PLE<uint8_t>;
