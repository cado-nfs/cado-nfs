#ifndef LEVEL4_PLE_INTERNAL_INL_HPP_
#define LEVEL4_PLE_INTERNAL_INL_HPP_

#include "cado_config.h"
#include "bblas_level4_ple_internal.hpp"
#include "bblas_simd.hpp"
#include "bblas_mat64.hpp"
#include "bblas_mat8.hpp"

#ifdef TIME_PLE
#include "timing.h"
struct timer_ple {
    double & t;
    timer_ple(double & t) : t(t) { t -= seconds(); }
    ~timer_ple() { t += seconds(); }
};
#define TIMER_PLE(t) timer_ple dummy(t)
#else
#define TIMER_PLE(t)    /**/
#endif

template<typename matrix>
int PLE<matrix>::find_pivot(unsigned int bi, unsigned int bj, unsigned int i, unsigned int j) const/*{{{*/
{
    TIMER_PLE(t_find_pivot);
    U mask = U(1) << j;
    for( ; bi < m ; bi++) {
        matrix & Y = X[bi * n + bj];
        for( ; i < B ; i++) {
            if (Y[i] & mask)
                return bi * B + i;
        }
        i = 0;
    }
    return -1;
}/*}}}*/

template<typename matrix>
void PLE<matrix>::propagate_pivot(unsigned int bi, unsigned int bj, unsigned int i, unsigned int j) const/*{{{*/
{
    TIMER_PLE(t_propagate_pivot);
    /* pivot row ii=bi*B+i to all rows below, but only for bits that are
     * right after column jj=bj*B+j.
     *
     *
     * This is done _only for the current column block bj !!!_
     */

    if (j == B-1) return;
    U mask = U(1) << j;
    ASSERT_ALWAYS(X[bi * n + bj][i] & mask);
    U c = X[bi * n + bj][i] & ((-mask) << 1);
    i++;
    if (i == B) { i = 0 ; bi++; }
    for( ; bi < m ; bi++) {
        matrix & Y = X[bi * n + bj];
        for( ; i < B ; i++) {
            Y[i] ^= c & -((Y[i] & mask) != 0);
        }
        i = 0;
    }
}/*}}}*/

/* These two specializations are very significant improvements.  */
#ifdef HAVE_AVX2
template<>
void PLE<mat64>::propagate_pivot(unsigned int bi, unsigned int bj, unsigned int i, unsigned int j) const/*{{{*/
{
    TIMER_PLE(t_propagate_pivot);
    typedef mat64 matrix;
    /* pivot row ii=bi*B+i to all rows below, but only for bits that are
     * right after column jj=bj*B+j.
     *
     *
     * This is done _only for the current column block bj !!!_
     */

    if (j == B-1) return;
    U mask = U(1) << j;
    ASSERT_ALWAYS(X[bi * n + bj][i] & mask);
    U c = X[bi * n + bj][i] & ((-mask) << 1);
    i++;
    constexpr const unsigned int SIMD = 4;
    for( ; (i & (SIMD-1)) && bi < m ; i++) {
        matrix & Y = X[bi * n + bj];
        Y[i] ^= c & -((Y[i] & mask) != 0);
    }
    if (i == B) { i = 0 ; bi++; }
    unsigned int iw = i / SIMD;
    __m256i mmask = _mm256_set1_epi64x(mask);
    __m256i cc = _mm256_set1_epi64x(c);
    for( ; bi < m ; bi++) {
        matrix & Y = X[bi * n + bj];
        __m256i *Yw = (__m256i *) Y.data();
        for( ; iw < B / SIMD ; iw++) {
            __m256i select = _mm256_and_si256(Yw[iw], mmask);
            select = _mm256_cmpeq_epi64(select, mmask);
            Yw[iw] = _mm256_xor_si256(Yw[iw], _mm256_and_si256(cc, select));
        }
        iw = 0;
    }
}/*}}}*/
#elif defined(HAVE_SSE41)
template<>
void PLE<mat64>::propagate_pivot(unsigned int bi, unsigned int bj, unsigned int i, unsigned int j) const/*{{{*/
{
    TIMER_PLE(t_propagate_pivot);
    typedef mat64 matrix;
    /* pivot row ii=bi*B+i to all rows below, but only for bits that are
     * right after column jj=bj*B+j.
     *
     *
     * This is done _only for the current column block bj !!!_
     */

    if (j == B-1) return;
    U mask = U(1) << j;
    ASSERT_ALWAYS(X[bi * n + bj][i] & mask);
    U c = X[bi * n + bj][i] & ((-mask) << 1);
    i++;
    constexpr const unsigned int SIMD = 2;
    for( ; (i & (SIMD-1)) && bi < m ; i++) {
        matrix & Y = X[bi * n + bj];
        Y[i] ^= c & -((Y[i] & mask) != 0);
    }
    if (i == B) { i = 0 ; bi++; }
    unsigned int iw = i / SIMD;
    __m128i mmask = _cado_mm_set1_epi64(mask);
    __m128i cc = _cado_mm_set1_epi64(c);
    for( ; bi < m ; bi++) {
        matrix & Y = X[bi * n + bj];
        __m128i *Yw = (__m128i *) Y.data();
        for( ; iw < B / SIMD ; iw++) {
            __m128i select = _mm_and_si128(Yw[iw], mmask);
            select = _mm_cmpeq_epi64(select, mmask);
            Yw[iw] = _mm_xor_si128(Yw[iw], _mm_and_si128(cc, select));
        }
        iw = 0;
    }
}/*}}}*/
#endif

template<>
void PLE<mat8>::propagate_pivot(unsigned int bi, unsigned int bj, unsigned int i, unsigned int j) const/*{{{*/
{
    TIMER_PLE(t_propagate_pivot);
    typedef mat8 matrix;
    /* pivot row ii=bi*B+i to all rows below, but only for bits that are
     * right after column jj=bj*B+j.
     *
     *
     * This is done _only for the current column block bj !!!_
     */

    if (j == B-1) return;
    U mask = U(1) << j;
    ASSERT_ALWAYS(X[bi * n + bj][i] & mask);
    U c = X[bi * n + bj][i] & ((-mask) << 1);
    i++;
    if (i == B) { i = 0 ; bi++; }
    __m64 mmask = _mm_set1_pi8(mask);
    __m64 cc = _mm_set1_pi8(c);
    if (i && bi < m) {
        matrix & Y = X[bi * n + bj];
        __m64 & Yw = * (__m64 *) Y.data();
        __m64 select = _mm_and_si64(Yw, mmask);
        select = _mm_cmpeq_pi8(select, mmask);
        /* Make sure we don't touch rows before i ! */
        select = _mm_and_si64(select, _mm_cvtsi64_m64(-(uint64_t(1) << (8*i))));
        Yw = _mm_xor_si64(Yw, _mm_and_si64(cc, select));
        bi++;
    }
    for( ; bi < m ; bi++) {
        matrix & Y = X[bi * n + bj];
        __m64 & Yw = * (__m64 *) Y.data();
        __m64 select = _mm_and_si64(Yw, mmask);
        select = _mm_cmpeq_pi8(select, mmask);
        Yw = _mm_xor_si64(Yw, _mm_and_si64(cc, select));
    }
}/*}}}*/

template<typename matrix>
void PLE<matrix>::propagate_row_permutations(unsigned int ii1, unsigned int bj0, std::vector<unsigned int>::const_iterator q0, std::vector<unsigned int>::const_iterator q1) const/*{{{*/
{
    TIMER_PLE(t_propagate_permutation);
    /* This propagates the pending permutations outside the current block
     * column.
     * Permutations are given by the range [q0..q1). More precisely,
     * the starting index is given by ii0 = ii1 - (q1 - q0), and the
     * image of index ii is given by
     *          *(q0 + ii - ii0) = *(q1 + ii - ii1)
     */
    for(unsigned bj = 0 ; bj < n ; bj++) {
        if (bj == bj0) continue;
        unsigned int ii = ii1 - (q1 - q0);
        for(auto q = q0 ; q != q1 ; q++, ii++) {
            if (ii == *q) continue;
            unsigned int bi = ii / B;
            unsigned int i  = ii & (B - 1);
            unsigned int pbi = *q / B;
            unsigned int pi  = *q & (B - 1);
            matrix & Y  = X[ bi * n + bj];
            matrix & pY = X[pbi * n + bj];
            U c = Y[i] ^ pY[pi];
            Y[i]   ^= c;
            pY[pi] ^= c;
        }
    }
}/*}}}*/

template<typename matrix>
void PLE<matrix>::move_L_fragments(unsigned int yii0, std::vector<unsigned int> const & Q) const/*{{{*/
{
    TIMER_PLE(t_move_l_fragments);
    /* This function receives the (yii0,yii0) coordinate of the first
     * item in the unit lower triangular part that we haven't completed
     * yet (or that we aren't sure that we have completed -- maybe we
     * have because it's already in place).
     * This function is called just before
     * a block is finished (either row block or column block). Either way,
     * we want to put the multipliers in the right positions. All multiplier
     * columns referenced by the list Q, with logical leading
     * coordinates (yyi0+x, Q[x]) must move to column
     * yyi0+x so that they are part of the unit lower triangular
     * matrix that we expect to find eventually.
     * Note 1: the multiplier column is _under_ the leading coordinate
     * (yyi0+x, Q[x]): only coordinates (yyi0+x+1, Q[x]) and downwards
     * are moved).
     * Note 2: by construction, (yii0+x)/B is constant as x runs through
     * [0..Q.size()-1], and so is Q[x]/B.
     *
     * This function assumes that the matrix coefficients that receive
     * coefficients from moved columns are set to zero beforehand, and
     * ensures that matrix entries from moved columns are set to zero
     * after the move (if they don't receive new column entries).
     */
    unsigned int ybi = yii0 / B;
    unsigned int yi  = yii0 & (B-1);
    unsigned int k = Q.size();
    unsigned int yii = yii0;
    ASSERT(!Q.empty());
    /* move the L fragments */
    unsigned int zbj = Q[0] / B;
    for(unsigned int dy = 0, dz ; dy < k ; dy+=dz, yi+=dz, yii+=dz) {
        unsigned int zjj = Q[dy];
        ASSERT(zjj / B == zbj);
        unsigned int zj  = zjj & (B-1);
        dz = 1;
        for( ; dy + dz < k && (Q[dy+dz] == Q[dy] + dz) ; dz++);
        if (yii == zjj) continue;
        /* We'll move a block of dz columns together. We first need to
         * "warm up" and transfer the upper unit lower triangular part as
         * we slowly enlarge the mask. */
        U mask = 1;
        for(unsigned int warmup = 1 ; warmup < dz ; warmup++) {
            /* we have yi + warmup = yii0 + dy + warmup
             *                     < yii0 + dy + dz
             *                     <= yii0 + dy + dz - 1
             *                     <= yii0 + k - 1
             *                     <= yii1 - 1
             * which on all accounts should still be in the same block.
             */
            ASSERT((yii + warmup) / B == ybi);
            U c = (X[ybi * n + zbj][yi + warmup] >> zj) & mask;
            X[ybi * n + ybi][yi + warmup] ^= (c << yi);
            X[ybi * n + zbj][yi + warmup] ^= (c << zj);
            /* is there a simpler way ? mask += mask + 1 ? */
            mask |= mask << 1;
        }
        /* we now move the part that comes _after_ the leading unit lower
         * triangular block.
         */
        unsigned int tbi = ybi;
        unsigned int ti = yi + dz;
        if (ti == B) { tbi++; ti = 0; }
        for( ; tbi < m ; tbi++) {
            for( ; ti < B ; ti++) {
                U c = (X[tbi * n + zbj][ti] >> zj) & mask;
                X[tbi * n + ybi][ti] ^= (c << yi);
                X[tbi * n + zbj][ti] ^= (c << zj);
            }
            ti = 0;
        }
    }
}/*}}}*/

#ifdef HAVE_AVX2
template<>
void PLE<mat64>::move_L_fragments(unsigned int yii0, std::vector<unsigned int> const & Q) const/*{{{*/
{
    TIMER_PLE(t_move_l_fragments);
    unsigned int ybi = yii0 / B;
    unsigned int yi  = yii0 & (B-1);
    unsigned int k = Q.size();
    unsigned int yii = yii0;
    ASSERT(!Q.empty());
    unsigned int zbj = Q[0] / B;
    for(unsigned int dy = 0, dz ; dy < k ; dy+=dz, yi+=dz, yii+=dz) {
        unsigned int zjj = Q[dy];
        ASSERT(zjj / B == zbj);
        unsigned int zj  = zjj & (B-1);
        dz = 1;
        for( ; dy + dz < k && (Q[dy+dz] == Q[dy] + dz) ; dz++);
        if (yii == zjj) continue;
        U mask = 1;
        for(unsigned int warmup = 1 ; warmup < dz ; warmup++) {
            /* we have yi + warmup = yii0 + dy + warmup
             *                     < yii0 + dy + dz
             *                     <= yii0 + dy + dz - 1
             *                     <= yii0 + k - 1
             *                     <= yii1 - 1
             * which on all accounts should still be in the same block.
             */
            ASSERT((yii + warmup) / B == ybi);
            U c = (X[ybi * n + zbj][yi + warmup] >> zj) & mask;
            X[ybi * n + ybi][yi + warmup] ^= (c << yi);
            X[ybi * n + zbj][yi + warmup] ^= (c << zj);
            /* is there a simpler way ? mask += mask + 1 ? */
            mask |= mask << 1;
        }
        unsigned int tbi = ybi;
        unsigned int ti = yi + dz;
        constexpr const unsigned int SIMD = 4;
        for ( ; (ti & (SIMD - 1)) && tbi < m ; ti++) {
            U c = (X[tbi * n + zbj][ti] >> zj) & mask;
            X[tbi * n + ybi][ti] ^= (c << yi);
            X[tbi * n + zbj][ti] ^= (c << zj);
        }
        if (ti == B) { tbi++; ti = 0; }
        __m256i mmask = _mm256_set1_epi64x(mask);
        unsigned int tiw = ti / SIMD;
        __m128i zjw = _mm_set1_epi64x(zj); // only lower 64 used
        __m128i yiw = _mm_set1_epi64x(yi); // only lower 64 used
        for( ; tbi < m ; tbi++) {
            __m256i *Xzjw = (__m256i *) X[tbi * n + zbj].data();
            __m256i *Xyiw = (__m256i *) X[tbi * n + ybi].data();
            for( ; tiw < B / SIMD ; tiw++) {
                __m256i & Xjw = Xzjw[tiw];
                __m256i & Xiw = Xyiw[tiw];
                __m256i cc = _mm256_and_si256(_mm256_srl_epi64(Xjw, zjw), mmask);
                Xiw = _mm256_xor_si256(Xiw, _mm256_sll_epi64(cc, yiw));
                Xjw = _mm256_xor_si256(Xjw, _mm256_sll_epi64(cc, zjw));
            }
            tiw = 0;
        }
    }
}/*}}}*/
#elif defined(HAVE_SSE41)
template<>
void PLE<mat64>::move_L_fragments(unsigned int yii0, std::vector<unsigned int> const & Q) const/*{{{*/
{
    TIMER_PLE(t_move_l_fragments);
    unsigned int ybi = yii0 / B;
    unsigned int yi  = yii0 & (B-1);
    unsigned int k = Q.size();
    unsigned int yii = yii0;
    ASSERT(!Q.empty());
    unsigned int zbj = Q[0] / B;
    for(unsigned int dy = 0, dz ; dy < k ; dy+=dz, yi+=dz, yii+=dz) {
        unsigned int zjj = Q[dy];
        ASSERT(zjj / B == zbj);
        unsigned int zj  = zjj & (B-1);
        dz = 1;
        for( ; dy + dz < k && (Q[dy+dz] == Q[dy] + dz) ; dz++);
        if (yii == zjj) continue;
        U mask = 1;
        for(unsigned int warmup = 1 ; warmup < dz ; warmup++) {
            /* we have yi + warmup = yii0 + dy + warmup
             *                     < yii0 + dy + dz
             *                     <= yii0 + dy + dz - 1
             *                     <= yii0 + k - 1
             *                     <= yii1 - 1
             * which on all accounts should still be in the same block.
             */
            ASSERT((yii + warmup) / B == ybi);
            U c = (X[ybi * n + zbj][yi + warmup] >> zj) & mask;
            X[ybi * n + ybi][yi + warmup] ^= (c << yi);
            X[ybi * n + zbj][yi + warmup] ^= (c << zj);
            /* is there a simpler way ? mask += mask + 1 ? */
            mask |= mask << 1;
        }
        unsigned int tbi = ybi;
        unsigned int ti = yi + dz;
        constexpr const unsigned int SIMD = 2;
        for ( ; (ti & (SIMD - 1)) && tbi < m ; ti++) {
            U c = (X[tbi * n + zbj][ti] >> zj) & mask;
            X[tbi * n + ybi][ti] ^= (c << yi);
            X[tbi * n + zbj][ti] ^= (c << zj);
        }
        if (ti == B) { tbi++; ti = 0; }
        __m128i mmask = _cado_mm_set1_epi64(mask);
        unsigned int tiw = ti / SIMD;
        __m128i zjw = _cado_mm_set1_epi64(zj); // only lower 64 used
        __m128i yiw = _cado_mm_set1_epi64(yi); // only lower 64 used
        for( ; tbi < m ; tbi++) {
            __m128i *Xzjw = (__m128i *) X[tbi * n + zbj].data();
            __m128i *Xyiw = (__m128i *) X[tbi * n + ybi].data();
            for( ; tiw < B / SIMD ; tiw++) {
                __m128i & Xjw = Xzjw[tiw];
                __m128i & Xiw = Xyiw[tiw];
                __m128i cc = _mm_and_si128(_mm_srl_epi64(Xjw, zjw), mmask);
                Xiw = _mm_xor_si128(Xiw, _mm_sll_epi64(cc, yiw));
                Xjw = _mm_xor_si128(Xjw, _mm_sll_epi64(cc, zjw));
            }
            tiw = 0;
        }
    }
}/*}}}*/
#endif

template<typename matrix>
void PLE<matrix>::trsm(unsigned int bi,/*{{{*/
        unsigned int bj,
        unsigned int yi0,
        unsigned int yi1) const
{
    TIMER_PLE(t_trsm);
    /* trsm is fairly trivial */
    for(unsigned int s = bj + 1 ; s < n ; s++) {
        matrix::trsm(X[bi * n + bi], X[bi * n + s], yi0, yi1);
    }
}/*}}}*/

template<typename matrix>
void PLE<matrix>::sub(unsigned int bi,/*{{{*/
        unsigned int bj,
        unsigned int yi0,
        unsigned int yi1,
        unsigned int ii) const
{
    TIMER_PLE(t_sub);
    unsigned int sbi = ii / B;
    unsigned int si  = ii & (B-1);
    for( ; sbi < m ; sbi++) {
        // i
        for(unsigned int sbj = bj + 1 ; sbj < n ; sbj++) {
            /* multiply columns [yi0..yi1-1] of block (sbi,bi) by
             * rows [yi0..yi1-1] of block (bi,sbj), and add that to
             * block (sbi,sbj)
             *
             * yes, we really mean block (sbi,bi). The indices [yi0..yi1)
             * are _really_ relative to that block.
             */
            matrix::addmul(X[sbi * n + sbj],
                    X[sbi * n + bi],
                    X[bi * n + sbj],
                    si, B, yi0, yi1);
        }
        si = 0;
    }
}/*}}}*/

template<typename matrix>
void PLE<matrix>::debug_stuff::apply_permutations(std::vector<unsigned int>::const_iterator p0, std::vector<unsigned int>::const_iterator p1)/*{{{*/
{
    /* apply the permutations to Xcc */
    for(unsigned int xii = 0 ; xii < (unsigned int) (p1 - p0) ; xii++) {
        unsigned int pii = p0[xii];
        if (xii == pii) continue;
        unsigned int xbi = xii / B;
        unsigned int xi = xii % B;
        unsigned int pbi = pii / B;
        unsigned int pi = pii % B;
        for(unsigned int bj = 0 ; bj < n ; bj++) {
            matrix & xY = X_target[xbi * n + bj];
            matrix & pY = X_target[pbi * n + bj];
            U c = xY[xi] ^ pY[pi];
            xY[xi] ^= c;
            pY[pi] ^= c;
        }
    }
}/*}}}*/

/* extract below the diagonal, only up to rank rr. The X field is
 * modified. */
template<typename matrix>
typename matrix::vector_type PLE<matrix>::debug_stuff::get_LL(unsigned int rr)/*{{{*/
{
    typename matrix::vector_type LL((m)*(m), 0);
    for(unsigned int bi = 0 ; bi < m ; bi++) {
        unsigned int bj = 0;
        for( ; bj <= bi && bj < iceildiv(rr, B) ; bj++) {
            for(unsigned int i = 0 ; i < B ; i++) {
                U c = X[bi*(n)+bj][i];
                unsigned int z = rr - bj * B;
                if (bj == bi && i < z)
                    z = i;
                if (z < B)
                    c &= (U(1) << z) - 1;
                LL[bi*(m)+bj][i] = c;
            }
        }
    }
    /* clear the blocks that we have just taken */
    for(unsigned int bi = 0 ; bi < m ; bi++) {
        for(unsigned int bj = 0 ; bj < n && bj < m ; bj++) {
            matrix::add(X[bi*n+bj], X[bi*n+bj], LL[bi*m+bj]);
        }
    }

    /* add implicit identity to L */
    for(unsigned int bi = 0 ; bi < m ; bi++) {
        for(unsigned int i = 0 ; i < B ; i++) {
            LL[bi*(m)+bi][i] ^= U(1) << i;
        }
    }

    return LL;
}/*}}}*/

template<typename matrix>
typename matrix::vector_type PLE<matrix>::debug_stuff::get_UU(unsigned int rr)/*{{{*/
{
    /* extract above the diagonal, only up to rank rr */
    typename matrix::vector_type UU((m)*(n), 0);
    for(unsigned int bi = 0 ; bi < m && bi < iceildiv(rr, B); bi++) {
        if (bi < n) {
            for(unsigned int i = 0 ; i < std::min(B, rr - bi * B) ; i++) {
                U c = X[bi*(n)+bi][i];
                c &= -(U(1) << i);
                UU[bi*(n)+bi][i] = c;
            }
        }
        for(unsigned int bj = bi + 1 ; bj < n ; bj++) {
            UU[bi*(n)+bj] = X[bi*(n)+bj];
            for(unsigned int i = rr - bi * B ; i < B ; i++) {
                UU[bi*(n)+bj][i] = 0;
            }
        }
    } 
    /* Finally, clear everything in Xcc that we haven't taken yet *//*{{{*/
    for(unsigned int bi = 0 ; bi < m ; bi++) {
        for(unsigned int bj = 0 ; bj < n ; bj++) {
            matrix::add(X[bi*n+bj], X[bi*n+bj], UU[bi*n+bj]);
        }
    }/*}}}*/
    return UU;
}/*}}}*/

template<typename matrix>
bool PLE<matrix>::debug_stuff::complete_check(typename matrix::vector_type const & LL, typename matrix::vector_type const & UU)/*{{{*/
{
    /* check that LL*UU + (remaining block in X) is equal to X_target */

    for(unsigned int bi = 0 ; bi < m ; bi++) {
        for(unsigned int bj = 0 ; bj < n ; bj++) {
            matrix C = X[bi*(n)+bj];
            for(unsigned int bk = 0 ; bk < m ; bk++) {
                matrix::addmul(C, LL[bi*(m)+bk], UU[bk*(n)+bj]);
            }
            ASSERT_ALWAYS(X_target[bi*(n)+bj] == C);
            if (X_target[bi*(n)+bj] != C) return false;
        }
    }
    return true;
}/*}}}*/

template<typename matrix>
bool PLE<matrix>::debug_stuff::check(matrix const * X0, std::vector<unsigned int>::const_iterator p0, unsigned int ii)/*{{{*/
{
    start_check(X0);
    apply_permutations(p0, p0 + ii);
    auto LL = get_LL(ii);
    auto UU = get_UU(ii);
    return complete_check(LL, UU);
}/*}}}*/

template<typename matrix>
std::vector<unsigned int> PLE<matrix>::operator()(debug_stuff * D)/*{{{*/
{
    std::vector<unsigned int> Lcols_pending;
    std::vector<unsigned int> pivs;
    size_t pos_q0 = 0;

#ifdef TIME_PLE
    ncalls++;
#endif
    TIMER_PLE(t_total);

    unsigned int ii = 0;
    for(unsigned int jj = 0 ; jj < n * B ; jj++) {
        if (ii >= m * B) break;
        unsigned int bi = ii / B;
        unsigned int i  = ii & (B-1);
        unsigned int bj = jj / B;
        unsigned int j  = jj & (B-1);

        int piv_ii = find_pivot(bi, bj, i, j);

        if (piv_ii >= 0) {
            unsigned int piv_bi = piv_ii / B;
            unsigned int piv_i  = piv_ii & (B-1);

            Lcols_pending.push_back(jj);

            /* Do the swap, locally.
             * In an MPI context, we would do the swap only on the
             * current block column. Here I'm not sure that it makes a
             * lot of sense to defer the swaps. Maybe a bit for locality.
             */
            pivs.push_back(piv_ii);

            if ((unsigned int) piv_ii != ii) {
#ifndef ACT_RIGHT_AWAY
                unsigned int s = bj;
#else
                for(unsigned s = 0 ; s < n ; s++)
#endif
                {
                    matrix & Y = X[bi * n + s];
                    matrix & piv_Y = X[piv_bi * n + s];
                    U c = Y[i] ^ piv_Y[piv_i];
                    Y[i]   ^= c;
                    piv_Y[piv_i] ^= c;
                }
            }

            /* A non-binary version would scale the multiplier column at
             * this point. Here we don't need to do that.
             */

            propagate_pivot(bi, bj, i, j);
        }

        /* Note that we have a column of rank-1 multipliers below
         * (ii,jj). But we would like to have it below (ii,ii) instead.
         * We'll keep it like this as long as we're still processing the
         * current block, and move all entries in one go when we reach
         * the end of the block.
         */

        bool finishing_block = (i == B-1) || (j == B-1);

        // for debugging. We don't do this by default because it slows
        // down debugging, and changes the course of the computation
        // anyway. But enabling this makes it possible to check the loop
        // invariant in simpler cases.
        // if (D) finishing_block = true;

        /* Time to increase ii and jj */
        ii += (piv_ii >= 0);

        ASSERT_ALWAYS(pos_q0 + Lcols_pending.size() == ii);

        if (!finishing_block) continue;
        if (Lcols_pending.empty()) {
            ASSERT_ALWAYS(pos_q0 == pivs.size());
            continue;
        }

#ifndef ACT_RIGHT_AWAY
        /* Note that we haven't increased bj yet, and that's on purpose */
        propagate_row_permutations(ii, bj, pivs.begin() + pos_q0, pivs.end());
#endif

        unsigned int yii0 = pos_q0;
        unsigned int yii1 = ii;
        unsigned int yi0 = yii0 & (B-1);
        unsigned int yi1 = yi0 + (yii1 - yii0);

        move_L_fragments(yii0, Lcols_pending);

        pos_q0 = pivs.size();
        Lcols_pending.clear();

        trsm(bi, bj, yi0, yi1);

        sub(bi, bj, yi0, yi1, ii);

        if (D) ASSERT_ALWAYS(D->check(X, pivs.begin(), ii));
    }
    ASSERT_ALWAYS(Lcols_pending.empty());
    return pivs;
}/*}}}*/

#ifdef TIME_PLE
template<typename matrix>
void PLE<matrix>::print_and_flush_stats()
{
    printf("PLE stats over %lu calls:\n", ncalls);
    printf("t_find_pivot: %g s (%.1f%%)\n", t_find_pivot / ncalls, 100.0 * t_find_pivot / t_total);
    printf("t_propagate_pivot: %g s (%.1f%%)\n", t_propagate_pivot / ncalls, 100.0 * t_propagate_pivot / t_total);
    printf("t_propagate_permutation: %g s (%.1f%%)\n", t_propagate_permutation / ncalls, 100.0 * t_propagate_permutation / t_total);
    printf("t_move_l_fragments: %g s (%.1f%%)\n", t_move_l_fragments / ncalls, 100.0 * t_move_l_fragments / t_total);
    printf("t_trsm: %g s (%.1f%%)\n", t_trsm / ncalls, 100.0 * t_trsm / t_total);
    printf("t_sub: %g s (%.1f%%)\n", t_sub / ncalls, 100.0 * t_sub / t_total);
    printf("t_total: %g s\n", t_total / ncalls);
    t_find_pivot = 0;
    t_propagate_pivot = 0;
    t_propagate_permutation = 0;
    t_move_l_fragments = 0;
    t_trsm = 0;
    t_sub = 0;
    t_total = 0;
    ncalls = 0;
}
#endif

#endif	/* LEVEL4_PLE_INTERNAL_INL_HPP_ */
