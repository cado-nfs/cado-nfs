#ifndef CADO_BBLAS_LEVEL4_PLE_INTERNAL_INL_HPP
#define CADO_BBLAS_LEVEL4_PLE_INTERNAL_INL_HPP

// IWYU pragma: private, include "bblas.hpp"
// IWYU pragma: friend ".*/bblas.*"

#include "cado_config.h"

#include <cstdint>

#include <algorithm>

#include "bblas_level4_ple_internal.hpp"
#include "bblas_simd.hpp"
#include "bblas_mat64.hpp"
#include "bblas_mat8.hpp"
#ifdef TIME_PLE
#include "timing.h"
#endif

#ifdef TIME_PLE
struct timer_ple {
    double & t;
    timer_ple(double & t) : t(t) { t -= seconds(); }
    ~timer_ple() { t += seconds(); }
};
#define TIMER_PLE(t) timer_ple dummy(t)
#else
#define TIMER_PLE(t)    /**/
#endif

template<typename T>
PLE<T>::PLE(bpack_view<T> b) : bpack_view<T>(b), weights(std::vector<unsigned int>(b.nrows(), 0))
{
    prio_to_data.reserve(b.nrows());
    data_to_prio.reserve(b.nrows());
    for(unsigned int ii = 0 ; ii < b.nrows() ; ii++) {
        data_to_prio.push_back(ii);
        prio_to_data.push_back(ii);
    }
}

template<typename T>
PLE<T>::PLE(bpack_view<T> b, std::vector<unsigned int> d) : bpack_view<T>(b), weights(d)
{
    ASSERT_ALWAYS(d.size() == b.nrows());
    /* We need to create a priority list. We'll try to maximize the
     * number of rows whose priority order matches their position.
     */
    struct prio_cmp {
        std::vector<unsigned int> const & v;
        bool use_row_index;
        prio_cmp(std::vector<unsigned int> const & v, bool use_row_index)
            : v(v)
            , use_row_index(use_row_index)
        {}
        bool operator()(unsigned int a, unsigned int b) const {
            int r;
            r = (v[b] < v[a]) - (v[a] < v[b]);
            if (r || !use_row_index) return r < 0;
            return a < b;
        }
    };

    prio_to_data.reserve(d.size());
    data_to_prio.assign(d.size(), UINT_MAX);
    /* First sort all rows by weight */
    for(unsigned int ii = 0 ; ii < d.size() ; ii++)
        prio_to_data.push_back(ii);
    std::ranges::sort(prio_to_data, prio_cmp(weights, true));
    std::vector<unsigned int> sort_again;
    for(unsigned int ii = 0 ; ii < d.size() ; ii++) {
        if (weights[prio_to_data[ii]] == weights[ii]) {
            prio_to_data[ii] = ii;
        } else {
            sort_again.push_back(ii);
            prio_to_data[ii] = UINT_MAX;
        }
    }
    std::ranges::sort(sort_again, prio_cmp(weights, true));
    auto it = sort_again.begin();
    for(unsigned int ii = 0 ; ii < d.size() ; ii++) {
        if (prio_to_data[ii] == UINT_MAX)
            prio_to_data[ii] = *it++;
        data_to_prio[prio_to_data[ii]] = ii;
    }
    ASSERT_ALWAYS(std::ranges::is_sorted(prio_to_data, prio_cmp(weights, false)));
}

template<typename T>
int PLE<T>::find_pivot(unsigned int bi, unsigned int bj, unsigned int i, unsigned int j)/*{{{*/
{
    TIMER_PLE(t_find_pivot);
    U mask = U(1) << j;
    /* This implementation makes use of the fact that expected
     * discrepancy between the min and max weight is constant. Therefore
     * we strive to do as few permutations as we can.
     */
    std::vector<unsigned int> steps;
    unsigned int ii0 = bi * B + i;
    unsigned int w = weights[prio_to_data[ii0]];
    unsigned int ii = ii0;

    for( ; ii < mblocks * B ; ii++) {
        unsigned int zii = prio_to_data[ii];
        if (cell(zii / B, bj)[zii % B] & mask)
            break;
        if (weights[zii] != w) {
            steps.push_back(ii);
            w = weights[zii];
        }
    }
    if (ii == mblocks * B)
        return -1;
    
    /* We'll swap row zii and row ii0, for sure (this might entail a
     * change even if ii = ii0 !). The swap doesn't happen right now, but
     * the changes to the prio_to_data list do, se we do them ahead of
     * time.
     *
     * The post-condition is that
     *  row ii0 is a pivot
     *  prio_to_data[ii0] = data_to_prio[ii0] = ii0
     *  prio_to_data[ii] >= ii0 for all ii >= ii0
     *  data_to_prio[ii] >= ii0 for all ii >= ii0
     *  data_to_prio[prio_to_data[ii]] == ii for all ii >= ii0
     *  weights[prio_to_data[ii0+1:end]] is sorted 
     *
     */

    unsigned int ii1  = data_to_prio[ii0];
    unsigned int zii  = prio_to_data[ii];

    /* There are actually two operations.
     * 1. change the physical position of two rows (rows zii and ii0)
     */

    prio_to_data[ii]  = ii0;
    prio_to_data[ii1] = zii;
    data_to_prio[ii0] = ii;
    data_to_prio[zii] = ii1;
    std::swap(weights[ii0], weights[zii]);

    unsigned int zii0 = prio_to_data[ii0];
    /*
     * 2. change the priority position of the rows (indices ii and ii0)
     *
     * We do this in several steps, since it might be a cycle.
     */
    prio_to_data[ii0] = ii0;
    prio_to_data[ii] = UINT_MAX;
    data_to_prio[ii0] = ii0;
    data_to_prio[zii0] = UINT_MAX;

    /*
     * This leaves open the question of whether we resolve the situation
     * by prio_to_data[ii] = zii0 and data_to_prio[zii0] = ii ; it's 
     * slightly more complicated.
     *
     * When we entered this function, w=weights[zii0=prio_to_data[ii0]] was
     * the minimum of all weights. But if we schedule this row zii0 at
     * position ii in the priority list, it may well be too light, and
     * could have to go earlier.
     *
     * This is where we remember all the increasing values of the weight.
     * Each of them entails an ordering change in the priority list.
     */

    unsigned int pos_hole = ii;
    /* make sure now that weights[prio[ii0+1:end]] is sorted */
    for(auto it = steps.rbegin() ; it != steps.rend() ; ++it) {
        unsigned int ii2 = *it;
        /* position ii2 in the priority list points to a row that is
         * heavier than w. Hence we want to place this row at position ii
         * in the priority least, and use position ii2 as the next
         * preferred position.
         */
        data_to_prio[prio_to_data[ii2]] = pos_hole;
        prio_to_data[pos_hole] = prio_to_data[ii2];
        pos_hole = ii2;
    }
    data_to_prio[zii0] = pos_hole;
    prio_to_data[pos_hole] = zii0;
#ifndef NDEBUG
    // check all post-conditions.
    ASSERT(cell(zii / B, bj)[zii % B] & mask);
    ASSERT(prio_to_data[ii0] == ii0);
    ASSERT(data_to_prio[ii0] == ii0);
    w = 0;
    for(ii = ii0 + 1 ; ii < mblocks * B ; ii++) {
        ASSERT(prio_to_data[ii] > ii0);
        ASSERT(data_to_prio[ii] > ii0);
        ASSERT(prio_to_data[data_to_prio[ii]] == ii);
        ASSERT(weights[prio_to_data[ii]] >= w);
        w = weights[prio_to_data[ii]];
    }
#endif
    /* return this only so that the actual permutation takes place */
    return zii;
}/*}}}*/

template<typename T>
void PLE<T>::propagate_pivot(unsigned int bi, unsigned int bj, unsigned int i, unsigned int j)/*{{{*/
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
    ASSERT_ALWAYS(cell(bi, bj)[i] & mask);
    U c = cell(bi, bj)[i] & ((-mask) << 1);
    i++;
    if (i == B) { i = 0 ; bi++; }
    for( ; bi < mblocks ; bi++) {
        bitmat<T> & Y = cell(bi, bj);
        for( ; i < B ; i++) {
            Y[i] ^= c & -((Y[i] & mask) != 0);
        }
        i = 0;
    }
}/*}}}*/

/* These two specializations are very significant improvements.  */
#ifdef HAVE_AVX2
template<>
void PLE<uint64_t>::propagate_pivot(unsigned int bi, unsigned int bj, unsigned int i, unsigned int j)/*{{{*/
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
    ASSERT_ALWAYS(cell(bi, bj)[i] & mask);
    U c = cell(bi, bj)[i] & ((-mask) << 1);
    i++;
    constexpr const unsigned int SIMD = 4;
    for( ; (i & (SIMD-1)) && bi < mblocks ; i++) {
        matrix & Y = cell(bi, bj);
        Y[i] ^= c & -((Y[i] & mask) != 0);
    }
    if (i == B) { i = 0 ; bi++; }
    unsigned int iw = i / SIMD;
    __m256i mmask = _mm256_set1_epi64x(mask);
    __m256i cc = _mm256_set1_epi64x(c);
    for( ; bi < mblocks ; bi++) {
        matrix & Y = cell(bi, bj);
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
void PLE<uint64_t>::propagate_pivot(unsigned int bi, unsigned int bj, unsigned int i, unsigned int j)/*{{{*/
{
    TIMER_PLE(t_propagate_pivot);
    typedef uint64_t T;
    /* pivot row ii=bi*B+i to all rows below, but only for bits that are
     * right after column jj=bj*B+j.
     *
     *
     * This is done _only for the current column block bj !!!_
     */

    if (j == B-1) return;
    U mask = U(1) << j;
    ASSERT_ALWAYS(cell(bi, bj)[i] & mask);
    U c = cell(bi, bj)[i] & ((-mask) << 1);
    i++;
    constexpr const unsigned int SIMD = 2;
    for( ; (i & (SIMD-1)) && bi < mblocks ; i++) {
        bitmat<T> & Y = cell(bi, bj);
        Y[i] ^= c & -((Y[i] & mask) != 0);
    }
    if (i == B) { i = 0 ; bi++; }
    unsigned int iw = i / SIMD;
    __m128i mmask = _cado_mm_set1_epi64(mask);
    __m128i cc = _cado_mm_set1_epi64(c);
    for( ; bi < mblocks ; bi++) {
        bitmat<T> & Y = cell(bi, bj);
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

#ifdef HAVE_MMX
template<>
void PLE<uint8_t>::propagate_pivot(unsigned int bi, unsigned int bj, unsigned int i, unsigned int j)/*{{{*/
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
    ASSERT_ALWAYS(cell(bi, bj)[i] & mask);
    U c = cell(bi, bj)[i] & ((-mask) << 1);
    i++;
    if (i == B) { i = 0 ; bi++; }
    __m64 mmask = _mm_set1_pi8(mask);
    __m64 cc = _mm_set1_pi8(c);
    if (i && bi < mblocks) {
        matrix & Y = cell(bi, bj);
        __m64 & Yw = * (__m64 *) Y.data();
        __m64 select = _mm_and_si64(Yw, mmask);
        select = _mm_cmpeq_pi8(select, mmask);
        /* Make sure we don't touch rows before i ! */
        /* Note that the deal about _mm_cvtsi64_m64 versus
         * _mm_cvtsi64x_m64 is ultimately very boring. The former is
         * exposed by gcc only on x86_64. The latter is a microsoft
         * intrinsic. Both boil down to a simple cast. Let's just do the
         * cast ourselves, and be done with it.
         */
        select = _mm_and_si64(select, (__m64) (-(uint64_t(1) << (8*i))));
        Yw = _mm_xor_si64(Yw, _mm_and_si64(cc, select));
        bi++;
    }
    for( ; bi < mblocks ; bi++) {
        matrix & Y = cell(bi, bj);
        __m64 & Yw = * (__m64 *) Y.data();
        __m64 select = _mm_and_si64(Yw, mmask);
        select = _mm_cmpeq_pi8(select, mmask);
        Yw = _mm_xor_si64(Yw, _mm_and_si64(cc, select));
    }
}/*}}}*/
#endif

template<typename T>
void PLE<T>::propagate_row_permutations(unsigned int ii1, unsigned int bj0, std::vector<unsigned int>::const_iterator q0, std::vector<unsigned int>::const_iterator q1)/*{{{*/
{
    TIMER_PLE(t_propagate_permutation);
    /* This propagates the pending permutations outside the current block
     * column.
     * Permutations are given by the range [q0..q1). More precisely,
     * the starting index is given by ii0 = ii1 - (q1 - q0), and the
     * image of index ii is given by
     *          *(q0 + ii - ii0) = *(q1 + ii - ii1)
     */
    for(unsigned bj = 0 ; bj < nblocks ; bj++) {
        if (bj == bj0) continue;
        unsigned int ii = ii1 - (q1 - q0);
        for(auto q = q0 ; q != q1 ; q++, ii++) {
            if (ii == *q) continue;
            unsigned int bi = ii / B;
            unsigned int i  = ii & (B - 1);
            unsigned int pbi = *q / B;
            unsigned int pi  = *q & (B - 1);
            bitmat<T> & Y  = cell(bi, bj);
            bitmat<T> & pY = cell(pbi, bj);
            U c = Y[i] ^ pY[pi];
            Y[i]   ^= c;
            pY[pi] ^= c;
        }
    }
}/*}}}*/

template<typename T>
void PLE<T>::move_L_fragments(unsigned int yii0, std::vector<unsigned int> const & Q)/*{{{*/
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
     * bitmat<T> that we expect to find eventually.
     * Note 1: the multiplier column is _under_ the leading coordinate
     * (yyi0+x, Q[x]): only coordinates (yyi0+x+1, Q[x]) and downwards
     * are moved).
     * Note 2: by construction, (yii0+x)/B is constant as x runs through
     * [0..Q.size()-1], and so is Q[x]/B.
     *
     * This function assumes that the bitmat<T> coefficients that receive
     * coefficients from moved columns are set to zero beforehand, and
     * ensures that bitmat<T> entries from moved columns are set to zero
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
            U c = (cell(ybi, zbj)[yi + warmup] >> zj) & mask;
            cell(ybi, ybi)[yi + warmup] ^= (c << yi);
            cell(ybi, zbj)[yi + warmup] ^= (c << zj);
            /* is there a simpler way ? mask += mask + 1 ? */
            mask |= mask << 1;
        }
        /* we now move the part that comes _after_ the leading unit lower
         * triangular block.
         */
        unsigned int tbi = ybi;
        unsigned int ti = yi + dz;
        if (ti == B) { tbi++; ti = 0; }
        for( ; tbi < mblocks ; tbi++) {
            for( ; ti < B ; ti++) {
                U c = (cell(tbi, zbj)[ti] >> zj) & mask;
                cell(tbi, ybi)[ti] ^= (c << yi);
                cell(tbi, zbj)[ti] ^= (c << zj);
            }
            ti = 0;
        }
    }
}/*}}}*/

#ifdef HAVE_AVX2
template<>
void PLE<uint64_t>::move_L_fragments(unsigned int yii0, std::vector<unsigned int> const & Q)/*{{{*/
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
            U c = (cell(ybi, zbj)[yi + warmup] >> zj) & mask;
            cell(ybi, ybi)[yi + warmup] ^= (c << yi);
            cell(ybi, zbj)[yi + warmup] ^= (c << zj);
            /* is there a simpler way ? mask += mask + 1 ? */
            mask |= mask << 1;
        }
        unsigned int tbi = ybi;
        unsigned int ti = yi + dz;
        constexpr const unsigned int SIMD = 4;
        for ( ; (ti & (SIMD - 1)) && tbi < mblocks ; ti++) {
            U c = (cell(tbi, zbj)[ti] >> zj) & mask;
            cell(tbi, ybi)[ti] ^= (c << yi);
            cell(tbi, zbj)[ti] ^= (c << zj);
        }
        if (ti == B) { tbi++; ti = 0; }
        __m256i mmask = _mm256_set1_epi64x(mask);
        unsigned int tiw = ti / SIMD;
        __m128i zjw = _mm_set1_epi64x(zj); // only lower 64 used
        __m128i yiw = _mm_set1_epi64x(yi); // only lower 64 used
        for( ; tbi < mblocks ; tbi++) {
            __m256i *Xzjw = (__m256i *) cell(tbi, zbj).data();
            __m256i *Xyiw = (__m256i *) cell(tbi, ybi).data();
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
void PLE<uint64_t>::move_L_fragments(unsigned int yii0, std::vector<unsigned int> const & Q)/*{{{*/
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
            U c = (cell(ybi, zbj)[yi + warmup] >> zj) & mask;
            cell(ybi, ybi)[yi + warmup] ^= (c << yi);
            cell(ybi, zbj)[yi + warmup] ^= (c << zj);
            /* is there a simpler way ? mask += mask + 1 ? */
            mask |= mask << 1;
        }
        unsigned int tbi = ybi;
        unsigned int ti = yi + dz;
        constexpr const unsigned int SIMD = 2;
        for ( ; (ti & (SIMD - 1)) && tbi < mblocks ; ti++) {
            U c = (cell(tbi, zbj)[ti] >> zj) & mask;
            cell(tbi, ybi)[ti] ^= (c << yi);
            cell(tbi, zbj)[ti] ^= (c << zj);
        }
        if (ti == B) { tbi++; ti = 0; }
        __m128i mmask = _cado_mm_set1_epi64(mask);
        unsigned int tiw = ti / SIMD;
        __m128i zjw = _cado_mm_set1_epi64(zj); // only lower 64 used
        __m128i yiw = _cado_mm_set1_epi64(yi); // only lower 64 used
        for( ; tbi < mblocks ; tbi++) {
            __m128i *Xzjw = (__m128i *) cell(tbi, zbj).data();
            __m128i *Xyiw = (__m128i *) cell(tbi, ybi).data();
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

template<typename T>
void PLE<T>::trsm(unsigned int bi,/*{{{*/
        unsigned int bj,
        unsigned int yi0,
        unsigned int yi1)
{
    TIMER_PLE(t_trsm);
    /* trsm is fairly trivial */
    for(unsigned int s = bj + 1 ; s < nblocks ; s++) {
        bitmat<T>::trsm(cell(bi, bi), cell(bi, s), yi0, yi1);
    }
}/*}}}*/

template<typename T>
void PLE<T>::sub(unsigned int bi,/*{{{*/
        unsigned int bj,
        unsigned int yi0,
        unsigned int yi1,
        unsigned int ii)
{
    TIMER_PLE(t_sub);
    unsigned int sbi = ii / B;
    unsigned int si  = ii & (B-1);
    for( ; sbi < mblocks ; sbi++) {
        // i
        for(unsigned int sbj = bj + 1 ; sbj < nblocks ; sbj++) {
            /* multiply columns [yi0..yi1-1] of block (sbi,bi) by
             * rows [yi0..yi1-1] of block (bi,sbj), and add that to
             * block (sbi,sbj)
             *
             * yes, we really mean block (sbi,bi). The indices [yi0..yi1)
             * are _really_ relative to that block.
             */
            bitmat<T>::addmul(cell(sbi, sbj),
                    cell(sbi, bi),
                    cell(bi, sbj),
                    si, B, yi0, yi1);
        }
        si = 0;
    }
}/*}}}*/

template<typename T>
void PLE<T>::debug_stuff::apply_permutations(std::vector<unsigned int>::const_iterator p0, std::vector<unsigned int>::const_iterator p1)/*{{{*/
{
    /* apply the permutations to Xcc */
    for(unsigned int xii = 0 ; xii < (unsigned int) (p1 - p0) ; xii++) {
        unsigned int pii = p0[xii];
        if (xii == pii) continue;
        unsigned int xbi = xii / B;
        unsigned int xi = xii % B;
        unsigned int pbi = pii / B;
        unsigned int pi = pii % B;
        for(unsigned int bj = 0 ; bj < nblocks ; bj++) {
            bitmat<T> & xY = target.cell(xbi, bj);
            bitmat<T> & pY = target.cell(pbi, bj);
            U c = xY[xi] ^ pY[pi];
            xY[xi] ^= c;
            pY[pi] ^= c;
        }
    }
}/*}}}*/

/* extract below the diagonal, only up to rank rr. The X field is
 * modified. */
template<typename T>
bpack<T> PLE<T>::debug_stuff::get_LL(unsigned int rr)/*{{{*/
{
    bpack<T> LL(nrows(), nrows());
    for(unsigned int bi = 0 ; bi < mblocks ; bi++) {
        unsigned int bj = 0;
        for( ; bj <= bi && bj < iceildiv(rr, B) ; bj++) {
            for(unsigned int i = 0 ; i < B ; i++) {
                U c = cell(bi, bj)[i];
                unsigned int z = rr - bj * B;
                if (bj == bi && i < z)
                    z = i;
                if (z < B)
                    c &= (U(1) << z) - 1;
                LL.cell(bi, bj)[i] = c;
            }
        }
    }
    /* clear the blocks that we have just taken */
    for(unsigned int bi = 0 ; bi < mblocks ; bi++) {
        for(unsigned int bj = 0 ; bj < nblocks && bj < mblocks ; bj++) {
            bitmat<T>::add(cell(bi, bj), cell(bi, bj), LL.cell(bi, bj));
        }
    }

    /* add implicit identity to L */
    for(unsigned int bi = 0 ; bi < mblocks ; bi++) {
        for(unsigned int i = 0 ; i < B ; i++) {
            LL.cell(bi, bi)[i] ^= U(1) << i;
        }
    }

    return LL;
}/*}}}*/

template<typename T>
bpack<T> PLE<T>::debug_stuff::get_UU(unsigned int rr)/*{{{*/
{
    /* extract above the diagonal, only up to rank rr */
    bpack<T> UU(nrows(), ncols());
    for(unsigned int bi = 0 ; bi < mblocks && bi < iceildiv(rr, B); bi++) {
        if (bi < nblocks) {
            for(unsigned int i = 0 ; i < std::min(B, rr - bi * B) ; i++) {
                U c = cell(bi, bi)[i];
                c &= -(U(1) << i);
                UU.cell(bi, bi)[i] = c;
            }
        }
        for(unsigned int bj = bi + 1 ; bj < nblocks ; bj++) {
            UU.cell(bi, bj) = cell(bi, bj);
            for(unsigned int i = rr - bi * B ; i < B ; i++) {
                UU.cell(bi, bj)[i] = 0;
            }
        }
    } 
    /* Finally, clear everything in Xcc that we haven't taken yet *//*{{{*/
    for(unsigned int bi = 0 ; bi < mblocks ; bi++) {
        for(unsigned int bj = 0 ; bj < nblocks ; bj++) {
            bitmat<T>::add(cell(bi, bj), cell(bi, bj), UU.cell(bi, bj));
        }
    }/*}}}*/
    return UU;
}/*}}}*/

template<typename T>
bool PLE<T>::debug_stuff::complete_check(bpack<T> const & LL, bpack<T> const & UU)/*{{{*/
{
    /* check that LL*UU + (remaining block in X) is equal to X_target */

    for(unsigned int bi = 0 ; bi < mblocks ; bi++) {
        for(unsigned int bj = 0 ; bj < nblocks ; bj++) {
            bitmat<T> C = cell(bi, bj);
            for(unsigned int bk = 0 ; bk < mblocks ; bk++) {
                bitmat<T>::addmul(C, LL.cell(bi, bk), UU.cell(bk, bj));
            }
            ASSERT_ALWAYS(target.cell(bi, bj) == C);
            if (target.cell(bi, bj) != C) return false;
        }
    }
    return true;
}/*}}}*/

template<typename T>
std::vector<unsigned int> PLE<T>::operator()(debug_stuff * D)/*{{{*/
{
    std::vector<unsigned int> Lcols_pending;
    std::vector<unsigned int> pivs;
    size_t pos_q0 = 0;

#ifdef TIME_PLE
    ncalls++;
#endif
    TIMER_PLE(t_total);

    unsigned int ii = 0;
    for(unsigned int jj = 0 ; jj < nblocks * B ; jj++) {
        if (ii >= mblocks * B) break;
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
                for(unsigned s = 0 ; s < nblocks ; s++)
#endif
                {
                    bitmat<T> & Y = cell(bi, s);
                    bitmat<T> & piv_Y = cell(piv_bi, s);
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

        if (D) ASSERT_ALWAYS(D->check(const_view(), pivs.begin(), ii));
    }
    ASSERT_ALWAYS(Lcols_pending.empty());
    return pivs;
}/*}}}*/

#ifdef TIME_PLE
template<typename T>
void PLE<T>::print_and_flush_stats()
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

#endif	/* BBLAS_LEVEL4_PLE_INTERNAL_INL_HPP_ */
