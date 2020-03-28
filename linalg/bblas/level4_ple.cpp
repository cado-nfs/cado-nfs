#include "cado.h"
#include "level4.hpp"
#include "level4_ple_internal.hpp"

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

int binary_blas_PLE(unsigned int * p, mat64 * X, unsigned int m, unsigned int n);

int PLE::find_pivot(unsigned int bi, unsigned int bj, unsigned int i, unsigned int j) const/*{{{*/
{
    const unsigned int B = 64;
    mat64 & Y = X[bi * n + bj];
    uint64_t mask = UINT64_C(1) << j;
    for( ; bi < m ; bi++) {
        for( ; i < B ; i++) {
            if (Y[i] & mask)
                return bi * B + i;
        }
        i = 0;
    }
    return -1;
}/*}}}*/

void PLE::propagate_pivot(unsigned int bi, unsigned int bj, unsigned int i, unsigned int j) const/*{{{*/
{
    /* row ii to all rows below, but only for bits that are
     * right after the pivot */

    const unsigned int B = 64;
    mat64 & Y = X[bi * n + bj];
    if (j == B-1) return;
    uint64_t mask = UINT64_C(1) << j;
    uint64_t c = Y[i] & ((-mask) << 1);
    i++;
    if (i == B) { i = 0 ; bi++; }
    for( ; bi < m ; bi++) {
        for( ; i < B ; i++) {
            Y[i] ^= c & -((Y[i] & mask) != 0);
        }
    }
}/*}}}*/

void PLE::propagate_permutations(unsigned int ii0, unsigned int bj0, unsigned int const * q0, unsigned int const * q1) const/*{{{*/
{
    const unsigned int B = 64;
    for(unsigned bj = 0 ; bj < n ; bj++) {
        if (bj == bj0) continue;
        unsigned int ii = ii0 - (q1 - q0);
        for(unsigned int const * q = q0 ; q != q1 ; q++, ii++) {
            if (ii == *q0) continue;
            unsigned int bi = ii / B;
            unsigned int i  = ii & (B - 1);
            unsigned int pbi = *q / B;
            unsigned int pi  = *q & (B - 1);
            mat64 & Y  = X[ bi * n + bj];
            mat64 & pY = X[pbi * n + bj];
            uint64_t c = Y[i] ^ pY[pi];
            Y[i]   ^= c;
            pY[pi] ^= c;
        }
    }
}/*}}}*/

void PLE::move_L_fragments(unsigned int yii0, std::vector<unsigned int> const & Q) const/*{{{*/
{
    /* This function receives the (yii0,yii0) coordinate of the first
     * item in the unit lower triangular part that we haven't completed
     * yet (or that we aren't sure that we have completed -- maybe we
     * have because it's already in place).
     * This function is called just before
     * a block is finished (either row block or column block). Either way,
     * we want to put the multipliers in the right positions. All multiplier
     * columns referenced by the list Q, with logical leading
     * coordinates (yyi0+x, Q[x]) must move to column
     * yyi0+x so that thaye are part of the unit lower triangular
     * matrix that we expect to find eventually.
     * Note that by construction, (yii0+x)/B is constant as x runs
     * through [0..Q.size()-1], and so is Q[x]/B.
     */
    const unsigned int B = 64;
    unsigned int ybi = yii0 / B;
    unsigned int yi  = yii0 & (B-1);
    unsigned int k = Q.size();
    unsigned int yii = yii0;
    ASSERT(!Q.size());
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
        uint64_t mask = 1;
        for(unsigned int warmup = 1 ; warmup < dz ; warmup++) {
            /* we have yi + warmup = yii0 + dy + warmup
             *                     < yii0 + dy + dz
             *                     <= yii0 + dy + dz - 1
             *                     <= yii0 + k - 1
             *                     <= yii1 - 1
             * which on all accounts should still be in the same block.
             */
            ASSERT((yi + warmup) / B == ybi);
            uint64_t c = (X[ybi * n + zbj][yi + warmup] >> zj) & mask;
            X[ybi * n + ybi][yi + warmup] ^= (c << yi);
            X[ybi * n + zbj][yi + warmup] ^= (c << zj);
            /* is there a simpler way ? mask += mask + 1 ? */
            mask |= mask << 1;
        }
        unsigned int tbi = ybi;
        unsigned int ti = yi + dz;
        if (ti == B) { tbi++; ti = 0; }
        for( ; tbi < m ; tbi++) {
            for( ; ti < B ; ti++) {
                uint64_t c = (X[tbi * n + zbj][ti] >> zj) & mask;
                X[tbi * n + tbi][ti] ^= (c << ti);
                X[tbi * n + zbj][ti] ^= (c << zj);
            }
        }
    }
}/*}}}*/

void PLE::trsm(unsigned int bi,/*{{{*/
        unsigned int bj,
        unsigned int yi0,
        unsigned int yi1) const
{
    /* trsm is fairly trivial */
    for(unsigned int s = bj + 1 ; s < n ; s++) {
        trsm64_general(X[bi * n + bi], X[bi * n + s], yi0, yi1);
    }
}/*}}}*/

void PLE::sub(unsigned int bi,/*{{{*/
        unsigned int bj,
        unsigned int yi0,
        unsigned int yi1,
        unsigned int ii) const
{
    const unsigned int B = 64;
    unsigned int sbi = ii / B;
    unsigned int si  = ii & (B-1);
    for( ; sbi < m ; sbi++) {
        // i
        for(unsigned int sbj = bj + 1 ; sbj < n ; sbj++) {
            /* multiply columns [yi0..yi1-1] of block (sbi,bj) by
             * rows [yi0..yi1-1] of block (bi,sbj), and add that to
             * block (sbi,sbj)
             */
            addmul_6464_6464_fragment_lookup4(X[sbi * n + sbj],
                    X[sbi * n + bj],
                    X[bi * n + sbj],
                    si, 64, yi0, yi1);
        }
        si = 0;
    }
}/*}}}*/

int PLE::operator()(unsigned int * p0)/*{{{*/
{
    std::vector<unsigned int> Lcols_pending;
    unsigned int * p = p0;
    unsigned int const * q0 = p;
    const unsigned int B = 64;

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
            *p++ = piv_ii;

            if ((unsigned int) piv_ii != ii) {
#ifndef ACT_RIGHT_AWAY
                unsigned int s = bj;
#else
                for(unsigned s = 0 ; s < n ; s++)
#endif
                {
                    mat64 & Y = X[bi * n + s];
                    mat64 & piv_Y = X[piv_bi * n + s];
                    uint64_t c = Y[i] ^ piv_Y[piv_i];
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

        /* Time to increase ii and jj */
        ii += (piv_ii >= 0);

        ASSERT_ALWAYS((q0 - p) + Lcols_pending.size() == ii);

        if (!finishing_block) continue;

#ifndef ACT_RIGHT_AWAY
        /* Note that we haven't increased bj yet, and that's on purpose */
        propagate_permutations(ii, bj, q0, p);
#endif

        unsigned int yii0 = q0 - p;
        unsigned int yii1 = ii;
        unsigned int yi0 = yii0 & (B-1);
        unsigned int yi1 = yii1 & (B-1);

        move_L_fragments(yii0, Lcols_pending);

        q0 = p;
        Lcols_pending.clear();

        trsm(bi, bj, yi0, yi1);

        sub(bi, bj, yi0, yi1, ii);
    }
    return p - p0;
}/*}}}*/

int binary_blas_PLE(unsigned int * p, mat64 * X, unsigned int m, unsigned int n)
{
    return PLE(X, m, n)(p);
}

