#include "cado.h" // IWYU pragma: keep

#include <cstdint>
#include <cstring>

#include "macros.h"
#include "blocklanczos_extraction.hpp"

/* given a 64 by 64 symmetric matrix A (given as 64 uint64_t's), and an input
 * bitmap S given in the form of one uint64_t, compute B and T so that
 * the following hold:
 *
 *    (i)    B = T * B = B * T
 *    (ii)   B * A * T = T
 *    (iii)  rank(B) = rank(A)
 *    (iv)   (1-S)*T = 1-S
 *
 * where S is identified with the diagonal matrix whose entries are given by S.
 *
 * In other words, this inverts (ii) a maximal minor (iii) in the matrix A,
 * selecting in priority (iv) for defining this minor the row indices
 * which are zero in S.
 *
 * The matrix A is clobbered.
 *
 * This routine does some allocation on the stack. Speed is not critical.
 */
uint64_t extraction_step(uint64_t * B, uint64_t * A, uint64_t S)
{
    int order[64], reorder[64];
    uint64_t B0[64];
    uint64_t T = 0;
    /* convert to a priority list, in a "save trees" style.  */
    for(int o=64,z=0,i=64;i-->0;) order[S&(uint64_t(1)<<i)?--o:z++]=i;
    for(int i = 0 ; i < 64 ; i++) B0[i] = uint64_t(1)<<i;
    for(int i = 0 ; i < 64 ; i++) reorder[i]=-1;
    for(int j = 0 ; j < 64 ; j++) {
        int const oj = order[j];
        uint64_t const mj = uint64_t(1)<<oj;
        int p = -1;
        for(int i = 0 ; i < 64 ; i++) {
            int const oi = order[i];
            uint64_t const mi = uint64_t(1)<<oi;
            if (T & mi) continue;
            if (A[oi] & mj) {
                p = i;
                break;
            }
        }
        if (p < 0) continue;

        int const op = order[p];
        /* Of course it's important to use indices op and oj here ! */
        reorder[op] = oj;
        uint64_t const mp = uint64_t(1) << op;
        /* We have a pivot, great. */
        ASSERT_ALWAYS(!(T & mp));
        T |= mp;
        /* add row op to all rows except itself */
        for(int i = 0 ; i < 64 ; i++) {
            if (i == p) continue;
            int const oi = order[i];
            uint64_t const x = ~-!(A[oi] & mj);
            B0[oi] ^= B0[op] & x;
            A[oi] ^= A[op] & x;
        }
    }
    /* Now at this point, we have some more work to do.
     *
     * A*T is now almost the identity matrix -- its square is diagonal
     * with zeros and ones, so A*T is just an involution. The reorder[]
     * array just has this involution.
     *
     * B is such that B*original_A = new_A.
     *
     * The matrix we want to return is new_A*T*b*T. We use reorder[] to
     * actually copy this into A.
     */
    memset(B, 0, sizeof(B0));
    for(int i = 0 ; i < 64 ; i++) {
        if (reorder[i] >= 0)
            B[reorder[i]]=B0[i]&T;
    }
    return T;
}


