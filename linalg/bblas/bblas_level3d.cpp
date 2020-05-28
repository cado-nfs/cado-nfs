#include "cado.h" // IWYU pragma: keep
#include <cstdint>         // for uint64_t
#include <algorithm>        // for min
#include "bblas_mat64.hpp"
#include "bblas_level3d.hpp"
#include "macros.h"         // for ASSERT

/**********************************************************************/
/* level 3d: solution of triangular linear systems
 *      trsm64
 *      trsm64_general
 */

/* implemented here:
 *    - trsm64: gen two 64*64 bit matrices L and U, with L unit lower
 *      triangular, replace U by L^-1*U.
 *    - trsm64_general: same, but apply only the square submatrix of L
 *      whose diagonal indices are in the given integer interval.
 *
 * Note that these functions do _not_ check that L is indeed unit lower
 * triangular. The coefficients at positions j>=i are simply not looked
 * up.
 */

void trsm64_general(mat64 const & L, mat64 & U, unsigned int n0, unsigned int n1)/*{{{*/
{
    // ASSERT(mat64_is_lowertriangular(L));
    // ASSERT(mat64_triangular_is_unit(L));
    ASSERT(n0 <= n1);
    if (n1 <= n0 + 1) return;
    /* need to determine the very first fragment before we can align */
    if (n0 % 4) {
        uint64_t c[8];
        unsigned int n0b = std::min(n0 + 4 - (n0 % 4), n1);
        uint64_t * uu = U.data() + n0;
        uint64_t const * ll= L.data() + n0;
        unsigned int d = n0b - n0;
        c[0] = 0;
        c[1] = uu[0];
        uint64_t m = 1;
        if (d >= 2) {
            c[2]=uu[1]^=c[(ll[1]>>n0)&1];
            c[3]=uu[1]^uu[0];
            m=3;
        }
        if (d >= 3) {
            c[4]=uu[2]^=c[(ll[2]>>n0)&3];
            c[5]=uu[2]^c[1];
            c[6]=uu[2]^c[2];
            c[7]=uu[2]^c[3];
            m=7;
        }
        for(unsigned int i = n0b ; i < n1 ; i++)
            U[i] ^= c[(L[i] >> n0)&m];
        n0 = n0b;
    }
    if (n1 <= n0 + 1) return;
    unsigned int b = n0;
    for(b = n0 ; b < n1 - (n1 % 4) ; b += 4) {
        uint64_t c[16];
        uint64_t * uu = U.data() + b;
        uint64_t const * ll = L.data() + b;
        c[0]=0;
        c[1]=uu[0];
        c[2]=uu[1]^=c[(ll[1]>>b)&1];
        c[3]=uu[1]^uu[0];
        c[4]=uu[2]^=c[(ll[2]>>b)&3];
        c[5]=uu[2]^c[1];
        c[6]=uu[2]^c[2];
        c[7]=uu[2]^c[3];
        c[8]=uu[3]^=c[(ll[3]>>b)&7];
        /* we might break here if b == 60 */
        c[9]=uu[3]^c[1];
        c[10]=uu[3]^c[2];
        c[11]=uu[3]^c[3];
        c[12]=uu[3]^c[4];
        c[13]=uu[3]^c[5];
        c[14]=uu[3]^c[6];
        c[15]=uu[3]^c[7];
        for(unsigned int i = b + 4 ; i < n1 ; i++)
            U[i] ^= c[(L[i] >> b)&15];
    }
    if (n1 % 4) {
        ASSERT(b < n1);
        unsigned int d = n1 % 4;
        ASSERT(n1 == b + d);
        if (d <= 1) return;
        uint64_t c[4];
        uint64_t * uu = U.data() + b;
        uint64_t const * ll = L.data() + b;
        c[0]=0;
        c[1]=uu[0];
        if (d >= 2) {
            c[2]=uu[1]^=c[(ll[1]>>b)&1];
            c[3]=uu[1]^uu[0];
        }
        if (d >= 3) {
                 uu[2]^=c[(ll[2]>>b)&3];
        }
    }
}/*}}}*/

void trsm64(mat64 const & L, mat64 & U)/*{{{*/
{
    trsm64_general(L, U, 0, 64);
}/*}}}*/

