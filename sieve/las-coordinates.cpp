#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdint>

#include <gmp.h>

#include "fb-types.hpp"
#include "gmp_aux.h"
#include "las-coordinates.hpp"
#include "las-config.hpp"
#include "las-qlattice.hpp"
#include "special-q.hpp"
#include "macros.h"


/*  Conversions between different representations for sieve locations:
 *
 * Coordinate systems for identifying an (a,b) pair.
 *
 * (a, b)  [short name "ab"] This is the one from textbooks. We always
 *         have b>=0 (only free relations have b=0, and we don't see them
 *         here). For siqs, we only have b=1.
 * (i, j)  [short name "ij"] las case: for a given special q, this is a point in
 *         the q-lattice. Given the lattice basis given by (a0 b0 a1 b1),
 *         this corresponds to the (a,b) pair equal to
 *         i*(a0,b0)+j*(a1,b1). By construction this should lead to one
 *         of the norms having doing.p as a factor.  i is within [-I/2,
 *         I/2[, and j is within [1, J[
 *         siqs case: (i, j) corresponds to the (a,b) pair equal to (rj+i*q, 1).
 *         For more detail on how rj is computed from j, see the method
 *         siqs_special_q_data::root_from_j.
 * (N, x)  [short name "Nx"] bucket number N, location x. N is within
 *         [0,nb_buckets[ and x within [0,bucket_region[ ; we have:
 *         N*bucket_region+x == (I/2+i)+j*I
 * X       [short name "X"]  N*bucket_region+x
 *
 * There are some change of coordinate functions in this file, and also
 * in the header file las-coordinates.hpp
 *
 * The naming scheme for the conversion functions is convert_FOO_to_BAR.
 * FOO or BAR are the "short names" mentioned above for the different
 * coordinate systems.
 *
 * In the single-coordinate system, X can exceed 32 bits (only if I>16,
 * in fact), while in the coordinate system (N,x), the coordinate x is
 * less than bucket_region, which is at most 2^16.
 */

void convert_X_to_ij(int & i, unsigned int & j, const uint64_t X, int logI)
{
    i = (X & ((1 << logI)-1)) - (1 << (logI - 1));
    j = X >> logI;
}

void convert_Nx_to_ij(int & i, unsigned int & j, const unsigned int N, const unsigned int x, int logI)
{
    convert_X_to_ij(i, j, convert_Nx_to_X(N, x), logI);
}

void convert_ij_to_X(uint64_t & x, int i, unsigned int j, int logI)
{
    x = (int64_t)i + (((uint64_t)j) << logI) + (uint64_t)(1 << (logI - 1));
}

void convert_ij_to_Nx(unsigned int & N, unsigned int & x, int i, unsigned int j, int logI)
{
    uint64_t xx;
    convert_ij_to_X(xx, i, j, logI);
    N = xx >> LOG_BUCKET_REGION;
    /* see CID 1453612; I _think_ that this assert should be enough */
    ASSERT_FOR_STATIC_ANALYZER(LOG_BUCKET_REGION <= 32);
    x = xx & (uint64_t)((1 << LOG_BUCKET_REGION) - 1);
}

#if 0 /* currently unused */
int convert_ab_to_X(unsigned int *x, const int64_t a, const uint64_t b, int logI, qlattice_basis const & Q)
{
    int i;
    unsigned int j;
    if (!Q.convert_ab_to_ij(i, j, a, b)) return 0;
    convert_ij_to_X(x, a, b, logI);
    return 1;
}
#endif

#if 0 /* currently unused */
int convert_ab_to_Nx(unsigned int * N, unsigned int *x, const int64_t a, const uint64_t b, int logI, qlattice_basis const & Q)
{
    int i;
    unsigned int j;
    if (!Q.convert_ab_to_ij(i, j, a, b)) return 0;
    convert_ij_to_Nx(N, x, a, b, logI);
    return 1;
}
#endif

