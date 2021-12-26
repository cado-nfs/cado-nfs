#include "cado.h" // IWYU pragma: keep
#include <cstdint>                      // for uint64_t, UINT64_C
#include "bblas_mat64.hpp"  // for mat64
#include "bblas_level4.hpp"

/* Computes e,mm such that mm=e*m is in row echelon form */
int
full_echelon_6464_imm(mat64& mm, mat64& e, mat64 const& m)
{
    mm = m;
    uint64_t mask = 1;
    // uint64_t cancelled_cols = 0;
    int r = 0;
    for (int j = 0; j < 64; j++, mask <<= 1)
        e[j] = mask;
    mask = 1;
    for (int j = 0; j < 64; j++, mask <<= 1) {
        int k = 0;
        uint64_t z = UINT64_C(1);
        uint64_t pr;
        for (k = 0; z; k++, z <<= 1) {
            pr = mm[k];
            if ((pr & mask) && !(pr & (mask - 1)))
                break;
        }
        if (!z)
            continue;
        z = mm[k];
        mm[k] = mm[j];
        mm[j] = z;
        z = e[k];
        e[k] = e[j];
        e[j] = z;
        r++;
        // cancelled_cols |= mask;
        uint64_t er = e[j];
        for (k = 0; k < 64; k++) {
            if (k == j)
                continue;
            uint64_t w = -((mm[k] & mask) != 0);
            mm[k] ^= pr & w;
            e[k] ^= er & w;
        }
    }
    return r;
}
