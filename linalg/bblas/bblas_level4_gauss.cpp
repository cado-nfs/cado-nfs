#include "cado.h" // IWYU pragma: keep

#include <cstdint>
#include <cstring>

#include <algorithm>
#include <array>

#include <gmp.h>

#include "bblas_gauss.h"
#include "bblas_level4.hpp"
#include "bblas_mat64.hpp" // for mat64

int gauss_6464_C(mat64 & mm, mat64 & e, mat64 const & m)
{
    mm = m;
    uint64_t * ee[64];
    for (int j = 0; j < 64; j++)
        ee[j] = &(e[j]);
    int const r = kernel((mp_limb_t *)mm.data(), (mp_limb_t **)ee, 64, 64,
                         64 / ULONG_BITS, 64 / ULONG_BITS);
    return r;
}

int gauss_6464_imm(mat64 & mm, mat64 & e, mat64 const & m)
{
    mm = m;
    uint64_t mask = 1;
    uint64_t taken = 0;
    // uint64_t cancelled_cols=0;
    int r = 0;
    for (int j = 0; j < 64; j++, mask <<= 1)
        e[j] = mask;
    mask = 1;
    for (int j = 0; j < 64; j++, mask <<= 1) {
        int k = 0;
        uint64_t z = 1;
        uint64_t pr;
        for (k = 0; z; k++, z <<= 1) {
            if (taken & z)
                continue;
            pr = mm[k];
            if (pr & mask)
                break;
        }

        if (!z)
            continue;
        taken |= z;
        r++;
        // cancelled_cols|=mask;
        uint64_t const er = e[k];
        for (k++; k < 64; k++) {
            uint64_t const w = -((mm[k] & mask) != 0);
            mm[k] ^= pr & w;
            e[k] ^= er & w;
        }
    }
    return r;
}

int gauss_128128_C(mat64 * m)
{
    std::array<mat64, 4> mm;
    std::ranges::copy(m, m + 4, std::begin(mm));
    int const r = kernel((mp_limb_t *) mm.front().data(),
            nullptr, 128, 128, 128 / ULONG_BITS, 128 / ULONG_BITS);
    return r;
}

#if 0
int gauss_128128_imm(uint64_t * m)
{
    mat64 mm[4] ATTRIBUTE((aligned(64)));
    uint64_t * pm = m;
    for(int j = 0 ; j < 64 ; j++, pm+=2) {
        mm[0][j] = pm[0];
        mm[1][j] = pm[1];
    }
    for(int j = 0 ; j < 64 ; j++, pm+=2) {
        mm[2][j] = pm[0];
        mm[3][j] = pm[1];
    }


    mat64 e;
    memcpy(mm,m,sizeof(mat64));
    uint64_t mask=1;
    uint64_t taken=0;
    uint64_t cancelled_cols=0;
    int r = 0;
    for(int j = 0 ; j < 64 ; j++, mask<<=1) e[j]=mask;
    mask=1;
    for(int j = 0 ; j < 64 ; j++, mask<<=1) {
        int k = 0;
        uint64_t z = UINT64_C(1);
        uint64_t pr;
        for(k = 0 ; z && !(((pr=mm[k])&mask) && !(taken&z)); k++, z<<=1) ;
        if (!z) continue;
        taken|=z;
        r++;
        cancelled_cols|=mask;
        uint64_t er = e[k];
        int k0=k;
#define TRIANGULAR_ONLY /* speeds up things by 20 to 25% */
#ifndef TRIANGULAR_ONLY
        k = 0;
#endif
        for( ; k < 64 ; k++) {
            if (k==k0) continue;
            uint64_t w = -((mm[k]&mask)!=0);
            mm[k]^=pr&w;
            e[k]^=er&w;
        }
    }
    return r;
}
#endif
