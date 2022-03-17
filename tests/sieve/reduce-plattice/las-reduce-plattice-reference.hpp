#ifndef LAS_REDUCE_PLATTICE_REFERENCE_HPP_
#define LAS_REDUCE_PLATTICE_REFERENCE_HPP_

#include "las-plattice.hpp"
#include "reduce-plattice/plattice-proxy.hpp"

/* XXX This code is known to be buggy. It does not pass
 * test_reduce_plattice.
 * (It's not used beyond testing. And even then, this code should simply
 * go away)
 */
int
reference (plattice_proxy *pli, const fbprime_t p, const fbroot_t r, uint32_t I)
{
    int32_t i0 = -((int32_t) p), i1 = (int32_t) r, j0 = 0, j1 = 1, k;
    const int32_t hI = (int32_t) I;
    const int32_t mhI = -hI;
    while ((i1 >= hI)) {
        k = i0 / i1; i0 %= i1; j0 -= k * j1;
        if ((i0 > mhI)) break;
        k = i1 / i0; i1 %= i0; j1 -= k * j0;
        /* We may conceivably unroll a bit more, or a bit less, here. Just
         * tuck in as many copies of the following block as you wish. */
        if (UNLIKELY(i1 < hI )) break;
        k = i0 / i1; i0 %= i1; j0 -= k * j1;
        if (UNLIKELY(i0 > mhI)) break;
        k = i1 / i0; i1 %= i0; j1 -= k * j0;
    }
    k = i1 - hI - i0;
    if (i1 > -i0) {
        if (UNLIKELY(!i0)) return 0;
        k /= i0; i1 -= k * i0; j1 -= k * j0;
    } else {
        if (UNLIKELY(!i1)) return 0;
        k /= i1; i0 += k * i1; j0 += k * j1;
    }
    pli->mi0 = - (int32_t) i0; pli->j0 = (uint32_t) j0; pli->i1 = (int32_t) i1; pli->j1 = (uint32_t) j1;
    return 1;
}


#endif	/* LAS_REDUCE_PLATTICE_REFERENCE_HPP_ */
