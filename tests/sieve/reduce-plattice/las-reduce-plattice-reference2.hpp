#ifndef LAS_REDUCE_PLATTICE_REFERENCE2_HPP_
#define LAS_REDUCE_PLATTICE_REFERENCE2_HPP_

#include "las-plattice.hpp"
#include "reduce-plattice/plattice-proxy.hpp"

/* XXX This code is known to be buggy. It does not pass
 * test_reduce_plattice.
 * (It's not used beyond testing. And even then, this code should simply
 * go away)
 */
int
reference2 (plattice_proxy *pli, const fbprime_t p, const fbroot_t r, const uint32_t I)
{
    const int32_t hI = (int32_t) I;
    int32_t i0 = - (int32_t) p, i1 = (int32_t) r, j0, j1;

#define RPA do {							\
    i0 += i1; j0 += j1;							\
    if (LIKELY(i0 + i1 * 4 > 0)) {					\
        int64_t c0 = i0, c1 = j0;					\
        c0 += i1; c1 += j1; if (LIKELY(c0 <= 0)) { i0 = c0; j0 = c1; }	\
        c0 += i1; c1 += j1; if (LIKELY(c0 <= 0)) { i0 = c0; j0 = c1; }	\
        c0 += i1; c1 += j1; if (LIKELY(c0 <= 0)) { i0 = c0; j0 = c1; }	\
    } else								\
    RPC;								\
} while (0)
#define RPB do {							\
    i1 += i0; j1 += j0;							\
    if (LIKELY(i1 + i0 * 4 < 0)) {					\
        int64_t c0 = i1, c1 = j1;						\
        c0 += i0; c1 += j0; if (LIKELY(c0 >= 0)) { i1 = c0; j1 = c1; }	\
        c0 += i0; c1 += j0; if (LIKELY(c0 >= 0)) { i1 = c0; j1 = c1; }	\
        c0 += i0; c1 += j0; if (LIKELY(c0 >= 0)) { i1 = c0; j1 = c1; }	\
    } else								\
    RPD;								\
} while (0)
#define RPC do {					\
    int64_t k = i0 / i1; i0 %= i1; j0 -= k * j1;	\
} while (0)
#define RPD do {					\
    int64_t k = i1 / i0; i1 %= i0; j1 -= k * j0;	\
} while (0)

    /* This code seems odd (looks after the i0 <= mhI loop),
       but gcc generates the fastest asm with it... */
    j0 = 0; j1 = 1;
    if (LIKELY(i1 >= hI)) {
        const int32_t mhI = -hI;
        RPC;
        while (UNLIKELY(i0 < -0X7FFFFFFF / 5)) {
            RPD;
            if (UNLIKELY(i1 < 0X7FFFFFFF / 5)) goto p15;
            RPC;
        }
        if (LIKELY(i0 <= mhI)) {
            do {
                RPB;
p15:
                if (UNLIKELY(i1 < hI)) break;
                RPA;
            } while (LIKELY(i0 <= mhI));
        }
    }

#undef RPA
#undef RPB
#undef RPC
#undef RPD

    int64_t k = i1 - hI - i0;
    if (i1 > -i0) {
        if (UNLIKELY(!i0)) return 0;
        k /= i0; i1 -= k * i0; j1 -= k * j0;
    } else {
        if (UNLIKELY(!i1)) return 0;
        k /= i1; i0 += k * i1; j0 += k * j1;
    }
    pli->mi0 = -(int32_t) i0; pli->j0 = (uint32_t) j0; pli->i1 = (int32_t) i1; pli->j1 = (uint32_t) j1;
    return 1;
}


#endif	/* LAS_REDUCE_PLATTICE_REFERENCE2_HPP_ */
