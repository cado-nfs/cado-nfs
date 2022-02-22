#ifndef LAS_REDUCE_PLATTICE_SIMD_HPP_
#define LAS_REDUCE_PLATTICE_SIMD_HPP_

#include <cstdint>

#include "reduce-plattice/plattice-proxy.hpp"

#include "reduce-plattice/simd-base.hpp"

#include "reduce-plattice/simd-avx2.hpp"        // empty if no support

#include "reduce-plattice/simd-avx512.hpp"        // empty if no support

/* We have to balance the probability of getting a large quotient with
 * the number of items that we're using simultaneously. The relative cost
 * of the full division versus the subtractive algorithm matters, too.
 * The thresholds in the specific instantiations have been determined
 * based on quick testing.
 */
template<typename T>
struct simd_div_threshold {
    static constexpr const int value = 4;
};

template<>
struct simd_div_threshold<simd_helper<uint32_t, 8>> {
    static constexpr const int value = 4;
};

#ifdef HAVE_AVX512F
template<>
struct simd_div_threshold<simd_helper<uint32_t, 16>> {
    static constexpr const int value = 5;
};
#endif

template<size_t N>
/* This does N instances of reduce_plattice in parallel */
void simd(plattice_proxy * pli, uint32_t I)
{
    typedef simd_helper<uint32_t, N> A;
    typedef typename A::mask mask;
    typedef typename A::type data;

    /* This is the main reduce_plattice loop */
    data zI = A::set1(I);
    /* This is just for fun, yes, we're doing N times the same
     * thing.
     */
    uint32_t explode[N] ATTR_ALIGNED(A::store_alignment);
    for(size_t j = 0 ; j < N ; j++) explode[j] = pli[j].mi0;
    data zmi0 = A::load(explode);
    for(size_t j = 0 ; j < N ; j++) explode[j] = pli[j].i1;
    data zi1 = A::load(explode);
    for(size_t j = 0 ; j < N ; j++) explode[j] = pli[j].j0;
    data zj0 = A::load(explode);
    for(size_t j = 0 ; j < N ; j++) explode[j] = pli[j].j1;
    data zj1 = A::load(explode);

    mask proceed = A::kor(
            A::cmpneq(zj0, A::setzero()),
            A::cmpneq(zj1, A::setzero()));

    mask flip;
    for(flip = A::zeromask() ; ; ) {
        /* as long as i1 >= I, proceed */
        proceed = A::mask_cmpge(proceed, zi1, zI);

        if (!A::mask2int(proceed)) break;

        mask toobig = A::mask_cmpge(proceed, zmi0, A::slli(zi1, simd_div_threshold<A>::value));
        if (UNLIKELY(A::mask2int(toobig))) {
            /* Some quotient is larger than 32. We must do a full
             * division.
             */
            data k = A::mask_div(A::setzero(), proceed, zmi0, zi1);
            /* note that A::mullo has a 10-cycle latency on
             * skylake... */
            zmi0 = A::sub(zmi0, A::mullo(k, zi1));
            zj0  = A::add(zj0,  A::mullo(k, zj1));
        } else {
            /* if mi0 >= i1, do a subtraction */
            mask subtract = A::mask_cmpge(proceed, zmi0, zi1);
            /* XXX in fact, it seems that subtract == proceed, right ? */
            zmi0 = A::mask_sub(zmi0, subtract, zmi0, zi1);
            zj0  = A::mask_add(zj0,  subtract, zj0,  zj1);
        }
        /* Any zmi0[j] which is now < zi1[j] deserves a swap */
        mask swap = A::cmplt(zmi0, zi1);
        data swapper;
        swapper = A::mask_bxor(A::setzero(), swap, zmi0, zi1);
        zmi0 = A::bxor(zmi0, swapper);
        zi1 = A::bxor(zi1, swapper);
        swapper = A::mask_bxor(A::setzero(), swap, zj0, zj1);
        zj0 = A::bxor(zj0, swapper);
        zj1 = A::bxor(zj1, swapper);
        flip = A::kxor(flip, swap);
    }

    proceed = A::kor(   A::cmpneq(zj0, A::setzero()),
            A::cmpneq(zj1, A::setzero()));

    mask haszero = A::mask_cmpeq(proceed, zi1, A::setzero());

    if (A::mask2int(haszero)) {
        /* This is exceptional. Explode back to single case */
        A::store(explode,  zmi0);
        for(size_t j = 0 ; j < N ; j++) pli[j].mi0 = explode[j];
        A::store(explode,  zi1);  
        for(size_t j = 0 ; j < N ; j++) pli[j].i1  = explode[j];
        A::store(explode,  zj0);  
        for(size_t j = 0 ; j < N ; j++) pli[j].j0  = explode[j];
        A::store(explode,  zj1);  
        for(size_t j = 0 ; j < N ; j++) pli[j].j1  = explode[j];

        int p = A::mask2int(proceed);

        for(size_t j = 0, m = A::mask2int(flip) ; j < N ; j++, m>>=1, p>>=1) {
            if (!(p & 1)) continue;
            if (pli[j].i1 == 0) {
                if (!(m&1))
                    pli[j].j0 = pli[j].j1 - pli[j].j0;
                pli[j].reduce_with_vertical_vector(I);
                continue;
            } else {
                int a = (pli[j].mi0 + pli[j].i1 - I) / pli[j].i1;
                pli[j].mi0 -= a * pli[j].i1;
                pli[j].j0  += a * pli[j].j1;
            }
            if (m&1) {
                std::swap(pli[j].mi0, pli[j].i1);
                std::swap(pli[j].j0, pli[j].j1);
            }
        }
    } else {
        /* proceed with the normal scenario */
        data sum = A::sub(A::add(zmi0, zi1), zI);
        mask toobig = A::cmpge(sum, A::slli(zi1, simd_div_threshold<A>::value));
        if (UNLIKELY(A::mask2int(toobig))) {
            data k = A::mask_div(A::setzero(), proceed, sum, zi1);
            /* note that the avx512 A::mullo has a 10-cycle latency on
             * skylake... */
            zmi0 = A::sub(zmi0, A::mullo(k, zi1));
            zj0  = A::add(zj0,  A::mullo(k, zj1));
        } else {
            for(mask q ; q = A::cmpge(sum, zi1), A::mask2int(q) ; ) {
                zmi0 = A::mask_sub(zmi0, q, zmi0, zi1);
                sum = A::mask_sub(sum, q, sum, zi1);
                zj0 = A::mask_add(zj0, q, zj0, zj1);
            }
        }
        /* based on fip, we should do a few swaps */
        data swapper;
        swapper = A::mask_bxor(A::setzero(), flip, zmi0, zi1);
        zmi0 = A::bxor(zmi0, swapper);
        zi1 = A::bxor(zi1, swapper);
        swapper = A::mask_bxor(A::setzero(), flip, zj0, zj1);
        zj0 = A::bxor(zj0, swapper);
        zj1 = A::bxor(zj1, swapper);

        A::store(explode,  zmi0);
        for(size_t j = 0 ; j < N ; j++) pli[j].mi0 = explode[j];
        A::store(explode,  zi1);  
        for(size_t j = 0 ; j < N ; j++) pli[j].i1  = explode[j];
        A::store(explode,  zj0);  
        for(size_t j = 0 ; j < N ; j++) pli[j].j0  = explode[j];
        A::store(explode,  zj1);  
        for(size_t j = 0 ; j < N ; j++) pli[j].j1  = explode[j];
    }
}


#endif	/* LAS_REDUCE_PLATTICE_SIMD_HPP_ */
