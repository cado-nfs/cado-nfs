#ifndef SIMD_AVX512_HPP_
#define SIMD_AVX512_HPP_

#include <cstdint>
#include <array>
#include <algorithm>

#ifdef HAVE_AVX512F

#include <x86intrin.h>

template<>
struct simd_helper<uint32_t, 16>
{
    typedef uint32_t T;
    static constexpr const size_t N = 16;
    static constexpr const size_t store_alignment = 64;
    typedef __m512i type;
    typedef __mmask16 mask;
    static inline mask zeromask() { return _mm512_int2mask(0); }
    static inline mask onemask() { return _mm512_int2mask(~0); }
    static inline type load(T const * p) { return _mm512_load_epi32(p); }
    static inline type mask_load(type const & src, mask m, T const * p) { return _mm512_mask_load_epi32(src, m, p); }
    static inline void store(T * p, type const a) { _mm512_store_epi32(p, a); }
    static inline type set1(T const x) { return _mm512_set1_epi32(x); }
    static inline int mask2int(mask x) { return _mm512_mask2int(x); }
    static inline mask cmpeq(type const a, type const b) { return _mm512_cmpeq_epu32_mask(a, b); }
    static inline mask cmpneq(type const a, type const b) { return _mm512_cmpneq_epu32_mask(a, b); }
    static inline mask cmpge(type const a, type const b) { return _mm512_cmpge_epu32_mask(a, b); }
    static inline mask mask_cmpge(mask const m, type const a, type const b) { return _mm512_mask_cmpge_epu32_mask(m, a, b); }
    static inline mask mask_cmpeq(mask const m, type const a, type const b) { return _mm512_mask_cmpeq_epu32_mask(m, a, b); }
    static inline mask mask_cmpneq(mask const m, type const a, type const b) { return _mm512_mask_cmpneq_epu32_mask(m, a, b); }
    static inline mask cmplt(type const a, type const b) { return _mm512_cmplt_epu32_mask(a, b); }
    static inline mask mask_cmplt(mask const m, type const a, type const b) { return _mm512_mask_cmplt_epu32_mask(m, a, b); }
    static inline type sub(type const a, type const b) { return _mm512_sub_epi32(a, b); }
    static inline type add(type const a, type const b) { return _mm512_add_epi32(a, b); }
    static inline type mask_sub(type const src, mask const m, type const a, type const b) { return _mm512_mask_sub_epi32(src, m, a, b); }
    static inline type mask_add(type const src, mask const m, type const a, type const b) { return _mm512_mask_add_epi32(src, m, a, b); }
    static inline type mask_bxor(type const src, mask const m, type const a, type const b) { return _mm512_mask_xor_epi32(src, m, a, b); }
    static inline type bxor(type const a, type const b) { return _mm512_xor_epi32(a, b); }
    static inline type mullo(type const a, type const b) { return _mm512_mullo_epi32(a, b); }
    static inline type setzero() { return _mm512_setzero_epi32(); }
    static inline type slli(type const a, unsigned int imm) { return _mm512_slli_epi32(a, imm); }
    static inline type div(type const & a, type const & b) {
        /* Unfortunately, _mm512_div_epu32 is a synthetic library
         * function, so we can't use it, lacking a proper library that can do
         * it. We have to cook our own...
         */
#if 0
        __m512i k = _mm512_div_epu32(a, b);
#else
        /* of course it's slow like hell */
        uint32_t explode_a[16] ATTR_ALIGNED(store_alignment);
        uint32_t explode_b[16] ATTR_ALIGNED(store_alignment);
        store(explode_a, a);
        store(explode_b, b);
        for(int i = 0 ; i < 16 ; i++)
            explode_a[i] = explode_a[i] / explode_b[i];
        return load(explode_a);
#endif
    }
    static inline type mask_div(type const & src, mask mm, type const & a, type const & b) {
        /* Unfortunately, _mm512_mask_div_epu32 is a synthetic library
         * function, so we can't use it, lacking a proper library that can do
         * it. We have to cook our own...
         */
#if 0
        __m512i k = _mm512_mask_div_epu32(src, mask, a, b);
#else
        /* of course it's slow like hell */
        uint32_t explode_a[16] ATTR_ALIGNED(store_alignment);
        uint32_t explode_b[16] ATTR_ALIGNED(store_alignment);
        store(explode_a, a);
        store(explode_b, b);
        for(int i = 0, m = mask2int(mm) ; m ; i++,m>>=1)
            if (m&1)
                explode_a[i] = explode_a[i] / explode_b[i];
        return mask_load(src, mm, explode_a);
#endif
    }

// intel documents both the _kXXX_mask16 and the _mm512_kXXX versions.
// gcc has both via the following #defines, but clang (at least clang7)
// seems to only have th eformer.
//
// #define _kand_mask16 _mm512_kand
// #define _kandn_mask16 _mm512_kandn
// #define _knot_mask16 _mm512_knot
// #define _kor_mask16 _mm512_kor
// #define _kxnor_mask16 _mm512_kxnor
// #define _kxor_mask16 _mm512_kxor
#if 0
    static inline mask knot(mask const a) { return _knot_mask16(a); }
    static inline mask kxor(mask const a, mask const b) { return _kxor_mask16(a, b); }
    static inline mask kand(mask const a, mask const b) { return _kand_mask16(a, b); }
    static inline mask kor(mask const a, mask const b) { return _kor_mask16(a, b); }
#endif

    static inline mask knot(mask const a) { return _mm512_knot(a); }
    static inline mask kxor(mask const a, mask const b) { return _mm512_kxor(a, b); }
    static inline mask kand(mask const a, mask const b) { return _mm512_kand(a, b); }
    static inline mask kor(mask const a, mask const b) { return _mm512_kor(a, b); }
};
#endif



#endif	/* SIMD_AVX512_HPP_ */
