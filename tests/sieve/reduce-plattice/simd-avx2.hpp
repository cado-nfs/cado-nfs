#ifndef SIMD_AVX2_HPP_
#define SIMD_AVX2_HPP_

#include <cstdint>
#include <array>
#include <algorithm>

#ifdef HAVE_AVX2

#include <x86intrin.h>

/* This one is a bit complicated because we want to take the occasion to
 * use the avx512 mask registers if we happen to have avx512f available.
 */
template<>
struct simd_helper<uint32_t, 8>
{
    typedef uint32_t T;
    static constexpr const size_t N = 8;
    static constexpr const size_t store_alignment = 32;
    typedef __m256i type;
#if defined(HAVE_AVX512F) && defined(HAVE_AVX512DQ) && defined(HAVE_AVX512VL)
    // _mm256_cmpeq_epu32_mask is AVX512VL
    // https://www.intel.com/content/www/us/en/develop/documentation/cpp-compiler-developer-guide-and-reference/top/compiler-reference/intrinsics/intrinsics-for-avx-512-additional-instructions/intrinsics-for-comparison-operations.html
    typedef __mmask8 mask;
    static inline mask zeromask() { return _cvtu32_mask8(0); }
    static inline mask onemask() { return _cvtu32_mask8(~0); }
    static inline int mask2int(mask x) { return _cvtmask8_u32(x); }
    static inline mask knot(mask const a) { return _knot_mask8(a); }
    static inline mask kxor(mask const a, mask const b) { return _kxor_mask8(a, b); }
    static inline mask kand(mask const a, mask const b) { return _kand_mask8(a, b); }
    static inline mask kor(mask const a, mask const b) { return _kor_mask8(a, b); }



    static inline mask cmpeq(type const a, type const b) { return _mm256_cmpeq_epu32_mask(a, b); }
    static inline mask cmpneq(type const a, type const b) { return _mm256_cmpneq_epu32_mask(a, b); }
    static inline mask cmpge(type const a, type const b) { return _mm256_cmpge_epu32_mask(a, b); }
    static inline mask cmplt(type const a, type const b) { return _mm256_cmplt_epu32_mask(a, b); }

    static inline mask mask_cmpge(mask const m, type const a, type const b) { return _mm256_mask_cmpge_epu32_mask(m, a, b); }
    static inline mask mask_cmpeq(mask const m, type const a, type const b) { return _mm256_mask_cmpeq_epu32_mask(m, a, b); }
    static inline mask mask_cmpneq(mask const m, type const a, type const b) { return _mm256_mask_cmpneq_epu32_mask(m, a, b); }
    static inline mask mask_cmplt(mask const m, type const a, type const b) { return _mm256_mask_cmplt_epu32_mask(m, a, b); }
    static inline type mask_sub(type const src, mask const m, type const a, type const b) { return _mm256_mask_sub_epi32(src, m, a, b); }
    static inline type mask_add(type const src, mask const m, type const a, type const b) { return _mm256_mask_add_epi32(src, m, a, b); }
    static inline type mask_bxor(type const src, mask const m, type const a, type const b) { return _mm256_mask_xor_epi32(src, m, a, b); }
#else
    typedef __m256i mask;
    static inline mask zeromask() { return _mm256_setzero_si256(); }
    static inline mask onemask() { return _mm256_set1_epi32(~0); }
    static inline int mask2int(mask x) { return _mm256_movemask_ps(_mm256_castsi256_ps(x)); }
    static inline mask knot(mask const a) { return kxor(a, onemask()); }
    static inline mask kxor(mask const a, mask const b) { return _mm256_xor_si256(a, b); }
    static inline mask kand(mask const a, mask const b) { return _mm256_and_si256(a, b); }
    static inline mask kor(mask const a, mask const b) { return _mm256_or_si256(a, b); }


    static inline mask cmpeq(type const a, type const b) { return _mm256_cmpeq_epi32(a, b); }
    static inline mask cmpneq(type const a, type const b) { return knot(_mm256_cmpeq_epi32(a, b)); }

    private:
    static inline type blendv(type const a, type const b, mask m)
    {
         return _mm256_castps_si256(
                 _mm256_blendv_ps(
                     _mm256_castsi256_ps(a),
                     _mm256_castsi256_ps(b),
                     _mm256_castsi256_ps(m)));
    }
    public:
    static inline mask cmpgt(type const a, type const b) {
        /* must convert signed comparison to unsigned */
        type opposite_sign = bxor(a, b);
        type m = _mm256_cmpgt_epi32(a,b);
        return blendv(zeromask(), onemask(), bxor(m, opposite_sign));
    }
    static inline mask cmpge(type const a, type const b) {
        /* must convert signed comparison to unsigned */
        type opposite_sign = bxor(a, b);
        /* compute b > a first */
        type m = _mm256_cmpgt_epi32(b, a);
        /* then a >= b  is  !(b > a) */
        return blendv(onemask(), zeromask(), bxor(m, opposite_sign));
    }
    static inline mask cmplt(type const a, type const b) { return knot(cmpge(a, b)); }
    static inline mask mask_cmpge(mask const m, type const a, type const b) { return kand(m, cmpge(a, b)); }
    static inline mask mask_cmpeq(mask const m, type const a, type const b) { return kand(m, cmpeq(a, b)); }
    static inline mask mask_cmpneq(mask const m, type const a, type const b) { return kand(m, cmpneq(a, b)); }
    static inline mask mask_cmplt(mask const m, type const a, type const b) { return kand(m, cmplt(a, b)); }
    static inline type mask_sub(type const src, mask const m, type const a, type const b) { return blendv(src, sub(a, b), m); }
    static inline type mask_add(type const src, mask const m, type const a, type const b) { return blendv(src, add(a, b), m); }
    static inline type mask_bxor(type const src, mask const m, type const a, type const b) { return blendv(src, bxor(a, b), m); }
#endif
    static inline type load(T const * p) { return _mm256_load_si256((type const *) p); }
    static inline void store(T * p, type const a) { _mm256_store_si256((type *) p, a); }
    static inline type set1(T const x) { return _mm256_set1_epi32(x); }

    static inline type sub(type const a, type const b) { return _mm256_sub_epi32(a, b); }
    static inline type add(type const a, type const b) { return _mm256_add_epi32(a, b); }
    static inline type bxor(type const a, type const b) { return _mm256_xor_si256(a, b); }
    static inline type mullo(type const a, type const b) { return _mm256_mullo_epi32(a, b); }
    static inline type setzero() { return _mm256_setzero_si256(); }
    static inline type slli(type const a, unsigned int imm) { return _mm256_slli_epi32(a, imm); }
    static inline type div(type const & a, type const & b) {
        /* Unfortunately, _mm256_div_epu32 is a synthetic library
         * function, so we can't use it, lacking a proper library that can do
         * it. We have to cook our own...
         */
#if 0
        __m256i k = _mm256_div_epu32(a, b);
#else
        /* of course it's slow like hell */
        uint32_t explode_a[8] ATTR_ALIGNED(store_alignment);
        uint32_t explode_b[8] ATTR_ALIGNED(store_alignment);
        store(explode_a, a);
        store(explode_b, b);
        for(int i = 0 ; i < 8 ; i++)
            explode_a[i] = explode_a[i] / explode_b[i];
        return load(explode_a);
#endif
    }
    static inline type mask_div(type const & src, mask mm, type const & a, type const & b) {
        /* Unfortunately, _mm256_mask_div_epu32 is a synthetic library
         * function, so we can't use it, lacking a proper library that can do
         * it. We have to cook our own...
         */
#if 0
        __m256i k = _mm256_mask_div_epu32(src, mask, a, b);
#else
        /* of course it's slow like hell */
        uint32_t explode_a[8] ATTR_ALIGNED(store_alignment);
        uint32_t explode_b[8] ATTR_ALIGNED(store_alignment);
        store(explode_a, a);
        store(explode_b, b);
        for(int i = 0, m = mask2int(mm) ; m ; i++,m>>=1)
            if (m&1)
                explode_a[i] = explode_a[i] / explode_b[i];
        // I'm not exactly convinced that there's a value in doing that and not
        // exploding src as well, while we're at it.
#if defined(HAVE_AVX512F) && defined(HAVE_AVX512DQ)
#ifdef HAVE_AVX512VL
        static_assert(std::is_same<__mmask8, mask>::value, "");
        return _mm256_mask_blend_epi32(mm, src, load(explode_a));
#else
        // the line below doesn't work: mm must be an immediate.
        // return _mm256_blend_epi32(src, load(explode_a), mask2int(mm));
        return blendv(src, load(explode_a), mm);
#endif
#else
        return blendv(src, load(explode_a), mm);
#endif
#endif
    }
};
#endif


#endif  /* SIMD_AVX2_HPP_ */
