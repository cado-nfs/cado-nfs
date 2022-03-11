#ifndef SIMD_BASE_HPP_
#define SIMD_BASE_HPP_

#include <cstdint>
#include <array>
#include <algorithm>

/* This simd class is functionally equivalent to the Intel simd
 * intrinsics, and can be replaced by explicit instantiations. For
 * debugging, it is nice to have the C++ code, of course !
 */
template<typename T, int N>
struct simd_helper {
    typedef std::array<T, N> type;
    static_assert(N <= 64, "masks cannot be larger than 64 bits");
    static constexpr const size_t store_alignment = sizeof(T);
    typedef uint64_t mask;
    static inline mask zeromask() { return 0; }
    static inline mask onemask() { return (mask(1) << N) - 1; }
    static inline type load(T const * p) { type r; std::copy_n(p, N, r.begin()); return r; }
    static inline type mask_load(type const & src, mask m, T const * p) {
        type r = src;
        for(size_t i = 0 ; i < N ; i++, m>>=1)
            if (m&1)
                r[i] = p[i];
        return r;
    }
    static inline void store(T * p, type const a) { std::copy_n(a.begin(), N, p); }
    static inline type set1(T const x) {
        type r;
        for(size_t i = 0 ; i < N ; i++) r[i] = x;
        return r;
    }
    static inline int mask2int(mask x) { return x; }
    static inline mask cmpneq(type const & a, type const & b) {
        mask m = 0;
        for(size_t i = 0 ; i < N ; i++) m |= mask(a[i] != b[i]) << i;
        return m;
    }
    static inline mask cmpeq(type const & a, type const & b) {
        mask m = 0;
        for(size_t i = 0 ; i < N ; i++) m |= mask(a[i] == b[i]) << i;
        return m;
    }
    static inline mask cmpge(type const & a, type const & b) {
        mask m = 0;
        for(size_t i = 0 ; i < N ; i++) m |= mask(a[i] >= b[i]) << i;
        return m;
    }
    static inline mask mask_cmpge(mask const m, type const a, type const b) { return m & cmpge(a, b); }
    static inline mask mask_cmpeq(mask const m, type const a, type const b) { return m & cmpeq(a, b); }
    static inline mask mask_cmpneq(mask const m, type const a, type const b) { return m & cmpneq(a, b); }
    static inline mask cmplt(type const & a, type const & b) {
        mask m = 0;
        for(size_t i = 0 ; i < N ; i++) m |= mask(a[i] < b[i]) << i;
        return m;
    }
    static inline mask mask_cmplt(mask const m, type const a, type const b) { return m & cmplt(a, b); }
    static inline type sub(type const & a, type const & b) {
        type r;
        for(size_t i = 0 ; i < N ; i++) r[i] = a[i] - b[i];
        return r;
    }
    static inline type add(type const & a, type const & b) {
        type r;
        for(size_t i = 0 ; i < N ; i++) r[i] = a[i] + b[i];
        return r;
    }
    static inline type mask_sub(type const & src, mask m, type const & a, type const & b) {
        type r = src;
        for(size_t i = 0 ; i < N ; i++, m>>=1) if (m&1) r[i] = a[i] - b[i];
        return r;
    }
    static inline type mask_add(type const & src, mask m, type const & a, type const & b) {
        type r = src;
        for(size_t i = 0 ; i < N ; i++, m>>=1) if (m&1) r[i] = a[i] + b[i];
        return r;
    }
    static inline type mask_bxor(type const & src, mask m, type const & a, type const & b) {
        type r = src;
        for(size_t i = 0 ; i < N ; i++, m>>=1) if (m&1) r[i] = a[i] ^ b[i];
        return r;
    }
    static inline type mullo(type const & a, type const & b) {
        type r;
        for(size_t i = 0 ; i < N ; i++) r[i] = a[i] * b[i];
        return r;
    }
    static inline type bxor(type const & a, type const & b) {
        type r;
        for(size_t i = 0 ; i < N ; i++) r[i] = a[i] ^ b[i];
        return r;
    }
    static inline type setzero() { return set1(0); }
    static inline type slli(type const a, unsigned int imm) {
        type r = a;
        for(size_t i = 0 ; i < N ; i++) r[i] <<= imm;
        return r;
    }
    static inline type div(type const & a, type const & b) {
        type r;
        for(size_t i = 0 ; i < N ; i++) r[i] = a[i] / b[i];
        return r;
    }
    static inline type mask_div(type const & src, mask m, type const & a, type const & b) {
        type r = src;
        for(size_t i = 0 ; i < N ; i++, m>>=1) if (m&1) r[i] = a[i] / b[i];
        return r;
    }
    static inline mask knot(mask const a) { return onemask() ^ a; }
    static inline mask kxor(mask const a, mask const b) { return a ^ b; }
    static inline mask kand(mask const a, mask const b) { return a & b; }
    static inline mask kor(mask const a, mask const b) { return a | b; }
};

#endif	/* SIMD_BASE_HPP_ */
