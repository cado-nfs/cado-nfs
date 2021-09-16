#include "cado.h"
#include <sys/types.h>
#include <unistd.h>
#include <stdint.h>
#include <vector>
#include <array>
#include <tuple>
#include <string>
#include <map>
#include <stdio.h>
#include "gcd.h"
#include "macros.h"
#include "getprime.h"
#include "fb-types.h"
#include "las-plattice.hpp"
#include "fmt/format.h"
#include <algorithm>
#if defined(HAVE_AVX512F) || defined(HAVE_AVX2) || defined(HAVE_AVX) || defined(HAVE_SSE41)
#include <x86intrin.h>
#endif

/* see plattice.sage */

/* The c code that is currently in production is at very least not
 * slower, and possibly even mildly faster than the old assembly or C
 * code. However the measurement is not easy. Note also that the new code
 * covers many more cases.
 *
 * This file contains several different implementations of
 * reduce_plattice, including experimental SIMD code. The fastest code is
 * the AVX-512 variant, but AVX-2 isn't bad. Note that the performance
 * depends a lot on the prime size, because that controls the loop
 * length. Both seem to uniformly beat the single data implementations.
 *
 * avx-2 machine  i5-7500 CPU @ 3.40GHz
 *
 * test-reduce-plattice -I 65536 -ntests 100000000 -B $((2**25)) -T -seed 1 -N
 *
 * production (two_legs): 100000000 tests in 3.8461s
 * simplistic: 100000000 tests in 3.8710s
 * using_64bit_mul: 100000000 tests in 3.8622s
 * swapping_loop: 100000000 tests in 4.0930s
 * swapping_loop2: 100000000 tests in 6.2752s
 * old_reference: 100000000 tests in 4.2740s [ BUGGY, MANY ERRORS REPORTED ]
 * reference2: 100000000 tests in 4.0639s [ BUGGY, MANY ERRORS REPORTED ]
 * reference2_asm: 100000000 tests in 3.9035s [ BUGGY, MANY ERRORS REPORTED ]
 * simd-avx2: 100000000 tests in 3.2203s
 *
 * test-reduce-plattice -I 65536 -ntests 100000000 -B $((2**32-1)) -T -seed 1 -N
 *
 * production (two_legs): 100000000 tests in 7.3322s
 * simplistic: 100000000 tests in 7.1406s
 * using_64bit_mul: 100000000 tests in 7.4497s
 * swapping_loop: 100000000 tests in 7.4919s
 * swapping_loop2: 100000000 tests in 11.9772s
 * old_reference: 100000000 tests in 4.4743s [ BUGGY, MANY ERRORS REPORTED ]
 * reference2: 100000000 tests in 4.3514s [ BUGGY, MANY ERRORS REPORTED ]
 * reference2_asm: 100000000 tests in 4.2292s [ BUGGY, MANY ERRORS REPORTED ]
 * simd-avx2: 100000000 tests in 5.1874s
 *
 * avx-512 machine Xeon(R) Gold 6130 CPU @ 2.10GHz
 *
 * production (two_legs): 100000000 tests in 7.6957s
 * simplistic: 100000000 tests in 7.5916s
 * using_64bit_mul: 100000000 tests in 7.9518s
 * swapping_loop: 100000000 tests in 8.6450s
 * swapping_loop2: 100000000 tests in 12.1884s
 * old_reference: 100000000 tests in 4.5391s [ BUGGY, MANY ERRORS REPORTED ]
 * reference2: 100000000 tests in 4.5022s [ BUGGY, MANY ERRORS REPORTED ]
 * reference2_asm: 100000000 tests in 4.4340s [ BUGGY, MANY ERRORS REPORTED ]
 * simd-avx512: 100000000 tests in 3.7655s
 *
 */

struct plattice : public plattice_info {
    /* use this default ctor only for the comparison with the reference
     * routines.
     */
    plattice () = default;

    using plattice_info::mi0;
    using plattice_info::j0;
    using plattice_info::i1;
    using plattice_info::j1;
    using plattice_info::check_post_conditions;

#ifndef NDEBUG
#define ASSERT_THROW(e, c) do { if (!(c)) throw (e)(#c); } while (0)
#else
#define ASSERT_THROW(e, c)
#endif
    struct error : public std::runtime_error {
        error(const char * s) : std::runtime_error(s) {}
    };
#define ASSERT_PLATTICE(c) ASSERT_THROW(error, c)

    /* The idea is tempting, but the cost of the full-width 64-bit
     * multiplication is really killing us.
     *
     * Also note that there are several nasty corner cases here and there
     * that force us to keep track of mi0 and i1 almost all the time
     * anyway.
     *
     * This code has many asserts that check that the 64-bit
     * representatives are kept in sync with the (-mi0,j0) and (i1,j1)
     * vectors. Only the debug code checks that. The non-debug code skips
     * this bookkeeping.
     */
    void using_64bit_mul(uint32_t I) {
        uint64_t Ix = uint64_t(I) << 32;
        constexpr uint64_t M = uint64_t(-1) << 32;
        constexpr uint64_t L MAYBE_UNUSED = uint32_t(-1);
        uint64_t ij0 = (uint64_t(-mi0) << 32) | j0;
        uint64_t ij1 = (uint64_t(i1) << 32) | j1;
        uint64_t t;
#define ASSERT_CONSISTENCY() do {					\
            ASSERT_PLATTICE(mi0 == -((uint32_t)(ij0>>32)));		\
            ASSERT_PLATTICE(i1 == ij1 >> 32);				\
            ASSERT_PLATTICE(j0 == (ij0 & L));				\
            ASSERT_PLATTICE(j1 == (ij1 & L));				\
} while (0)

        for( ;; ) {
            ASSERT_PLATTICE(j0 <= j1);
            ASSERT_CONSISTENCY();
            t = ij1 - Ix;
            i1 = ij1 >> 32;
            ASSERT_PLATTICE((i1 < I) == (t >= ij1));
            mi0 = -((uint32_t)(ij0>>32));
            if (t >= ij1) { // i1 < I) {
                j1 = ij1;
                if (!i1) {
                    j0 = ij1 - ij0;
                    lattice_with_vertical_vector(I);
                    return;
                }
                ASSERT_PLATTICE(mi0 + i1 >= I);

                int a = (mi0 + i1 - I) / i1;
                ij0 += a * ij1;
#ifndef NDEBUG
                mi0 -= a * i1;
                j0  += a * j1;
#endif
                ASSERT_CONSISTENCY();
                mi0 = -((uint32_t)(ij0>>32));
                j0 = ij0;
                return;
            }
            if (mi0 < i1 * 3) {
                {
                    ij0 += ij1;
                    mi0 -= i1; // see below
#ifndef NDEBUG
                    j0 += j1;
#endif
                }
                /* what we really want to test is mi0 >= i1. However the
                 * j coordinates come into play when comparing -ij0 with
                 * ij1. 
                 */
                if (mi0 >= i1) { // if (-ij0 >= ij1) {
                    ij0 += ij1;
#ifndef NDEBUG
                    mi0 -= i1; j0 += j1;
#endif
                }
            } else
            {
                int k = mi0 / i1;
                // int k = (-(ij0 & M)) / (ij1 & M);
                ASSERT_PLATTICE(k);
                ij0 += k * ij1;
#ifndef NDEBUG
                mi0 -= k * i1;
                j0 += k * j1;
#endif
                ASSERT_CONSISTENCY();
            }

            ASSERT_PLATTICE(j1 <= j0);

            ASSERT_CONSISTENCY();
            t = (ij0&M) + Ix;
            ASSERT_PLATTICE((mi0 < I) == ((-(ij0&M)) < Ix));
            ASSERT_PLATTICE((mi0 < I) == (t != 0 && t <= Ix));
            mi0 = -((uint32_t)(ij0>>32));
            i1 = ij1 >> 32;
            if (t != 0 && t <= Ix) { // mi0 < I) {
                if (!mi0) {
                    mi0 = i1;
                    j0 = ij1;
                    j1 = ij0;
                    i1 = 0;
                    lattice_with_vertical_vector(I);
                    return;
                }
                int a = (mi0 + i1 - I) / mi0;
                ij1 += a * ij0;
#ifndef NDEBUG
                i1 -= a * mi0;
                j1 += a * j0;
#endif
                ASSERT_CONSISTENCY();
                i1 = ij1 >> 32;
                j1 = ij1;
                j0 = ij0;
                return;
            }
            if (i1 < 3 * mi0) {
                {
                    ij1 += ij0;
                    i1 -= mi0; // see below
#ifndef NDEBUG
                    j1 += j0;
#endif
                }
                // same remark as above
                if (i1 >= mi0) { // if (ij1 >= -ij0) {
                    ij1 += ij0;
#ifndef NDEBUG
                    i1 -= mi0; j1 += j0;
#endif
                }
            } else
            {
                int k = i1 / mi0;
                ASSERT(k);
                ij1 += k * ij0;
#ifndef NDEBUG
                i1 -= k * mi0; j1 += k * j0;
#endif
                ASSERT_CONSISTENCY();
            }
        }
#undef ASSERT_CONSISTENCY
    }

    void simplistic(uint32_t I) {
        /* This is the main reduce_plattice loop */
        for( ;; ) {
            ASSERT(j0 <= j1);
            if (i1 < I) {
                if (i1 == 0) {
                    // Lo=matrix([ (mi0, j1-j0), (i1, j1)])
                    j0 = j1 - j0;
                    lattice_with_vertical_vector(I);
                    return;
                }
                ASSERT(mi0 + i1 >= I);
                int a = (mi0 + i1 - I) / i1;
                mi0 -= a * i1;
                j0  += a * j1;
                return;
            }
            {
                int k = mi0 / i1; mi0 -= k * i1; j0 += k * j1;
                ASSERT(k);
            }

            ASSERT(j1 <= j0);
            if (mi0 < I) {
                if (mi0 == 0) {
                    mi0 = i1;
                    i1 = j0 ; j0 = j1 ; j1 = i1;
                    i1 = 0;
                    lattice_with_vertical_vector(I);
                    return;
                }
                ASSERT(mi0 + i1 >= I);
                int a = (mi0 + i1 - I) / mi0;
                i1 -= a * mi0;
                j1 += a * j0;
                return;
            }
            {
                int k = i1 / mi0; i1 -= k * mi0; j1 += k * j0;
                ASSERT(k);
            }
        }
    }

    void swapping_loop(uint32_t I)
    {
        /* This is the main reduce_plattice loop */
        int flip;
        for(flip = 0 ; i1 >= I; flip ^= 1 ) {
            /* do partial unrolling for the frequent case where the
             * quotient is either 1 or 2.
             * this has a significant overall impact
             */
#if 0
            if (mi0 < i1 * 3) {
                { mi0 -= i1; j0 += j1; }
                if (mi0 >= i1) { mi0 -= i1; j0 += j1; }
            } else
#endif
            {
                int k = mi0 / i1; mi0 -= k * i1; j0 += k * j1;
            }
            std::swap(mi0, i1);
            std::swap(j0, j1);
        }
        /* an "UNLIKELY" macro here actually has an adverse
         * effect...  */
        if (i1 == 0) {
            if (!flip)
                // Lo=matrix([ (mi0, j1-j0), (i1, j1)])
                j0 = j1 - j0;
            lattice_with_vertical_vector(I);
        } else {
            int a = (mi0 + i1 - I) / i1;
            mi0 -= a * i1;
            j0  += a * j1;
        }
    }

    void instrumented_two_legs(uint32_t I, std::map<int, unsigned long> & T) {
        /* This is the main reduce_plattice loop */
        for( ;; ) {
            if (i1 < I) {
                /* an "UNLIKELY" macro here actually has an adverse
                 * effect...  */
                if (i1 == 0) {
                    // Lo=matrix([ (mi0, j1-j0), (i1, j1)])
                    j0 = j1 - j0;
                    lattice_with_vertical_vector(I);
                    return;
                }
                ASSERT(mi0 + i1 >= I);
                int a = (mi0 + i1 - I) / i1;
                mi0 -= a * i1;
                j0  += a * j1;
                return;
            }
            {
                int k = mi0 / i1; mi0 -= k * i1; j0 += k * j1;
                T[k]++;
            }
            if (mi0 < I) {
                /* an "UNLIKELY" macro here actually has an adverse
                 * effect...  */
                if (mi0 == 0) {
                    mi0 = i1;
                    i1 = j0 ; j0 = j1 ; j1 = i1;
                    i1 = 0;
                    lattice_with_vertical_vector(I);
                    return;
                }
                ASSERT(mi0 + i1 >= I);
                int a = (mi0 + i1 - I) / mi0;
                i1 -= a * mi0;
                j1 += a * j0;
                return;
            }
            {
                int k = i1 / mi0; i1 -= k * mi0; j1 += k * j0;
                T[k]++;
            }
        }
    }

    /* XXX
     * beware: this constructor takes I, but it shadows a constructor in
     * the production code which takes only logI !!!
     */
    plattice(const unsigned long q, const unsigned long r, bool proj, uint32_t I) : plattice_info()
    {
        initial_basis(q, r, proj);
        /* At this point, (mi0,j0) represents itself, i.e. a vector with
         * two positive coordinates.
         * Note that j0==0
         */
        // ASSERT_ALWAYS(check_pre_conditions(I));
        bool needs_special_treatment = (i1 == 0 || (j1 > 1 && mi0 < I));
        if (needs_special_treatment) {
            lattice_with_vertical_vector(I);
            return;
        }
        two_legs(I);
        // simplistic(I);
        // using_64bit_mul(I);
        // swapping_loop(I);
    }

    bool early(uint32_t I) {
        bool needs_special_treatment = (i1 == 0 || (j1 > 1 && mi0 < I));
        if (needs_special_treatment)
            lattice_with_vertical_vector(I);
        return needs_special_treatment;
    }

    public:
    using plattice_info::initial_basis;
    using plattice_info::two_legs;
    using plattice_info::lattice_with_vertical_vector;

    // friend void instrumented_two_legs(plattice *pli, const unsigned long q, const unsigned long r, bool proj, uint32_t I, std::map<int, unsigned long> & T);
};

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
        std::fill_n(r.begin(), N, x);
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

#ifdef HAVE_AVX2
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
#if defined(HAVE_AVX512F) && defined(HAVE_AVX512DQ)
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
#if defined(HAVE_AVX512F) && defined(HAVE_AVX512DQ)
#ifdef HAVE_AVX512VL
        static_assert(std::is_same<__mmask8, mask>::value);
        return _mm256_mask_blend_epi32(mm, src, load(explode_a));
#else
        return _mm256_blend_epi32(src, load(explode_a), mask2int(mm));
#endif
#else
        return blendv(src, load(explode_a), mm);
#endif
#endif
    }
};
#endif

#ifdef HAVE_AVX512F

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
void simd(plattice * pli, uint32_t I)
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
                pli[j].lattice_with_vertical_vector(I);
            } else {
                int a = (pli[j].mi0 + pli[j].i1 - I) / pli[j].i1;
                pli[j].mi0 -= a * pli[j].i1;
                pli[j].j0  += a * pli[j].j1;
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

void swapping_loop2(plattice & pli, uint32_t I)
{
    uint32_t mi0 = pli.mi0;
    uint32_t i1 = pli.i1;
    uint32_t j0 = pli.j0;
    uint32_t j1 = pli.j1;
    /* This is the main reduce_plattice loop */
    int flip;
    for(flip = 0 ; ; ) {
        int proceed = i1 >= I;
        if (!proceed) break;
        int toobig = proceed & (mi0 >= (i1 << 5));
        if (UNLIKELY(toobig)) {
            /* Some quotient is larger than 32. We must do a full
             * division.
             */
            uint32_t k = proceed ? (mi0 / i1) : 0;
            mi0 -= k * i1; j0 += k * j1;
        } else {
            int subtract = mi0 >= i1;
            mi0 = subtract ? (mi0 - i1) : mi0;
            j0  = subtract ? ( j0 + j1) : j0;
        }
        /* Any zmi0[j] which is now < zi1[j] deserves a swap */
        int swap = mi0 < i1;
        uint32_t swapper;
        swapper = swap ? (mi0 ^ i1) : 0; mi0 = mi0 ^ swapper; i1 = i1 ^ swapper;
        swapper = swap ?  (j0 ^ j1) : 0;  j0 =  j0 ^ swapper; j1 = j1 ^ swapper;
        flip = flip ^ swap;
    }
    /* an "UNLIKELY" macro here actually has an adverse
     * effect...  */
    if (i1 == 0) {
        if (!flip)
            // Lo=matrix([ (mi0, j1-j0), (i1, j1)])
            j0 = j1 - j0;
        pli.mi0 = mi0;
        pli.i1 = i1;
        pli.j0 = j0;
        pli.j1 = j1;
        pli.lattice_with_vertical_vector(I);
        return;
    } else {
        int a = (mi0 + i1 - I) / i1;
        mi0 -= a * i1;
        j0  += a * j1;
    }
    pli.mi0 = mi0;
    pli.i1 = i1;
    pli.j0 = j0;
    pli.j1 = j1;
}

int
reference (plattice *pli, const fbprime_t p, const fbroot_t r, uint32_t I)
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

int
reference2 (plattice *pli, const fbprime_t p, const fbroot_t r, const uint32_t I)
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

#if defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM) && !(defined(__APPLE_CC__) && defined(__llvm__) && __APPLE_CC__ == 5621) && !defined(__INTEL_COMPILER)
/* icc can't compile the avx512 code if we happen to enable this code.
 * https://community.intel.com/t5/Intel-C-Compiler/asm-callq-and-kand-mask8-intrinsic-generate-vkmovb-quot-no-such/m-p/1140906
 */
#define TEST_ASSEMBLY_CODE_DELETED_BY_5f258ce8b
#endif

#ifdef TEST_ASSEMBLY_CODE_DELETED_BY_5f258ce8b
    int
reference2_asm (plattice *pli, const fbprime_t p, const fbroot_t r, const uint32_t I)
{
    const int32_t hI = (int32_t) I;
    int32_t i0 = - (int32_t) p, i1 = (int32_t) r, j0, j1;

    /* Mac OS X 10.8 embarks a version of llvm which crashes on the code
     * below (could be that the constraints are exerting too much of the
     * compiler's behaviour).
     *
     * See tracker #16540
     */

#define RPA(LABEL)      \
    "addl %2, %0\n"						\
    "leal (%0,%2,4), %%edx\n"					\
    "addl %3, %1\n"						\
    "testl %%edx, %%edx\n"					\
    "movl %0, %%eax\n"						\
    "jng " LABEL "\n"						\
    "addl %2, %%eax\n"						\
    "leal (%1,%3,1), %%edx\n"					\
    "cmovngl %%eax, %0\n"						\
    "cmovngl %%edx, %1\n"						\
    "addl %3, %%edx\n"						\
    "addl %2, %%eax\n"						\
    "cmovngl %%edx, %1\n"						\
    "cmovngl %%eax, %0\n"						\
    "addl %3, %%edx\n"						\
    "addl %2, %%eax\n"						\
    "cmovngl %%edx, %1\n"						\
    "cmovngl %%eax, %0\n"
#define RPB(LABEL)      \
    "addl %0, %2\n"						\
    "leal (%2,%0,4), %%edx\n"					\
    "addl %1, %3\n"						\
    "testl %%edx, %%edx\n"					\
    "movl %2, %%eax\n"						\
    "jns " LABEL "\n"						\
    "addl %0, %%eax\n"						\
    "leal (%1,%3,1), %%edx\n"					\
    "cmovnsl %%eax, %2\n"						\
    "cmovnsl %%edx, %3\n"						\
    "addl %1, %%edx\n"						\
    "addl %0, %%eax\n"						\
    "cmovnsl %%edx, %3\n"						\
    "cmovnsl %%eax, %2\n"						\
    "addl %1, %%edx\n"						\
    "addl %0, %%eax\n"						\
    "cmovnsl %%edx, %3\n"						\
    "cmovnsl %%eax, %2\n"
#define RPC     \
    "cltd\n"							\
    "idivl %2\n"							\
    "imull %3, %%eax\n"						\
    "movl %%edx, %0\n"						\
    "subl %%eax, %1\n"
#define RPD "cltd\n"    \
    "idivl %0\n"							\
    "imull %1, %%eax\n"						\
    "movl %%edx, %2\n"						\
    "subl %%eax, %3\n"

    int32_t mhI;
    __asm__ __volatile__ (
            "xorl %1, %1\n"
            "cmpl %2, %5\n"
            "movl $0x1, %3\n"
            "jg 9f\n"
            "movl %5, %%eax\n"
            "negl %%eax\n"
            "movl %%eax, %4\n"
            "movl %0, %%eax\n"
            "cltd\n"
            "idivl %2\n"
            "subl %%eax, %1\n"
            "cmpl $0xe6666667, %%edx\n"
            "movl %%edx, %0\n"
            "jl 0f\n"
            ".balign 8\n"
            "1:\n"
            "cmpl %0, %4\n"
            "jl 9f\n"
            RPB("3f")
            "2:\n"
            "cmpl %2, %5\n"
            "jg 9f\n"
            RPA("4f")
            "jmp 1b\n"
            ".balign 8\n"
            "3:\n"
            RPD
            "jmp 2b\n"
            ".balign 8\n"
            "4:\n"
            RPC
            "jmp 1b\n"
            ".balign 8\n"
            "0:\n"
            "movl %2, %%eax\n"
            "cltd\n"
            "idivl %0\n"
            "imull %1, %%eax\n"
            "subl %%eax, %3\n"
            "cmpl $0x19999999, %%edx\n"
            "movl %%edx, %2\n"
            "jle 2b\n"
            "movl %0, %%eax\n"
            "cltd\n"
            "idivl %2\n"
            "imull %3, %%eax\n"
            "subl %%eax, %1\n"
            "cmpl $0xe6666667, %%edx\n"
            "movl %%edx, %0\n"
            "jge 1b\n"
            "jmp 0b\n"
            "9:\n"
            : "+&r"(i0), "=&r"(j0), "+&r"(i1), "=&r"(j1),
        "=&rm"(mhI) : "rm"(hI) : "%rax", "%rdx", "cc");

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
#endif

struct post_condition_error : public std::exception {};

struct a_test {
    unsigned long p;
    unsigned long q;
    unsigned long r;
    int k;
};

struct test_wrap {
    uint32_t I = 512;
    bool nocheck = false;
    bool timing = false;
    bool failed = false;
    std::vector<a_test> tests;
};

struct call_production {
    static constexpr const size_t batch_count = 1;
    static constexpr const bool has_known_bugs = false;
    static constexpr const bool old_interface = false;
    static constexpr const char * what = "production (two_legs)";
    static inline void call(plattice & L, uint32_t I, a_test = a_test()) {
        L.two_legs(I);
    }
};

struct call_simplistic {
    static constexpr const size_t batch_count = 1;
    static constexpr const char * what = "simplistic";
    static constexpr const bool has_known_bugs = false;
    static constexpr const bool old_interface = false;
    static inline void call(plattice & L, uint32_t I, a_test = a_test()) {
        L.simplistic(I);
    }
};


struct call_using_64bit_mul {
    static constexpr const size_t batch_count = 1;
    static constexpr const char * what = "using_64bit_mul";
    static constexpr const bool has_known_bugs = false;
    static constexpr const bool old_interface = false;
    static inline void call(plattice & L, uint32_t I, a_test = a_test()) {
        L.using_64bit_mul(I);
    }
};

struct call_swapping_loop {
    static constexpr const size_t batch_count = 1;
    static constexpr const char * what = "swapping_loop";
    static constexpr const bool has_known_bugs = false;
    static constexpr const bool old_interface = false;
    static inline void call(plattice & L, uint32_t I, a_test = a_test()) {
        L.swapping_loop(I);
    }
};

struct call_swapping_loop2 {
    static constexpr const size_t batch_count = 1;
    static constexpr const char * what = "swapping_loop2";
    static constexpr const bool has_known_bugs = false;
    static constexpr const bool old_interface = false;
    static inline void call(plattice & L, uint32_t I, a_test = a_test()) {
        swapping_loop2(L, I);
    }
};

struct call_old_reference {
    static constexpr const size_t batch_count = 1;
    static constexpr const char * what = "old_reference";
    static constexpr const bool has_known_bugs = true;
    static constexpr const bool old_interface = true;
    static inline void call(plattice & L, uint32_t I, a_test const & aa) {
        reference(&L, aa.q, aa.r, I);
    }
};

struct call_reference2 {
    static constexpr const size_t batch_count = 1;
    static constexpr const char * what = "reference2";
    static constexpr const bool has_known_bugs = true;
    static constexpr const bool old_interface = true;
    static inline void call(plattice & L, uint32_t I, a_test const & aa) {
        reference2(&L, aa.q, aa.r, I);
    }
};

#ifdef TEST_ASSEMBLY_CODE_DELETED_BY_5f258ce8b
struct call_reference2_asm {
    static constexpr const size_t batch_count = 1;
    static constexpr const char * what = "reference2_asm";
    static constexpr const bool has_known_bugs = true;
    static constexpr const bool old_interface = true;
    static inline void call(plattice & L, uint32_t I, a_test const & aa) {
        reference2_asm(&L, aa.q, aa.r, I);
    }
};
#endif

template<size_t N>
struct call_simd_base {
    static constexpr const size_t batch_count = N;
    static constexpr const bool has_known_bugs = false;
    static constexpr const bool old_interface = false;
    static inline void call(plattice * pL, uint32_t I) {
        simd<N>(pL, I);
    }
};

struct call_simd_avx512 : public call_simd_base<16> {
    static constexpr const char * what = "simd-avx512";
};

struct call_simd_avx2 : public call_simd_base<8> {
    static constexpr const char * what = "simd-avx2";
};

template<size_t N>
struct call_simd : public call_simd_base<N> {};
template<>
struct call_simd<4> : public call_simd_base<4> {
    static constexpr const char * what = "simd-synthetic<4> (sse-4.1, maybe)";
};
template<>
struct call_simd<2> : public call_simd_base<2> {
    static constexpr const char * what = "simd-synthetic<2>";
};

template<typename T>
inline 
typename std::enable_if<T::old_interface, unsigned long>::type
test_inner(plattice * L, test_wrap const & tw, a_test const * aa)
{
    T::call(*L, tw.I, *aa);
    return L->mi0;
}

template<typename T>
inline
typename std::enable_if<!T::old_interface && T::batch_count == 1, unsigned long>::type
test_inner(plattice * L, test_wrap const & tw, a_test const * aa)
{
    bool proj = aa->r > aa->q;
    L->initial_basis(aa->q, proj ? (aa->r - aa->q) : aa->r, proj);
    if (!L->early(tw.I))
        T::call(*L, tw.I);
    return L->mi0;
}

template<typename T>
inline
typename std::enable_if<!T::old_interface && T::batch_count != 1, unsigned long>::type
test_inner(plattice * L, test_wrap const & tw, a_test const * aa)
{
    constexpr size_t N = T::batch_count;
    size_t j = 0;
    bool normal = true;
    /* prepare the lattices for all calls. */
    for(j = 0 ; j < N ; j++) {
        // unsigned long p = tests[j].p;
        unsigned long q = aa[j].q;
        unsigned long r = aa[j].r;
        bool proj = r > q;
        L[j].initial_basis(q, proj ? (r-q) : r, proj);
        if (L[j].early(tw.I))
            normal = false;
    }
    if (normal) {
        T::call(L, tw.I);
    } else {
        for(j = 0 ; j < N ; j++)
            L[j].two_legs(tw.I);
    }
    return L[0].mi0;
}

template<typename T>
void test_correctness(test_wrap & tw)
{
    if (tw.nocheck) return;
    constexpr size_t N = T::batch_count;
    int nfailed = 0;
    std::vector<a_test> const & tests = tw.tests;
    std::string thiscode = T::what;
    if (T::has_known_bugs)
        thiscode += " (known buggy)";

    for(size_t i = 0 ; i < tests.size() ; i += N) {
        plattice L[N];
        size_t j = 0;
        const char * when = "pre";
        try {
            test_inner<T>(L, tw, &(tests[i]));
            when = "post";
            for(j = 0 ; j < N && i + j < tests.size() ; j++) {
                if (!L[j].check_post_conditions(tw.I))
                    throw post_condition_error();
            }
        } catch (plattice::error const & e) {
            a_test const & aa = tests[i+j];
            std::string c = "failed in-algorithm check";
            std::string t = fmt::format(
                    FMT_STRING("p^k={}^{} r={}"), aa.p, aa.k, aa.r);
            std::string msg = fmt::format(FMT_STRING("{}: {} check for {}\n"),
                    thiscode, c, t);
            fputs(msg.c_str(), stderr);
            if (!T::has_known_bugs) tw.failed = true;
        } catch (post_condition_error const & e) {
            if (nfailed < 16) {
                a_test const & aa = tests[i+j];
                std::string c = fmt::format("failed check ({})", when);
                std::string t = fmt::format(
                        FMT_STRING("p^k={}^{} r={}"), aa.p, aa.k, aa.r);
                std::string msg = fmt::format(
                        FMT_STRING("{}: {} check for {}\n"), thiscode, c, t);
                fputs(msg.c_str(), stderr);
                if (!T::has_known_bugs) tw.failed = true;
                if (++nfailed >= 16)
                    fprintf(stderr, "%s: stopped reporting errors, go fix your program\n", thiscode.c_str());
            }
        }
    }
}

template<typename T>
inline
typename std::enable_if<T::batch_count != 1, unsigned long>::type
test_speed(test_wrap & tw)
{
    test_correctness<T>(tw);
    if (!tw.timing) return 0;
    clock_t clk0 = clock();
    unsigned long dummy_local = 0;
    for(size_t i = 0 ; i < tw.tests.size() ; i += T::batch_count) {
        plattice L[T::batch_count];
        dummy_local += test_inner<T>(L, tw, &tw.tests[i]);
    }
    clock_t clk1 = clock();
    if (tw.timing) {
        printf("# %s: %zu tests in %.4fs\n",
                T::what, tw.tests.size(), ((double)(clk1-clk0))/CLOCKS_PER_SEC);
    }
    return dummy_local;
}

template<typename T>
inline
typename std::enable_if<T::batch_count == 1, unsigned long>::type
test_speed(test_wrap & tw)
{
    test_correctness<T>(tw);
    if (!tw.timing) return 0;
    clock_t clk0 = clock();
    unsigned long dummy_local = 0;
    for(a_test const & aa : tw.tests) {
        plattice L[T::batch_count];
        dummy_local += test_inner<T>(L, tw, &aa);
    }
    clock_t clk1 = clock();
    if (tw.timing) {
        printf("# %s: %zu tests in %.4fs\n",
                T::what, tw.tests.size(), ((double)(clk1-clk0))/CLOCKS_PER_SEC);
    }
    return dummy_local;
}

int main(int argc, char * argv[])
{
    setbuf(stdout, NULL);
    setbuf(stderr, NULL);

    test_wrap tw;

    uint32_t B = 1e6;
    size_t ntests = 1000;
    unsigned long seed = 0;
    bool quiet = false;

    for( ; argc > 1 ; argv++,argc--) {
        if (strcmp(argv[1], "-B") == 0) {
            argv++,argc--;
            B = atoi(argv[1]);
        } else if (strcmp(argv[1], "-I") == 0) {
            argv++,argc--;
            tw.I = atoi(argv[1]);
        } else if (strcmp(argv[1], "-ntests") == 0) {
            argv++,argc--;
            ntests = atoi(argv[1]);
        } else if (strcmp(argv[1], "-seed") == 0) {
            argv++,argc--;
            seed = atol(argv[1]);
        } else if (strcmp(argv[1], "-q") == 0) {
            quiet = true;
        } else if (strcmp(argv[1], "-N") == 0) {
            tw.nocheck = true;
        } else if (strcmp(argv[1], "-T") == 0) {
            tw.timing = true;
        }
    }

    std::vector<std::tuple<unsigned long, int>> prime_powers;
    /* count primes up to B, and add prime powers as well */
    prime_info pi;
    prime_info_init(pi);
    for (unsigned long p = 2 ; p < B ; p = getprime_mt (pi)) {
        int k = 1;
        unsigned long q = p;
        for( ; q < B ; q*=p, k++) {
            if (q < tw.I) continue;
            prime_powers.emplace_back(p, k);
        }
    }
    prime_info_clear(pi);

    std::map<int, std::map<std::string, int> > stats;

    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    if (!seed) {
        seed = getpid();
        if (!quiet) printf("seeding with pid = %lu\n", seed);
    }
    gmp_randseed_ui(rstate, seed);

    std::vector<a_test> & tests = tw.tests;

    /* interesting with I=0x20000 */
    tests.emplace_back(a_test { 5, 390625, 771795, 8 });

    /* we have sporadic failures with this one */
    tests.emplace_back(a_test { 1579, 1579, 1579, 1 });

    for( ; tw.tests.size() < ntests ; ) {
        unsigned long j = gmp_urandomm_ui(rstate, prime_powers.size());
        unsigned long p;
        int k;
        std::tie(p, k) = prime_powers[j];
        unsigned long q = 1;
        for(int s = k ; s-- ; q*=p) ;
        unsigned long r = gmp_urandomm_ui(rstate, q + q / p);
        bool proj = r > q;
        if (proj)
            r = q + p * (r - q);

        a_test aa { p, q, r, k };

        tests.emplace_back(aa);

        if (!quiet) {
            std::string desc = aa.r > aa.q ? "proj" : "affine";
            if (aa.p == 2) desc += "+even";
            stats[aa.k][desc]++;
        }
    }

    {
        /* verify that all tests are in range */
        auto jt = tests.begin();
        for(auto const & aa : tests) {
            plattice L;
            bool proj = aa.r > aa.q;
            L.initial_basis(aa.q, proj ? (aa.r - aa.q) : aa.r, proj);
            if (L.check_pre_conditions(tw.I)) {
                *jt++ = aa;
            } else {
                fprintf(stderr, "skipping out-of-range test for p^k=%lu^%d r=%lu\n", aa.p, aa.k, aa.r);
            }
        }
        tests.erase(jt, tests.end());
    }


    /* declare it as volatile so that the compiler doesn't outsmart us.
     */
    volatile unsigned long dummy = 0;

    dummy += test_speed<call_production>(tw);

    dummy += test_speed<call_simplistic>(tw);
    dummy += test_speed<call_using_64bit_mul>(tw);
    dummy += test_speed<call_swapping_loop>(tw);
    dummy += test_speed<call_swapping_loop2>(tw);


    if (tw.timing) {
        /* as far as correctness is concerned, we're not interested in
         * this code. We *know* that it is buggy. But when we report
         * timings, we want to know how it fares !
         */
        dummy += test_speed<call_old_reference>(tw);
        dummy += test_speed<call_reference2>(tw);
#ifdef TEST_ASSEMBLY_CODE_DELETED_BY_5f258ce8b
        dummy += test_speed<call_reference2_asm>(tw);
#endif
    }

#ifdef HAVE_AVX512F
    dummy += test_speed<call_simd_avx512>(tw);
#endif

#ifdef HAVE_AVX2
    dummy += test_speed<call_simd_avx2>(tw);
#endif

    dummy += test_speed<call_simd<4>>(tw);
    dummy += test_speed<call_simd<2>>(tw);

    if (!quiet) {
        /* Finish with that one, but here we don't care about the
         * correctness, we only want to do the frequency table of quotients.
         */
        std::map<int, unsigned long> T;
        for(auto const & aa : tests) {
            bool proj = aa.r > aa.q;
            plattice L;
            L.initial_basis(aa.q, proj ? (aa.r - aa.q) : aa.r, proj);
            if (!L.early(tw.I))
                L.instrumented_two_legs(tw.I, T);
        }
        for(auto & kv : stats) {
            int k = kv.first;
            printf("%d:", k);
            for(auto & dn : kv.second) {
                printf(" %s:%d", dn.first.c_str(), dn.second);
            }
            printf("\n");
        }
        unsigned long sumT = 0;
        for(auto const & x : T) {
            sumT += x.second;
        }
        unsigned long cumulative = 0;
        for(auto const & x : T) {
            printf("%d: %lu (%.1f%%) (%.1f%%)\n", x.first, x.second, 100.0 * x.second / sumT, 100.0 * (cumulative += x.second) / sumT);
            if (x.first >= 16) break;
        }
    }
    gmp_randclear(rstate);
    return tw.failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
