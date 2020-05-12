#include "cado.h" // IWYU pragma: keep

#include <stdlib.h>        // for abort
#include <cstddef>         // for size_t
#include <cstdint>         // for uint8_t, uint32_t, uint64_t
#ifdef HAVE_SSSE3
#include <emmintrin.h>
#include <tmmintrin.h>  // IWYU pragma: keep
#endif
#ifdef HAVE_AVX2
#include <immintrin.h>
#endif
#ifdef HAVE_ARM_NEON
#include <arm_neon.h>
#endif

#include "intrinsics.hpp"  // for adds, andnot, loadu, set1, _and, set0
#include "macros.h"        // for ATTRIBUTE, ATTRIBUTE_ARTIFICIAL, ASSERT
#include "ularith.h"       // for ularith_submod_ul_ul

#include "las-sieve2357.hpp"

/* Specialize adds() for use in sieve2357, assuming that sieving very small
   primes does not overflow an element */
template<>
inline uint32_t ATTRIBUTE((__always_inline__)) ATTRIBUTE_ARTIFICIAL
adds<uint32_t, uint8_t>(const uint32_t a, const uint32_t b)
{
  return a + b; /* WARNING, we assume no carry between elements here! */
}

template<>
inline uint64_t ATTRIBUTE((__always_inline__)) ATTRIBUTE_ARTIFICIAL
adds<uint64_t, uint8_t>(const uint64_t a, const uint64_t b)
{
  return a + b; /* WARNING, we assume no carry between elements here! */
}

static inline unsigned long ATTRIBUTE((__always_inline__)) ATTRIBUTE_ARTIFICIAL
modsub(const unsigned long a, const unsigned long b, const unsigned long m)
{
  unsigned long r;
  ularith_submod_ul_ul(&r, a, b, m);
  return r;
}

/* This "patterns_base" is only here to placate the C++ compiler: when it
 * instantiates, it really wants to see the specializations of the static
 * data members above the code that uses them. There are various ways to
 * do so, some uglier than others, and none that really looks nice. Here
 * is just one example.
 */
template<typename ELEMTYPE>
struct patterns_base {
    static const ELEMTYPE mask2[64];
    static const ELEMTYPE mask4[64];
    static const ELEMTYPE mask8[64];
    static const ELEMTYPE mask16[64];
    static const ELEMTYPE mask32[64];
    static const ELEMTYPE mask3[64];
    static const ELEMTYPE mask5[64];
    static const ELEMTYPE mask7[64];
};

static constexpr const uint8_t ff = ~(uint8_t)0;
template<> const uint8_t patterns_base<uint8_t>::mask2[64] =
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   ff, 0, ff, 0, ff, 0, ff, 0, ff, 0, ff, 0, ff, 0, ff, 0, ff, 0, ff, 0, ff, 0, ff, 0, ff, 0, ff, 0, ff, 0, ff, 0};
template<> const uint8_t patterns_base<uint8_t>::mask4[64] =
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   ff, 0, 0, 0, ff, 0, 0, 0, ff, 0, 0, 0, ff, 0, 0, 0, ff, 0, 0, 0, ff, 0, 0, 0, ff, 0, 0, 0, ff, 0, 0, 0};
template<> const uint8_t patterns_base<uint8_t>::mask8[64] =
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   ff, 0, 0, 0, 0, 0, 0, 0, ff, 0, 0, 0, 0, 0, 0, 0, ff, 0, 0, 0, 0, 0, 0, 0, ff, 0, 0, 0, 0, 0, 0, 0};
template<> const uint8_t patterns_base<uint8_t>::mask16[64] =
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   ff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
template<> const uint8_t patterns_base<uint8_t>::mask32[64] =
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   ff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
template<> const uint8_t patterns_base<uint8_t>::mask3[64] =
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   ff, 0, 0, ff, 0, 0, ff, 0, 0, ff, 0, 0, ff, 0, 0, ff, 0, 0, ff, 0, 0, ff, 0, 0, ff, 0, 0, ff, 0, 0, ff, 0};
template<> const uint8_t patterns_base<uint8_t>::mask5[64] =
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   ff, 0, 0, 0, 0, ff, 0, 0, 0, 0, ff, 0, 0, 0, 0, ff, 0, 0, 0, 0, ff, 0, 0, 0, 0, ff, 0, 0, 0, 0, ff, 0};
template<> const uint8_t patterns_base<uint8_t>::mask7[64] =
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   ff, 0, 0, 0, 0, 0, 0, ff, 0, 0, 0, 0, 0, 0, ff, 0, 0, 0, 0, 0, 0, ff, 0, 0, 0, 0, 0, 0, ff, 0, 0, 0};

template<typename ELEMTYPE>
class patterns : public patterns_base<ELEMTYPE> {
public:
  /* A demultiplexer that returns the correct mask array for a given "stride"
     value. */
  static inline const ELEMTYPE *get_mask(unsigned int stride) {
      switch (stride) {
          case 2: return patterns_base<ELEMTYPE>::mask2;
          case 4: return patterns_base<ELEMTYPE>::mask4;
          case 8: return patterns_base<ELEMTYPE>::mask8;
          case 16: return patterns_base<ELEMTYPE>::mask16;
          case 32: return patterns_base<ELEMTYPE>::mask32;
          case 3: return patterns_base<ELEMTYPE>::mask3;
          case 5: return patterns_base<ELEMTYPE>::mask5;
          case 7: return patterns_base<ELEMTYPE>::mask7;
          default: abort();
      }
  }

  /* Return a mask of stride "stride", shifted by "shift", via an unaligned
     memory read */
  template <typename SIMDTYPE>
  static inline SIMDTYPE
  get_shifted_mask(const unsigned int stride, const fbprime_t shift)
  {
      const ELEMTYPE * mask = get_mask(stride);
      const ELEMTYPE * p = &mask[32 - shift];
      return loadu<SIMDTYPE, ELEMTYPE>(p);
  }

  /* Return a sieving pattern of stride "stride", shifted by "shift", with
     value "elem" in locations where it hits */
  template <typename SIMDTYPE>
  static inline SIMDTYPE
  get_pattern(const unsigned int stride, const fbprime_t shift, const ELEMTYPE elem)
  {
      const SIMDTYPE shifted_mask = get_shifted_mask<SIMDTYPE>(stride, shift);
      return _and<SIMDTYPE, ELEMTYPE>(shifted_mask, set1<SIMDTYPE, ELEMTYPE>(elem));
  }
};

template <typename SIMDTYPE, typename ELEMTYPE>
inline void
sieve2357<SIMDTYPE, ELEMTYPE>::sieve_odd_prime(SIMDTYPE * const result, const ELEMTYPE logp,
    const fbprime_t stride, const fbprime_t idx, const SIMDTYPE even_mask)
{
    fbprime_t offset = idx;
    for (size_t i = 0; i < stride; i++) {
        SIMDTYPE pattern = patterns<ELEMTYPE>::template get_pattern<SIMDTYPE>(stride, offset, logp);
        pattern = andnot<SIMDTYPE, ELEMTYPE>(even_mask, pattern);
        result[i] = adds<SIMDTYPE, ELEMTYPE>(result[i], pattern);
        offset = modsub(offset, sizeof(SIMDTYPE) % stride, stride);
    }
}

template <typename SIMDTYPE, typename ELEMTYPE>
inline void
sieve2357<SIMDTYPE, ELEMTYPE>::big_loop_set(SIMDTYPE * __restrict__ sievearray, const SIMDTYPE * const __restrict__ pattern235,
    const SIMDTYPE * const __restrict__ pattern7, size_t l7)
{
  size_t im15 = 0;
  for (size_t i = l7; i > 0; i--) {
    *sievearray++ = adds<SIMDTYPE, ELEMTYPE>(pattern235[im15], pattern7[0]);
    *sievearray++ = adds<SIMDTYPE, ELEMTYPE>(pattern235[im15 + 1], pattern7[1]);
    *sievearray++ = adds<SIMDTYPE, ELEMTYPE>(pattern235[im15 + 2], pattern7[2]);
    *sievearray++ = adds<SIMDTYPE, ELEMTYPE>(pattern235[im15 + 3], pattern7[3]);
    *sievearray++ = adds<SIMDTYPE, ELEMTYPE>(pattern235[im15 + 4], pattern7[4]);
    *sievearray++ = adds<SIMDTYPE, ELEMTYPE>(pattern235[im15 + 5], pattern7[5]);
    *sievearray++ = adds<SIMDTYPE, ELEMTYPE>(pattern235[im15 + 6], pattern7[6]);
    if (im15 >= 8) {im15 -= 8;} else {im15 += 7;}
  }
}

template <typename SIMDTYPE, typename ELEMTYPE>
inline void
sieve2357<SIMDTYPE, ELEMTYPE>::big_loop_add(SIMDTYPE * __restrict__ sievearray, const SIMDTYPE * const __restrict__ pattern235,
    const SIMDTYPE * const __restrict__ pattern7, size_t l7)
{
  size_t im15 = 0;
  for (size_t i = l7; i > 0; i--) {
    SIMDTYPE t;
    t = adds<SIMDTYPE, ELEMTYPE>(pattern235[im15], pattern7[0]);
    sievearray[0] = adds<SIMDTYPE, ELEMTYPE>(t, sievearray[0]);
    t = adds<SIMDTYPE, ELEMTYPE>(pattern235[im15 + 1], pattern7[1]);
    sievearray[1] = adds<SIMDTYPE, ELEMTYPE>(t, sievearray[1]);
    t = adds<SIMDTYPE, ELEMTYPE>(pattern235[im15 + 2], pattern7[2]);
    sievearray[2] = adds<SIMDTYPE, ELEMTYPE>(t, sievearray[2]);
    t = adds<SIMDTYPE, ELEMTYPE>(pattern235[im15 + 3], pattern7[3]);
    sievearray[3] = adds<SIMDTYPE, ELEMTYPE>(t, sievearray[3]);
    t = adds<SIMDTYPE, ELEMTYPE>(pattern235[im15 + 4], pattern7[4]);
    sievearray[4] = adds<SIMDTYPE, ELEMTYPE>(t, sievearray[4]);
    t = adds<SIMDTYPE, ELEMTYPE>(pattern235[im15 + 5], pattern7[5]);
    sievearray[5] = adds<SIMDTYPE, ELEMTYPE>(t, sievearray[5]);
    t = adds<SIMDTYPE, ELEMTYPE>(pattern235[im15 + 6], pattern7[6]);
    sievearray[6] = adds<SIMDTYPE, ELEMTYPE>(t, sievearray[6]);
    sievearray += 7;
    if (im15 >= 8) {im15 -= 8;} else {im15 += 7;}
  }
}

#if 0 && defined(HAVE_AVX2) && defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
inline void
big_loop_set<__m256i, uint8_t>(__m256i * __restrict__ sievearray, const __m256i * __restrict__ pattern235,
         const __m256i * const __restrict__ pattern7, const size_t l7)
{
  __m256i r, pattern235_0 = pattern235[0],
             pattern235_1 = pattern235[1],
             pattern235_2 = pattern235[2],
             pattern235_3 = pattern235[3],
             pattern235_4 = pattern235[4],
             pattern235_5 = pattern235[5],
             pattern235_6 = pattern235[6];

  const __m256i * pattern235_plus_7 = pattern235 + 7;
  const __m256i * const pattern235_plus_8 = pattern235 + 8;

  if (l7 == 0)
      return;

#define LOOP_BLOCK(i) \
               "vpaddusb %[p235_"CADO_STRINGIZE(i)"], %[p7_"CADO_STRINGIZE(i)"], %[r] \n\t"	/* p015 */ \
               "vmovdqa "CADO_STRINGIZE(i)"*32(%[pattern235]), %[p235_"CADO_STRINGIZE(i)"] \n\t"/* p23  */ \
               "vmovdqa %[r], "CADO_STRINGIZE(i)"*32(%[sievearray]) \n\t"			/* p237 p4 */

  for (size_t i = l7 - 1; i > 0; i--) {
      size_t t;
      __asm__ volatile (
               LOOP_BLOCK(0)
               "cmp %[pp8], %[pattern235] \n\t"				// p0156 carry iff pp8 > pattern235
               LOOP_BLOCK(1)
               LOOP_BLOCK(2)
               LOOP_BLOCK(3)
               LOOP_BLOCK(4)
               LOOP_BLOCK(5)
               LOOP_BLOCK(6)
               "lea 7*32(%[pattern235]), %[t] \n\t"			// p15
               "lea -8*32(%[pattern235]), %[pattern235] \n\t"		// p15
               "cmovc %[t], %[pattern235] \n\t"				// p06
               "add $7*32, %[sievearray] \n\t"				// p0156
                                                                        // p6 for dec/jnz
               : [sievearray] "+r" (sievearray), [r] "=&x" (r), [t] "=&r" (t),
                 [pattern235] "+r" (pattern235_plus_7), [p235_0] "+x" (pattern235_0),
                 [p235_1] "+x" (pattern235_1), [p235_2] "+x" (pattern235_2),
                 [p235_3] "+x" (pattern235_3), [p235_4] "+x" (pattern235_4),
                 [p235_5] "+x" (pattern235_5), [p235_6] "+x" (pattern235_6)
               : [pp8] "r" (pattern235_plus_8), [p7_0] "x" (pattern7[0]),
                 [p7_1] "x" (pattern7[1]), [p7_2] "x" (pattern7[2]),
                 [p7_3] "x" (pattern7[3]), [p7_4] "x" (pattern7[4]),
                 [p7_5] "x" (pattern7[5]), [p7_6] "x" (pattern7[6])
               : "cc" //, "memory"
      );
  }
  *sievearray++ = _mm256_adds_epu8(pattern235_0, pattern7[0]);
  *sievearray++ = _mm256_adds_epu8(pattern235_1, pattern7[1]);
  *sievearray++ = _mm256_adds_epu8(pattern235_2, pattern7[2]);
  *sievearray++ = _mm256_adds_epu8(pattern235_3, pattern7[3]);
  *sievearray++ = _mm256_adds_epu8(pattern235_4, pattern7[4]);
  *sievearray++ = _mm256_adds_epu8(pattern235_5, pattern7[5]);
  *sievearray++ = _mm256_adds_epu8(pattern235_6, pattern7[6]);
}
#endif

template <typename SIMDTYPE, typename ELEMTYPE>
inline SIMDTYPE
sieve2357<SIMDTYPE, ELEMTYPE>::sieve2(const fbprime_t q, const fbprime_t idx, const uint8_t logp)
{
    const size_t N MAYBE_UNUSED = sizeof(SIMDTYPE) / sizeof(ELEMTYPE);
    ASSERT(q <= N);
    switch (q) {
        case 2: return patterns<ELEMTYPE>::template get_pattern<SIMDTYPE>(2, idx, logp);
        case 4: return patterns<ELEMTYPE>::template get_pattern<SIMDTYPE>(4, idx, logp);
        case 8: return patterns<ELEMTYPE>::template get_pattern<SIMDTYPE>(8, idx, logp);
        case 16: return patterns<ELEMTYPE>::template get_pattern<SIMDTYPE>(16, idx, logp);
        case 32: return patterns<ELEMTYPE>::template get_pattern<SIMDTYPE>(32, idx, logp);
        default: abort();
    }
}

/* If skip_mod_2 == 0: returns {0, 0, 0, ..., 0}
   If skip_mod_2 == 1: returns {0xff, 0, 0xff, 0, ..., 0xff, 0}
       (0xff at even indices, 0 at odd. Indices are 0-based.)
   If skip_mod_2 == 2: returns {0, 0xff, 0, 0xff, ..., 0, 0xff}
       (0 at even indices, 0xff at odd)
   This mask will be used in an andnot operation. In case 0, the andnot
   zeroes out nothing, in case 1 it zeroes out elements at even indices,
   and in case 2 it zeroes out elements at odd indices. */
template <typename SIMDTYPE, typename ELEMTYPE>
inline SIMDTYPE
sieve2357<SIMDTYPE, ELEMTYPE>::get_even_mask(const int skip_mod_2)
{
    const bool only_odd = (skip_mod_2 & 1) != 0;
    const bool only_even = (skip_mod_2 & 2) != 0;
    ASSERT (!(only_odd && only_even));
    const ELEMTYPE *mask2 = patterns<ELEMTYPE>::get_mask(2);
    const ELEMTYPE *even_mask;
    if (only_odd) {
        even_mask = mask2 + 32; /* ff, 0, ff, 0, ... */
    } else if (only_even) {
        even_mask = mask2 + 31; /* 0, ff, 0, ff, ... */
    } else {
        even_mask = mask2; /* 0, 0, 0, 0, ... */
    }
    return loadu<SIMDTYPE, ELEMTYPE>(even_mask);
}

template <typename SIMDTYPE, typename ELEMTYPE>
void
sieve2357<SIMDTYPE, ELEMTYPE>::sieve(SIMDTYPE * const sievearray, const size_t arraylen, const prime_t *primes,
      const int skip_mod_2, const int update_operation, where_am_I & w MAYBE_UNUSED)
{
  const SIMDTYPE zero = set0<SIMDTYPE, ELEMTYPE>();
  const SIMDTYPE even_mask = get_even_mask(skip_mod_2);
  const size_t N = sizeof(SIMDTYPE) / sizeof(ELEMTYPE);
  SIMDTYPE pattern2 = zero;

  /* Sieve projective primes with q=1 */
  for ( ; primes->q == 1 ; primes++) {
    pattern2 = adds<SIMDTYPE, ELEMTYPE>(pattern2, set1<SIMDTYPE, ELEMTYPE>(primes->logp));
  }

  /* Sieve powers of 2 */
  for ( ; primes->q != 0 && primes->q % 2 == 0 ; primes++) {
    pattern2 = adds<SIMDTYPE, ELEMTYPE>(pattern2, sieve2(primes->q, primes->idx, primes->logp));
  }

  pattern2 = andnot<SIMDTYPE, ELEMTYPE>(even_mask, pattern2);

  /* This triggers a warning on debian-8-arm32
   *
   *    /home/ci/jenkins/workspace/master/compile-debian-8-arm32/sieve/las-sieve2357.cpp:310:12: warning: requested alignment 16 is larger than 8 [-Wattributes]
   SIMDTYPE pattern23[3] = {pattern2, pattern2, pattern2};
   */
#ifdef HAVE_ALIGNAS
  alignas(sizeof(SIMDTYPE)) 
#endif
  SIMDTYPE pattern23[3] = {pattern2, pattern2, pattern2};
  for ( ; primes->q == 3 ; primes++) {
    sieve_odd_prime(pattern23, primes->logp, 3, primes->idx, even_mask);
  }

#ifdef HAVE_ALIGNAS
  alignas(sizeof(SIMDTYPE)) 
#endif
  SIMDTYPE pattern235[15 + 6] = {
    pattern23[0], pattern23[1], pattern23[2],
    pattern23[0], pattern23[1], pattern23[2],
    pattern23[0], pattern23[1], pattern23[2],
    pattern23[0], pattern23[1], pattern23[2],
    pattern23[0], pattern23[1], pattern23[2],
    pattern23[0], pattern23[1], pattern23[2],
    pattern23[0], pattern23[1], pattern23[2],
  };
#ifdef HAVE_ALIGNAS
  alignas(sizeof(SIMDTYPE)) 
#endif
  SIMDTYPE pattern5[5] = {zero, zero, zero, zero, zero};
  for ( ; primes->q == 5 ; primes++) {
    sieve_odd_prime(pattern5, primes->logp, 5, primes->idx, even_mask);
  }
 
  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < 5; j++) {
      pattern235[i*5 + j] = adds<SIMDTYPE, ELEMTYPE>(pattern235[i*5 + j], pattern5[j]);
    }
  }

  /* Pad the list at the end so that we don't have to check for mod 15
     wrap-around after each increment. */
  pattern235[15 + 0] = pattern235[0];
  pattern235[15 + 1] = pattern235[1];
  pattern235[15 + 2] = pattern235[2];
  pattern235[15 + 3] = pattern235[3];
  pattern235[15 + 4] = pattern235[4];
  pattern235[15 + 5] = pattern235[5];

#ifdef HAVE_ALIGNAS
  alignas(sizeof(SIMDTYPE))
#endif
  SIMDTYPE pattern7[7] = {zero, zero, zero, zero, zero, zero, zero};
  for ( ; primes->q == 7 ; primes++) {
    sieve_odd_prime(pattern7, primes->logp, 7, primes->idx, even_mask);
  }

  ASSERT_ALWAYS(primes->q == 0);

  const size_t l7 = arraylen / 7 / N;
  
  if (update_operation == update_set) {
    big_loop_set(sievearray, pattern235, pattern7, l7);
    for (size_t i = l7 * 7; i < arraylen / N; i++) {
      sievearray[i] = adds<SIMDTYPE, ELEMTYPE>(pattern235[i % 15], pattern7[i % 7]);
    }
  } else if (update_operation == update_add) {
    big_loop_add(sievearray, pattern235, pattern7, l7);
    for (size_t i = l7 * 7; i < arraylen / N; i++) {
      SIMDTYPE t;
      t = adds<SIMDTYPE, ELEMTYPE>(pattern235[i % 15], pattern7[i % 7]);
      sievearray[i] = adds<SIMDTYPE, ELEMTYPE>(t, sievearray[i]);
    }
  } else
    abort();
}

template class patterns<uint8_t>;

#if GNUC_VERSION_ATLEAST(6,1,0)
/* https://gcc.gnu.org/bugzilla/show_bug.cgi?id=69884 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
#endif
template class sieve2357<uint32_t, uint8_t>;
template class sieve2357<uint64_t, uint8_t>;
#ifdef HAVE_SSSE3
template class sieve2357<__m128i, uint8_t>;
#endif
#ifdef HAVE_AVX2
template class sieve2357<__m256i, uint8_t>;
#endif
#ifdef HAVE_ARM_NEON
template class sieve2357<uint8x16_t, uint8_t>;
#endif
#if GNUC_VERSION_ATLEAST(6,1,0)
#pragma GCC diagnostic pop
#endif
