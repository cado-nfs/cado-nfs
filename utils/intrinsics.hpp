#ifndef INTRINSICS_HPP_
#define INTRINSICS_HPP_

/* A header file that make SIMD intrinsics that perform the same operation
   available under the same name, templated by the SIMD word type and the
   element type. */

#include <limits>
#include <macros.h>

#ifdef HAVE_SSE
#include <xmmintrin.h>
#endif

#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif

#ifdef HAVE_SSE3
#include <pmmintrin.h>
#endif

#ifdef HAVE_SSSE3
#include <tmmintrin.h>
#endif

#ifdef HAVE_SSE41
#include <smmintrin.h>
#endif

#ifdef HAVE_SSE42
#include <nmmintrin.h>
#endif

#ifdef HAVE_AVX
/* Also covers AVX2 */
#include <immintrin.h>
#endif
#ifdef HAVE_ARM_NEON
#include <arm_neon.h>
#endif

template <typename SIMDTYPE, typename ELEMTYPE>
static SIMDTYPE set0();

template <typename SIMDTYPE, typename ELEMTYPE>
static SIMDTYPE set1(ELEMTYPE);

template <typename SIMDTYPE, typename ELEMTYPE>
static SIMDTYPE adds(SIMDTYPE, SIMDTYPE);

template <typename SIMDTYPE, typename ELEMTYPE>
static SIMDTYPE _and(SIMDTYPE, SIMDTYPE);

template <typename SIMDTYPE, typename ELEMTYPE>
static SIMDTYPE andnot(SIMDTYPE);

template <typename SIMDTYPE, typename ELEMTYPE>
static SIMDTYPE loadu(const ELEMTYPE *);

template<typename SIMDTYPE, typename ELEMTYPE>
inline SIMDTYPE ATTRIBUTE_ALWAYS_INLINE ATTRIBUTE_ARTIFICIAL
adds(const SIMDTYPE a, const SIMDTYPE b)
{
  const size_t N = sizeof(SIMDTYPE) / sizeof(ELEMTYPE);
  SIMDTYPE r;
  ELEMTYPE *ap = (ELEMTYPE *)&a;
  ELEMTYPE *bp = (ELEMTYPE *)&b;
  ELEMTYPE *rp = (ELEMTYPE *)&r;
  /* This is slow, but correct */
  for (size_t i = 0; i < N; i++) {
    if (std::numeric_limits<ELEMTYPE>::max() - ap[i] < bp[i])
      rp[i] = std::numeric_limits<ELEMTYPE>::max();
    else
      rp[i] = ap[i] + bp[i];
  }
  return r;
}

template<typename SIMDTYPE, typename ELEMTYPE>
inline SIMDTYPE ATTRIBUTE_ALWAYS_INLINE ATTRIBUTE_ARTIFICIAL
set0()
{
  return 0;
}

template<>
inline uint32_t ATTRIBUTE_ALWAYS_INLINE ATTRIBUTE_ARTIFICIAL
set1<uint32_t, uint8_t>(const uint8_t c)
{
  uint32_t r = (uint32_t) c;
  r += r << 8;
  r += r << 16;
  return r;
}

template<>
inline uint64_t ATTRIBUTE_ALWAYS_INLINE ATTRIBUTE_ARTIFICIAL
set1<uint64_t, uint8_t>(const uint8_t c)
{
  uint64_t r = (uint32_t) c;
  r += r << 8;
  r += r << 16;
  r += r << 32;
  return r;
}

template<typename SIMDTYPE, typename ELEMTYPE>
inline SIMDTYPE ATTRIBUTE_ALWAYS_INLINE ATTRIBUTE_ARTIFICIAL
_and(const SIMDTYPE a, const SIMDTYPE b)
{
  return a & b;
}

template<typename SIMDTYPE, typename ELEMTYPE>
inline SIMDTYPE ATTRIBUTE_ALWAYS_INLINE ATTRIBUTE_ARTIFICIAL
andnot(const SIMDTYPE a, const SIMDTYPE b)
{
  return (~a) & b;
}

template<typename SIMDTYPE, typename ELEMTYPE>
inline SIMDTYPE ATTRIBUTE_ALWAYS_INLINE ATTRIBUTE_ARTIFICIAL
loadu(const ELEMTYPE *p)
{
  return * (const SIMDTYPE *)p;
}

#ifdef HAVE_SSE2

template<>
inline __m128i ATTRIBUTE_ALWAYS_INLINE ATTRIBUTE_ARTIFICIAL
set0<__m128i, uint8_t>()
{
  return _mm_setzero_si128 ();
}

template<>
inline __m128i ATTRIBUTE_ALWAYS_INLINE ATTRIBUTE_ARTIFICIAL
set1<__m128i, uint8_t>(const uint8_t c)
{
  return _mm_set1_epi8(c);
}


template<>
inline __m128i ATTRIBUTE_ALWAYS_INLINE ATTRIBUTE_ARTIFICIAL
adds<__m128i, uint8_t>(const __m128i a, const __m128i b)
{
  return _mm_adds_epu8(a, b);
}

template<>
inline __m128i ATTRIBUTE_ALWAYS_INLINE ATTRIBUTE_ARTIFICIAL
_and<__m128i, uint8_t>(const __m128i a, const __m128i b)
{
  return _mm_and_si128(a, b);
}

template<>
inline __m128i ATTRIBUTE_ALWAYS_INLINE ATTRIBUTE_ARTIFICIAL
andnot<__m128i, uint8_t>(const __m128i a, const __m128i b)
{
  return _mm_andnot_si128(a, b);
}

template<>
inline __m128i ATTRIBUTE_ALWAYS_INLINE ATTRIBUTE_ARTIFICIAL
loadu<__m128i, uint8_t>(const uint8_t *p)
{
  return _mm_loadu_si128((const __m128i *) p);
}

#endif

#ifdef HAVE_AVX

template<>
inline __m256i ATTRIBUTE_ALWAYS_INLINE ATTRIBUTE_ARTIFICIAL
set0<__m256i, uint8_t>()
{
  return _mm256_setzero_si256 ();
}

template<>
inline __m256i ATTRIBUTE_ALWAYS_INLINE ATTRIBUTE_ARTIFICIAL
loadu<__m256i>(const uint8_t *p)
{
  return _mm256_loadu_si256((const __m256i *) p);
}

#endif

#ifdef HAVE_AVX2

template<>
inline __m256i ATTRIBUTE_ALWAYS_INLINE ATTRIBUTE_ARTIFICIAL
set1<__m256i, uint8_t>(const uint8_t c)
{
  __m128i t = _mm_cvtsi64_si128(c);
  return _mm256_broadcastb_epi8(t);
}

template<>
inline __m256i ATTRIBUTE_ALWAYS_INLINE ATTRIBUTE_ARTIFICIAL
adds<__m256i, uint8_t>(const __m256i a, const __m256i b)
{
  return _mm256_adds_epu8(a, b);
}

template<>
inline __m256i ATTRIBUTE_ALWAYS_INLINE ATTRIBUTE_ARTIFICIAL
_and<__m256i, uint8_t>(const __m256i a, const __m256i b)
{
  return _mm256_and_si256(a, b);
}

template<>
inline __m256i ATTRIBUTE_ALWAYS_INLINE ATTRIBUTE_ARTIFICIAL
andnot<__m256i, uint8_t>(const __m256i a, const __m256i b)
{
  return _mm256_andnot_si256(a, b);
}

#endif

#ifdef HAVE_ARM_NEON

template<>
inline uint8x16_t ATTRIBUTE_ALWAYS_INLINE ATTRIBUTE_ARTIFICIAL
set0<uint8x16_t, uint8_t>()
{
  return vdupq_n_u8(0);
}

template<>
inline uint8x16_t ATTRIBUTE_ALWAYS_INLINE ATTRIBUTE_ARTIFICIAL
set1<uint8x16_t, uint8_t>(const uint8_t c)
{
   return vdupq_n_u8(c);
}

template<>
inline uint8x16_t ATTRIBUTE_ALWAYS_INLINE ATTRIBUTE_ARTIFICIAL
adds<uint8x16_t, uint8_t>(const uint8x16_t a, const uint8x16_t b)
{
  return vqaddq_u8(a, b);
}

template<>
inline uint8x16_t ATTRIBUTE_ALWAYS_INLINE ATTRIBUTE_ARTIFICIAL
_and<uint8x16_t, uint8_t>(const uint8x16_t a, const uint8x16_t b)
{
  return vandq_u8(a, b);
}

template<>
inline uint8x16_t ATTRIBUTE_ALWAYS_INLINE ATTRIBUTE_ARTIFICIAL
andnot<uint8x16_t, uint8_t>(const uint8x16_t a, const uint8x16_t b)
{
  return vandq_u8(vmvnq_u8(a), b);
}

template<>
inline uint8x16_t ATTRIBUTE_ALWAYS_INLINE ATTRIBUTE_ARTIFICIAL
loadu<uint8x16_t>(const uint8_t *p)
{
  return vld1q_u8(p);
}
#endif

#endif
