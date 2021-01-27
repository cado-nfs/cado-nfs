#ifndef BBLAS_SIMD_HPP_
#define BBLAS_SIMD_HPP_

// IWYU pragma: private, include "bblas.hpp"
// IWYU pragma: friend ".*/bblas.*"

#include "bblas.hpp"

#ifdef HAVE_MMX
#include <mmintrin.h>
#endif

#if defined(HAVE_SSE2) && ULONG_BITS == 64
#include <emmintrin.h>
/*  helper macros for sse-2. Copied from gf2x */
/* {{{ _mm_cvtsi64_m64 is not consistent across compiler versions... */
#if defined(__GNUC__) && __GNUC__ == 4 &&__GNUC_MINOR__ == 1
#define _cado_mm_cvtsi64_m64(u) _mm_cvtsi64x_m64((u))
#else
#define _cado_mm_cvtsi64_m64(u) _mm_cvtsi64_m64((u))
#endif
/* }}} */
/* {{{ _cado_mm_setr_epi64 _m128i from 2 int64_t's */
#define _cado_mm_setr_epi64(lo, hi)                      		\
    _mm_setr_epi64(                                      		\
            _cado_mm_cvtsi64_m64((int64_t) (lo)),       		\
            _cado_mm_cvtsi64_m64((int64_t) (hi))        		\
        )
/* }}} */
/* {{{ _cado_mm_set1_epi64 _m128i from 1 int64_t's */
#define _cado_mm_set1_epi64(u) _mm_set1_epi64( _cado_mm_cvtsi64_m64((int64_t) (u)))
/* }}} */
/* {{{ _cado_mm_setr_epi64_c _m128i from 2 int64_t CONSTANTS (and try to get suffix right) */
#define _cado_mm_setr_epi64_c(lo, hi)                    		\
    _mm_setr_epi64(                                      		\
            _cado_mm_cvtsi64_m64(INT64_C(lo)),          		\
            _cado_mm_cvtsi64_m64(INT64_C(hi))           		\
        )
/* }}} */
/* {{{ _cado_mm_set1_epi64_c _m128i from 1 int64_t CONSTANT (and try to get suffix right) */
#define _cado_mm_set1_epi64_c(u) _mm_set1_epi64( _cado_mm_cvtsi64_m64(INT64_C(u)))
/* }}} */
/* {{{ same for 32-bits (which, for some, have SSE-2) */
#define _cado_mm_setr_epi32(a0, a1, a2, a3)				\
    _mm_setr_epi32(                                      		\
            (int32_t) (a0),						\
            (int32_t) (a1),						\
            (int32_t) (a2),						\
            (int32_t) (a3)						\
            )
#define _cado_mm_set1_epi32(u) _mm_set1_epi32( (int32_t) (u))
#define _cado_mm_setr_epi32_c(a0, a1, a2, a3)				\
    _mm_setr_epi32(                                      		\
            (INT32_C(a0)),          					\
            (INT32_C(a1)),           					\
            (INT32_C(a2)),          					\
            (INT32_C(a3))           					\
        )
#define _cado_mm_set1_epi32_c(u) _mm_set1_epi32(INT32_C(u))
/* }}} */
/*  */
#endif
#ifdef HAVE_SSE41
#include <smmintrin.h>
#endif
#ifdef HAVE_AVX2
#include <immintrin.h>
#endif


#endif	/* BBLAS_SIMD_HPP_ */
