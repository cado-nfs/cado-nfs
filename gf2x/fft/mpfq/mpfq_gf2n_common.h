#ifndef MPFQ_GF2N_COMMON_H_
#define MPFQ_GF2N_COMMON_H_

#include "mpfq.h"

#ifdef __cplusplus
extern "C" {
#endif

// type for small field char 2
typedef struct {
    int io_type;
} mpfq_2_field_struct;

typedef mpfq_2_field_struct mpfq_2_field[1];
typedef mpfq_2_field_struct * mpfq_2_dst_field;

#if GNUC_VERSION(4,3,0) || GNUC_VERSION(4,3,1)
#warning "Your GCC version is buggy. Binary fields may fail randomly"
#endif

/* include some helpers for SSE intrinsics */
#ifdef  __SSE2__        /* maybe a configure-time check would be better */
#include <emmintrin.h>
#include <stdint.h>
#if defined(__GNUC__) && __GNUC__ == 4 &&__GNUC_MINOR__ == 1
#define _mpfq_mm_cvtsi64_m64(u) _mm_cvtsi64x_m64((u))
#else
#define _mpfq_mm_cvtsi64_m64(u) _mm_cvtsi64_m64((u))
#endif
/* _m128i from 2 int64_t's */
#define _mpfq_mm_setr_epi64(lo, hi)                                     \
    _mm_setr_epi64(                                                     \
            _mpfq_mm_cvtsi64_m64((int64_t) (lo)),                       \
            _mpfq_mm_cvtsi64_m64((int64_t) (hi))                        \
        )
/* _m128i from 1 int64_t's */
#define _mpfq_mm_set1_epi64(u) _mm_set1_epi64( _mpfq_mm_cvtsi64_m64((int64_t) (u)))
/* _m128i from 2 int64_t CONSTANTS (and try to get suffix right) */
#define _mpfq_mm_setr_epi64_c(lo, hi)                                   \
    _mm_setr_epi64(                                                     \
            _mpfq_mm_cvtsi64_m64(INT64_C(lo)),                          \
            _mpfq_mm_cvtsi64_m64(INT64_C(hi))                           \
        )
/* _m128i from 1 int64_t CONSTANT (and try to get suffix right) */
#define _mpfq_mm_set1_epi64_c(u) _mm_set1_epi64( _mpfq_mm_cvtsi64_m64(INT64_C(u)))
/* and same for 32-bits (which, for some, have SSE-2) */
#define _mpfq_mm_setr_epi32(a0, a1, a2, a3)                             \
    _mm_setr_epi32(                                                     \
            (int32_t) (a0),                                             \
            (int32_t) (a1),                                             \
            (int32_t) (a2),                                             \
            (int32_t) (a3)                                              \
            )
#define _mpfq_mm_set1_epi32(u) _mm_set1_epi32( (int32_t) (u))
#define _mpfq_mm_setr_epi32_c(a0, a1, a2, a3)                           \
    _mm_setr_epi32(                                                     \
            (INT32_C(a0)),                                              \
            (INT32_C(a1)),                                              \
            (INT32_C(a2)),                                              \
            (INT32_C(a3))                                               \
        )
#define _mpfq_mm_set1_epi32_c(u) _mm_set1_epi32(INT32_C(u))
#endif  /* __SSE2__ */

#ifdef __cplusplus
}
#endif

#endif	/* MPFQ_GF2N_COMMON_H_ */
