#ifndef UTILS_MPFR_AUXX_HPP_
#define UTILS_MPFR_AUXX_HPP_

#include <cstdint>
#include <cmath>

#include <gmp.h>
#include <mpfr.h>

#include "mpfr_aux.h"
#include "utils_cxx.hpp"

/* A C++ wrapper around functions in mpfr_aux.h. The functions here delegate
 * to the proper C function depending on the argument type. We keep them
 * out of global name space to avoid obscuring errors.
 *
 * Note that we use the cado_mpfr_ prefix (inside the mpfr_auxx
 * namespace!) because most mpfr names are macros.
 */

namespace mpfr_auxx
{

#define MPFR_AUXX_DEFINE_FUNC_aXr_INT(OP)                               \
    static inline void cado_mpfr_##OP(mpfr_ptr a, mpfr_srcptr b,        \
                                      mpfr_rnd_t rnd)                   \
    {                                                                   \
        mpfr_##OP(a, b, rnd);                                           \
    }                                                                   \
                                                                        \
    template <typename T>                                               \
    static inline void                                                  \
    cado_mpfr_##OP(mpfr_ptr a, T const b, mpfr_rnd_t rnd)               \
    requires cado::converts_via<T, long>                                \
    {                                                                   \
        mpfr_##OP##_si(a, b, rnd);                                      \
    }                                                                   \
                                                                        \
    template <typename T>                                               \
    static inline void                                                  \
    cado_mpfr_##OP(mpfr_ptr a, T const b, mpfr_rnd_t rnd)               \
    requires cado::converts_via<T, unsigned long>                       \
    {                                                                   \
        mpfr_##OP##_ui(a, b, rnd);                                      \
    }                                                                   \
                                                                        \
    template <typename T>                                               \
    static inline void                                                  \
    cado_mpfr_##OP(mpfr_ptr a, T const b, mpfr_rnd_t rnd)               \
        requires(!cado::converts_via<T, long> &&                        \
                  cado::converts_via<T, int64_t>)                       \
    {                                                                   \
        mpfr_##OP##_int64(a, b, rnd);                                   \
    }                                                                   \
                                                                        \
    template <typename T>                                               \
    static inline void                                                  \
    cado_mpfr_##OP(mpfr_ptr a, T const b, mpfr_rnd_t rnd)               \
        requires(!cado::converts_via<T, unsigned long> &&               \
                  cado::converts_via<T, uint64_t>)                      \
    {                                                                   \
        mpfr_##OP##_uint64(a, b, rnd);                                  \
    }                                                                   \
                                                                        \
    static inline void                                                  \
    cado_mpfr_##OP(mpfr_ptr a, mpz_srcptr b, mpfr_rnd_t rnd)            \
    {                                                                   \
        mpfr_##OP##_z(a, b, rnd);                                       \
    }

#define MPFR_AUXX_DEFINE_FUNC_aXr_FP_META(OP, TYP, SUF)                 \
    static inline void                                                  \
    cado_mpfr_##OP(mpfr_ptr a, TYP const b, mpfr_rnd_t rnd)             \
    {                                                                   \
        mpfr_##OP##SUF(a, b, rnd);                                      \
    }

    /* same as mpfr, except that we preserve the sign bit of NAN */
#define MPFR_AUXX_DEFINE_FUNC_aXr_FP_META_KEEP_NAN_SIGN(OP, TYP, SUF)   \
    static inline void                                                  \
    cado_mpfr_##OP(mpfr_ptr a, TYP const b, mpfr_rnd_t rnd)             \
    {                                                                   \
        mpfr_##OP##SUF(a, b, rnd);                                      \
        if (std::isnan(b))                                              \
            mpfr_setsign(a, a, std::signbit(b), MPFR_RNDN);             \
    }

#define MPFR_AUXX_DEFINE_FUNC_aXr(OP)                                   \
        MPFR_AUXX_DEFINE_FUNC_aXr_INT(OP)                               \
        MPFR_AUXX_DEFINE_FUNC_aXr_FP_META_KEEP_NAN_SIGN(OP, float, _d)  \
        MPFR_AUXX_DEFINE_FUNC_aXr_FP_META_KEEP_NAN_SIGN(OP, double, _d) \
        MPFR_AUXX_DEFINE_FUNC_aXr_FP_META_KEEP_NAN_SIGN(OP, long double, _ld)

MPFR_AUXX_DEFINE_FUNC_aXr(init_set)
MPFR_AUXX_DEFINE_FUNC_aXr(set)


/*****************************************************************/
#define MPFR_AUXX_DEFINE_FUNC_tX_INT(OP)				\
    static inline int cado_mpfr_##OP(mpfr_srcptr a, mpfr_srcptr b)	\
    {									\
        return mpfr_##OP(a, b);						\
    }									\
    template <typename T>						\
    static inline int							\
    cado_mpfr_##OP(mpfr_srcptr a, T const b)				\
        requires cado::converts_via<T, long>				\
    {									\
        return mpfr_##OP##_si(a, b);					\
    }									\
    template <typename T>						\
    static inline int							\
    cado_mpfr_##OP(mpfr_srcptr a, T const b)				\
        requires cado::converts_via<T, unsigned long>			\
    {									\
        return mpfr_##OP##_ui(a, b);					\
    }									\
    template <typename T>						\
    static inline int							\
    cado_mpfr_##OP(mpfr_srcptr a, T const b)				\
        requires(!cado::converts_via<T, long> &&			\
                  cado::converts_via<T, int64_t>)			\
    {									\
        return mpfr_##OP##_int64(a, b);					\
    }									\
    template <typename T>						\
    static inline int							\
    cado_mpfr_##OP(mpfr_srcptr a, T const b)				\
        requires(!cado::converts_via<T, unsigned long> &&		\
                  cado::converts_via<T, uint64_t>)			\
    {									\
        return mpfr_##OP##_uint64(a, b);				\
    }

#define MPFR_AUXX_DEFINE_FUNC_tX_FP_META(OP, TYP, SUF)			\
    static inline int							\
    cado_mpfr_##OP(mpfr_srcptr a, TYP b)				\
    {									\
        return mpfr_##OP##SUF(a, b);					\
    }

#define MPFR_AUXX_DEFINE_FUNC_tX(OP)                                    \
    MPFR_AUXX_DEFINE_FUNC_tX_INT(OP)                                    \
    MPFR_AUXX_DEFINE_FUNC_tX_FP_META(OP, float, _d)			\
    MPFR_AUXX_DEFINE_FUNC_tX_FP_META(OP, double, _d)			\
    MPFR_AUXX_DEFINE_FUNC_tX_FP_META(OP, long double, _ld)

MPFR_AUXX_DEFINE_FUNC_tX(cmp)

static inline int
cado_mpfr_cmp(mpfr_srcptr a, mpz_srcptr b)
{
    return mpfr_cmp_z(a, b);
}

/*****************************************************************/

#define MPFR_AUXX_DEFINE_FUNC_abXr_INT(OP)                              \
    static inline int cado_mpfr_##OP(mpfr_ptr a,                        \
                                     mpfr_srcptr b, mpfr_srcptr c,      \
                                     mpfr_rnd_t rnd)                    \
    {                                                                   \
        return ::mpfr_##OP(a, b, c, rnd);                               \
    }                                                                   \
    template <typename T>                                               \
    static inline int                                                   \
    cado_mpfr_##OP(mpfr_ptr a, mpfr_srcptr b, const T c,                \
                   mpfr_rnd_t rnd)                                      \
    requires cado::converts_via<T, unsigned long>                       \
    {                                                                   \
        return mpfr_##OP##_ui(a, b, c, rnd);                            \
    }                                                                   \
    template <typename T>                                               \
    static inline int cado_mpfr_##OP(mpfr_ptr a, mpfr_srcptr b,         \
                                     const T c, mpfr_rnd_t rnd)         \
    requires(!cado::converts_via<T, unsigned long> &&                   \
              cado::converts_via<T, uint64_t>)                          \
    {                                                                   \
        return mpfr_##OP##_uint64(a, b, c, rnd);                        \
    }                                                                   \
    template <typename T>                                               \
    static inline int                                                   \
    cado_mpfr_##OP(mpfr_ptr a, mpfr_srcptr b, const T c, mpfr_rnd_t rnd)\
    requires cado::converts_via<T, long>                             \
    {                                                                   \
        return mpfr_##OP##_si(a, b, c, rnd);                            \
    }                                                                   \
    template <typename T>                                               \
    static inline int cado_mpfr_##OP(mpfr_ptr a, mpfr_srcptr b,         \
                                     const T c, mpfr_rnd_t rnd)         \
    requires(!cado::converts_via<T, long> &&                            \
              cado::converts_via<T, int64_t>)                           \
    {                                                                   \
        return mpfr_##OP##_int64(a, b, c, rnd);                         \
    }
#define MPFR_AUXX_DEFINE_FUNC_abXr_FP_META(OP, TYP, SUF)                \
    static inline int                                                   \
    cado_mpfr_##OP(mpfr_ptr a, mpfr_srcptr b, const TYP c,            \
                   mpfr_rnd_t rnd)                                      \
    {                                                                   \
        return mpfr_##OP##SUF(a, b, c, rnd);                             \
    }

#define MPFR_AUXX_DEFINE_FUNC_aXbr_INT(OP, FIXUP)                           \
    template <typename T>                                               \
    static inline int                                                   \
    cado_mpfr_##OP(mpfr_ptr a, const T b, mpfr_srcptr c,                \
                   mpfr_rnd_t rnd)                                      \
    requires cado::converts_via<T, unsigned long>                       \
    {                                                                   \
        int r = mpfr_##OP##_ui(a, c, b, rnd);                           \
        FIXUP;                                                          \
        return r;                                                       \
    }                                                                   \
    template <typename T>                                               \
    static inline int cado_mpfr_##OP(mpfr_ptr a, const T b,             \
                                     mpfr_srcptr c, mpfr_rnd_t rnd)     \
    requires(!cado::converts_via<T, unsigned long> &&                   \
              cado::converts_via<T, uint64_t>)                          \
    {                                                                   \
        int r = mpfr_##OP##_uint64(a, c, b, rnd);                       \
        FIXUP;                                                          \
        return r;                                                       \
    }                                                                   \
    template <typename T>                                               \
    static inline int                                                   \
    cado_mpfr_##OP(mpfr_ptr a, const T b, mpfr_srcptr c,                \
                   mpfr_rnd_t rnd)                                      \
    requires cado::integral_fits_v<T, long>                             \
    {                                                                   \
        int r = mpfr_##OP##_si(a, c, b, rnd);                           \
        FIXUP;                                                          \
        return r;                                                       \
    }                                                                   \
    template <typename T>                                               \
    static inline int cado_mpfr_##OP(mpfr_ptr a, const T b,             \
                                     mpfr_srcptr c, mpfr_rnd_t rnd)     \
    requires(!cado::converts_via<T, long> &&                            \
              cado::converts_via<T, int64_t>)                           \
    {                                                                   \
        int r = mpfr_##OP##_int64(a, c, b, rnd);                        \
        FIXUP;                                                          \
        return r;                                                       \
    }

#define MPFR_AUXX_DEFINE_FUNC_aXbr_FP_META(OP, FIXUP, TYP, SUF)         \
    static inline int                                                   \
    cado_mpfr_##OP(mpfr_ptr a, const TYP b, mpfr_srcptr c,              \
                   mpfr_rnd_t rnd)                                      \
    {                                                                   \
        int r = mpfr_##OP##SUF(a, c, b, rnd);                           \
        FIXUP;                                                          \
        return r;                                                       \
    }

#define MPFR_AUXX_DEFINE_FUNC_abXr(OP)                                  \
    MPFR_AUXX_DEFINE_FUNC_abXr_INT(OP)                                  \
    MPFR_AUXX_DEFINE_FUNC_abXr_FP_META(OP, float, _d)                   \
    MPFR_AUXX_DEFINE_FUNC_abXr_FP_META(OP, double, _d)                   \
    MPFR_AUXX_DEFINE_FUNC_abXr_FP_META(OP, long double, _ld)

#define MPFR_AUXX_DEFINE_FUNC_aXbr(OP, FIXUP)                           \
    MPFR_AUXX_DEFINE_FUNC_aXbr_INT(OP, FIXUP)                           \
    MPFR_AUXX_DEFINE_FUNC_aXbr_FP_META(OP, FIXUP, float, _d)            \
    MPFR_AUXX_DEFINE_FUNC_aXbr_FP_META(OP, FIXUP, double, _d)           \
    MPFR_AUXX_DEFINE_FUNC_aXbr_FP_META(OP, FIXUP, long double, _ld)

MPFR_AUXX_DEFINE_FUNC_abXr(add)
MPFR_AUXX_DEFINE_FUNC_abXr(sub)
MPFR_AUXX_DEFINE_FUNC_abXr(mul)
MPFR_AUXX_DEFINE_FUNC_abXr(addmul)
MPFR_AUXX_DEFINE_FUNC_abXr(submul)
MPFR_AUXX_DEFINE_FUNC_abXr(div)
MPFR_AUXX_DEFINE_FUNC_abXr(remainder)

/* Add these for convenience only */
MPFR_AUXX_DEFINE_FUNC_aXbr(sub, mpfr_neg(a, a, rnd); r = -r)
MPFR_AUXX_DEFINE_FUNC_aXbr(div, mpfr_pow_si(a, a, -1, rnd); r = -r)
MPFR_AUXX_DEFINE_FUNC_aXbr(add, /* no fixup for commutative op */)
MPFR_AUXX_DEFINE_FUNC_aXbr(mul, /* no fixup for commutative op */)

} /* namespace mpfr_auxx */

#endif /* UTILS_MPFR_AUXX_HPP_ */
