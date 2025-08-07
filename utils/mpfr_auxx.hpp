#ifndef UTILS_MPFR_AUXX_HPP_
#define UTILS_MPFR_AUXX_HPP_

#include <cstdint>

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

#define MPFR_AUXX_DEFINE_FUNC2(DTYPE, OP)                               \
    static inline void cado_mpfr_##OP(DTYPE a, mpfr_srcptr b,           \
                                      mpfr_rnd_t rnd)                   \
    {                                                                   \
        mpfr_##OP(a, b, rnd);                                           \
    }                                                                   \
                                                                        \
    template <typename T>                                               \
    static inline void                                                  \
    cado_mpfr_##OP(DTYPE a, T const b, mpfr_rnd_t rnd)                  \
    requires cado::converts_via<T, long>                                \
    {                                                                   \
        mpfr_##OP##_si(a, b, rnd);                                      \
    }                                                                   \
                                                                        \
    template <typename T>                                               \
    static inline void                                                  \
    cado_mpfr_##OP(DTYPE a, T const b, mpfr_rnd_t rnd)                  \
    requires cado::converts_via<T, unsigned long>                       \
    {                                                                   \
        mpfr_##OP##_ui(a, b, rnd);                                      \
    }                                                                   \
                                                                        \
    template <typename T>                                               \
    static inline void                                                  \
    cado_mpfr_##OP(DTYPE a, T const b, mpfr_rnd_t rnd)                  \
        requires(!cado::converts_via<T, long> &&                        \
                  cado::converts_via<T, int64_t>)                       \
    {                                                                   \
        mpfr_##OP##_int64(a, b, rnd);                                   \
    }                                                                   \
                                                                        \
    template <typename T>                                               \
    static inline void                                                  \
    cado_mpfr_##OP(DTYPE a, T const b, mpfr_rnd_t rnd)                  \
        requires(!cado::converts_via<T, unsigned long> &&               \
                  cado::converts_via<T, uint64_t>)                      \
    {                                                                   \
        mpfr_##OP##_uint64(a, b, rnd);                                  \
    }

MPFR_AUXX_DEFINE_FUNC2(mpfr_ptr, init_set)
MPFR_AUXX_DEFINE_FUNC2(mpfr_ptr, set)

/*****************************************************************/
static inline int cado_mpfr_cmp(mpfr_srcptr a, mpfr_srcptr b)
{
    return mpfr_cmp(a, b);
}

template <typename T>
static inline int
cado_mpfr_cmp(mpfr_srcptr a, T const b)
    requires cado::converts_via<T, long>
{
    return mpfr_cmp_si(a, b);
}

template <typename T>
static inline int
cado_mpfr_cmp(mpfr_srcptr a, T const b)
    requires cado::converts_via<T, unsigned long>
{
    return mpfr_cmp_ui(a, b);
}

template <typename T>
static inline int
cado_mpfr_cmp(mpfr_srcptr a, T const b)
    requires(!cado::converts_via<T, long> &&
              cado::converts_via<T, int64_t>)
{
    return mpfr_cmp_int64(a, b);
}

template <typename T>
static inline int
cado_mpfr_cmp(mpfr_srcptr a, T const b)
    requires(!cado::converts_via<T, unsigned long> &&
              cado::converts_via<T, uint64_t>)
{
    return mpfr_cmp_uint64(a, b);
}

static inline int
cado_mpfr_cmp(mpfr_srcptr a, float b)
{
    /* no mpfr_cmp_* for float */
    return mpfr_cmp_d(a, b);
}

static inline int
cado_mpfr_cmp(mpfr_srcptr a, double b)
{
    return mpfr_cmp_d(a, b);
}

static inline int
cado_mpfr_cmp(mpfr_srcptr a, long double b)
{
    return mpfr_cmp_ld(a, b);
}

/*****************************************************************/

#define MPFR_AUXX_DEFINE_FUNC3(OP)                                      \
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

#define MPFR_AUXX_DEFINE_FUNC3_REFLEX(OP, FIXUP)                        \
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

MPFR_AUXX_DEFINE_FUNC3(add)
MPFR_AUXX_DEFINE_FUNC3(sub)
MPFR_AUXX_DEFINE_FUNC3(mul)
MPFR_AUXX_DEFINE_FUNC3(addmul)
MPFR_AUXX_DEFINE_FUNC3(submul)
MPFR_AUXX_DEFINE_FUNC3(div)
MPFR_AUXX_DEFINE_FUNC3(remainder)

/* Add these for convenience only */
MPFR_AUXX_DEFINE_FUNC3_REFLEX(sub, mpfr_neg(a, a, rnd); r = -r)
MPFR_AUXX_DEFINE_FUNC3_REFLEX(div, mpfr_pow_si(a, a, -1, rnd); r = -r)
MPFR_AUXX_DEFINE_FUNC3_REFLEX(add, /* no fixup for commutative op */)
MPFR_AUXX_DEFINE_FUNC3_REFLEX(mul, /* no fixup for commutative op */)

} /* namespace mpfr_auxx */

#endif /* UTILS_MPFR_AUXX_HPP_ */
