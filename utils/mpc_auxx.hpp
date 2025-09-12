#ifndef UTILS_MPC_AUXX_HPP_
#define UTILS_MPC_AUXX_HPP_

#include <cstdint>
#include <complex>

#include <mpc.h>

#include "mpc_aux.h"
#include "utils_cxx.hpp"

/* A C++ wrapper around functions in mpc_aux.h. The functions here delegate
 * to the proper C function depending on the argument type. We keep them
 * out of global name space to avoid obscuring errors.
 *
 * Note that we use the cado_mpc_ prefix (inside the mpc_auxx
 * namespace!) because most mpc names are macros.
 */

namespace mpc_auxx
{
#define MPC_AUXX_DEFINE_FUNC2(DTYPE, OP)                                \
    static inline void cado_mpc_##OP(DTYPE a, mpc_srcptr b,             \
                                      mpc_rnd_t rnd)                    \
    {                                                                   \
        mpc_##OP(a, b, rnd);                                            \
    }                                                                   \
                                                                        \
    template <typename T>                                               \
    static inline void                                                  \
    cado_mpc_##OP(DTYPE a, T const b, mpc_rnd_t rnd)                    \
    requires cado::converts_via<T, long>                                \
    {                                                                   \
        mpc_##OP##_si(a, b, rnd);                                       \
    }                                                                   \
                                                                        \
    template <typename T>                                               \
    static inline void                                                  \
    cado_mpc_##OP(DTYPE a, T const b, mpc_rnd_t rnd)                    \
    requires cado::converts_via<T, unsigned long>                       \
    {                                                                   \
        mpc_##OP##_ui(a, b, rnd);                                       \
    }                                                                   \
                                                                        \
    template <typename T>                                               \
    static inline void                                                  \
    cado_mpc_##OP(DTYPE a, T const b, mpc_rnd_t rnd)                    \
        requires(!cado::converts_via<T, long> &&                        \
                  cado::converts_via<T, int64_t>)                       \
    {                                                                   \
        mpc_##OP##_int64(a, b, rnd);                                    \
    }                                                                   \
                                                                        \
    template <typename T>                                               \
    static inline void                                                  \
    cado_mpc_##OP(DTYPE a, T const b, mpc_rnd_t rnd)                    \
        requires(!cado::converts_via<T, unsigned long> &&               \
                  cado::converts_via<T, uint64_t>)                      \
    {                                                                   \
        mpc_##OP##_uint64(a, b, rnd);                                   \
    }                                                                   \
                                                                        \
    static inline void                                                  \
    cado_mpc_##OP(DTYPE a, float const b, mpc_rnd_t rnd)                \
    {                                                                   \
        /* use _d for floats. mpc's _f functions are for mpf ! */       \
        mpc_##OP##_d(a, b, rnd);                                        \
    }                                                                   \
                                                                        \
    static inline void                                                  \
    cado_mpc_##OP(DTYPE a, double const b, mpc_rnd_t rnd)               \
    {                                                                   \
        mpc_##OP##_d(a, b, rnd);                                        \
    }                                                                   \
                                                                        \
    static inline void                                                  \
    cado_mpc_##OP(DTYPE a, long double const b, mpc_rnd_t rnd)          \
    {                                                                   \
        mpc_##OP##_ld(a, b, rnd);                                       \
    }                                                                   \
                                                                        \
    static inline void                                                  \
    cado_mpc_##OP(DTYPE a, std::complex<float> const b, mpc_rnd_t rnd)  \
    {                                                                   \
        /* use _d for floats. mpc's _f functions are for mpf ! */       \
        mpc_##OP##_d_d(a, b.real(), b.imag(), rnd);                     \
    }                                                                   \
                                                                        \
    static inline void                                                  \
    cado_mpc_##OP(DTYPE a, std::complex<double> const b, mpc_rnd_t rnd) \
    {                                                                   \
        mpc_##OP##_d_d(a, b.real(), b.imag(), rnd);                     \
    }                                                                   \
                                                                        \
    static inline void                                                  \
    cado_mpc_##OP(DTYPE a, std::complex<long double> const b, mpc_rnd_t rnd) \
    {                                                                   \
        mpc_##OP##_ld_ld(a, b.real(), b.imag(), rnd);                   \
    }

/* This is missing in mpc.h */
static inline void mpc_init_set(mpc_ptr a, mpc_srcptr b, mpc_rnd_t rnd)
{
    mpc_init2(a, mpfr_get_default_prec());
    mpc_set(a, b, rnd);
}

MPC_AUXX_DEFINE_FUNC2(mpc_ptr, init_set)
MPC_AUXX_DEFINE_FUNC2(mpc_ptr, set)


/*****************************************************************/
static inline int cado_mpc_cmp(mpc_srcptr a, mpc_srcptr b)
{
    return mpc_cmp(a, b);
}

template <typename T>
static inline int
cado_mpc_cmp(mpc_srcptr a, T const b)
    requires cado::converts_via<T, long>
{
    return mpc_cmp_si(a, b);
}

template <typename T>
static inline int
cado_mpc_cmp(mpc_srcptr a, T const b)
    requires cado::converts_via<T, unsigned long>
{
    return mpc_cmp_ui(a, b);
}

template <typename T>
static inline int
cado_mpc_cmp(mpc_srcptr a, T const b)
    requires(!cado::converts_via<T, long> &&
              cado::converts_via<T, int64_t>)
{
    return mpc_cmp_int64(a, b);
}

template <typename T>
static inline int
cado_mpc_cmp(mpc_srcptr a, T const b)
    requires(!cado::converts_via<T, unsigned long> &&
              cado::converts_via<T, uint64_t>)
{
    return mpc_cmp_uint64(a, b);
}

static inline int
cado_mpc_cmp(mpc_srcptr a, float b)
{
    /* no mpc_cmp_* for float */
    return mpc_cmp_d(a, b);
}

static inline int
cado_mpc_cmp(mpc_srcptr a, double b)
{
    return mpc_cmp_d(a, b);
}

static inline int
cado_mpc_cmp(mpc_srcptr a, long double b)
{
    return mpc_cmp_ld(a, b);
}

static inline int
cado_mpc_cmp(mpc_srcptr a, std::complex<float> b)
{
    /* no mpc_cmp_* for float */
    return mpc_cmp_d_d(a, b.real(), b.imag());
}

static inline int
cado_mpc_cmp(mpc_srcptr a, std::complex<double> b)
{
    return mpc_cmp_d_d(a, b.real(), b.imag());
}

static inline int
cado_mpc_cmp(mpc_srcptr a, std::complex<long double> b)
{
    return mpc_cmp_ld_ld(a, b.real(), b.imag());
}

/*****************************************************************/

#define MPC_AUXX_DEFINE_FUNC3(OP)                                       \
    static inline int cado_mpc_##OP(mpc_ptr a,                          \
                                    mpc_srcptr b, mpc_srcptr c,         \
                                    mpc_rnd_t rnd)                      \
    {                                                                   \
        return ::mpc_##OP(a, b, c, rnd);                                \
    }                                                                   \
    template <typename T>                                               \
    static inline int                                                   \
    cado_mpc_##OP(mpc_ptr a, mpc_srcptr b, const T c, mpc_rnd_t rnd)    \
    requires cado::converts_via<T, unsigned long>                       \
    {                                                                   \
        return mpc_##OP##_ui(a, b, c, rnd);                             \
    }                                                                   \
    template <typename T>                                               \
    static inline int cado_mpc_##OP(mpc_ptr a, mpc_srcptr b,            \
                                    const T c, mpc_rnd_t rnd)           \
    requires(!cado::converts_via<T, unsigned long> &&                   \
              cado::converts_via<T, uint64_t>)                          \
    {                                                                   \
        return mpc_##OP##_uint64(a, b, c, rnd);                         \
    }                                                                   \
    template <typename T>                                               \
    static inline int                                                   \
    cado_mpc_##OP(mpc_ptr a, mpc_srcptr b, const T c, mpc_rnd_t rnd)    \
    requires cado::converts_via<T, long>                                \
    {                                                                   \
        return mpc_##OP##_si(a, b, c, rnd);                             \
    }                                                                   \
    template <typename T>                                               \
    static inline int cado_mpc_##OP(mpc_ptr a, mpc_srcptr b,            \
                                     const T c, mpc_rnd_t rnd)          \
    requires(!cado::converts_via<T, long> &&                            \
              cado::converts_via<T, int64_t>)                           \
    {                                                                   \
        return mpc_##OP##_int64(a, b, c, rnd);                          \
    }                                                                   \
                                                                        \
    static inline void                                                  \
    cado_mpc_##OP(mpc_ptr a, mpc_srcptr b, float c, mpc_rnd_t rnd)      \
    {                                                                   \
        /* use _d for floats. mpc's _f functions are for mpf ! */       \
        mpc_##OP##_d(a, b, c, rnd);                                     \
    }                                                                   \
                                                                        \
    static inline void                                                  \
    cado_mpc_##OP(mpc_ptr a, mpc_srcptr b, double c, mpc_rnd_t rnd)     \
    {                                                                   \
        mpc_##OP##_d(a, b, c, rnd);                                     \
    }                                                                   \
                                                                        \
    static inline void                                                  \
    cado_mpc_##OP(mpc_ptr a, mpc_srcptr b, long double c, mpc_rnd_t rnd)\
    {                                                                   \
        mpc_##OP##_ld(a, b, c, rnd);                                    \
    }                                                                   \
                                                                        \
    static inline void                                                  \
    cado_mpc_##OP(mpc_ptr a, mpc_srcptr b, std::complex<float> c, mpc_rnd_t rnd)          \
    {                                                                   \
        /* use _d for floats. mpc's _f functions are for mpf ! */       \
        mpc_##OP##_d_d(a, b, c.real(), c.imag(), rnd);                  \
    }                                                                   \
                                                                        \
    static inline void                                                  \
    cado_mpc_##OP(mpc_ptr a, mpc_srcptr b, std::complex<double> c, mpc_rnd_t rnd)          \
    {                                                                   \
        mpc_##OP##_d_d(a, b, c.real(), c.imag(), rnd);                  \
    }                                                                   \
                                                                        \
    static inline void                                                  \
    cado_mpc_##OP(mpc_ptr a, mpc_srcptr b, std::complex<long double> c, mpc_rnd_t rnd)          \
    {                                                                   \
        mpc_##OP##_ld_ld(a, b, c.real(), c.imag(), rnd);                \
    }

#define MPC_AUXX_DEFINE_FUNC3_REFLEX(OP, FIXUP)                         \
    template <typename T>                                               \
    static inline int                                                   \
    cado_mpc_##OP(mpc_ptr a, const T b, mpc_srcptr c,                   \
                   mpc_rnd_t rnd)                                       \
    requires cado::converts_via<T, unsigned long>                       \
    {                                                                   \
        int r = mpc_##OP##_ui(a, c, b, rnd);                            \
        FIXUP;                                                          \
        return r;                                                       \
    }                                                                   \
    template <typename T>                                               \
    static inline int cado_mpc_##OP(mpc_ptr a, const T b,               \
                                     mpc_srcptr c, mpc_rnd_t rnd)       \
    requires(!cado::converts_via<T, unsigned long> &&                   \
              cado::converts_via<T, uint64_t>)                          \
    {                                                                   \
        int r = mpc_##OP##_uint64(a, c, b, rnd);                        \
        FIXUP;                                                          \
        return r;                                                       \
    }                                                                   \
    template <typename T>                                               \
    static inline int                                                   \
    cado_mpc_##OP(mpc_ptr a, const T b, mpc_srcptr c,                   \
                   mpc_rnd_t rnd)                                       \
    requires cado::integral_fits_v<T, long>                             \
    {                                                                   \
        int r = mpc_##OP##_si(a, c, b, rnd);                            \
        FIXUP;                                                          \
        return r;                                                       \
    }                                                                   \
    template <typename T>                                               \
    static inline int cado_mpc_##OP(mpc_ptr a, const T b,               \
                                     mpc_srcptr c, mpc_rnd_t rnd)       \
    requires(!cado::converts_via<T, long> &&                            \
              cado::converts_via<T, int64_t>)                           \
    {                                                                   \
        int r = mpc_##OP##_int64(a, c, b, rnd);                         \
        FIXUP;                                                          \
        return r;                                                       \
    }

MPC_AUXX_DEFINE_FUNC3(add)
MPC_AUXX_DEFINE_FUNC3(sub)
MPC_AUXX_DEFINE_FUNC3(mul)
MPC_AUXX_DEFINE_FUNC3(addmul)
MPC_AUXX_DEFINE_FUNC3(submul)
MPC_AUXX_DEFINE_FUNC3(div)
MPC_AUXX_DEFINE_FUNC3(remainder)

/* Add these for convenience only */
MPC_AUXX_DEFINE_FUNC3_REFLEX(sub, mpc_neg(a, a, rnd); r = -r)
MPC_AUXX_DEFINE_FUNC3_REFLEX(add, /* no fixup for commutative op */)
MPC_AUXX_DEFINE_FUNC3_REFLEX(mul, /* no fixup for commutative op */)

} /* namespace mpc_auxx */


#endif	/* UTILS_MPC_AUXX_HPP_ */
