#ifndef UTILS_MPFR_AUXX_HPP_
#define UTILS_MPFR_AUXX_HPP_

#include <cstdint>

#include <type_traits>

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

/* NOLINTBEGIN(bugprone-macro-parentheses) */
#define MPFR_AUXXX_EXPOSE_TEMPLATE_T_RETTYPE_IF_T_FITS(T, ret, fit)            \
    template <typename T>                                                      \
    static inline                                                              \
        typename std::enable_if<integral_fits_v<T, fit>, ret>::type

/* This is not accepted by gcc-9.2.0, however it seems correct as far as
 * I can tell. See bug #21817
 *
 */
#define ALT_MPFR_AUXXX_EXPOSE_TEMPLATE_T_RETTYPE_IF_T_FITS(T, ret, fit)        \
    template <typename T, integral_fits_t<T, fit> = false> static inline ret
/* NOLINTEND(bugprone-macro-parentheses) */

/*****************************************************************/
static inline int cado_mpfr_set(mpfr_ptr a, mpfr_srcptr b, mpfr_rnd_t rnd)
{
    return mpfr_set(a, b, rnd);
}

MPFR_AUXXX_EXPOSE_TEMPLATE_T_RETTYPE_IF_T_FITS(T, int, long)
cado_mpfr_set(mpfr_ptr a, T const b, mpfr_rnd_t rnd)
{
    return mpfr_set_si(a, b, rnd);
}

MPFR_AUXXX_EXPOSE_TEMPLATE_T_RETTYPE_IF_T_FITS(T, int, unsigned long)
cado_mpfr_set(mpfr_ptr a, T const b, mpfr_rnd_t rnd)
{
    return mpfr_set_ui(a, b, rnd);
}

inline int cado_mpfr_set(mpfr_ptr a, int64_t const b, mpfr_rnd_t rnd)
{
    return mpfr_set_int64(a, b, rnd);
}

inline int cado_mpfr_set(mpfr_ptr a, uint64_t const b, mpfr_rnd_t rnd)
{
    return mpfr_set_uint64(a, b, rnd);
}

/*****************************************************************/
static inline void cado_mpfr_init_set(mpfr_ptr a, mpfr_srcptr b, mpfr_rnd_t rnd)
{
    mpfr_init_set(a, b, rnd);
}

MPFR_AUXXX_EXPOSE_TEMPLATE_T_RETTYPE_IF_T_FITS(T, void, long)
cado_mpfr_init_set(mpfr_ptr a, T const b, mpfr_rnd_t rnd)
{
    mpfr_init_set_si(a, b, rnd);
}

MPFR_AUXXX_EXPOSE_TEMPLATE_T_RETTYPE_IF_T_FITS(T, void, unsigned long)
cado_mpfr_init_set(mpfr_ptr a, T const b, mpfr_rnd_t rnd)
{
    mpfr_init_set_ui(a, b, rnd);
}

inline void cado_mpfr_init_set(mpfr_ptr a, int64_t const b, mpfr_rnd_t rnd)
{
    mpfr_init_set_int64(a, b, rnd);
}

inline void cado_mpfr_init_set(mpfr_ptr a, uint64_t const b, mpfr_rnd_t rnd)
{
    mpfr_init_set_uint64(a, b, rnd);
}

/*****************************************************************/
static inline int cado_mpfr_cmp(mpfr_srcptr a, mpfr_srcptr b)
{
    return mpfr_cmp(a, b);
}

MPFR_AUXXX_EXPOSE_TEMPLATE_T_RETTYPE_IF_T_FITS(T, int, long)
cado_mpfr_cmp(mpfr_srcptr a, T const b)
{
    return mpfr_cmp_si(a, b);
}

MPFR_AUXXX_EXPOSE_TEMPLATE_T_RETTYPE_IF_T_FITS(T, int, unsigned long)
cado_mpfr_cmp(mpfr_srcptr a, T const b)
{
    return mpfr_cmp_ui(a, b);
}

inline int cado_mpfr_cmp(mpfr_srcptr a, int64_t const b)
{
    return mpfr_cmp_int64(a, b);
}

inline int cado_mpfr_cmp(mpfr_srcptr a, uint64_t const b)
{
    return mpfr_cmp_uint64(a, b);
}
/*****************************************************************/

#define MPFR_AUXX_DEFINE_FUNC3(OP)                                             \
    static inline int cado_mpfr_##OP(mpfr_ptr a, mpfr_srcptr b, mpfr_srcptr c, \
                                     mpfr_rnd_t rnd)                           \
    {                                                                          \
        return ::mpfr_##OP(a, b, c, rnd);                                      \
    }                                                                          \
    MPFR_AUXXX_EXPOSE_TEMPLATE_T_RETTYPE_IF_T_FITS(T, int, unsigned long)      \
    mpfr_##OP(mpfr_ptr a, mpfr_srcptr b, const T c, mpfr_rnd_t rnd)            \
    {                                                                          \
        return mpfr_##OP##_ui(a, b, c, rnd);                                   \
    }                                                                          \
    static inline int cado_mpfr_##OP(mpfr_ptr a, mpfr_srcptr b,                \
                                     const uint64_t c, mpfr_rnd_t rnd)         \
    {                                                                          \
        return mpfr_##OP##_uint64(a, b, c, rnd);                               \
    }                                                                          \
    MPFR_AUXXX_EXPOSE_TEMPLATE_T_RETTYPE_IF_T_FITS(T, int, long)               \
    mpfr_##OP(mpfr_ptr a, mpfr_srcptr b, const T c, mpfr_rnd_t rnd)            \
    {                                                                          \
        return mpfr_##OP##_si(a, b, c, rnd);                                   \
    }                                                                          \
    static inline int cado_mpfr_##OP(mpfr_ptr a, mpfr_srcptr b,                \
                                     const int64_t c, mpfr_rnd_t rnd)          \
    {                                                                          \
        return mpfr_##OP##_int64(a, b, c, rnd);                                \
    }

MPFR_AUXX_DEFINE_FUNC3(add)
MPFR_AUXX_DEFINE_FUNC3(sub)
MPFR_AUXX_DEFINE_FUNC3(mul)
MPFR_AUXX_DEFINE_FUNC3(addmul)
MPFR_AUXX_DEFINE_FUNC3(submul)
MPFR_AUXX_DEFINE_FUNC3(div)
MPFR_AUXX_DEFINE_FUNC3(remainder)

/* Add these for convenience only */
MPFR_AUXXX_EXPOSE_TEMPLATE_T_RETTYPE_IF_T_FITS(T, void, unsigned long)
cado_mpfr_sub(mpfr_ptr a, T b, mpfr_srcptr c, mpfr_rnd_t rnd)
{
    mpfr_sub_ui(a, c, b, rnd);
    mpfr_neg(a, a, rnd);
}

static inline void cado_mpfr_sub(mpfr_ptr a, uint64_t b, mpfr_srcptr c,
                                 mpfr_rnd_t rnd)
{
    mpfr_sub_uint64(a, c, b, rnd);
    mpfr_neg(a, a, rnd);
}

MPFR_AUXXX_EXPOSE_TEMPLATE_T_RETTYPE_IF_T_FITS(T, void, unsigned long)
cado_mpfr_add(mpfr_ptr a, T b, mpfr_srcptr c, mpfr_rnd_t rnd)
{
    mpfr_add_ui(a, c, b, rnd);
}

static inline void cado_mpfr_add(mpfr_ptr a, uint64_t b, mpfr_srcptr c,
                                 mpfr_rnd_t rnd)
{
    mpfr_add_uint64(a, c, b, rnd);
}

} /* namespace mpfr_auxx */

#endif /* UTILS_MPFR_AUXX_HPP_ */
