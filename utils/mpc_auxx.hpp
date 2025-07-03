#ifndef UTILS_MPC_AUXX_HPP_
#define UTILS_MPC_AUXX_HPP_

#include <cstdint>

#include <type_traits>

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

/* NOLINTBEGIN(bugprone-macro-parentheses) */
#define MPC_AUXXX_EXPOSE_TEMPLATE_T_RETTYPE_IF_T_FITS(T, ret, fit)            \
    template <typename T>                                                      \
    static inline                                                              \
        typename std::enable_if<integral_fits_v<T, fit>, ret>::type

/* This is not accepted by gcc-9.2.0, however it seems correct as far as
 * I can tell. See bug #21817
 *
 */
#define ALT_MPC_AUXXX_EXPOSE_TEMPLATE_T_RETTYPE_IF_T_FITS(T, ret, fit)        \
    template <typename T, integral_fits_t<T, fit> = false> static inline ret
/* NOLINTEND(bugprone-macro-parentheses) */

/*****************************************************************/
static inline int cado_mpc_set(mpc_ptr a, mpc_srcptr b, mpc_rnd_t rnd)
{
    return mpc_set(a, b, rnd);
}

MPC_AUXXX_EXPOSE_TEMPLATE_T_RETTYPE_IF_T_FITS(T, int, long)
cado_mpc_set(mpc_ptr a, T const b, mpc_rnd_t rnd)
{
    return mpc_set_si(a, b, rnd);
}

MPC_AUXXX_EXPOSE_TEMPLATE_T_RETTYPE_IF_T_FITS(T, int, unsigned long)
cado_mpc_set(mpc_ptr a, T const b, mpc_rnd_t rnd)
{
    return mpc_set_ui(a, b, rnd);
}

inline int cado_mpc_set(mpc_ptr a, int64_t const b, mpc_rnd_t rnd)
{
    return mpc_set_int64(a, b, rnd);
}

inline int cado_mpc_set(mpc_ptr a, uint64_t const b, mpc_rnd_t rnd)
{
    return mpc_set_uint64(a, b, rnd);
}

/*****************************************************************/
static inline void cado_mpc_init_set(mpc_ptr a, mpc_srcptr b, mpc_rnd_t rnd)
{
    mpc_init2(a, mpfr_get_default_prec());
    mpc_set(a, b, rnd);
}

MPC_AUXXX_EXPOSE_TEMPLATE_T_RETTYPE_IF_T_FITS(T, void, long)
cado_mpc_init_set(mpc_ptr a, T const b, mpc_rnd_t rnd)
{
    mpc_init_set_si(a, b, rnd);
}

MPC_AUXXX_EXPOSE_TEMPLATE_T_RETTYPE_IF_T_FITS(T, void, unsigned long)
cado_mpc_init_set(mpc_ptr a, T const b, mpc_rnd_t rnd)
{
    mpc_init_set_ui(a, b, rnd);
}

inline void cado_mpc_init_set(mpc_ptr a, int64_t const b, mpc_rnd_t rnd)
{
    mpc_init_set_int64(a, b, rnd);
}

inline void cado_mpc_init_set(mpc_ptr a, uint64_t const b, mpc_rnd_t rnd)
{
    mpc_init_set_uint64(a, b, rnd);
}

/*****************************************************************/
static inline int cado_mpc_cmp(mpc_srcptr a, mpc_srcptr b)
{
    return mpc_cmp(a, b);
}

MPC_AUXXX_EXPOSE_TEMPLATE_T_RETTYPE_IF_T_FITS(T, int, long)
cado_mpc_cmp(mpc_srcptr a, T const b)
{
    return mpc_cmp_si(a, b);
}

MPC_AUXXX_EXPOSE_TEMPLATE_T_RETTYPE_IF_T_FITS(T, int, unsigned long)
cado_mpc_cmp(mpc_srcptr a, T const b)
{
    return mpc_cmp_ui(a, b);
}

inline int cado_mpc_cmp(mpc_srcptr a, int64_t const b)
{
    return mpc_cmp_int64(a, b);
}

inline int cado_mpc_cmp(mpc_srcptr a, uint64_t const b)
{
    return mpc_cmp_uint64(a, b);
}
/*****************************************************************/

#define MPC_AUXX_DEFINE_FUNC3(OP)                                             \
    static inline int cado_mpc_##OP(mpc_ptr a, mpc_srcptr b, mpc_srcptr c, \
                                     mpc_rnd_t rnd)                           \
    {                                                                          \
        return ::mpc_##OP(a, b, c, rnd);                                      \
    }                                                                          \
    MPC_AUXXX_EXPOSE_TEMPLATE_T_RETTYPE_IF_T_FITS(T, int, unsigned long)      \
    mpc_##OP(mpc_ptr a, mpc_srcptr b, const T c, mpc_rnd_t rnd)            \
    {                                                                          \
        return mpc_##OP##_ui(a, b, c, rnd);                                   \
    }                                                                          \
    static inline int cado_mpc_##OP(mpc_ptr a, mpc_srcptr b,                \
                                     const uint64_t c, mpc_rnd_t rnd)         \
    {                                                                          \
        return mpc_##OP##_uint64(a, b, c, rnd);                               \
    }                                                                          \
    MPC_AUXXX_EXPOSE_TEMPLATE_T_RETTYPE_IF_T_FITS(T, int, long)               \
    mpc_##OP(mpc_ptr a, mpc_srcptr b, const T c, mpc_rnd_t rnd)            \
    {                                                                          \
        return mpc_##OP##_si(a, b, c, rnd);                                   \
    }                                                                          \
    static inline int cado_mpc_##OP(mpc_ptr a, mpc_srcptr b,                \
                                     const int64_t c, mpc_rnd_t rnd)          \
    {                                                                          \
        return mpc_##OP##_int64(a, b, c, rnd);                                \
    }

MPC_AUXX_DEFINE_FUNC3(add)
MPC_AUXX_DEFINE_FUNC3(sub)
MPC_AUXX_DEFINE_FUNC3(mul)
MPC_AUXX_DEFINE_FUNC3(addmul)
MPC_AUXX_DEFINE_FUNC3(submul)
MPC_AUXX_DEFINE_FUNC3(div)
MPC_AUXX_DEFINE_FUNC3(remainder)

/* Add these for convenience only */
MPC_AUXXX_EXPOSE_TEMPLATE_T_RETTYPE_IF_T_FITS(T, void, unsigned long)
cado_mpc_sub(mpc_ptr a, T b, mpc_srcptr c, mpc_rnd_t rnd)
{
    mpc_sub_ui(a, c, b, rnd);
    mpc_neg(a, a, rnd);
}

static inline void cado_mpc_sub(mpc_ptr a, uint64_t b, mpc_srcptr c,
                                 mpc_rnd_t rnd)
{
    mpc_sub_uint64(a, c, b, rnd);
    mpc_neg(a, a, rnd);
}

MPC_AUXXX_EXPOSE_TEMPLATE_T_RETTYPE_IF_T_FITS(T, void, unsigned long)
cado_mpc_add(mpc_ptr a, T b, mpc_srcptr c, mpc_rnd_t rnd)
{
    mpc_add_ui(a, c, b, rnd);
}

static inline void cado_mpc_add(mpc_ptr a, uint64_t b, mpc_srcptr c,
                                 mpc_rnd_t rnd)
{
    mpc_add_uint64(a, c, b, rnd);
}

} /* namespace mpc_auxx */


#endif	/* UTILS_MPC_AUXX_HPP_ */
