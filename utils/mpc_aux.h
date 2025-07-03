#ifndef UTILS_MPC_AUX_H_
#define UTILS_MPC_AUX_H_

#include "cado_config.h" // for ULONG_BITS

#include <limits.h>
#include <stdint.h>

#include "macros.h"
#include <mpc.h>

#ifdef __cplusplus
extern "C" {
#endif

/* {{{ mpc_set_{int64,uin64} */
#if ULONG_BITS < 64
/* Set z to q. On 32-bit machines, we can't use mpc_set_ui! */
static inline int mpc_set_uint64(mpc_ptr z, uint64_t q, mpc_rnd_t rnd)
{
    if (q <= ULONG_MAX)
        return mpc_set_ui(z, (unsigned long)q, rnd);
    else {
        ASSERT_ALWAYS(sizeof(unsigned long) == 4);
        int const r0 = mpc_set_ui(z, (unsigned long)(q >> 32), rnd);
        mpc_mul_2exp(z, z, 32, rnd);
        int const r1 = mpc_add_ui(z, z, (unsigned long)(q & 4294967295UL), rnd);
        return r0 ? r0 : r1;
    }
}

static inline void mpc_set_int64(mpc_ptr z, int64_t q, mpc_rnd_t rnd)
{
    if (LONG_MIN <= q && q <= LONG_MAX)
        return mpc_set_si(z, (long)q, rnd);
    else if (q >= 0)
        return mpc_set_uint64(z, (uint64_t)q, rnd);
    else {
        int const r = mpc_set_uint64(z, -(uint64_t)q, rnd);
        mpc_neg(z, z, rnd);
        return -r;
    }
}
#else
static inline int mpc_set_uint64(mpc_ptr a, uint64_t const b, mpc_rnd_t rnd)
{
    return mpc_set_ui(a, b, rnd);
}
static inline int mpc_set_int64(mpc_ptr a, int64_t const b, mpc_rnd_t rnd)
{
    return mpc_set_si(a, b, rnd);
}
#endif
/* }}} */

/* Some of these functions are useful additions to the mpc types, in the
 * context of interaction with other types
 */
static inline void mpc_init_set_si(mpc_ptr z, long x, mpc_rnd_t rnd)
{
    mpc_init2(z, mpfr_get_default_prec());
    mpc_set_si(z, x, rnd);
}

static inline void mpc_init_set_ui(mpc_ptr z, unsigned long x, mpc_rnd_t rnd)
{
    mpc_init2(z, mpfr_get_default_prec());
    mpc_set_ui(z, x, rnd);
}
static inline int mpc_cmp_ui(mpc_srcptr a, unsigned long c)
{
    if (c <= LONG_MAX) {
        return mpc_cmp_si_si(a, long(c), 0);
    } else {
        return mpfr_cmp_ui(mpc_realref(a), c) == 0 && mpfr_cmp_ui(mpc_imagref(a), 0) == 0;
    }
}
static inline int mpc_sub_si(mpc_ptr a, mpc_srcptr b, long c, mpc_rnd_t rnd)
{
    int rr = mpfr_set(mpc_imagref(a), mpc_imagref(b), MPC_RND_IM(rnd));
    int ri = mpfr_sub_si(mpc_realref(a), mpc_realref(b), c, MPC_RND_RE(rnd));
    return MPC_INEX(rr, ri);
}
static inline int mpc_div_si(mpc_ptr a, mpc_srcptr b, long c, mpc_rnd_t rnd)
{
    int rr = mpfr_div_si(mpc_imagref(a), mpc_imagref(b), c, MPC_RND_IM(rnd));
    int ri = mpfr_div_si(mpc_realref(a), mpc_realref(b), c, MPC_RND_RE(rnd));
    return MPC_INEX(rr, ri);
}


/* mpc_init_set_{int64,uint64} {{{ */
/* for consistency with the other mpc_init_* ctor functions, the
 * mpc_t is initalized with mpc_get_default_prec()
 */
#if ULONG_BITS < 64
static inline void mpc_init_set_uint64(mpc_ptr z, uint64_t x, mpc_rnd_t rnd)
{
    mpc_init2(z, mpc_get_default_prec());
    mpc_set_uint64(z, x, rnd);
}

static inline void mpc_init_set_int64(mpc_ptr z, int64_t x, mpc_rnd_t rnd)
{
    mpc_init2(z, mpc_get_default_prec());
    mpc_set_int64(z, x, rnd);
}

#else
static inline void mpc_init_set_uint64(mpc_ptr a, const uint64_t b,
                                       mpc_rnd_t rnd)
{
    mpc_init_set_ui(a, b, rnd);
}
static inline void mpc_init_set_int64(mpc_ptr a, int64_t const b, mpc_rnd_t rnd)
{
    mpc_init_set_si(a, b, rnd);
}
#endif
// }}}

#if ULONG_BITS < 64
#define MPC_AUX_DEFINE_COMPARISON(OP)                                          \
    static inline int mpc_##OP##_uint64(mpc_srcptr a, uint64_t c)              \
    {                                                                          \
        if (c <= ULONG_MAX)                                                    \
            return mpc_##OP##_ui(a, (unsigned long)c);                         \
        mpc_t cc;                                                              \
        mpc_init2(cc, 64);                                                     \
        mpc_set_uint64(cc, c, MPC_RNDNN);                                      \
        int const r = mpc_##OP(a, cc);                                         \
        mpc_clear(cc);                                                         \
        return r;                                                              \
    }                                                                          \
    static inline int mpc_##OP##_int64(mpc_srcptr a, int64_t c)                \
    {                                                                          \
        if (LONG_MIN <= c && c <= LONG_MAX)                                    \
            return mpc_##OP##_ui(a, (long)c);                                  \
        mpc_t cc;                                                              \
        mpc_init2(cc, 64);                                                     \
        mpc_set_int64(cc, c, MPC_RNDNN);                                       \
        int const r = mpc_##OP(a, cc);                                         \
        mpc_clear(cc);                                                         \
        return r;                                                              \
    }
#define MPC_AUX_DEFINE_FUNC3(OP)                                               \
    static inline int mpc_##OP##_uint64(mpc_ptr a, mpc_srcptr b, uint64_t c,   \
                                        mpc_rnd_t rnd)                         \
    {                                                                          \
        if (c <= ULONG_MAX)                                                    \
            return mpc_##OP##_ui(a, b, (unsigned long)c, rnd);                 \
        mpc_t cc;                                                              \
        mpc_init2(cc, 64);                                                     \
        mpc_set_uint64(cc, c, MPC_RNDNN);                                      \
        int const r = mpc_##OP(a, b, cc, rnd);                                 \
        mpc_clear(cc);                                                         \
        return r;                                                              \
    }                                                                          \
    static inline int mpc_##OP##_int64(mpc_ptr a, mpc_srcptr b, int64_t c,     \
                                       mpc_rnd_t rnd)                          \
    {                                                                          \
        if (LONG_MIN <= c && c <= LONG_MAX)                                    \
            return mpc_##OP##_ui(a, b, (long)c, rnd);                          \
        mpc_t cc;                                                              \
        mpc_init2(cc, 64);                                                     \
        mpc_set_int64(cc, c, MPC_RNDNN);                                       \
        int const r = mpc_##OP(a, b, cc, rnd);                                 \
        mpc_clear(cc);                                                         \
        return r;                                                              \
    }
#else
#define MPC_AUX_DEFINE_COMPARISON(OP)                                          \
    static inline int mpc_##OP##_uint64(mpc_srcptr a, uint64_t c)              \
    {                                                                          \
        return mpc_##OP##_ui(a, c);                                            \
    }                                                                          \
    static inline int mpc_##OP##_int64(mpc_srcptr a, int64_t c)                \
    {                                                                          \
        return mpc_##OP##_si(a, c);                                            \
    }
#define MPC_AUX_DEFINE_FUNC3(OP)                                               \
    static inline int mpc_##OP##_uint64(mpc_ptr a, mpc_srcptr b, uint64_t c,   \
                                        mpc_rnd_t rnd)                         \
    {                                                                          \
        return mpc_##OP##_ui(a, b, c, rnd);                                    \
    }                                                                          \
    static inline int mpc_##OP##_int64(mpc_ptr a, mpc_srcptr b, int64_t c,     \
                                       mpc_rnd_t rnd)                          \
    {                                                                          \
        return mpc_##OP##_si(a, b, c, rnd);                                    \
    }
#endif

/* {{{ mpc_{add,sub}mul_{ui,si} */

/* mpc doesn't have any addmul functions, but it does have FMA's */

static inline int mpc_addmul(mpc_ptr a, mpc_srcptr b, mpc_srcptr c,
                             mpc_rnd_t rnd)
{
    return mpc_fma(a, b, c, a, rnd);
}

static inline int mpc_submul(mpc_ptr a, mpc_srcptr b, mpc_srcptr c,
                             mpc_rnd_t rnd)
{
    mpc_neg(a, a, rnd);
    int const r = mpc_fma(a, b, c, a, rnd);
    mpc_neg(a, a, rnd);
    return -r;
}

static inline int mpc_addmul_ui(mpc_ptr a, mpc_srcptr b, unsigned long const c,
                                mpc_rnd_t rnd)
{
    mpc_t cc;
    mpc_init2(cc, ULONG_BITS);
    mpc_set_ui(cc, c, MPC_RNDNN);
    int const r = mpc_addmul(a, b, cc, rnd);
    mpc_clear(cc);
    return r;
}

static inline int mpc_submul_ui(mpc_ptr a, mpc_srcptr b, unsigned long const c,
                                mpc_rnd_t rnd)
{
    mpc_t cc;
    mpc_init2(cc, ULONG_BITS);
    mpc_set_ui(cc, c, MPC_RNDNN);
    int const r = mpc_submul(a, b, cc, rnd);
    mpc_clear(cc);
    return r;
}

static inline int mpc_addmul_si(mpc_ptr a, mpc_srcptr b, long const c,
                                mpc_rnd_t rnd)
{
    mpc_t cc;
    mpc_init2(cc, ULONG_BITS);
    mpc_set_si(cc, c, MPC_RNDNN);
    int const r = mpc_addmul(a, b, cc, rnd);
    mpc_clear(cc);
    return r;
}

static inline int mpc_submul_si(mpc_ptr a, mpc_srcptr b, long const c,
                                mpc_rnd_t rnd)
{
    mpc_t cc;
    mpc_init2(cc, ULONG_BITS);
    mpc_set_si(cc, c, MPC_RNDNN);
    int const r = mpc_submul(a, b, cc, rnd);
    mpc_clear(cc);
    return r;
}

/* }}} */

/* {{{ mpc_{div,remainder}_{ui,si} */

static inline int mpc_remainder(mpc_ptr a, mpc_srcptr b, mpc_srcptr c,
                                mpc_rnd_t rnd)
{
    mpfr_prec_t pr, pi;
    mpc_get_prec2(&pr, &pi, a);
    mpc_t qq;
    mpc_init2(qq, MAX(pr, pi));
    mpc_div(qq, b, c, rnd);
    mpfr_rint_round(mpc_realref(qq), mpc_realref(qq), MPC_RND_RE(rnd));
    mpfr_rint_round(mpc_imagref(qq), mpc_imagref(qq), MPC_RND_IM(rnd));
    /* now compute b-q*c */
    mpc_neg(qq, qq, rnd);
    int const r = mpc_fma(a, qq, c, b, rnd);
    mpc_clear(qq);
    return r;
}

static inline int mpc_remainder_ui(mpc_ptr a, mpc_srcptr b,
                                   unsigned long const c, mpc_rnd_t rnd)
{
    mpc_t cc;
    mpc_init2(cc, ULONG_BITS);
    mpc_set_ui(cc, c, MPC_RNDNN);
    int const r = mpc_remainder(a, b, cc, rnd);
    mpc_clear(cc);
    return r;
}

static inline int mpc_remainder_si(mpc_ptr a, mpc_srcptr b, long const c,
                                   mpc_rnd_t rnd)
{
    mpc_t cc;
    mpc_init2(cc, ULONG_BITS);
    mpc_set_si(cc, c, MPC_RNDNN);
    int const r = mpc_remainder(a, b, cc, rnd);
    mpc_clear(cc);
    return r;
}
/* }}} */

/*  mpc_{set,cmp,add,sub,mul,addmul,submul,div,remainder}_{int64,uin64} */
MPC_AUX_DEFINE_COMPARISON(cmp)
MPC_AUX_DEFINE_FUNC3(add)
MPC_AUX_DEFINE_FUNC3(sub)
MPC_AUX_DEFINE_FUNC3(mul)
MPC_AUX_DEFINE_FUNC3(addmul)
MPC_AUX_DEFINE_FUNC3(submul)
MPC_AUX_DEFINE_FUNC3(div)
MPC_AUX_DEFINE_FUNC3(remainder)

#ifdef __cplusplus
}
#endif

#endif /* UTILS_MPC_AUX_H_ */
