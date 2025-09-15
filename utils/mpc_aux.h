#ifndef UTILS_MPC_AUX_H_
#define UTILS_MPC_AUX_H_

#include "cado_config.h" // for ULONG_BITS

#include <limits.h>
#include <stdint.h>

#include <mpc.h>

#include "macros.h"
#include "mpfr_aux.h"

/* we add here some functions that are, in our opinion, missing from the
 * mpc interface. In many cases the code is about right, but sometimes we
 * fall short of providing a correctly rounded implementation, which
 * would need some extra work
 */

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

/* comparison functions must follow the same convention as mpc_cmp ! */
static inline int mpc_cmp_ui(mpc_srcptr a, unsigned long c)
{
    if (c <= LONG_MAX) {
        return mpc_cmp_si_si(a, long(c), 0);
    } else {
        int const rr = mpfr_cmp_ui(mpc_realref(a), c);
        int const ri = mpfr_cmp_ui(mpc_imagref(a), 0);
        return MPC_INEX(rr, ri);
    }
}
static inline int mpc_cmp_d(mpc_srcptr a, double c)
{
    int const rr = mpfr_cmp_d(mpc_realref(a), c);
    int const ri = mpfr_cmp_ui(mpc_imagref(a), 0);
    return MPC_INEX(rr, ri);
}
static inline int mpc_cmp_ld(mpc_srcptr a, long double c)
{
    int const rr = mpfr_cmp_ld(mpc_realref(a), c);
    int const ri = mpfr_cmp_ui(mpc_imagref(a), 0);
    return MPC_INEX(rr, ri);
}
static inline int mpc_cmp_fr(mpc_srcptr a, mpfr_srcptr c)
{
    int const rr = mpfr_cmp(mpc_realref(a), c);
    int const ri = mpfr_cmp_ui(mpc_imagref(a), 0);
    return MPC_INEX(rr, ri);
}
static inline int mpc_cmp_d_d(mpc_srcptr a, double cr, double ci)
{
    int const rr = mpfr_cmp_d(mpc_realref(a), cr);
    int const ri = mpfr_cmp_d(mpc_imagref(a), ci);
    return MPC_INEX(rr, ri);
}
static inline int mpc_cmp_ld_ld(mpc_srcptr a, long double cr, long double ci)
{
    int const rr = mpfr_cmp_ld(mpc_realref(a), cr);
    int const ri = mpfr_cmp_ld(mpc_imagref(a), ci);
    return MPC_INEX(rr, ri);
}

static inline int mpc_sub_si(mpc_ptr a, mpc_srcptr b, long c, mpc_rnd_t rnd)
{
    int const rr = mpfr_set(mpc_imagref(a), mpc_imagref(b), MPC_RND_RE(rnd));
    int const ri = mpfr_sub_si(mpc_realref(a), mpc_realref(b), c, MPC_RND_IM(rnd));
    return MPC_INEX(rr, ri);
}
static inline int mpc_div_si(mpc_ptr a, mpc_srcptr b, long c, mpc_rnd_t rnd)
{
    int const rr = mpfr_div_si(mpc_imagref(a), mpc_imagref(b), c, MPC_RND_RE(rnd));
    int const ri = mpfr_div_si(mpc_realref(a), mpc_realref(b), c, MPC_RND_IM(rnd));
    return MPC_INEX(rr, ri);
}


/* mpc_init_set_{int64,uint64} {{{ */
/* for consistency with the other mpc_init_* ctor functions, the
 * mpc_t is initalized with mpfr_get_default_prec()
 */
#if ULONG_BITS < 64
static inline void mpc_init_set_uint64(mpc_ptr z, uint64_t x, mpc_rnd_t rnd)
{
    mpc_init2(z, mpfr_get_default_prec());
    mpc_set_uint64(z, x, rnd);
}

static inline void mpc_init_set_int64(mpc_ptr z, int64_t x, mpc_rnd_t rnd)
{
    mpc_init2(z, mpfr_get_default_prec());
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

// {{{ mpc_init_set_{d,ld,d_d,ld_ld}
static inline void mpc_init_set_d(mpc_ptr z, double x, mpc_rnd_t rnd)
{
    mpc_init2(z, mpfr_get_default_prec());
    mpc_set_d(z, x, rnd);
}

static inline void mpc_init_set_ld(mpc_ptr z, long double x, mpc_rnd_t rnd)
{
    mpc_init2(z, mpfr_get_default_prec());
    mpc_set_ld(z, x, rnd);
}

static inline void mpc_init_set_d_d(mpc_ptr z, double x, double y, mpc_rnd_t rnd)
{
    mpc_init2(z, mpfr_get_default_prec());
    mpc_set_d_d(z, x, y, rnd);
}

static inline void mpc_init_set_ld_ld(mpc_ptr z, long double x, long double y, mpc_rnd_t rnd)
{
    mpc_init2(z, mpfr_get_default_prec());
    mpc_set_ld_ld(z, x, y, rnd);
}
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

static inline int mpc_fms(mpc_ptr res, mpc_srcptr x, mpc_srcptr y, mpc_srcptr z,
                             mpc_rnd_t rnd)
{
    mpc_t mz;
    mpfr_prec_t pr, pi;
    mpc_get_prec2(&pr, &pi, z);
    mpc_init3(mz, pr, pi);
    int const r = mpc_fma(res, x, y, mz, rnd);
    mpc_neg(res, res, rnd);
    mpc_clear(mz);
    return r;
}

static inline int mpc_submul(mpc_ptr a, mpc_srcptr b, mpc_srcptr c,
                             mpc_rnd_t rnd)
{
    /* This one isn't copied from mpfr_aux.h, since mpc does not have an
     * fms operation */
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

#define MPC_AUX_FP_OP1(OP, TYP, SUF)                                         \
    static inline int mpc_##OP##SUF(mpc_ptr z, mpc_srcptr x,                 \
                                   TYP a, mpc_rnd_t rnd)                     \
    {                                                                        \
        int const rr = mpfr_##OP##SUF(mpc_realref(z), mpc_realref(x),        \
                                     a, MPC_RND_RE(rnd));                    \
        int const ri = mpfr_set(mpc_imagref(z), mpc_imagref(x),              \
                                MPC_RND_IM(rnd));                            \
        return MPC_INEX(rr, ri);                                             \
    }                                                                        \
    static inline int mpc_##OP##SUF##SUF(mpc_ptr z, mpc_srcptr x,            \
                                   TYP ar, TYP ai, mpc_rnd_t rnd)            \
    {                                                                        \
        int const rr = mpfr_##OP##SUF(mpc_realref(z), mpc_realref(x),        \
                                     ar, MPC_RND_RE(rnd));                   \
        int const ri = mpfr_##OP##SUF(mpc_imagref(z), mpc_imagref(x),        \
                                     ai, MPC_RND_IM(rnd));                   \
        return MPC_INEX(rr, ri);                                             \
    }

#define MPC_AUX_MUL_OP1(OP, TYP, SUF1, SUF2)                                 \
    static inline int mpc_##OP##SUF1(mpc_ptr z, mpc_srcptr x,                \
                                   TYP a, mpc_rnd_t rnd)                     \
    {                                                                        \
        int const rr = mpfr_##OP##SUF2(mpc_realref(z), mpc_realref(x),       \
                                     a, MPC_RND_RE(rnd));                    \
        int const ri = mpfr_##OP##SUF2(mpc_imagref(z), mpc_imagref(x),       \
                                     a, MPC_RND_IM(rnd));                    \
        return MPC_INEX(rr, ri);                                             \
    }

#define MPC_AUX_MUL_OP2(OP, TYP, SUF, combiner)                              \
    static inline int mpc_##OP##SUF##SUF(mpc_ptr z, mpc_srcptr x,            \
                                   TYP ar, TYP ai, mpc_rnd_t rnd)            \
    {                                                                        \
        /* zr = xr*ar - xi*ai                                                \
         * zi = xr*ai + xi*ar                                                \
         * we need two temps because of possible overlaps.                   \
         */                                                                  \
        mpc_t t;                                                             \
        mpfr_prec_t pr,pi;                                                   \
        mpc_get_prec2(&pr, &pi, z);                                          \
        mpc_init3(t, pr, pi);                                                \
        mpfr_mul##SUF(mpc_realref(t), mpc_realref(x), ar, MPC_RND_RE(rnd));  \
        mpfr_mul##SUF(mpc_imagref(t), mpc_realref(x), ai, MPC_RND_IM(rnd));  \
        int const rr = mpfr_submul##SUF(mpc_realref(t),                      \
                                        mpc_imagref(x), ai, MPC_RND_RE(rnd));\
        int const ri = mpfr_addmul##SUF(mpc_imagref(t),                      \
                                        mpc_imagref(x), ar, MPC_RND_IM(rnd));\
        combiner;                                                            \
        mpc_clear(t);                                                        \
        return MPC_INEX(rr, ri);                                             \
    }
#define MPC_AUX_MUL_OP3(OP, TYP, SUF)                                   \
        MPC_AUX_MUL_OP2(OP, TYP, SUF, mpc_swap(z, t))                   \
        MPC_AUX_MUL_OP2(add##OP, TYP, SUF, mpc_add(z, z, t, rnd))       \
        MPC_AUX_MUL_OP2(sub##OP, TYP, SUF, mpc_sub(z, z, t, rnd))

#define MPC_AUX_DIV_OP(OP, TYP, SUF)                                         \
    static inline int mpc_##OP##SUF##SUF(mpc_ptr z, mpc_srcptr x,            \
                                        TYP ar, TYP ai, mpc_rnd_t rnd)       \
    {                                                                        \
        /* zr = xr*ar - xi*ai                                                \
         * zi = xr*ai + xi*ar                                                \
         * we need two temps because of possible overlaps.                   \
         */                                                                  \
        mpc_t t;                                                             \
        mpfr_prec_t pr,pi;                                                   \
        mpc_get_prec2(&pr, &pi, z);                                          \
        mpc_init3(t, pr, pi);                                                \
        mpc_set##SUF##SUF(t, ar, ai, rnd);                                   \
        int const r = mpc_##OP(z, x, t, rnd);                                \
        mpc_clear(t);                                                        \
        return r;                                                            \
    }

MPC_AUX_FP_OP1(add, double, _d)
MPC_AUX_FP_OP1(add, long double, _ld)
MPC_AUX_FP_OP1(sub, double, _d)
MPC_AUX_FP_OP1(sub, long double, _ld)
/* multiplication by a real number */
MPC_AUX_MUL_OP1(mul, double, _d, _d)
MPC_AUX_MUL_OP1(mul, long double, _ld, _ld)
MPC_AUX_MUL_OP1(addmul, double, _d, _d)
MPC_AUX_MUL_OP1(addmul, long double, _ld, _ld)
MPC_AUX_MUL_OP1(addmul, mpfr_srcptr, _fr,)
MPC_AUX_MUL_OP1(submul, double, _d, _d)
MPC_AUX_MUL_OP1(submul, long double, _ld, _ld)
MPC_AUX_MUL_OP1(submul, mpfr_srcptr, _fr,)
/* multiplication by a complex number */
MPC_AUX_MUL_OP3(mul, double, _d)
MPC_AUX_MUL_OP3(mul, long double, _ld)
/* division by a real number */
MPC_AUX_MUL_OP1(div, double, _d, _d)
MPC_AUX_MUL_OP1(div, long double, _ld, _ld)
MPC_AUX_MUL_OP1(remainder, double, _d, _d)
MPC_AUX_MUL_OP1(remainder, long double, _ld, _ld)
MPC_AUX_MUL_OP1(remainder, mpfr_srcptr, _fr,)
/* division by a complex number is much less pretty! */
MPC_AUX_DIV_OP(div, double, _d)
MPC_AUX_DIV_OP(div, long double, _ld)
MPC_AUX_DIV_OP(remainder, double, _d)
MPC_AUX_DIV_OP(remainder, long double, _ld)

#ifdef __cplusplus
}
#endif

#endif /* UTILS_MPC_AUX_H_ */
