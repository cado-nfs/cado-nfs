#ifndef UTILS_MPFR_AUX_H_
#define UTILS_MPFR_AUX_H_

#include "cado_config.h" // for ULONG_BITS

#include <limits.h>
#include <stdint.h>

#include <mpfr.h>

#ifdef __cplusplus
extern "C" {
#endif

/* {{{ mpfr_set_{int64,uin64} */
#if ULONG_BITS < 64
/* Set z to q. On 32-bit machines, we can't use mpfr_set_ui! */
static inline int mpfr_set_uint64(mpfr_ptr z, uint64_t q, mpfr_rnd_t rnd)
{
    if (q <= ULONG_MAX)
        return mpfr_set_ui(z, (unsigned long)q, rnd);
    else {
        ASSERT_ALWAYS(sizeof(unsigned long) == 4);
        int const r0 = mpfr_set_ui(z, (unsigned long)(q >> 32), rnd);
        mpfr_mul_2exp(z, z, 32, rnd);
        int const r1 =
            mpfr_add_ui(z, z, (unsigned long)(q & 4294967295UL), rnd);
        return r0 ? r0 : r1;
    }
}

static inline void mpfr_set_int64(mpfr_ptr z, int64_t q, mpfr_rnd_t rnd)
{
    if (LONG_MIN <= q && q <= LONG_MAX)
        return mpfr_set_si(z, (long)q, rnd);
    else if (q >= 0)
        return mpfr_set_uint64(z, (uint64_t)q, rnd);
    else {
        int const r = mpfr_set_uint64(z, -(uint64_t)q, rnd);
        mpfr_neg(z, z, rnd);
        return -r;
    }
}
#else
static inline int mpfr_set_uint64(mpfr_ptr a, uint64_t const b, mpfr_rnd_t rnd)
{
    return mpfr_set_ui(a, b, rnd);
}
static inline int mpfr_set_int64(mpfr_ptr a, int64_t const b, mpfr_rnd_t rnd)
{
    return mpfr_set_si(a, b, rnd);
}
#endif
/* }}} */

/* mpfr_init_set_{int64,uint64} {{{ */
/* for consistency with the other mpfr_init_* ctor functions, the
 * mpfr_t is initalized with mpfr_get_default_prec()
 */
#if ULONG_BITS < 64
static inline void mpfr_init_set_uint64(mpfr_ptr z, uint64_t x, mpfr_rnd_t rnd)
{
    mpfr_init2(z, mpfr_get_default_prec());
    mpfr_set_uint64(z, x, rnd);
}

static inline void mpfr_init_set_int64(mpfr_ptr z, int64_t x, mpc_rnd_t rnd)
{
    mpfr_init2(z, mpfr_get_default_prec());
    mpfr_set_int64(z, x, rnd);
}

#else
static inline void mpfr_init_set_uint64(mpfr_ptr a, const uint64_t b,
                                        mpfr_rnd_t rnd)
{
    mpfr_init_set_ui(a, b, rnd);
}
static inline void mpfr_init_set_int64(mpfr_ptr a, int64_t const b,
                                       mpfr_rnd_t rnd)
{
    mpfr_init_set_si(a, b, rnd);
}
#endif
// }}}

#if ULONG_BITS < 64
#define MPFR_AUX_DEFINE_COMPARISON(OP)                                         \
    static inline int mpfr_##OP##_uint64(mpfr_srcptr a, uint64_t c)            \
    {                                                                          \
        if (c <= ULONG_MAX)                                                    \
            return mpfr_##OP##_ui(a, (unsigned long)c);                        \
        mpfr_t cc;                                                             \
        mpfr_init2(cc, 64);                                                    \
        mpfr_set_uint64(cc, c, MPFR_RNDN);                                     \
        int const r = mpfr_##OP(a, cc);                                        \
        mpfr_clear(cc);                                                        \
        return r;                                                              \
    }                                                                          \
    static inline int mpfr_##OP##_int64(mpfr_srcptr a, int64_t c)              \
    {                                                                          \
        if (LONG_MIN <= c && c <= LONG_MAX)                                    \
            return mpfr_##OP##_ui(a, (long)c);                                 \
        mpfr_t cc;                                                             \
        mpfr_init2(cc, 64);                                                    \
        mpfr_set_int64(cc, c, MPFR_RNDN);                                      \
        int const r = mpfr_##OP(a, cc);                                        \
        mpfr_clear(cc);                                                        \
        return r;                                                              \
    }
#define MPFR_AUX_DEFINE_FUNC3(OP)                                              \
    static inline int mpfr_##OP##_uint64(mpfr_ptr a, mpfr_srcptr b,            \
                                         uint64_t c, mpfr_rnd_t rnd)           \
    {                                                                          \
        if (c <= ULONG_MAX)                                                    \
            return mpfr_##OP##_ui(a, b, (unsigned long)c, rnd);                \
        mpfr_t cc;                                                             \
        mpfr_init2(cc, 64);                                                    \
        mpfr_set_uint64(cc, c, MPFR_RNDN);                                     \
        int const r = mpfr_##OP(a, b, cc, rnd);                                \
        mpfr_clear(cc);                                                        \
        return r;                                                              \
    }                                                                          \
    static inline int mpfr_##OP##_int64(mpfr_ptr a, mpfr_srcptr b, int64_t c,  \
                                        mpfr_rnd_t rnd)                        \
    {                                                                          \
        if (LONG_MIN <= c && c <= LONG_MAX)                                    \
            return mpfr_##OP##_ui(a, b, (long)c, rnd);                         \
        mpfr_t cc;                                                             \
        mpfr_init2(cc, 64);                                                    \
        mpfr_set_int64(cc, c, MPFR_RNDN);                                      \
        int const r = mpfr_##OP(a, b, cc, rnd);                                \
        mpfr_clear(cc);                                                        \
        return r;                                                              \
    }
#else
#define MPFR_AUX_DEFINE_COMPARISON(OP)                                         \
    static inline int mpfr_##OP##_uint64(mpfr_srcptr a, uint64_t c)            \
    {                                                                          \
        return mpfr_##OP##_ui(a, c);                                           \
    }                                                                          \
    static inline int mpfr_##OP##_int64(mpfr_srcptr a, int64_t c)              \
    {                                                                          \
        return mpfr_##OP##_si(a, c);                                           \
    }
#define MPFR_AUX_DEFINE_FUNC3(OP)                                              \
    static inline int mpfr_##OP##_uint64(mpfr_ptr a, mpfr_srcptr b,            \
                                         uint64_t c, mpfr_rnd_t rnd)           \
    {                                                                          \
        return mpfr_##OP##_ui(a, b, c, rnd);                                   \
    }                                                                          \
    static inline int mpfr_##OP##_int64(mpfr_ptr a, mpfr_srcptr b, int64_t c,  \
                                        mpfr_rnd_t rnd)                        \
    {                                                                          \
        return mpfr_##OP##_si(a, b, c, rnd);                                   \
    }
#endif

/* {{{ mpfr_{add,sub}mul_{ui,si} */

/* mpfr doesn't have any addmul functions, but it does have FMA's */

static inline int mpfr_addmul(mpfr_ptr a, mpfr_srcptr b, mpfr_srcptr c,
                              mpfr_rnd_t rnd)
{
    return mpfr_fma(a, b, c, a, rnd);
}

static inline int mpfr_submul(mpfr_ptr a, mpfr_srcptr b, mpfr_srcptr c,
                              mpfr_rnd_t rnd)
{
    int const r = mpfr_fms(a, b, c, a, rnd);
    mpfr_neg(a, a, rnd);
    return -r;
}

static inline int mpfr_addmul_ui(mpfr_ptr a, mpfr_srcptr b,
                                 unsigned long const c, mpfr_rnd_t rnd)
{
    mpfr_t cc;
    mpfr_init2(cc, ULONG_BITS);
    mpfr_set_ui(cc, c, MPFR_RNDN);
    int const r = mpfr_addmul(a, b, cc, rnd);
    mpfr_clear(cc);
    return r;
}

static inline int mpfr_submul_ui(mpfr_ptr a, mpfr_srcptr b,
                                 unsigned long const c, mpfr_rnd_t rnd)
{
    mpfr_t cc;
    mpfr_init2(cc, ULONG_BITS);
    mpfr_set_ui(cc, c, MPFR_RNDN);
    int const r = mpfr_submul(a, b, cc, rnd);
    mpfr_clear(cc);
    return r;
}

static inline int mpfr_addmul_si(mpfr_ptr a, mpfr_srcptr b, long const c,
                                 mpfr_rnd_t rnd)
{
    mpfr_t cc;
    mpfr_init2(cc, ULONG_BITS);
    mpfr_set_si(cc, c, MPFR_RNDN);
    int const r = mpfr_addmul(a, b, cc, rnd);
    mpfr_clear(cc);
    return r;
}

static inline int mpfr_submul_si(mpfr_ptr a, mpfr_srcptr b, long const c,
                                 mpfr_rnd_t rnd)
{
    mpfr_t cc;
    mpfr_init2(cc, ULONG_BITS);
    mpfr_set_si(cc, c, MPFR_RNDN);
    int const r = mpfr_submul(a, b, cc, rnd);
    mpfr_clear(cc);
    return r;
}

static inline int mpfr_addmul_d(mpfr_ptr a, mpfr_srcptr b, double c,
                              mpfr_rnd_t rnd)
{
    mpfr_t cc;
    /* TODO: should we just write 53 here ? */
    mpfr_init2(cc, mpfr_get_prec(a));
    mpfr_set_d(cc, c, MPFR_RNDN);
    int const r = mpfr_addmul(a, b, cc, rnd);
    mpfr_clear(cc);
    return r;
}

static inline int mpfr_submul_d(mpfr_ptr a, mpfr_srcptr b, double c,
                              mpfr_rnd_t rnd)
{
    mpfr_t cc;
    /* TODO: should we just write 53 here ? */
    mpfr_init2(cc, mpfr_get_prec(a));
    mpfr_set_d(cc, c, MPFR_RNDN);
    int const r = mpfr_submul(a, b, cc, rnd);
    mpfr_clear(cc);
    return r;
}

static inline int mpfr_addmul_ld(mpfr_ptr a, mpfr_srcptr b, long double c,
                              mpfr_rnd_t rnd)
{
    mpfr_t cc;
    /* TODO: should we just write 80 here ? */
    mpfr_init2(cc, mpfr_get_prec(a));
    mpfr_set_ld(cc, c, MPFR_RNDN);
    int const r = mpfr_addmul(a, b, cc, rnd);
    mpfr_clear(cc);
    return r;
}

static inline int mpfr_submul_ld(mpfr_ptr a, mpfr_srcptr b, long double c,
                              mpfr_rnd_t rnd)
{
    mpfr_t cc;
    /* TODO: should we just write 80 here ? */
    mpfr_init2(cc, mpfr_get_prec(a));
    mpfr_set_ld(cc, c, MPFR_RNDN);
    int const r = mpfr_submul(a, b, cc, rnd);
    mpfr_clear(cc);
    return r;
}

/* }}} */

/* {{{ mpfr_{div,remainder}_{ui,si} */

static inline int mpfr_remainder_ui(mpfr_ptr a, mpfr_srcptr b,
                                    unsigned long const c, mpfr_rnd_t rnd)
{
    mpfr_t cc;
    mpfr_init2(cc, ULONG_BITS);
    mpfr_set_ui(cc, c, MPFR_RNDN);
    int const r = mpfr_remainder(a, b, cc, rnd);
    mpfr_clear(cc);
    return r;
}

static inline int mpfr_remainder_si(mpfr_ptr a, mpfr_srcptr b, long const c,
                                    mpfr_rnd_t rnd)
{
    mpfr_t cc;
    mpfr_init2(cc, ULONG_BITS);
    mpfr_set_si(cc, c, MPFR_RNDN);
    int const r = mpfr_remainder(a, b, cc, rnd);
    mpfr_clear(cc);
    return r;
}
/* }}} */

/*  mpfr_{set,cmp,add,sub,mul,addmul,submul,div,remainder}_{int64,uin64} */
MPFR_AUX_DEFINE_COMPARISON(cmp)
MPFR_AUX_DEFINE_FUNC3(add)
MPFR_AUX_DEFINE_FUNC3(sub)
MPFR_AUX_DEFINE_FUNC3(mul)
MPFR_AUX_DEFINE_FUNC3(addmul)
MPFR_AUX_DEFINE_FUNC3(submul)
MPFR_AUX_DEFINE_FUNC3(div)
MPFR_AUX_DEFINE_FUNC3(remainder)

#define MPFR_AUX_DEFINE_D_LD_FUNCTION(OP, TYP, SUF)			\
    static inline int mpfr_##OP##SUF(mpfr_ptr a, mpfr_srcptr b,		\
                                     long double c, mpfr_rnd_t rnd)	\
    {									\
        mpfr_t r;							\
        mpfr_init2(r, mpfr_get_prec(a));				\
        mpfr_set##SUF(r, c, rnd);					\
        int const res = mpfr_##OP(a, b, r, rnd);			\
        mpfr_clear(r);							\
        return res;							\
    }

MPFR_AUX_DEFINE_D_LD_FUNCTION(add, long double, _ld)
MPFR_AUX_DEFINE_D_LD_FUNCTION(sub, long double, _ld)
MPFR_AUX_DEFINE_D_LD_FUNCTION(mul, long double, _ld)
MPFR_AUX_DEFINE_D_LD_FUNCTION(div, long double, _ld)
MPFR_AUX_DEFINE_D_LD_FUNCTION(remainder, double, _d)
MPFR_AUX_DEFINE_D_LD_FUNCTION(remainder, long double, _ld)

#ifdef __cplusplus
}
#endif

#endif /* UTILS_MPFR_AUX_H_ */
