#ifndef MATMUL_FACADE_H_
#define MATMUL_FACADE_H_

/* This file is included by all matmul implementations so as to enable
 * the macros MATMUL_NAME and MATMUL_T.
 *
 * It does not make sense to have it included from the higher-level code
 * (in fact, it will error out).
 */
#include "matmul.h"

/* Include all matmul implementations here */

#if !defined(MM_MPFQ_LAYER)
#error "Please compile this file with the MM_MPFQ_LAYER macro defined"
#endif

/* Now some cpp glue which sets up the different options */
#define MATMUL_NAME_INNER(abase, impl, func) matmul_ ## abase ## _ ## impl ## _ ## func
#define MATMUL_NAME_INNER0(abase, impl, func) MATMUL_NAME_INNER(abase, impl, func)
#define MATMUL_NAME(func) MATMUL_NAME_INNER0(MM_MPFQ_LAYER, MM_IMPL, func)

#define MATMUL_T(func) matmul_ ## func ## _t

#define REBIND_F(mm, func) \
        mm->bind->func = (MATMUL_T(func)) & MATMUL_NAME(func)

#define SET_IMPL_INNER(mm, mpfq_, impl_) do {  \
    mm->bind->impl = # impl_ ;                  \
} while (0)

#define SET_IMPL_INNER0(mm, mpfq_, impl_) SET_IMPL_INNER(mm, mpfq_, impl_)

#define SET_IMPL(mm) SET_IMPL_INNER0(mm, MM_MPFQ_LAYER, MM_IMPL)

#define REBIND_ALL(mm) do {				\
        REBIND_F(mm, build_cache);		        \
        REBIND_F(mm, reload_cache);			\
        REBIND_F(mm, save_cache);		        \
        REBIND_F(mm, mul);			        \
        REBIND_F(mm, report);				\
        REBIND_F(mm, clear);				\
        REBIND_F(mm, init);				\
        REBIND_F(mm, auxv);				\
        REBIND_F(mm, aux);			       	\
        SET_IMPL(mm);                                   \
    } while (0)

#ifdef __cplusplus
extern "C" {
#endif

/* those are defined on a per-impl basis, in the impl specific source
 * files (most often it's called matmul-MM_IMPL.c , or .cpp) */
extern matmul_ptr MATMUL_NAME(init)(void *, param_list pl, int);
extern void MATMUL_NAME(build_cache)(matmul_ptr, uint32_t *);
extern int MATMUL_NAME(reload_cache)(matmul_ptr);
extern void MATMUL_NAME(save_cache)(matmul_ptr);
extern void MATMUL_NAME(mul)(matmul_ptr, void *, const void *, int);
extern void MATMUL_NAME(report)(matmul_ptr, double scale);
extern void MATMUL_NAME(clear)(matmul_ptr mm);
extern void MATMUL_NAME(auxv)(matmul_ptr mm, int op, va_list ap);
extern void MATMUL_NAME(aux)(matmul_ptr mm, int op, ...);

/* this one is defined collectively for all impls in matmul_facade.c */
extern void MATMUL_NAME(rebind)(matmul_ptr mm);

#ifdef __cplusplus
}
#endif

#endif	/* MATMUL_FACADE_H_ */
