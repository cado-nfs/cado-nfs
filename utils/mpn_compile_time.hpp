#ifndef MPN_COMPILE_TIME_HPP_
#define MPN_COMPILE_TIME_HPP_

#include <stddef.h>
#include <gmp.h>
#include <type_traits>
#include "gmp-hacks.h"

/* This header provides overloads for some basic mpn functions such as:
 *      mpn_copyi
 *      mpn_copyd
 *      mpn_zero
 *      mpn_zero_p
 *      mpn_cmp
 *      mpn_cmp_ui (well, it doesn't exist but we provide it)
 *
 * in such a way that the size_t argument can be an
 * std::integral_constant or a runtime_constant. In the former case,
 * simple header-only loops are emitted, so that the compiler can avoid
 * the overhead of calling the gmp function.
 */
namespace mpn_compile_time {

    namespace details {
        template<typename X>
            struct mpn_helper;

        template<>
            struct mpn_helper<size_t> {
                /* On the N>0 requirement, the GMP doc says:
                 *   A common requirement for all functions is that each
                 *   source area needs at least one limb.  No size argument
                 *   may be zero.
                 */
                static inline void copyi(mp_limb_t * dst, const mp_limb_t * src, size_t N) {
                    ASSERT(N);
                    ::mpn_copyi(dst, src, N);
                }
                static inline void copyd(mp_limb_t * dst, const mp_limb_t * src, size_t N) {
                    ASSERT(N);
                    ::mpn_copyd(dst, src, N);
                }
                static inline void zero(mp_limb_t * dst, size_t N) {
                    ASSERT(N);
                    ::mpn_zero(dst, N);
                }
                static inline int zero_p(mp_limb_t const * x, size_t N) {
                    ASSERT(N);
                    return ::mpn_zero_p(x, N);
                }
                static inline int cmp(mp_limb_t const * a, mp_limb_t const * b, size_t N) {
                    ASSERT(N);
                    return ::mpn_cmp(a, b, N);
                }
                static inline int cmp_ui(mp_limb_t const * x, unsigned long a, mp_size_t N) {
                    ASSERT(N);
                    int r = (a < x[0]) - (x[0] < a);
                    if (r || N==1) return r;
                    return !zero_p(x + 1, N - 1);
                }
                static inline void set_ui(mp_limb_t * x, unsigned long a, mp_size_t N) {
                    ASSERT(N);
                    x[0] = a;
                    if (N==1) return;
                    ::mpn_zero(x + 1, N - 1);
                }
                static inline void SET_MPZ(mp_limb_t * DST, size_t NLIMBS, mpz_srcptr SRC) {
                    ASSERT(NLIMBS);
                    return ::MPN_SET_MPZ(DST, NLIMBS, SRC);
                }
            };

        /* TODO: avoid unrolling when N is a large integral constant */
        template<size_t N>
            struct mpn_helper<std::integral_constant<size_t, N>> {
                typedef std::integral_constant<size_t, N> SIZE_T;
                static inline void copyi(mp_limb_t * dst, const mp_limb_t * src, SIZE_T)
                {
                    /* This is a constant loop */
                    for(size_t i = 0 ; i < N ; i++)
                        dst[i] = src[i];
                }
                static inline void copyd(mp_limb_t * dst, const mp_limb_t * src, SIZE_T)
                {
                    /* This is a constant loop */
                    for(size_t i = N ; i-- ; )
                        dst[i] = src[i];
                }
                static inline void zero(mp_limb_t * dst, SIZE_T)
                {
                    /* This is a constant loop */
                    for(size_t i = 0 ; i < N ; i++)
                        dst[i] = 0;
                }
                static inline int zero_p(mp_limb_t const * x, SIZE_T) {
                    /* This is a constant loop */
                    for(size_t i = 0 ; i < N ; i++)
                        if (x[i]) return 0;
                    return 1;
                }
                static inline int cmp(mp_limb_t const * a, mp_limb_t const * b, SIZE_T) {
                    /* This is a constant loop */
                    for(size_t i = N ; i-- ; ) {
                        int r = (a[i] < b[i]) - (b[i] < a[i]);
                        if (r) return r;
                    }
                    return 0;
                }
                static inline int cmp_ui(mp_limb_t const * x, unsigned long a, SIZE_T) {
                    int r = (a < x[0]) - (x[0] < a);
                    if (r > 0) return r;
                    /* This is a constant loop */
                    for(size_t i = 1 ; i < N ; i++)
                        if (x[i]) return 1;
                    return r;
                }
                static inline void set_ui(mp_limb_t * x, unsigned long a, SIZE_T) {
                    x[0] = a;
                    ::mpn_zero(x + 1, N - 1);
                }
                static inline void SET_MPZ(mp_limb_t * DST, SIZE_T, mpz_srcptr SRC) {
                    /* Use the runtime version */
                    return ::MPN_SET_MPZ(DST, N, SRC);
                }
            };
    }

    template<typename X>
        inline void mpn_copyi(mp_limb_t * dst, const mp_limb_t * src, X N) {
            return details::mpn_helper<X>::copyi(dst, src, N);
        }
    template<typename X>
        inline void mpn_copyd(mp_limb_t * dst, const mp_limb_t * src, X N) {
            return details::mpn_helper<X>::copyd(dst, src, N);
        }
    template<typename X>
        inline void mpn_zero(mp_limb_t * dst, X N) {
            return details::mpn_helper<X>::zero(dst, N);
        }
    template<typename X>
        inline int mpn_zero_p(mp_limb_t const * x, X N) {
            return details::mpn_helper<X>::zero_p(x, N);
        }
    template<typename X>
        inline int mpn_cmp(mp_limb_t const * a, mp_limb_t const * b, X N) {
            return details::mpn_helper<X>::cmp(a, b, N);
        }
    template<typename X>
        inline int mpn_cmp_ui(mp_limb_t const * x, unsigned long a, X N) {
            return details::mpn_helper<X>::cmp_ui(x, a, N);
        }
    template<typename X>
        inline void mpn_set_ui(mp_limb_t * x, unsigned long a, X N) {
            details::mpn_helper<X>::set_ui(x, a, N);
        }
    template<typename X>
        inline void MPN_SET_MPZ(mp_limb_t * dst, X N, mpz_srcptr src)
        {
            details::mpn_helper<X>::SET_MPZ(dst, N, src);
        }
}


#endif	/* MPN_COMPILE_TIME_HPP_ */
