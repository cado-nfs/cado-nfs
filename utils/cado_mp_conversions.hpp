#ifndef UTILS_CADO_MP_CONVERSIONS_HPP_
#define UTILS_CADO_MP_CONVERSIONS_HPP_

#include "cado_config.h"        // IWYU pragma: keep

#include <limits>
#include <complex>
#include <algorithm>

#include <gmp.h>

#include "cxx_mpz.hpp"
#include "cado_type_traits.hpp"
#include "cado_math_aux.hpp"

#ifdef HAVE_MPFR
#include <mpfr.h>
#include "cxx_mpfr.hpp"
#endif

#ifdef HAVE_MPC
#include <mpc.h>
#include "cxx_mpc.hpp"
#endif


namespace cado_math_aux {
    /* for floating point types, std::numeric_limits<T>::digits counts
     * the implicit bit as well (when there is one -- IEEE 80-bit
     * extended precision doesn't have one).
     *
     * E.g. on sysv 64-bit abi (meaning standard linux):
     * static_assert(std::numeric_limits<double>::digits==53, "AA");
     * static_assert(std::numeric_limits<long double>::digits==64, "AA");
     */

    template<typename T>
        static inline cxx_mpz mpz_from(T c)
        requires std::is_floating_point_v<T>
        {
            /* This converts to an mpz integer with unit accuracy (of
             * course digits below the unit are lost). Rounding is
             * towards zero.
             *
             * It seems that we don't need to fiddle with the rounding
             * mode here, the default truncation behavior of type casts
             * (C99 § 6.3.1.4.1) is what we need.
             */
            if (c == 0) return 0;
            int e;
            T x = cado_math_aux::frexp(c, &e);
            /* x == c * 2^-e */
            /* x is in (-1,0.5], [0.5,1) */

            constexpr int mantissa_bits = std::numeric_limits<T>::digits;
            constexpr int ui_bits = std::numeric_limits<unsigned long>::digits;
            /* First convert the complete mantissa to a mantissa_bits-bit
             * signed integer. This is a constant number of conversion
             * operations, decided at compile time (and in most cases
             * it's actually just one operation, except for quadruple
             * precision long doubles on 64-bit platforms, or doubles and
             * beyond on 32-bit platforms).
             */
            cxx_mpz z = 0;

            int sx = sgn(x);

            x *= sx;
            /* now x is in [0.5, 1) -- we could possibly swallow one
             * extra bit but that would be ridiculous gains since the
             * number of iterations below would likely be unchanged */

            for(int b = 0 ; b < mantissa_bits ; b += ui_bits) {
                /* at this point we have
                 * sx * c * 2^(-e+b) = z + (x) with:
                 *
                 *  - sx in {-1,0,1}
                 *  - x in [0,1)
                 *
                 */

                const int ell = std::min(mantissa_bits - b, ui_bits);
                const T x_ell = cado_math_aux::ldexp(x, ell);
                const auto xr = (unsigned long) x_ell;
                /* x*(2^ell) is in (0, 2^ell). We must
                 * truncate it (not round to nearest) if we want an
                 * unsigned long,
                 * since ell <= ui_bits.
                 *
                 * Furthermore, we have
                 * eps = x*(2^ell)-xr in [0,1)
                 */

                /* sx*c*2^(-e+b) = z + (xr+eps)/2^ell
                 * sx*c*2^(-e+b+ell) = z*2^ell + xr + eps
                 */
                cxx_mpz zi = xr;

                mpz_mul_2exp(z, z, ell);
                mpz_add(z, z, zi);

                /* prepare next iteration. I have small fears that there
                 * could be corner cases that allow the compiler to
                 * outsmart us. */
                x = x_ell - T(xr);
            }

            if (e > mantissa_bits) {
                mpz_mul_2exp(z, z, e - mantissa_bits);
            } else {
                mpz_tdiv_q_2exp(z, z, mantissa_bits - e);
            }
            if (sx < 0)
                mpz_neg(z, z);
            return z;
        }

    template<typename T>
    static inline cxx_mpz mpz_from(T const & c)
    requires std::is_same_v<T, cxx_mpz>
    { return c; }

#ifdef HAVE_MPFR
    template<typename T>
    static inline cxx_mpz mpz_from(T const & c)
    requires std::is_same_v<T, cxx_mpfr>
    {
        cxx_mpz res;
        mpfr_get_z(res, c, MPFR_RNDN);
        return res;
    }
#endif

    template<typename T>
    static inline void exact_form(cxx_mpz & m, int & e, T x)
    requires cado_math_aux::is_real_v<T>
    {
        e = 0;
        /* frexp returns x in (-1,-0.5] u [0.5,1). We want an integer, so
         * we want to scale this up */
        T mantissa = cado_math_aux::frexp(x, &e);
        constexpr int d = std::numeric_limits<T>::digits;
        m = mpz_from<T>(cado_math_aux::ldexp(mantissa, d-1));
        e -= (d-1);
    }

    template<typename T>
    static inline void exact_form(cxx_mpz & m, int & e, T const & x)
    requires cado_math_aux::is_integral_v<T>
    {
        m = x;
        e = 0;
    }
    template<typename T>
        T ulp(T r)
        requires std::is_floating_point_v<T>
        {
        /*
         * Let next(r) be the smallest floating point
         * number strictly larger than r.
         *
         * ilog(b) computes I such that |r| * 2^-I is in [1,2)
         *
         * we want ulp(r) = next(|r|)-|r|
         * we also have ulp(r)/2^I = next(|r|/2^I)-|r|
         * because |r| is in the same binade as 1,
         *  ulp(r)/2^I == next(1)-1 = epsilon().
         */
        return cado_math_aux::ldexp(std::numeric_limits<T>::epsilon(), std::ilogb(r));
    }

    template<typename T>
        struct converter_from_mpz;

    /* This is the same as mpz_get_ld, but should get more mantissa bits
     * correct if we are to use it with quad precision floating points.
     */
    template<typename T>
        requires std::is_floating_point_v<T>
        struct converter_from_mpz<T> {
            T operator()(mpz_srcptr z) const {
                T ld = 0;
                cxx_mpz zr = z;

                int b = mpz_sizeinbase(z, 2);
                int e = b - std::numeric_limits<double>::max_exponent;
                if (e > 0) {
                    /* First scale the input number down. There are still enough
                     * bits remaining to fill the mantissa, of course.
                     */
                    mpz_tdiv_q_2exp(zr, zr, e);
                } else {
                    e = 0;
                }

                constexpr int M = std::numeric_limits<double>::digits;
                constexpr int LM = std::numeric_limits<T>::digits;
                cxx_mpz t;
                for(int b = 0 ; b + M < LM ; b += M) {
                    double d = mpz_get_d (zr);
                    mpz_set_d (t, d);
                    mpz_sub (zr, zr, t);
                    ld += d;
                }
                ld += mpz_get_d (zr);
                ld = cado_math_aux::ldexp(ld, e);
                return ld;
            }
        };

    template<typename T>
        requires cado::converts_via<T, int64_t>
        struct converter_from_mpz<T> {
            T operator()(mpz_srcptr z) const {
                return mpz_get_int64(z);
            }
        };

    template<typename T>
        requires cado::converts_via<T, uint64_t>
        struct converter_from_mpz<T> {
            T operator()(mpz_srcptr z) const {
                return mpz_get_uint64(z);
            }
        };

    template<typename T>
        struct converter_from_mpz<std::complex<T>> {
            std::complex<T> operator()(mpz_srcptr z) const {
                return { converter_from_mpz<T>()(z), 0 };
            }
        };

    template<>
        struct converter_from_mpz<cxx_mpz> {
            cxx_mpz operator()(mpz_srcptr z) const {
                return { z };
            }
        };

    
#ifdef HAVE_MPFR
    template<>
        struct converter_from_mpz<cxx_mpfr> {
            cxx_mpfr operator() (mpz_srcptr z) const
            {
                cxx_mpfr res;
                mpfr_set_z(res, z, MPFR_RNDN);
                return res;
            }
        };
#endif

#ifdef HAVE_MPC
    template<>
        struct converter_from_mpz<cxx_mpc> {
            cxx_mpc operator() (mpz_srcptr z) const
            {
                cxx_mpc res;
                mpc_set_z(res, z, MPFR_RNDN);
                return res;
            }
        };
#endif

    template<typename T>
        static inline T mpz_get(mpz_srcptr z) {
            return converter_from_mpz<T>()(z);
        }

} /* namespace cado_math_aux */


#endif	/* UTILS_CADO_MP_CONVERSIONS_HPP_ */
