#ifndef UTILS_CADO_MATH_AUX_HPP_
#define UTILS_CADO_MATH_AUX_HPP_

#include <cstdint>
#include <cmath>
#include <cstddef>
#include <limits>

#include <gmp.h>
#include "cxx_mpz.hpp"
#include "macros.h"
#include "gmp_aux.h"

namespace cado_math_aux
{
    /*
    template<typename T> T abs(T x);
    template<> inline double abs<double>(double x) { return fabs(x); }
    template<> inline float abs<float>(float x) { return fabsf(x); }
    template<> inline long double abs<long double>(long double x) { return fabsl(x); }
    */

    template<typename T, int n>
        struct multiply_by_poweroftwo {
            constexpr T operator()(T x) {
                return multiply_by_poweroftwo<T, n-1>()(2*x);
            }
        };

    template<typename T>
        struct multiply_by_poweroftwo<T, 0> {
            constexpr T operator()(T x) {
                return x;
            }
        };

    /* This compares two floating point values with a relative error margin
    */
    template<typename T>
        static bool approx_eq_relative(
                const T d1,
                const T d2,
                const T err_margin)
        {
            using namespace cado_math_aux;
            return abs(d1) * (1. - err_margin) <= abs(d2) && abs(d2) <= abs(d1) * (1. + err_margin);
        }


    TEMPLATE_ENABLED_ON_TEMPLATE_ARG(
            typename T,
            !std::numeric_limits<T>::is_integer)
        static bool equal_within_ulps(T x, T y, std::size_t n)
        {
            const T m = std::min(std::fabs(x), std::fabs(y));
            const int exp = m < std::numeric_limits<T>::min()
                ? std::numeric_limits<T>::min_exponent - 1
                : std::ilogb(m);
            return std::fabs(x - y) <= n * std::ldexp(std::numeric_limits<T>::epsilon(), exp);
        }

    TEMPLATE_ENABLED_ON_TEMPLATE_ARG(
            typename T,
            !std::numeric_limits<T>::is_integer)
        static int accurate_bits(T reference, T computed)
        {
            ASSERT_ALWAYS(reference != 0);
            T c = std::fabs((computed-reference)/reference);
            return c == 0 ? INT_MAX : -std::ilogb(c);
        }

    template<typename T>
    void exact_form(cxx_mpz & m, int & e, T x)
    {
        int xe;
        T mantissa = std::frexp(x, &xe);
        constexpr int d = std::numeric_limits<T>::digits;
        mpz_set_d(m, multiply_by_poweroftwo<T, d-1>(mantissa));
        e -= (d-1);
    }

    template<typename T>
    T mpz_get(mpz_srcptr x);
    template<> inline float mpz_get<float>(mpz_srcptr x) { return mpz_get_d(x); }
    template<> inline double mpz_get<double>(mpz_srcptr x) { return mpz_get_d(x); }
    /* This one is in gmp_aux.h */
    template<> inline long double mpz_get<long double>(mpz_srcptr x) { return mpz_get_ld(x); }

    template<typename T> void do_not_outsmart_me(T &) {}
#if defined(__i386)
    template<> void do_not_outsmart_me<double>(double & x) {
        volatile double mx = x; x = mx;
    }
#endif

    template <typename T> int sgn(T val) {
        return (T(0) < val) - (val < T(0));
    }

    template<typename T>
        cxx_mpz mpz_from(T c)
        {
            /* This converts to an mpz integer with unit accuracy (of
             * course digits below the unit are lost)
             */
            if (c == 0) return 0;
            int e;
            T x = std::frexp(c, &e);
            /* x == c * 2^e */
            x -= sgn(c) * T(0.5);
            /* x + sgn(c)/2 == c * 2^e */
            /* x is in (-0.5, 0.5) */
            /* y is in (-2^(digits-1), 2^(digits-1)), which is always
             * castable to an int64_t */
            T y = std::ldexp(x, std::numeric_limits<T>::digits);
            /* y / 2^digits == c * 2^e - sgn(c)/2 */
            /* y == c * 2^(e-digits) - sgn(c)*2^(digits-1) */
            static_assert(std::numeric_limits<T>::digits <= 64);
            cxx_mpz z;
            mpz_set_si(z, sgn(c));
            mpz_mul_2exp(z, z, std::numeric_limits<T>::digits-1);
            /* z = sgn(c)*2^(digits-1) */
            mpz_add_int64(z, z, int64_t(y));
            if (e > std::numeric_limits<T>::digits) {
                mpz_mul_2exp(z, z, e - std::numeric_limits<T>::digits);
            } else {
                mpz_div_2exp(z, z, std::numeric_limits<T>::digits - e);
            }
            return z;
        }
}


#endif	/* UTILS_CADO_MATH_AUX_HPP_ */
