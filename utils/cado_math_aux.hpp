#ifndef CADO_UTILS_CADO_MATH_AUX_HPP
#define CADO_UTILS_CADO_MATH_AUX_HPP

#include <cmath>
#include <cfenv>
#include <cstddef>
#include <cstdint>
#include <climits>

#include <algorithm>
#include <array>
#include <limits>
#include <type_traits>

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
            constexpr T operator()(T x) const {
                return multiply_by_poweroftwo<T, n-1>()(2*x);
            }
        };

    template<typename T>
        struct multiply_by_poweroftwo<T, 0> {
            constexpr T operator()(T x) const {
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

    template<typename T> static inline void do_not_outsmart_me(T &) {}
#if defined(__i386)
    template<> inline void do_not_outsmart_me<double>(double & x) {
        volatile double mx = x; x = mx;
    }
#endif

    template <typename T> int sgn(T val) {
        return (T(0) < val) - (val < T(0));
    }

    struct temporary_round_mode {
        int saved;
        explicit temporary_round_mode(int mode)
            : saved(fegetround())
        {
            fesetround(mode);
        }
        temporary_round_mode(temporary_round_mode const &) = delete;
        temporary_round_mode(temporary_round_mode &&) = delete;
        temporary_round_mode & operator=(temporary_round_mode const &) = delete;
        temporary_round_mode & operator=(temporary_round_mode &&) = delete;
        ~temporary_round_mode() { fesetround(saved); }
    };

    static inline bool rounding_towards_zero_works()
    {
        /* This runtime check can detect dysfunctional math environments.
         * valgrind is one of them, unfortunately.
         */
        temporary_round_mode dummy(FE_TOWARDZERO);
        const double d0 = (uint64_t(1) << 52) + 1;
        const double d1 = std::ldexp(d0, 53);
        /* We can't initialize it as d = d0 + d1 because that wouldn't
         * necessarily obey the rounding mode
         */
        volatile double d = d0;
        d += d1;
        return d == d1;
    }

    static inline bool valgrind_long_double_hopeless()
    {
        /* Another one. This time, the result is even more useless. Note
         * that in fact, I'm not even sure that frexp _works_ under
         * valgrind. ld prints as 0x8p+16381 and frexp returns 0. Both
         * are clearly bogus.
         */
        volatile double d0 = 1.79769313486231570815e+308; // 0x1.fffffffffffffp+1023;
        volatile double d1 = 1.9958403095347196e+292;     // 0x1.fffffffffffffp+970;
        volatile long double ld = d0;
        ld += d1;
        int e;
        std::frexp(ld, &e);
        return e != 1025;
    }


    /* for floating point types, std::numeric_limits<T>::digits counts
     * the implicit bit as well (when there is one -- IEEE 80-bit
     * extended precision doesn't have one).
     *
     * E.g. on sysv 64-bit abi (meaning standard linux):
     * static_assert(std::numeric_limits<double>::digits==53, "AA");
     * static_assert(std::numeric_limits<long double>::digits==64, "AA");
     */

    template<typename T>
        cxx_mpz mpz_from(T c)
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
            T x = std::frexp(c, &e);
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
                const T x_ell = std::ldexp(x, ell);
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
    void exact_form(cxx_mpz & m, int & e, T x)
    {
        int xe = 0;
        T mantissa = std::frexp(x, &xe);
        constexpr int d = std::numeric_limits<T>::digits;
        m = mpz_from<T>(std::ldexp(mantissa, d-1));
        e -= (d-1);
    }

    template<typename T>
    typename std::enable_if<std::is_floating_point<T>::value, T>::type
    ulp(T r) {
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
        return std::ldexp(std::numeric_limits<T>::epsilon(), std::ilogb(r));
    }

    /* This is the same as mpz_get_ld, but should get more mantissa bits
     * correct if we are to use it with quad precision floating points.
     */
    template<typename T>
    typename std::enable_if<std::is_floating_point<T>::value, T>::type
    mpz_get (mpz_srcptr z)
    {
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
        ld = std::ldexp(ld, e);
        return ld;
    }



    template<typename T>
        class constant_time_square_root {
            static constexpr T mid(T a, T b) { return (a+b)/2; }
            static constexpr T above(T m, T x) { return m*m > x; }
            static constexpr T recurse(T a, T b, T m, T x)
            { return above(m, x) ? sqrt(a, m, x) : sqrt(m, b, x); }
            static constexpr T sqrt(T a, T b, T m, T x)
            { return m == a ? a : recurse(a, b, m, x); }
            static constexpr T sqrt(T a, T b, T x)
            { return sqrt(a, b, mid(a, b), x); }
            public:
            static constexpr typename std::enable_if<std::is_integral<T>::value, T>::type sqrt(T x)
            { return sqrt(T(0), (T(1) << (std::numeric_limits<T>::digits/2)), x); }
        };

    template<typename T>
        static constexpr typename std::enable_if<std::is_integral<T>::value, T>::type constant_sqrt(T x) {
            return constant_time_square_root<T>::sqrt(x);
        }

    /* some fun with compile time code */
    /* {{{ compile time evaluation of 2^k mod n */
    template<int k, int n> struct pow2_mod {
        struct even {
            static constexpr int value() {
                typedef pow2_mod<k/2, n> T;
                return (T::value() * T::value()) % n;
            }
        };
        struct odd {
            static constexpr int value() { return (even::value() + even::value()) % n; }
        };
        static constexpr int value() {
            return (k % 2) ? odd::value() : even::value();
        };
    };
    template<int n> struct pow2_mod<0, n>
    {
        static constexpr int value() { return 1; }
    };
    /* }}} */
    /* {{{ compile time bezout relation */
    template<int a, int b,
        int u0 = 1, int v0 = 0, int g0 = a,
        int u1 = 0, int v1 = 1, int g1 = b>
            struct bezout_relation;
    template<int a, int b,
        int u0, int v0, int g0,
        int u1, int v1, int g1>
            struct bezout_relation : public bezout_relation<
                                     a, b,
                                     u1, v1, g1, 
                                     u0 - (g0/g1) * u1, v0 - (g0 / g1) * v1, g0 % g1>
    {};
    template<int a, int b,
        int u0, int v0, int g0,
        int u1, int v1>
            struct bezout_relation<a, b, u0, v0, g0, u1, v1, 0> {
                static constexpr int gcd() { return g0; }
                static constexpr int u() { return u0; }
                static constexpr int v() { return v0; }
            };
    /* }}} */
    /* {{{ compile-time map on an fixed-size array */
    template<int N, template <int, int> class F, int left=N, int... Rest>
        struct array_map : public array_map<N, F, left - 1, F<N, left - 1>::value(), Rest...> {
        };

    template<int N, template <int, int> class F, int... Rest>
        struct array_map<N, F, 0, Rest...> {
            static constexpr
                std::array<int, N>
                value { Rest... };
        };

    template<int N, template<int, int> class F, int... Rest>
        constexpr std::array<int, N> array_map<N, F, 0, Rest...>::value; // c++11
    /* }}} */
    /* {{{ compute -1/i mod n at compile time */
    template<int n, int k> struct minus_inverse_mod_n {
        typedef bezout_relation<n, n-k> B;
        static constexpr int value() { return B::gcd() == 1 ? (B::v() < 0 ? n + B::v() : B::v()) : 0; }
    };
    /* }}} */
    /* {{{ inverse modulo powers of two */
    template<int k, typename T, T current = (3*k)^2, int c = 5, bool done = c >= std::numeric_limits<T>::digits>
        struct invmod;

    template<int k, typename T, T current, int c>
        struct invmod<k, T, current, c, false>
        : public invmod<k, T, current * (2 - current * k), c*2>
        {
            static_assert(k & 1, "k must be odd");
            static_assert(!std::is_signed<T>::value, "T must be an unsigned type");
        };
    template<int k, typename T, T current, int c>
        struct invmod<k, T, current, c, true> {
            static constexpr T value() { return current; }
        };
    /* }}} */

    /* for any type that defines the <<= and += operators, provide code
     * that computes a constant multiple of a given value (of any
     * acceptable input type for +=).
     * Optionally, a chooser template (such as at_most<8>) can be passed
     * in order to control when to fall back on multiplication at
     * runtime. If this evaluates to false, then a *= operator on the
     * destination type is emitted.
     */

    namespace details {
        template<unsigned int n, bool is_odd = n & 1>
            struct mul_c_impl;

        template<unsigned int n>
            struct mul_c_impl<n, false> {
                typedef mul_c_impl<n/2> super;
                template<typename T, typename U>
                    static void mul(T & t, U const & a) {
                        super::mul(t, a);
                        t <<= 1;
                    }
                static constexpr int number_of_shifts() {
                    return 1 + super::number_of_shifts();
                }
                static constexpr int number_of_additions() {
                    return super::number_of_shifts();
                }
            };
        template<unsigned int n>
            struct mul_c_impl<n, true> {
                typedef mul_c_impl<n/2> super;
                template<typename T, typename U>
                    static void mul(T & t, U const & a) {
                        super::mul(t, a);
                        t <<= 1;
                        t += a;
                    }
                static constexpr int number_of_shifts() {
                    return 1 + super::number_of_shifts();
                }
                static constexpr int number_of_additions() {
                    return 1 + super::number_of_shifts();
                }
            };

        template<>
            struct mul_c_impl<1> {
                template<typename T, typename U>
                    static void mul(T & t, U const & a) {
                        t = a;
                    }
                static constexpr int number_of_shifts() {
                    return 0;
                }
                static constexpr int number_of_additions() {
                    return 0;
                }
            };

        template<>
            struct mul_c_impl<0> {
                template<typename T, typename U>
                    static void mul(T & t, U const &) {
                        t = 0;
                    }
                static constexpr int number_of_shifts() {
                    return 0;
                }
                static constexpr int number_of_additions() {
                    return 0;
                }
            };

        template<unsigned int n, bool is_small>
            struct mul_c_impl_choice
            : public mul_c_impl<n> {
            };

        template<unsigned int n>
            struct mul_c_impl_choice<n, false> {
                template<typename T, typename U>
                    static void mul(T & t, U const & a) {
                        t = a;
                        t *= n;
                    }
            };


        /* addmul is different because it doesn't use double-and-add.
         * It's not very useful, so I'm commenting it out.
         */
#if 0
        template<unsigned int n>
            struct addmul_c_impl {
                template<typename T, typename U>
                    static void addmul(T & t, U const & a) {
                        addmul_c_impl<n-1>::addmul(t, a);
                        t += a;
                    }
            };
        template<>
            struct addmul_c_impl<0> {
                template<typename T, typename U>
                    static void addmul(T &, U const &) {
                    }
            };

        template<unsigned int n, bool is_very_small, typename chooser_mul>
            struct addmul_c_impl_choice
            : public addmul_c_impl<n> {};

        template<unsigned int n, typename chooser_mul>
            struct addmul_c_impl_choice<n, false, chooser_mul> {
                template<typename T, typename U>
                    static void addmul(T & t, U const & a) {
                        T t2;
                        mul_c<n, chooser_mul>(t2, a);
                        t += t2;
                    }
            };
#endif
    }

    /* Three chooser types. Evaluation to true means "use an addition
     * chain".
     */
    template<unsigned int max>
        struct at_most {
            template<unsigned int n> struct choose : public
                                                     std::integral_constant<bool, (n <= max)> {
                                                     };
        };

    struct always {
        template<unsigned int n> struct choose : public
                                                 std::true_type {
                                                 };
    };

    struct never {
        template<unsigned int n> struct choose : public
                                                 std::false_type {
                                                 };
    };

    template<unsigned int n, typename chooser = at_most<8>, typename T, typename U>
        void mul_c(T & t, U const & a, chooser const & = chooser{})
        {
            details::mul_c_impl_choice<n, chooser::template choose<n>::value>::mul(t, a);
        }
#if 0
    /* see above. the code is fine, but I'm yet to be convinced that it's
     * useful
     */
    template<unsigned int n, typename chooser = at_most<2>, typename chooser_mul, typename T, typename U>
        void addmul_c(T & t, U const & a, chooser const & = {}, chooser_mul const & = {})
        {
            details::addmul_c_impl_choice<n, chooser::template choose<n>::value, chooser_mul>::addmul(t, a);
        }
#endif
}


#endif	/* CADO_UTILS_CADO_MATH_AUX_HPP */
