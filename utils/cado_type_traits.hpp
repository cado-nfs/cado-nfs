#ifndef UTILS_COERCION_HPP_
#define UTILS_COERCION_HPP_

#include "cado_config.h"        // IWYU pragma: keep

#include <complex>
#include <type_traits>
#include <array>

#include "cxx_mpz.hpp"
#ifdef HAVE_MPFR
#include "cxx_mpfr.hpp"
#endif
#ifdef HAVE_MPC
#include "cxx_mpc.hpp"
#endif

namespace cado {
    struct nop_function {
        template<typename... Args>
            void operator()(Args&& ...) const { }
    };

    template<typename A>
    concept is_std_array_v =
           requires { typename A::value_type; }
        && requires { std::tuple_size_v<A>; }
        && std::is_same_v<A, std::array<typename A::value_type, std::tuple_size_v<A>>>;

    static_assert(is_std_array_v<std::array<int, 3>>);
    static_assert(!is_std_array_v<int[3]>);

    template<typename A>
    concept is_bounded_array_or_std_array_v =
        (std::is_bounded_array_v<A> || is_std_array_v<A>);

    static_assert(is_bounded_array_or_std_array_v<std::array<int, 3>>);
    static_assert(is_bounded_array_or_std_array_v<int[3]>);


    template<typename T> struct make_signed;
    template<typename T> struct make_unsigned;

    template<typename T> requires requires { typename std::make_signed<T>; }
    struct make_signed<T> : public std::make_signed<T> {};
    template<typename T> requires requires { typename std::make_unsigned<T>; }
    struct make_unsigned<T> : public std::make_unsigned<T> {};

    template<typename T> using make_signed_t = make_signed<T>::type;
    template<typename T> using make_unsigned_t = make_unsigned<T>::type;

    /* because this header file is always included after cxx_mpz, let's
     * move the definitions here.
     */
    template<> struct make_signed<cxx_mpz> { using type = cxx_mpz; };
    template<> struct make_unsigned<cxx_mpz> { using type = cxx_mpz; };
    template<> struct make_signed<cxx_mpq> { using type = cxx_mpq; };
    template<> struct make_unsigned<cxx_mpq> { using type = cxx_mpq; };
} /* namespace cado */


namespace cado_math_aux {

namespace is_coercible_details {
    template<typename T, typename U> struct impl_strict;
    template<typename T, typename U> struct impl_loose;

    template<typename From, typename To> struct impl_strict : public std::false_type {};
    // { static constexpr bool value = false; };

    template<> struct impl_strict<int, long> : public std::true_type {};

    template<> struct impl_strict<int, cxx_mpz> : public std::true_type {};
    template<> struct impl_strict<long, cxx_mpz> : public std::true_type {};

    template<> struct impl_strict<int, float> : public std::true_type {};
    template<> struct impl_strict<long, float> : public std::true_type {};
    template<> struct impl_strict<cxx_mpz, float> : public std::true_type {};

    template<> struct impl_strict<int, double> : public std::true_type {};
    template<> struct impl_strict<long, double> : public std::true_type {};
    template<> struct impl_strict<cxx_mpz, double> : public std::true_type {};
    template<> struct impl_strict<float, double> : public std::true_type {};

    template<> struct impl_strict<int, long double> : public std::true_type {};
    template<> struct impl_strict<long, long double> : public std::true_type {};
    template<> struct impl_strict<cxx_mpz, long double> : public std::true_type {};
    template<> struct impl_strict<float, long double> : public std::true_type {};
    template<> struct impl_strict<double, long double> : public std::true_type {};

#ifdef HAVE_MPFR
    template<> struct impl_strict<int, cxx_mpfr> : public std::true_type {};
    template<> struct impl_strict<long, cxx_mpfr> : public std::true_type {};
    template<> struct impl_strict<cxx_mpz, cxx_mpfr> : public std::true_type {};
    template<> struct impl_strict<float, cxx_mpfr> : public std::true_type {};
    template<> struct impl_strict<double, cxx_mpfr> : public std::true_type {};
    template<> struct impl_strict<long double, cxx_mpfr> : public std::true_type {};
#endif
    template<typename U> struct impl_strict<int, std::complex<U>> : public std::true_type {};
    template<typename U> struct impl_strict<long, std::complex<U>> : public std::true_type {};
    template<typename T, typename U>
    requires impl_loose<T, U>::value
    struct impl_strict<T, std::complex<U>> : public std::true_type {};
#ifdef HAVE_MPC
    template<> struct impl_strict<int, cxx_mpc> : public std::true_type {};
    template<> struct impl_strict<long, cxx_mpc> : public std::true_type {};
    template<> struct impl_strict<cxx_mpz, cxx_mpc> : public std::true_type {};
    template<> struct impl_strict<float, cxx_mpc> : public std::true_type {};
    template<> struct impl_strict<double, cxx_mpc> : public std::true_type {};
    template<> struct impl_strict<cxx_mpfr, cxx_mpc> : public std::true_type {};
#endif

    template<typename From, typename To> struct impl_loose : public impl_strict<From, To> {};

    template<typename T> struct impl_loose<T, T> : public std::true_type {};

} /* namespace is_coercible_details */

template<typename From, typename To>
using is_coercible = is_coercible_details::impl_loose<From, To>;

template<typename From, typename To>
using is_strictly_coercible = is_coercible_details::impl_strict<From, To>;

template<typename From, typename To>
inline constexpr bool is_coercible_v = is_coercible<From, To>::value;

template<typename From, typename To>
inline constexpr bool is_strictly_coercible_v = is_strictly_coercible<From, To>::value;

template<typename From, typename To>
using is_coercible_t = std::enable_if_t<is_coercible_v<From, To>, bool>;

template<typename From, typename To>
using is_strictly_coercible_t = std::enable_if_t<is_strictly_coercible_v<From, To>, bool>;


template<typename T> struct is_integral : public std::is_integral<T> {};

template<> struct is_integral<cxx_mpz> : public std::true_type {};
template<typename X>
inline constexpr bool is_integral_v = is_integral<X>::value;

/* A type is real if it makes sense to do Newton iterations and it's
 * totally ordered.
 */
template<typename T> struct is_real : public std::false_type {};
template<> struct is_real<float> : public std::true_type {};
template<> struct is_real<double> : public std::true_type {};
template<> struct is_real<long double> : public std::true_type {};
#ifdef HAVE_MPFR
template<> struct is_real<cxx_mpfr> : public std::true_type {};
#endif
template<typename X>
inline constexpr bool is_real_v = is_real<X>::value;
template<typename X>
using is_real_t = std::enable_if_t<is_real_v<X>, bool>;

/* A complex type requires algorithms such as Jenkins-Traub for
 * rootfinding.
 */
template<typename T> struct is_complex : public std::false_type {};
template<> struct is_complex<std::complex<float>> : public std::true_type {};
template<> struct is_complex<std::complex<double>> : public std::true_type {};
template<> struct is_complex<std::complex<long double>> : public std::true_type {};
#ifdef HAVE_MPC
template<> struct is_complex<cxx_mpc> : public std::true_type {};
#endif
template<typename X>
inline constexpr bool is_complex_v = is_complex<X>::value;
template<typename X>
using is_complex_t = std::enable_if_t<is_complex_v<X>, bool>;

static_assert(is_coercible_v<int, long>);
static_assert(!is_coercible_v<double, long>);
static_assert(is_strictly_coercible_v<int, long>);
static_assert(!is_strictly_coercible_v<int, int>);
static_assert(is_strictly_coercible_v<cxx_mpz, double>);
static_assert(is_strictly_coercible_v<cxx_mpz, std::complex<double>>);
static_assert(is_strictly_coercible_v<double, std::complex<double>>);

} /* namespace cado_math_aux */

#endif	/* UTILS_COERCION_HPP_ */
