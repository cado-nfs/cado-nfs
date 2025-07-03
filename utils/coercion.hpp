#ifndef UTILS_COERCION_HPP_
#define UTILS_COERCION_HPP_

#include "cado_config.h"

#include <type_traits>

#include "cxx_mpz.hpp"
#ifdef HAVE_MPFR
#include "cxx_mpfr.hpp"
#endif
#ifdef HAVE_MPC
#include "cxx_mpc.hpp"
#endif

namespace is_coercible_details {
    template<typename From, typename To> struct impl : public std::false_type {};

    template<> struct impl<int, long> : public std::true_type {};

    template<> struct impl<int, cxx_mpz> : public std::true_type {};
    template<> struct impl<long, cxx_mpz> : public std::true_type {};

    template<> struct impl<int, float> : public std::true_type {};
    template<> struct impl<long, float> : public std::true_type {};
    template<> struct impl<cxx_mpz, float> : public std::true_type {};

    template<> struct impl<int, double> : public std::true_type {};
    template<> struct impl<long, double> : public std::true_type {};
    template<> struct impl<cxx_mpz, double> : public std::true_type {};
    template<> struct impl<float, double> : public std::true_type {};

#ifdef HAVE_MPFR
    template<> struct impl<int, cxx_mpfr> : public std::true_type {};
    template<> struct impl<long, cxx_mpfr> : public std::true_type {};
    template<> struct impl<cxx_mpz, cxx_mpfr> : public std::true_type {};
    template<> struct impl<float, cxx_mpfr> : public std::true_type {};
    template<> struct impl<double, cxx_mpfr> : public std::true_type {};
#endif
#ifdef HAVE_MPC
    template<> struct impl<int, cxx_mpc> : public std::true_type {};
    template<> struct impl<long, cxx_mpc> : public std::true_type {};
    template<> struct impl<cxx_mpz, cxx_mpc> : public std::true_type {};
    template<> struct impl<float, cxx_mpc> : public std::true_type {};
    template<> struct impl<double, cxx_mpc> : public std::true_type {};
    template<> struct impl<cxx_mpfr, cxx_mpc> : public std::true_type {};
#endif



} /* namespace is_coercible_details */

template<typename From, typename To>
using is_coercible = typename is_coercible_details::impl<From, To>::value;



#endif	/* UTILS_COERCION_HPP_ */
