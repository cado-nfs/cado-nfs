#ifndef UTILS_CADO_CONSTANTS_HPP_
#define UTILS_CADO_CONSTANTS_HPP_

#include <cmath>

#include <complex>

#include <numbers>
#include "number_context.hpp"

#ifdef HAVE_MPFR
#include <mpfr.h>
#include "cxx_mpfr.hpp"
#endif

#ifdef HAVE_MPC
#include <mpc.h>
#include "cxx_mpc.hpp"
#endif

namespace cado_math_aux {
    /* the cado_math_aux mathematical constants are not constexpr because
     * we want to be able to instantiate them with variable precision
     * types. (they remain constexpr for pod types though)
     */
    namespace details {
        template<typename T>
            struct pi_impl {
                static constexpr T value(cado::number_context<T> = {}) { return std::numbers::pi_v<T>; }
            };
        template<typename T>
            struct pi_impl<std::complex<T>> {
                static constexpr std::complex<T> value(cado::number_context<std::complex<T>> = {}) { return std::numbers::pi_v<T>; }
            };
#ifdef HAVE_MPFR
        template<>
            struct pi_impl<cxx_mpfr> {
                static cxx_mpfr value(cado::number_context<cxx_mpfr> const & tr = {}) {
                    cxx_mpfr res = tr(0);
                    mpfr_const_pi(res, MPFR_RNDN);
                    return res;
                }
            };
#endif
#ifdef HAVE_MPC
        template<>
            struct pi_impl<cxx_mpc> {
                static cxx_mpc value(cado::number_context<cxx_mpc> const & tr = {}) {
                    cxx_mpc res = tr(0);
                    mpfr_const_pi(mpc_realref(res), MPFR_RNDN);
                    return res;
                }
            };
#endif
    } /* namespace details */
    template<typename T>
    inline T pi_v(cado::number_context<T> const & tr = {}) {
        return details::pi_impl<T>::value(tr);
    }
} /* namespace cado_math_aux */

#endif	/* UTILS_CADO_CONSTANTS_HPP_ */
