#ifndef UTILS_NUMBER_CONTEXT_HPP_
#define UTILS_NUMBER_CONTEXT_HPP_

#include "cado_config.h"        // IWYU pragma: keep

#include <cctype>
#include <cstddef>

/* we have many different number types: int, long, double, float, but
 * also cxx_mpz or cxx_mpfr. In most cases, there is no ambiguity as to
 * what "0" means.
 *
 * cxx_mpfr is a different beast, though: we need some auxiliary
 * infomation (the precision) in order to properly interpret a given
 * value.
 */

#include <limits>
#include <string>
#include <type_traits>
#include <complex>

#include <gmp.h>
#ifdef HAVE_MPFR
#include <mpfr.h>
#endif
#ifdef HAVE_MPC
#include <mpc.h>
#endif

#include "number_literal.hpp"
#include "cado_math_aux.hpp"
#include "cxx_mpz.hpp"
#ifdef HAVE_MPFR
#include "cxx_mpfr.hpp"
#endif
#ifdef HAVE_MPC
#include "cxx_mpc.hpp"
#endif
#include "cado_parsing_base.hpp"

namespace cado {
    template<typename T> struct number_context;
    template<typename T> struct number_context {
        number_context() = default;
        template<typename U>
        explicit number_context(U const &) {}
        T operator()(std::string const&) const;
        T operator()(number_literal const & N) const
            requires std::is_integral_v<T>
        {
            if (N.has_point || N.has_exponent || N.is_imaginary)
                throw parse_error();
            return operator()(N.full);
        }
        T operator()(number_literal const & N) const
            requires(!std::is_integral_v<T>)
        {
            return operator()(N.full);
        }
        template<typename U>
        T operator()(U x) const
        requires (std::is_integral_v<U> || std::is_same_v<T, U>)
        { return x; }

        T operator()(mpz_srcptr x) const { return cado_math_aux::mpz_get<T>(x); }
        T operator()(cxx_mpz const & x) const { return cado_math_aux::mpz_get<T>(x); }
    };

    template<> struct number_context<cxx_mpz> {
        number_context() = default;
        template<typename U>
        explicit number_context(U const &) {}
        cxx_mpz operator()(std::string const & x) const {
            cxx_mpz z;
            mpz_set_str(z, x.c_str(), 0);
            return z;
        }
        cxx_mpz operator()(number_literal const & N) const {
            if (N.has_point || N.has_exponent)
                throw parse_error();
            return operator()(N.full);
        }
        template<typename U>
        cxx_mpz operator()(U x) const
        requires std::is_integral_v<U>
        { return x; }

        /* these are just as trivial as the above, but we need overloads
         * because they're not integral types */
        cxx_mpz operator()(mpz_srcptr x) const { return x; }
        cxx_mpz operator()(cxx_mpz const & x) const { return x; }
    };
#ifdef HAVE_MPFR
    template<> struct number_context<cxx_mpfr> {
        mpfr_prec_t prec = mpfr_get_default_prec();
        number_context() = default;
        explicit number_context(mpfr_prec_t prec)
            : prec(prec)
        {}

        template<typename U>
        explicit number_context(number_context<U> const &)
        requires (std::is_floating_point_v<U> || std::is_integral_v<U>)
        : prec(std::numeric_limits<U>::digits)
        {}

        explicit number_context(cxx_mpfr const & x)
            : prec(x.prec())
        {}
        cxx_mpfr operator()(std::string const & s) const {
            cxx_mpfr res;
            mpfr_set_prec(res, prec);
            const int r = mpfr_set_str(res, s.c_str(), 0, MPFR_RNDN);
            if (r != 0)
                throw parse_error();
            if (mpfr_nan_p(res.x) && s.starts_with('-'))
                mpfr_setsign(res, res, 1, MPFR_RNDN);

            return res;
        }
        cxx_mpfr operator()(number_literal const & N) const {
            return operator()(N.full);
        }
        template<typename U>
        cxx_mpfr operator()(U x) const
        requires std::is_integral_v<U>
        {
            cxx_mpfr y;
            mpfr_set_prec(y, prec);
            mpfr_auxx::cado_mpfr_set(y, x, MPFR_RNDN);
            return y;
        }
        cxx_mpfr operator()(cxx_mpfr const & x) const
        {
            cxx_mpfr y;
            mpfr_set_prec(y, prec);
            mpfr_auxx::cado_mpfr_set(y, x, MPFR_RNDN);
            return y;
        }
        cxx_mpfr operator()(mpz_srcptr x) const
        {
            cxx_mpfr y;
            mpfr_set_prec(y, prec);
            mpfr_auxx::cado_mpfr_set(y, x, MPFR_RNDN);
            return y;
        }
        cxx_mpfr operator()(cxx_mpz const & x) const
        {
            cxx_mpfr y;
            mpfr_set_prec(y, prec);
            mpfr_auxx::cado_mpfr_set(y, x, MPFR_RNDN);
            return y;
        }
    };
#endif

    template<typename T> struct number_context<std::complex<T>> {
        number_context() = default;
        template<typename U>
        explicit number_context(U const &) {}
        std::complex<T> operator()(std::string const & s) const;
        std::complex<T> operator()(number_literal const & N) const
        {
            T h(number_context<T>()(N.strip_i()));
            if (N.is_imaginary)
                return { 0, h };
            else
                return { h, 0 };
        }
        template<typename U>
        std::complex<T> operator()(U x) const
        requires (std::is_integral_v<U> || std::is_same_v<std::complex<T>, U>)
        { return x; }

        std::complex<T> operator()(mpz_srcptr x) const { return cado_math_aux::mpz_get<std::complex<T>>(x); }
        std::complex<T> operator()(cxx_mpz const & x) const { return cado_math_aux::mpz_get<std::complex<T>>(x); }
    };

#ifdef HAVE_MPC
    /* it's essentially the same thing. If we want to extend, it might
     * make sense to specify whether we want L2 or Linf approximations,
     * or if we want skewed precision for real and imaginary part, but
     * that's about it.
     */
    template<> struct number_context<cxx_mpc> : public number_context<cxx_mpfr> {
        using number_context<cxx_mpfr>::number_context;

        explicit number_context(cxx_mpc const & x)
            : number_context<cxx_mpfr>(x.prec())
        {}
        template<typename U>
        cxx_mpc operator()(U const & u) const {
            return cxx_mpc(number_context<cxx_mpfr>::operator()(u));
        }
        cxx_mpc operator()(std::string const & s) const;
        cxx_mpc operator()(number_literal const & N) const {
            cxx_mpc res;
            mpc_set_prec(res, prec);
            cxx_mpfr h(number_context<cxx_mpfr>::operator()(N.strip_i()));
            if (N.is_imaginary) {
                mpfr_swap(mpc_imagref(res), h);
                mpfr_set_ui(mpc_realref(res), 0, MPFR_RNDN);
            } else {
                mpfr_swap(mpc_realref(res), h);
                mpfr_set_ui(mpc_imagref(res), 0, MPFR_RNDN);
            }
            return res;
        }
        cxx_mpc operator()(cxx_mpc const & u) const {
            cxx_mpc res;
            mpc_set_prec(res, prec);
            mpc_set(res, u, MPC_RNDNN);
            return res;
        }
    };
#endif

#define CADO_NUMBER_CONTEXT_CALL_STRING(T, code)			\
    template<> inline T                                                 \
    cado::number_context<T>::operator()(std::string const & s) const	\
    {									\
        size_t pos;							\
        auto ret = code;						\
        if (pos != s.size())						\
            throw parse_error();					\
        return ret;							\
    }

#define CADO_NUMBER_CONTEXT_DEFINE_EXPLICIT_INSTANTIATION(T)            \
    template struct number_context<T>;

#define CADO_NUMBER_CONTEXT_DO(T, code)                                 \
    CADO_NUMBER_CONTEXT_CALL_STRING(T, code)                            \
    CADO_NUMBER_CONTEXT_DEFINE_EXPLICIT_INSTANTIATION(T)

#define CADO_NUMBER_CONTEXT_PROXY(T, U)                                 \
    CADO_NUMBER_CONTEXT_CALL_STRING(T, number_context<U>()(s))          \
    CADO_NUMBER_CONTEXT_DEFINE_EXPLICIT_INSTANTIATION(T)


CADO_NUMBER_CONTEXT_DO(int, std::stoi(s, &pos))
CADO_NUMBER_CONTEXT_DO(long, std::stol(s, &pos))
CADO_NUMBER_CONTEXT_DO(unsigned int, std::stoul(s, &pos))
CADO_NUMBER_CONTEXT_DO(unsigned long, std::stoul(s, &pos))

#ifndef LONG_LONG_IS_EXACTLY_LONG
CADO_NUMBER_CONTEXT_DO(long long, std::stoll(s, &pos))
#endif
#ifndef UNSIGNED_LONG_LONG_IS_EXACTLY_UNSIGNED_LONG
CADO_NUMBER_CONTEXT_DO(unsigned long long, std::stoull(s, &pos))
#endif

#ifdef UINT64_T_IS_EXACTLY_UNSIGNED_LONG
#elif defined(UINT64_T_IS_EXACTLY_UNSIGNED_LONG_LONG)
#elif defined(UINT64_T_IS_EXACTLY_UNSIGNED)
#elif defined(UINT64_T_IS_COMPATIBLE_WITH_UNSIGNED_LONG)
CADO_NUMBER_CONTEXT_PROXY(uint64_t, unsigned long)
#elif defined(UINT64_T_IS_COMPATIBLE_WITH_UNSIGNED_LONG_LONG)
CADO_NUMBER_CONTEXT_PROXY(uint64_t, unsigned long long)
#elif defined(UINT64_T_IS_COMPATIBLE_WITH_UNSIGNED)
CADO_NUMBER_CONTEXT_PROXY(uint64_t, unsigned)
#else
#error "We need to know how to create uint64_t from basic types"
#endif

#ifdef INT64_T_IS_EXACTLY_LONG
#elif defined(INT64_T_IS_EXACTLY_LONG_LONG)
#elif defined(INT64_T_IS_EXACTLY_INT
#elif defined(INT64_T_IS_COMPATIBLE_WITH_LONG)
CADO_NUMBER_CONTEXT_PROXY(int64_t, long)
#elif defined(INT64_T_IS_COMPATIBLE_WITH_LONG_LONG)
CADO_NUMBER_CONTEXT_PROXY(int64_t, long long)
#elif defined(INT64_T_IS_COMPATIBLE_WITH_INT
CADO_NUMBER_CONTEXT_PROXY(int64_t, int)
#else
#error "We need to know how to create int64_t from basic types"
#endif

#ifdef UINT32_T_IS_EXACTLY_UNSIGNED_LONG
#elif defined(UINT32_T_IS_EXACTLY_UNSIGNED)
#elif defined(UINT32_T_IS_COMPATIBLE_WITH_UNSIGNED_LONG)
CADO_NUMBER_CONTEXT_PROXY(uint32_t, unsigned long)
#elif defined(UINT32_T_IS_COMPATIBLE_WITH_UNSIGNED)
CADO_NUMBER_CONTEXT_PROXY(uint32_t, unsigned)
#else
#error "We need to know how to create uint32_t from basic types"
#endif

#ifdef INT32_T_IS_EXACTLY_LONG
#elif defined(INT32_T_IS_EXACTLY_INT)
#elif defined(INT32_T_IS_COMPATIBLE_WITH_LONG)
CADO_NUMBER_CONTEXT_PROXY(int32_t, long)
#elif defined(INT32_T_IS_COMPATIBLE_WITH_INT)
CADO_NUMBER_CONTEXT_PROXY(int32_t, int)
#else
#error "We need to know how to create int32_t from basic types"
#endif

CADO_NUMBER_CONTEXT_DO(float, std::stof(s, &pos))
CADO_NUMBER_CONTEXT_DO(double, std::stod(s, &pos))
CADO_NUMBER_CONTEXT_DO(long double, std::stold(s, &pos))

/* we don't have a specialization of operator(), so we mustn't
 * declare one. What we have is an implicit instantiation, which happens
 * in the cpp file.
 */
extern template struct number_context<std::complex<float>>;
extern template struct number_context<std::complex<double>>;
extern template struct number_context<std::complex<long double>>;

} /* namespace cado */


#endif	/* UTILS_NUMBER_CONTEXT_HPP_ */
