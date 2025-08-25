#ifndef UTILS_NUMBER_CONTEXT_HPP_
#define UTILS_NUMBER_CONTEXT_HPP_

#include "cado_config.h"        // IWYU pragma: keep

/* we have many different number types: int, long, double, float, but
 * also cxx_mpz or cxx_mpfr. In most cases, there is no ambiguity as to
 * what "0" means.
 *
 * cxx_mpfr is a different beast, though: we need some auxiliary
 * infomation (the precision) in order to properly interpret a given
 * value.
 */

#include <string>
#include <type_traits>

#include <gmp.h>

#include "number_literal.hpp"
#include "cxx_mpz.hpp"
#ifdef HAVE_MPFR
#include "cxx_mpfr.hpp"
#endif

namespace cado {
    template<typename T> struct number_context {
        number_context() = default;
        template<typename U>
        explicit number_context(U const &) {}
        T operator()(std::string const&) const;
        T operator()(number_literal const & N) const
            requires std::is_integral_v<T>
        {
            if (N.has_point || N.has_exponent)
                throw number_literal::parse_error();
            return operator()(N.full);
        }
        T operator()(number_literal const & N) const
            requires(!std::is_integral_v<T>)
        {
            return operator()(N.full);
        }
        template<typename U>
        T interpret_integral(U x) const { return T(x); }
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
                throw number_literal::parse_error();
            return operator()(N.full);
        }
        template<typename U>
        cxx_mpz interpret_integral(U x) const { return cxx_mpz(x); }
    };
#ifdef HAVE_MPFR
    template<> struct number_context<cxx_mpfr> {
        mpfr_prec_t prec = mpfr_get_default_prec();
        number_context() = default;
        explicit number_context(mpfr_prec_t prec)
            : prec(prec)
        {}
        explicit number_context(cxx_mpfr const & x)
            : prec(x.prec())
        {}
        cxx_mpfr operator()(std::string const & s) const {
            cxx_mpfr res;
            mpfr_set_prec(res, prec);
            const int r = mpfr_set_str(res, s.c_str(), 0, MPFR_RNDN);
            if (r != 0)
                throw number_literal::parse_error();
            if (mpfr_nan_p(res.x) && s.starts_with('-'))
                mpfr_setsign(res, res, 1, MPFR_RNDN);

            return res;
        }
        cxx_mpfr operator()(number_literal const & N) const {
            return operator()(N.full);
        }
        template<typename U>
        cxx_mpfr interpret_integral(U x) const {
            cxx_mpfr y;
            mpfr_set_prec(y, prec);
            mpfr_auxx::cado_mpfr_set(y, x, MPFR_RNDN);
            return y;
        }
    };
#endif


    template<>
    inline int number_context<int>::operator()(std::string const & s) const
    {
        size_t pos;
        auto ret = std::stoi(s, &pos);
        if (pos != s.size())
            throw number_literal::parse_error();
        return ret;
    }
    template<>
    inline unsigned int number_context<unsigned int>::operator()(std::string const & s) const
    {
        /* stoul will do */
        size_t pos;
        auto ret = std::stoul(s, &pos);
        if (pos != s.size())
            throw number_literal::parse_error();
        return ret;
    }
    template<>
    inline long number_context<long>::operator()(std::string const & s) const
    {
        size_t pos;
        auto ret = std::stol(s, &pos);
        if (pos != s.size())
            throw number_literal::parse_error();
        return ret;
    }
    template<>
    inline unsigned long number_context<unsigned long>::operator()(std::string const & s) const
    {
        size_t pos;
        auto ret = std::stoul(s, &pos);
        if (pos != s.size())
            throw number_literal::parse_error();
        return ret;
    }
#ifndef LONG_LONG_IS_EXACTLY_LONG
    template<>
    inline long long number_context<long long>::operator()(std::string const & s) const
    {
        size_t pos;
        auto ret = std::stoll(s, &pos);
        if (pos != s.size())
            throw number_literal::parse_error();
        return ret;
    }
#endif
#ifndef UNSIGNED_LONG_LONG_IS_EXACTLY_UNSIGNED_LONG
    template<>
    inline unsigned long long number_context<unsigned long long>::operator()(std::string const & s) const
    {
        size_t pos;
        auto ret = std::stoull(s, &pos);
        if (pos != s.size())
            throw number_literal::parse_error();
        return ret;
    }
#endif

#ifdef UINT64_T_IS_EXACTLY_UNSIGNED_LONG
#elif defined(UINT64_T_IS_EXACTLY_UNSIGNED_LONG_LONG)
#elif defined(UINT64_T_IS_EXACTLY_UNSIGNED)
#elif defined(UINT64_T_IS_COMPATIBLE_WITH_UNSIGNED_LONG)
    template<>
    inline uint64_t number_context<uint64_t>::operator()(std::string const & s) const
    {
        return number_context<unsigned long>()(s);
    }
#elif defined(UINT64_T_IS_COMPATIBLE_WITH_UNSIGNED_LONG_LONG)
    template<>
    inline uint64_t number_context<uint64_t>::operator()(std::string const & s) const
    {
        return number_context<unsigned long long>()(s);
    }
#elif defined(UINT64_T_IS_COMPATIBLE_WITH_UNSIGNED)
    template<>
    inline uint64_t number_context<uint64_t>::operator()(std::string const & s) const
    {
        return number_context<unsigned int>()(s);
    }
#else
#error "We need to know how to create uint64_t from basic types"
#endif

#ifdef INT64_T_IS_EXACTLY_LONG
#elif defined(INT64_T_IS_EXACTLY_LONG_LONG)
#elif defined(INT64_T_IS_EXACTLY_INT
#elif defined(INT64_T_IS_COMPATIBLE_WITH_LONG)
    template<>
    inline int64_t number_context<int64_t>::operator()(std::string const & s) const
    {
        return number_context<long>()(s);
    }
#elif defined(INT64_T_IS_COMPATIBLE_WITH_LONG_LONG)
    template<>
    inline int64_t number_context<int64_t>::operator()(std::string const & s) const
    {
        return number_context<long long>()(s);
    }
#elif defined(INT64_T_IS_COMPATIBLE_WITH_INT)
    template<>
    inline int64_t number_context<int64_t>::operator()(std::string const & s) const
    {
        return number_context<int>()(s);
    }
#else
#error "We need to know how to create int64_t from basic types"
#endif

#ifdef UINT32_T_IS_EXACTLY_UNSIGNED_LONG
#elif defined(UINT32_T_IS_EXACTLY_UNSIGNED)
#elif defined(UINT32_T_IS_COMPATIBLE_WITH_UNSIGNED_LONG)
    template<>
    inline uint32_t number_context<uint32_t>::operator()(std::string const & s) const
    {
        return number_context<unsigned long>()(s);
    }
#elif defined(UINT32_T_IS_COMPATIBLE_WITH_UNSIGNED)
    template<>
    inline uint32_t number_context<uint32_t>::operator()(std::string const & s) const
    {
        return number_context<unsigned int>()(s);
    }
#else
#error "We need to know how to create uint32_t from basic types"
#endif

#ifdef INT32_T_IS_EXACTLY_LONG
#elif defined(INT32_T_IS_EXACTLY_INT)
#elif defined(INT32_T_IS_COMPATIBLE_WITH_LONG)
    template<>
    inline int32_t number_context<int32_t>::operator()(std::string const & s) const
    {
        return number_context<long>()(s);
    }
#elif defined(INT32_T_IS_COMPATIBLE_WITH_INT)
    template<>
    inline int32_t number_context<int32_t>::operator()(std::string const & s) const
    {
        return number_context<int>()(s);
    }
#else
#error "We need to know how to create int32_t from basic types"
#endif

    template<>
    inline float number_context<float>::operator()(std::string const & s) const
    {
        size_t pos;
        auto ret = std::stof(s, &pos);
        if (pos != s.size())
            throw number_literal::parse_error();
        return ret;
    }
    template<>
    inline double number_context<double>::operator()(std::string const & s) const
    {
        size_t pos;
        auto ret = std::stod(s, &pos);
        if (pos != s.size())
            throw number_literal::parse_error();
        return ret;
    }
    template<>
    inline long double number_context<long double>::operator()(std::string const & s) const
    {
        size_t pos;
        auto ret = std::stold(s, &pos);
        if (pos != s.size())
            throw number_literal::parse_error();
        return ret;
    }
} /* namespace cado */


#endif	/* UTILS_NUMBER_CONTEXT_HPP_ */
