#ifndef UTILS_CXX_MPFR_HPP_
#define UTILS_CXX_MPFR_HPP_

/* This is a thin c++ wrapper around mpfr. It is much less successful and
 * useful than the cxx_mpz wrapper, because of the conventions that are
 * in place with mpfr types: those imply that the rounding precision is
 * taken from the destination type, which of course makes sense, but
 * unfortunately it doesn't play well with the want we want to work with
 * c++ objects here.
 *
 * For this same reason, most overloads are killed on purpose
 *
 * So this cxx_mpfr is really _only_ an object holder.
 *
 * The rounding mode is MPFR_RNDN.
 *
 * The precision must be set by the caller
 */
#include <cstdint>
#include <cstdlib>

#include <limits>
#include <ostream>
#include <type_traits>

#include <mpfr.h>
#include "fmt/ostream.h"
#include "fmt/base.h"

#include "is_non_narrowing_conversion.hpp"
#include "macros.h"

struct cxx_mpfr {
public:
    mpfr_t x;
    // NOLINTBEGIN(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
    cxx_mpfr() { mpfr_init(x); }

    /* not even these are clearly correct
     *
    template <typename T, typename std::enable_if<
        std::is_integral<T>::value &&
        std::is_signed<T>::value &&
        cado_math_aux::is_non_narrowing_conversion<T, int64_t>::value,
        int>::type = 0 >
    // NOLINTNEXTLINE(hicpp-explicit-conversions,google-explicit-constructor)
    cxx_mpfr (const T & rhs) {
        mpfr_init2(x, std::numeric_limits<T>::digits);
        mpfr_set(x, int64_t(rhs), MPFR_RNDN);
    }
    template <typename T, typename std::enable_if<
        std::is_integral<T>::value &&
        std::is_signed<T>::value &&
        cado_math_aux::is_non_narrowing_conversion<T, int64_t>::value,
        int>::type = 0 >
    cxx_mpfr & operator=(const T a) {
        // set x to _exactly a
        mpfr_set_prec(x, std::numeric_limits<T>::digits);
        mpfr_set(x, int64_t(a), MPFR_RNDN);
        return *this;
    }
    template <typename T, typename std::enable_if<
        std::is_integral<T>::value &&
        !std::is_signed<T>::value &&
        cado_math_aux::is_non_narrowing_conversion<T, uint64_t>::value,
        int>::type = 0 >
    // NOLINTNEXTLINE(hicpp-explicit-conversions,google-explicit-constructor)
    cxx_mpfr (const T & rhs) {
        mpfr_init2(x, std::numeric_limits<T>::digits);
        mpfr_set(x, uint64_t(rhs), MPFR_RNDN);
    }
    template <typename T, typename std::enable_if<
        std::is_integral<T>::value &&
        !std::is_signed<T>::value &&
        cado_math_aux::is_non_narrowing_conversion<T, uint64_t>::value,
        int>::type = 0 >
    cxx_mpfr & operator=(const T a) {
        mpfr_set_prec(x, std::numeric_limits<T>::digits);
        mpfr_set(x, uint64_t(a), MPFR_RNDN);
        return *this;
    }
    */

    ~cxx_mpfr() { mpfr_clear(x); }
    // NOLINTNEXTLINE(hicpp-explicit-conversions,google-explicit-constructor)
    cxx_mpfr(mpfr_srcptr a) {
        mpfr_init2(x, mpfr_get_prec(a));
        mpfr_set(x, a, MPFR_RNDN);
    }
    cxx_mpfr(cxx_mpfr const & o) : cxx_mpfr(o.x)
    {}
    cxx_mpfr & operator=(cxx_mpfr const & o) {
        if (&o != this) {
            mpfr_set_prec(x, mpfr_get_prec(o.x));
            mpfr_set(x, o.x, MPFR_RNDN);
        }
        return *this;
    }
    // NOLINTEND(cppcoreguidelines-pro-type-member-init,hicpp-member-init)

#if __cplusplus >= 201103L
    cxx_mpfr(cxx_mpfr && o) noexcept
        : cxx_mpfr()
    {
        mpfr_swap(x, o.x);
    }
    cxx_mpfr& operator=(cxx_mpfr && o) noexcept {
        if (&o != this)
            mpfr_swap(x, o.x);
        return *this;
    }
#endif
    // NOLINTBEGIN(hicpp-explicit-conversions,google-explicit-constructor)
    operator mpfr_ptr() { return x; }
    operator mpfr_srcptr() const { return x; }
    /* it is very important to have the conversion to bool, otherwise the
     * implicit conversion to mpfr_ptr wins!
     */
    explicit operator bool() { return !mpfr_zero_p(x); }
    explicit operator bool() const { return !mpfr_zero_p(x); }
    // NOLINTEND(hicpp-explicit-conversions,google-explicit-constructor)
    mpfr_ptr operator->() { return x; }
    mpfr_srcptr operator->() const { return x; }

    /* tricky one. we'd like to live without it.
    explicit operator uint64_t() const {return mpfr_get_uint64(x);}
    */
};

#if GNUC_VERSION_ATLEAST(4,3,0)
extern void mpfr_init(cxx_mpfr & pl) __attribute__((error("mpfr_init must not be called on a mpfr reference -- it is the caller's business (via a ctor)")));
extern void mpfr_init2(cxx_mpfr & pl, mpfr_prec_t) __attribute__((error("mpfr_init must not be called on a mpfr reference -- it is the caller's business (via a ctor)")));
extern void mpfr_clear(cxx_mpfr & pl) __attribute__((error("mpfr_clear must not be called on a mpfr reference -- it is the caller's business (via a dtor)")));
#endif

inline std::ostream& operator<<(std::ostream& os, cxx_mpfr const& x) {
    char * s = nullptr;
    mpfr_asprintf(&s, "%Rf", x.x);
    os << s;
    free(s);
    return os;
}
/*
inline std::istream& operator>>(std::istream& is, cxx_mpfr & x) {
    std::string s;
    if (!(is >> s))
        return s;
    choke me;
    can use neither mpfr_set_str nor mpfr_strtofr here...
    return is >> (mpfr_ptr) x;
}
*/

namespace fmt {
    template <> struct formatter<cxx_mpfr>: ostream_formatter {};
}


#endif	/* UTILS_CXX_MPFR_HPP_ */
