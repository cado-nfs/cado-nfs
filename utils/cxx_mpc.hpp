#ifndef UTILS_CXX_MPC_HPP_
#define UTILS_CXX_MPC_HPP_

#include <cstdint>
#include <cstdlib>

#include <limits>
#include <ostream>
#include <type_traits>

#include "fmt/base.h"
#include "fmt/ostream.h"
#include <mpc.h>

#include "is_non_narrowing_conversion.hpp"
#include "macros.h"

#include "mpc_auxx.hpp"

struct cxx_mpc {
  public:
    mpc_t x;
    // NOLINTBEGIN(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
    cxx_mpc() { mpc_init2(x, mpfr_get_default_prec()); }

    /* calling the mpc macros on cxx_mpc objects leads to diagnostics
     * errors with clang unless we first cast them to mpc_srcptr, which
     * a bit burdensome.
     */
    mpfr_prec_t prec() const { return mpc_get_prec(x); }

    template <
        typename T,
        typename std::enable_if<
            std::is_integral<T>::value && std::is_signed<T>::value &&
                cado_math_aux::is_non_narrowing_conversion<T, int64_t>::value,
            int>::type = 0>
    // NOLINTNEXTLINE(hicpp-explicit-conversions,google-explicit-constructor)
    cxx_mpc(T const & rhs)
    {
        mpc_init2(x, std::numeric_limits<T>::digits);
        mpc_auxx::cado_mpc_set(x, int64_t(rhs), MPFR_RNDN);
    }
    template <
        typename T,
        typename std::enable_if<
            std::is_integral<T>::value && std::is_signed<T>::value &&
                cado_math_aux::is_non_narrowing_conversion<T, int64_t>::value,
            int>::type = 0>
    cxx_mpc & operator=(T const a)
    {
        mpc_set_prec(x, std::numeric_limits<T>::digits);
        mpc_auxx::cado_mpc_set(x, int64_t(a), MPFR_RNDN);
        return *this;
    }
    template <
        typename T,
        typename std::enable_if<
            std::is_integral<T>::value && !std::is_signed<T>::value &&
                cado_math_aux::is_non_narrowing_conversion<T, uint64_t>::value,
            int>::type = 0>
    // NOLINTNEXTLINE(hicpp-explicit-conversions,google-explicit-constructor)
    cxx_mpc(T const & rhs)
    {
        mpc_init2(x, std::numeric_limits<T>::digits);
        mpc_auxx::cado_mpc_set(x, uint64_t(rhs), MPFR_RNDN);
    }
    template <
        typename T,
        typename std::enable_if<
            std::is_integral<T>::value && !std::is_signed<T>::value &&
                cado_math_aux::is_non_narrowing_conversion<T, uint64_t>::value,
            int>::type = 0>
    cxx_mpc & operator=(T const a)
    {
        mpc_set_prec(x, std::numeric_limits<T>::digits);
        mpc_auxx::cado_mpc_set(x, uint64_t(a), MPFR_RNDN);
        return *this;
    }

    ~cxx_mpc() { mpc_clear(x); }
    // NOLINTNEXTLINE(hicpp-explicit-conversions,google-explicit-constructor)
    cxx_mpc(mpc_srcptr a)
    {
        mpc_init2(x, mpc_get_prec(a));
        mpc_set(x, a, MPFR_RNDN);
    }
    cxx_mpc(cxx_mpc const & o)
        : cxx_mpc(o.x)
    {
    }
    cxx_mpc & operator=(cxx_mpc const & o)
    {
        if (&o != this) {
            mpc_set_prec(x, mpc_get_prec(o.x));
            mpc_set(x, o.x, MPFR_RNDN);
        }
        return *this;
    }
    // NOLINTEND(cppcoreguidelines-pro-type-member-init,hicpp-member-init)

#if __cplusplus >= 201103L
    cxx_mpc(cxx_mpc && o) noexcept
        : cxx_mpc()
    {
        mpc_swap(x, o.x);
    }
    cxx_mpc & operator=(cxx_mpc && o) noexcept
    {
        if (&o != this)
            mpc_swap(x, o.x);
        return *this;
    }
#endif
    // NOLINTBEGIN(hicpp-explicit-conversions,google-explicit-constructor)
    operator mpc_ptr() { return x; }
    operator mpc_srcptr() const { return x; }
    /* it is very important to have the conversion to bool, otherwise the
     * implicit conversion to mpc_ptr wins!
     */
    explicit operator bool() { return mpc_cmp_si_si(x, 0, 0) != 0; }
    explicit operator bool() const { return mpc_cmp_si_si(x, 0, 0) != 0; }
    // NOLINTEND(hicpp-explicit-conversions,google-explicit-constructor)
    mpc_ptr operator->() { return x; }
    mpc_srcptr operator->() const { return x; }

    /* tricky one. we'd like to live without it.
    explicit operator uint64_t() const {return mpc_get_uint64(x);}
    */
};

#if GNUC_VERSION_ATLEAST(4, 3, 0)
extern void mpc_init(cxx_mpc & pl)
    __attribute__((error("mpc_init must not be called on a mpc reference -- "
                         "it is the caller's business (via a ctor)")));
extern void mpc_init2(cxx_mpc & pl, mpfr_prec_t)
    __attribute__((error("mpc_init must not be called on a mpc reference -- "
                         "it is the caller's business (via a ctor)")));
extern void mpc_clear(cxx_mpc & pl)
    __attribute__((error("mpc_clear must not be called on a mpc reference -- "
                         "it is the caller's business (via a dtor)")));
#endif

inline std::ostream & operator<<(std::ostream & os, cxx_mpc const & x)
{
    char * s = nullptr;
    mpfr_asprintf(&s, "%Rf+i*%Rf", mpc_realref(x.x), mpc_imagref(x.x));
    os << s;
    free(s);
    return os;
}
/*
inline std::istream& operator>>(std::istream& is, cxx_mpc & x) {
    std::string s;
    if (!(is >> s))
        return s;
    choke me;
    can use neither mpc_set_str nor mpc_strtofr here...
    return is >> (mpc_ptr) x;
}
*/

namespace fmt
{
template <> struct formatter<cxx_mpc> : ostream_formatter {
};
} // namespace fmt

/* Now here's a layer we're not particularly happy with */

/* NOLINTBEGIN(bugprone-macro-parentheses) */
#define CXX_MPC_DEFINE_CMP(OP)                                                \
    inline bool operator OP(cxx_mpc const & a, cxx_mpc const & b)            \
    {                                                                          \
        return mpc_cmp(a, b) OP 0;                                            \
    }                                                                          \
    inline bool operator OP(mpc_srcptr a, cxx_mpc const & b)                 \
    {                                                                          \
        return mpc_cmp(a, b) OP 0;                                            \
    }                                                                          \
    inline bool operator OP(cxx_mpc const & a, mpc_srcptr b)                 \
    {                                                                          \
        return mpc_cmp(a, b) OP 0;                                            \
    }                                                                          \
    template <typename T,                                                      \
              std::enable_if_t<std::is_integral<T>::value, int> = 0>           \
    inline bool operator OP(cxx_mpc const & a, const T b)                     \
    {                                                                          \
        return mpc_auxx::cado_mpc_cmp(a, b) OP 0;                            \
    }                                                                          \
    template <typename T,                                                      \
              std::enable_if_t<std::is_integral<T>::value, int> = 0>           \
    inline bool operator OP(const T a, cxx_mpc const & b)                     \
    {                                                                          \
        return 0 OP mpc_auxx::cado_mpc_cmp(b, a);                            \
    }
/* NOLINTEND(bugprone-macro-parentheses) */

CXX_MPC_DEFINE_CMP(==)
CXX_MPC_DEFINE_CMP(!=)
CXX_MPC_DEFINE_CMP(<)
CXX_MPC_DEFINE_CMP(>)
CXX_MPC_DEFINE_CMP(<=)
CXX_MPC_DEFINE_CMP(>=)

#if __cplusplus >= 202002L
inline bool operator<=>(cxx_mpc const & a, cxx_mpc const & b)
{
    return mpc_auxx::cado_mpc_cmp(b, a);
}
#endif

#define CXX_MPC_DEFINE_TERNARY(OP, TEXTOP)                                    \
    inline cxx_mpc operator OP(cxx_mpc const & a, cxx_mpc const & b)        \
    {                                                                          \
        cxx_mpc r;                                                            \
        mpc_set_prec(r, mpc_get_prec(mpc_srcptr(a)));                       \
        mpc_##TEXTOP(r, a, b, MPFR_RNDN);                                     \
        return r;                                                              \
    }                                                                          \
    template <typename T,                                                      \
              std::enable_if_t<std::is_integral<T>::value, int> = 0>           \
    inline cxx_mpc operator OP(cxx_mpc const & a, T const b)                 \
    {                                                                          \
        cxx_mpc r;                                                            \
        mpc_set_prec(r, mpc_get_prec(mpc_srcptr(a)));                       \
        mpc_auxx::cado_mpc_##TEXTOP(r, a, b, MPFR_RNDN);                     \
        return r;                                                              \
    }                                                                          \
    template <typename T,                                                      \
              std::enable_if_t<std::is_integral<T>::value, int> = 0>           \
    inline cxx_mpc operator OP(T const a, cxx_mpc const & b)                 \
    {                                                                          \
        cxx_mpc r;                                                            \
        mpc_set_prec(r, mpc_get_prec(mpc_srcptr(b)));                       \
        mpc_auxx::cado_mpc_##TEXTOP(r, b, a, MPFR_RNDN);                     \
        return r;                                                              \
    }                                                                          \
    inline cxx_mpc & operator OP##=(cxx_mpc & a, cxx_mpc const & b)	\
    {									\
        mpc_##TEXTOP(a, a, b, MPFR_RNDN);				\
        return a;							\
    }									\
    template <typename T, std::enable_if_t<std::is_integral<T>::value, int> = 0>        \
    inline cxx_mpc & operator OP##=(cxx_mpc & a, T const b)		\
    {									\
        mpc_auxx::cado_mpc_##TEXTOP(a, a, b, MPFR_RNDN);		\
        return a;							\
    }

CXX_MPC_DEFINE_TERNARY(+, add)
CXX_MPC_DEFINE_TERNARY(-, sub)
CXX_MPC_DEFINE_TERNARY(*, mul)
CXX_MPC_DEFINE_TERNARY(/, div)
CXX_MPC_DEFINE_TERNARY(%, remainder)

inline cxx_mpc operator-(cxx_mpc const & a)
{
    cxx_mpc r;
    mpfr_set_prec(mpc_realref(r), mpfr_get_prec(mpfr_srcptr(mpc_realref(a))));
    mpfr_set_prec(mpc_imagref(r), mpfr_get_prec(mpfr_srcptr(mpc_realref(a))));
    mpc_neg(r, a, MPFR_RNDN);
    return r;
}

inline cxx_mpc & operator<<=(cxx_mpc & a, unsigned long const s)
{
    mpfr_mul_2exp(mpc_realref(a), mpc_realref(a), s, MPFR_RNDN);
    mpfr_mul_2exp(mpc_imagref(a), mpc_imagref(a), s, MPFR_RNDN);
    return a;
}
inline cxx_mpc operator<<(cxx_mpc const & a, unsigned long const s)
{
    cxx_mpc r {a};
    return r <<= s;
}

inline cxx_mpc & operator>>=(cxx_mpc & a, unsigned long const s)
{
    mpfr_div_2exp(mpc_realref(a), mpc_realref(a), s, MPFR_RNDN);
    mpfr_div_2exp(mpc_imagref(a), mpc_imagref(a), s, MPFR_RNDN);
    return a;
}
inline cxx_mpc operator>>(cxx_mpc const & a, unsigned long const s)
{
    cxx_mpc r {a};
    return r >>= s;
}

#endif	/* UTILS_CXX_MPC_HPP_ */
