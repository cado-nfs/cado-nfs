#ifndef UTILS_EXTRA_COMPLEX_OVERLOADS_HPP_
#define UTILS_EXTRA_COMPLEX_OVERLOADS_HPP_

/*
 * The different overloads that are defined on std::complex cannot be
 * used for mixed integer/complex arithmetic. This provides some band
 * aids.
 *
 * https://en.cppreference.com/w/cpp/numeric/complex/operator_arith3.html
 *
 * In addition, we also provide comparison operators, which of course are
 * highly debatable. We do this only with the intent of sorting data.
 * Even then, it's absolutely conceivable that adding these tests achieve
 * nothing (behold the fact that these tests are not at all robust to
 * inifinitesimal alterations of the operands).
 */
#include <complex>
#include <compare>

template<typename T>
static inline std::partial_ordering operator<=>(std::complex<T> const & a, std::complex<T> const & b)
{
    std::partial_ordering const c = a.imag() <=> b.imag();
    if (c == std::partial_ordering::unordered) return c;
    std::partial_ordering const r = a.real() <=> b.real();
    if (c != 0) return c;
    return r;
}
template <typename T, typename U>
static inline std::partial_ordering operator<=>(std::complex<T> const & a, const U b)
    requires(std::is_integral_v<U> || std::is_floating_point_v<U>)
{
    std::partial_ordering const c = a.imag() <=> 0;
    if (c == std::partial_ordering::unordered) return c;
    std::partial_ordering const r = a.real() <=> b;
    if (c != 0) return c;
    return r;
}
template <typename T, typename U>
static inline bool operator==(std::complex<T> const & a, const U b)
    requires(std::is_integral_v<U> || std::is_floating_point_v<U>)
{
    std::partial_ordering const c = a.imag() <=> 0;
    if (c != 0) return false;
    std::partial_ordering const r = a.real() <=> b;
    return r == 0;
}

#if 0
constexpr int ordering_as_int(std::partial_ordering cmp) noexcept {
    return (cmp < 0) ? -1 : ((cmp == 0) ? 0 : 1);
}
template <typename T, typename U>
static inline int operator<=>(const U a, std::complex<T> const & b)
    requires(std::is_integral_v<U> || std::is_floating_point_v<U>)
{
    int const c = -ordering_as_int(b.imag() <=> 0);
    int const r = -ordering_as_int(b.real() <=> a);
    return (c << 1) + r;
}

/* NOLINTBEGIN(bugprone-macro-parentheses) */
#define CADO_STD_COMPLEX_DEFINE_CMP_SYM(OP)                             \
    template<typename T>                                                \
    inline bool operator OP(std::complex<T> const & a,                  \
                            std::complex<T> const & b)                  \
    {                                                                   \
        return ((a) <=> (b)) OP 0;                                      \
    }
#define CADO_STD_COMPLEX_DEFINE_CMP_ASYM(OP)                            \
    template <typename T, typename U>                                   \
    inline bool operator OP(std::complex<T> const & a, const U b)       \
    requires(std::is_integral_v<U> || std::is_floating_point_v<U>)      \
    {                                                                   \
        return ((a) <=> (b)) OP 0;                                      \
    }                                                                   \
    template <typename T, typename U>                                   \
    inline bool operator OP(const U a, std::complex<T> const & b)       \
    requires(std::is_integral_v<U> || std::is_floating_point_v<U>)      \
    {                                                                   \
        return 0 OP ((b) <=> (a));                                      \
    }
#define CADO_STD_COMPLEX_DEFINE_CMP(OP)                                 \
    CADO_STD_COMPLEX_DEFINE_CMP_SYM(OP)                                 \
    CADO_STD_COMPLEX_DEFINE_CMP_ASYM(OP)                                \
/* NOLINTEND(bugprone-macro-parentheses) */

CADO_STD_COMPLEX_DEFINE_CMP_ASYM(==)
CADO_STD_COMPLEX_DEFINE_CMP_ASYM(!=)
CADO_STD_COMPLEX_DEFINE_CMP(<)
CADO_STD_COMPLEX_DEFINE_CMP(>)
CADO_STD_COMPLEX_DEFINE_CMP(<=)
CADO_STD_COMPLEX_DEFINE_CMP(>=)
#endif

#define CADO_STD_COMPLEX_MIXED_OPERATOR_OVERLOAD(OP)			\
    template <typename T, typename U>					\
    static inline std::complex<T> OP(					\
            std::complex<T> const & a,					\
            const U b)							\
        requires std::is_integral_v<U>					\
    {									\
        return OP(a, T(b));						\
    }									\
    template <typename T, typename U>					\
    static inline std::complex<T> OP(				        \
            const U a,							\
            std::complex<T> const & b)					\
        requires std::is_integral_v<U>					\
    {									\
        return OP(T(a), b);						\
    }

CADO_STD_COMPLEX_MIXED_OPERATOR_OVERLOAD(operator+)
CADO_STD_COMPLEX_MIXED_OPERATOR_OVERLOAD(operator-)
CADO_STD_COMPLEX_MIXED_OPERATOR_OVERLOAD(operator*)
CADO_STD_COMPLEX_MIXED_OPERATOR_OVERLOAD(operator/)


#endif	/* UTILS_EXTRA_COMPLEX_OVERLOADS_HPP_ */
