#ifndef UTILS_CXX_MPFR_HPP_
#define UTILS_CXX_MPFR_HPP_

/* This is a thin c++ wrapper around mpfr.
 *
 * It is enumbered by the conventions in place with mpfr types: those
 * imply that the rounding precision is taken from the destination type,
 * which of course makes sense, but unfortunately doesn't play well
 * with the way we want to work with c++ objects here.
 *
 * So we're on shaky grounds. As far as using cxx_mpfr's is concerned, we
 * will assume that if we insist on using overloads, it means that:
 *  - rounding is always MPF_RNDN
 *  - the precision is the preision of the first *INPUT* (or
 *  input-output) operand of cxx_mpfr type.
 *
 * Finer grain usage is still possible with the mpfr_ functions. Note
 * however that many mpfr functions are macros, which leads to all sorts
 * of warnings.
 */
#include <cstdint>
#include <cstdlib>

#include <limits>
#include <istream>
#include <ostream>
#include <type_traits>

#include "fmt/base.h"
#include <mpfr.h>

#include "utils_cxx.hpp"
#include "macros.h"
#include "mpfr_auxx.hpp"

struct cxx_mpfr {
  public:
    mpfr_t x;
    // NOLINTBEGIN(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
    cxx_mpfr() { mpfr_init(x); }

    /* calling the mpfr macros on cxx_mpfr objects leads to diagnostics
     * errors with clang unless we first cast them to mpfr_srcptr, which
     * a bit burdensome.
     */
    mpfr_prec_t prec() const { return mpfr_get_prec(x); }
    int sgn() const { return mpfr_sgn(x); }

    /* mpfr code determines the working precision from the _target_ type.
     * This doesn't work in a context where we want the result to be
     * returned.
     *
     * In most cases, we do away with this by taking the precision of the
     * first cxx_mpfr argument as defining the precision of the result.
     * If this is not what is desired, it's always possible to use
     * combinations of mpfr_set_prec and the classical multi-operand
     * functions of mpfr, like mpfr_add or mpfr_mul
     *
     * Construction from immediate integers, as we have below, is of
     * course annoying because there's no precision to speak of, beyond
     * the precision that is inherent to T, or possibly the global,
     * default precision.
     *
     * There's also the "originating precision" in the case of operator=,
     * meaning the precision of the object being assigned to.  We rule
     * that out because that is just not how c++ works, and would cause
     * all sorts of trouble: Writing "a = x;" should not depend on the
     * prior value of a, period.
     *
     * A consequence of the latter is that whatever can be done for
     * ctor's can also be done for assignments, and probably should.
     *
     * We can do two things.
     *  - stick to what mpfr itself does.  All the mpfr_init_set
     *    functions use the global default precision, so we can use that.
     *  - set the precision to the number of bits that are necessary to
     *    hold the source type.
     *
     * The former introduces a dependency on a global variable.which is
     * bad, and also leaves us with no convenient way to set a cxx_mpfr to the
     * _exact_ value of an integer! The only thing we can do is
     * mpfr_set_prec followed by mpfr_set_ui.
     *
     * The apparent "consistency" bonus of the first approach is, in
     * fact, limited, since the mpfr_init_* functions cannot be called on
     * cxx_mpfr types.
     *
     * It seems that in both cases, we're going to need mpfr_prec_round
     * quite often anyway. For this reason, the second approach, which at
     * least gives the opportunity to represent integers exactly, is
     * preferred.
     */

    template <typename T>
    // NOLINTNEXTLINE(hicpp-explicit-conversions,google-explicit-constructor)
    cxx_mpfr(T const & rhs)
    requires cado::converts_via<T, int64_t>
    {
        mpfr_init2(x, std::numeric_limits<T>::digits);
        mpfr_auxx::cado_mpfr_set(x, int64_t(rhs), MPFR_RNDN);
    }
    template <typename T>
    // NOLINTNEXTLINE(hicpp-explicit-conversions,google-explicit-constructor)
    cxx_mpfr(T const & rhs)
    requires cado::converts_via<T, uint64_t>
    {
        mpfr_init2(x, std::numeric_limits<T>::digits);
        mpfr_auxx::cado_mpfr_set(x, uint64_t(rhs), MPFR_RNDN);
    }

    template <typename T>
    cxx_mpfr & operator=(T const a)
    requires cado::converts_via<T, int64_t>
    {
        mpfr_set_prec(x, std::numeric_limits<T>::digits);
        mpfr_auxx::cado_mpfr_set(x, int64_t(a), MPFR_RNDN);
        return *this;
    }
    template <typename T>
    cxx_mpfr & operator=(T const a)
    requires cado::converts_via<T, uint64_t>
    {
        mpfr_set_prec(x, std::numeric_limits<T>::digits);
        mpfr_auxx::cado_mpfr_set(x, uint64_t(a), MPFR_RNDN);
        return *this;
    }

    explicit cxx_mpfr(double rhs)
    {
        mpfr_init2(x, std::numeric_limits<decltype(rhs)>::digits);
        mpfr_set_d(x, rhs, MPFR_RNDN);
    }

    explicit cxx_mpfr(long double rhs)
    {
        mpfr_init2(x, std::numeric_limits<decltype(rhs)>::digits);
        mpfr_set_ld(x, rhs, MPFR_RNDN);
    }

    cxx_mpfr & operator=(double a)
    {
        mpfr_set_prec(x, std::numeric_limits<decltype(a)>::digits);
        mpfr_set_d(x, a, MPFR_RNDN);
        return *this;
    }
    cxx_mpfr & operator=(long double a)
    {
        mpfr_set_prec(x, std::numeric_limits<decltype(a)>::digits);
        mpfr_set_ld(x, a, MPFR_RNDN);
        return *this;
    }

    ~cxx_mpfr() { mpfr_clear(x); }
    // NOLINTNEXTLINE(hicpp-explicit-conversions,google-explicit-constructor)
    cxx_mpfr(mpfr_srcptr a)
    {
        mpfr_init2(x, mpfr_get_prec(a));
        mpfr_set(x, a, MPFR_RNDN);
    }
    cxx_mpfr(cxx_mpfr const & o)
        : cxx_mpfr(o.x)
    {
    }
    cxx_mpfr & operator=(cxx_mpfr const & o)
    {
        if (&o != this) {
            mpfr_set_prec(x, mpfr_get_prec(o.x));
            mpfr_set(x, o.x, MPFR_RNDN);
        }
        return *this;
    }
    // NOLINTEND(cppcoreguidelines-pro-type-member-init,hicpp-member-init)

    cxx_mpfr(cxx_mpfr && o) noexcept
        : cxx_mpfr()
    {
        mpfr_swap(x, o.x);
    }
    cxx_mpfr & operator=(cxx_mpfr && o) noexcept
    {
        if (&o != this)
            mpfr_swap(x, o.x);
        return *this;
    }
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
    struct input_with_precision {
        cxx_mpfr & x;
        mpfr_prec_t p;
    };
};

#if GNUC_VERSION_ATLEAST(4, 3, 0)
extern void mpfr_init(cxx_mpfr & pl)
    __attribute__((error("mpfr_init must not be called on a mpfr reference -- "
                         "it is the caller's business (via a ctor)")));
extern void mpfr_init2(cxx_mpfr & pl, mpfr_prec_t)
    __attribute__((error("mpfr_init must not be called on a mpfr reference -- "
                         "it is the caller's business (via a ctor)")));
extern void mpfr_clear(cxx_mpfr & pl)
    __attribute__((error("mpfr_clear must not be called on a mpfr reference -- "
                         "it is the caller's business (via a dtor)")));
#endif

std::ostream & operator<<(std::ostream & os, cxx_mpfr const & x);

extern std::istream & operator>>(std::istream & is, cxx_mpfr::input_with_precision xp);

inline std::istream& operator>>(std::istream& is, cxx_mpfr & x)
{
    return is >> cxx_mpfr::input_with_precision { .x = x, .p = mpfr_get_default_prec() };
}

namespace fmt
{
    namespace detail {
        FMT_TYPE_CONSTANT(cxx_mpfr, double_type);
    }

    template <>
        /* reimplement native_formatter<double>, from base.h. We can't
         * inherit since the specs_ field is private there.
         */
        struct formatter<cxx_mpfr>
        {
            static constexpr fmt::detail::type TYPE = fmt::detail::type::double_type;
            using Char = char;
            private:
            fmt::detail::dynamic_format_specs<Char> specs_;
            public:
            // using nonlocking = void;
            FMT_CONSTEXPR auto parse(parse_context<Char>& ctx) -> const Char*
            {
                if (ctx.begin() == ctx.end() || *ctx.begin() == '}') 
                    return ctx.begin();
                return parse_format_specs(ctx.begin(), ctx.end(), specs_, ctx, TYPE);
            }

            auto format(cxx_mpfr const & x, format_context& ctx) const
                -> format_context::iterator;
        };
} // namespace fmt

inline int operator<=>(cxx_mpfr const & a, cxx_mpfr const & b)
{
    return mpfr_auxx::cado_mpfr_cmp(a, b);
}
inline int operator<=>(mpfr_srcptr a, cxx_mpfr const & b)
{
    return mpfr_auxx::cado_mpfr_cmp(a, b);
}
inline int operator<=>(cxx_mpfr const & a, mpfr_srcptr b)
{
    return mpfr_auxx::cado_mpfr_cmp(a, b);
}
template <typename T>
static inline int operator<=>(cxx_mpfr const & a, const T b)
    requires(std::is_integral_v<T> || std::is_floating_point_v<T>)
{
    return mpfr_auxx::cado_mpfr_cmp(a, b);
}
template <typename T>
static inline int operator<=>(const T a, cxx_mpfr const & b)
    requires(std::is_integral_v<T> || std::is_floating_point_v<T>)
{
    return -mpfr_auxx::cado_mpfr_cmp(b, a);
}

/* NOLINTBEGIN(bugprone-macro-parentheses) */
#define CXX_MPFR_DEFINE_CMP(OP)                                         \
    inline bool operator OP(cxx_mpfr const & a, cxx_mpfr const & b)     \
    {                                                                   \
        return ((a) <=> (b)) OP 0;                                      \
    }                                                                   \
    inline bool operator OP(mpfr_srcptr a, cxx_mpfr const & b)          \
    {                                                                   \
        return ((a) <=> (b)) OP 0;                                      \
    }                                                                   \
    inline bool operator OP(cxx_mpfr const & a, mpfr_srcptr b)          \
    {                                                                   \
        return ((a) <=> (b)) OP 0;                                      \
    }                                                                   \
    template <typename T>                                               \
    inline bool operator OP(cxx_mpfr const & a, const T b)              \
    requires(std::is_integral_v<T> || std::is_floating_point_v<T>)      \
    {                                                                   \
        return ((a) <=> (b)) OP 0;                                      \
    }                                                                   \
    template <typename T>                                               \
    inline bool operator OP(const T a, cxx_mpfr const & b)              \
    requires(std::is_integral_v<T> || std::is_floating_point_v<T>)      \
    {                                                                   \
        return 0 OP ((b) <=> (a));                                      \
    }
/* NOLINTEND(bugprone-macro-parentheses) */

CXX_MPFR_DEFINE_CMP(==)
CXX_MPFR_DEFINE_CMP(!=)
CXX_MPFR_DEFINE_CMP(<)
CXX_MPFR_DEFINE_CMP(>)
CXX_MPFR_DEFINE_CMP(<=)
CXX_MPFR_DEFINE_CMP(>=)

#define CXX_MPFR_DEFINE_TERNARY(OP, TEXTOP)				\
    inline cxx_mpfr operator OP(cxx_mpfr const & a, cxx_mpfr const & b)	\
    {									\
        cxx_mpfr r;							\
        mpfr_set_prec(r, mpfr_get_prec(mpfr_srcptr(a)));		\
        mpfr_auxx::cado_mpfr_##TEXTOP(r, a, b, MPFR_RNDN);		\
        return r;							\
    }									\
    template <typename T>						\
    inline cxx_mpfr operator OP(cxx_mpfr const & a, T const b)		\
    requires std::is_integral_v<T>					\
    {									\
        cxx_mpfr r;							\
        mpfr_set_prec(r, mpfr_get_prec(mpfr_srcptr(a)));		\
        mpfr_auxx::cado_mpfr_##TEXTOP(r, a, b, MPFR_RNDN);		\
        return r;							\
    }									\
    template <typename T>						\
    inline cxx_mpfr operator OP(T const a, cxx_mpfr const & b)		\
    requires std::is_integral_v<T>					\
    {									\
        cxx_mpfr r;							\
        mpfr_set_prec(r, mpfr_get_prec(mpfr_srcptr(b)));		\
        mpfr_auxx::cado_mpfr_##TEXTOP(r, a, b, MPFR_RNDN);		\
        return r;							\
    }									\
    inline cxx_mpfr & operator OP##=(cxx_mpfr & a, cxx_mpfr const & b)	\
    {									\
        mpfr_auxx::cado_mpfr_##TEXTOP(a, a, b, MPFR_RNDN);		\
        return a;							\
    }									\
    template <typename T>						\
    inline cxx_mpfr & operator OP##=(cxx_mpfr & a, T const b)		\
    requires std::is_integral_v<T>					\
    {									\
        mpfr_auxx::cado_mpfr_##TEXTOP(a, a, b, MPFR_RNDN);		\
        return a;							\
    }

CXX_MPFR_DEFINE_TERNARY(+, add)
CXX_MPFR_DEFINE_TERNARY(-, sub)
CXX_MPFR_DEFINE_TERNARY(*, mul)
CXX_MPFR_DEFINE_TERNARY(/, div)
CXX_MPFR_DEFINE_TERNARY(%, remainder)

inline cxx_mpfr operator-(cxx_mpfr const & a)
{
    cxx_mpfr r;
    mpfr_set_prec(r, mpfr_get_prec(mpfr_srcptr(a)));
    mpfr_neg(r, a, MPFR_RNDN);
    return r;
}

inline cxx_mpfr & operator<<=(cxx_mpfr & a, unsigned long const s)
{
    mpfr_mul_2exp(a, a, s, MPFR_RNDN);
    return a;
}
inline cxx_mpfr operator<<(cxx_mpfr const & a, unsigned long const s)
{
    cxx_mpfr r {a};
    mpfr_mul_2exp(r, r, s, MPFR_RNDN);
    return r;
}

inline cxx_mpfr & operator>>=(cxx_mpfr & a, unsigned long const s)
{
    mpfr_div_2exp(a, a, s, MPFR_RNDN);
    return a;
}
inline cxx_mpfr operator>>(cxx_mpfr const & a, unsigned long const s)
{
    cxx_mpfr r {a};
    mpfr_div_2exp(r, r, s, MPFR_RNDN);
    return r;
}

#endif /* UTILS_CXX_MPFR_HPP_ */
