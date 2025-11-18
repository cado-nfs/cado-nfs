#ifndef CADO_UTILS_CADO_MATH_AUX_HPP
#define CADO_UTILS_CADO_MATH_AUX_HPP

#include "cado_config.h"

#include <cmath>
#include <cfenv>
#include <cstdint>
#include <climits>

#include <type_traits>
#include <complex>

#include <gmp.h>

#include "cxx_mpz.hpp"
#include "utils_cxx.hpp"
#include "gmp_aux.h"
#include "macros.h"

#ifdef HAVE_MPFR
#include <mpfr.h>
#include "cxx_mpfr.hpp"
#endif

#ifdef HAVE_MPC
#include <mpc.h>
#include "cxx_mpc.hpp"
#include "mpc_aux.h"
#endif

namespace cado_math_aux
{

#if 0
    /* {{{ POC: wrap around std::pow */
    template<typename T, typename = void>
        constexpr bool has_pow = false;
    template<typename T>
        constexpr bool has_pow<T,
                std::void_t<decltype(std::pow(std::declval<T>(),1))>
                  > = true;

    static_assert(has_pow<float>);
    static_assert(has_pow<double>);

    template<typename T>
        T pow(T const & x, int e)
        requires has_pow<T>
        {
            return std::pow(x, e);
        }
    template<typename T>
        T pow(T const & x, T const & e)
        requires has_pow<T>
        {
            return std::pow(x, e);
        }

#ifdef HAVE_MPFR
    static_assert(!has_pow<cxx_mpfr>);
    inline cxx_mpfr pow(cxx_mpfr const & x, int e)
    {
        cxx_mpfr res;
        mpfr_set_prec(res, x.prec());
        mpfr_pow_si(res, x, e, MPFR_RNDN);
        return res;
    }
    inline cxx_mpfr pow(cxx_mpfr const & x, cxx_mpfr const & e)
    {
        cxx_mpfr res;
        mpfr_set_prec(res, x.prec());
        mpfr_pow(res, x, e, MPFR_RNDN);
        return res;
    }
#endif
    /* }}} */
#endif

    /* {{{ simple wrappers around std::pow, + cxx_mpfr extensions */
    /* see comment about ldexpf and friends
    inline float pow(float x, float e) { return std::pow(x, e); }
    inline double pow(double x, double e) { return std::pow(x, e); }
    inline long double pow(long double x, long double e) { return std::pow(x, e); }
    */
    using std::pow;
    // inline float pow(float x, long e) { return std::pow(x, e); }
    // inline double pow(double x, int e) { return std::pow(x, e); }
    // inline long double pow(long double x, int e) { return std::pow(x, e); }
#ifdef HAVE_MPFR
    inline cxx_mpfr pow(cxx_mpfr const & x, cxx_mpfr const & e)
    {
        cxx_mpfr res;
        mpfr_set_prec(res, x.prec());
        mpfr_pow(res, x, e, MPFR_RNDN);
        return res;
    }
    inline cxx_mpfr pow(cxx_mpfr const & x, int e)
    {
        cxx_mpfr res;
        mpfr_set_prec(res, x.prec());
        mpfr_pow_si(res, x, e, MPFR_RNDN);
        return res;
    }
#endif
#ifdef HAVE_MPC
    inline cxx_mpc pow(cxx_mpc const & x, cxx_mpc const & e)
    {
        cxx_mpc res;
        mpc_set_prec(res, x.prec());
        mpc_pow(res, x, e, MPC_RNDNN);
        return res;
    }
    inline cxx_mpc pow(cxx_mpc const & x, int e)
    {
        cxx_mpc res;
        mpc_set_prec(res, x.prec());
        mpc_pow_si(res, x, e, MPC_RNDNN);
        return res;
    }
#endif
    inline cxx_mpz pow(cxx_mpz const & x, int e)
    {
        ASSERT_ALWAYS(e >= 0);
        cxx_mpz res;
        mpz_pow_ui(res, x, e);
        return res;
    }
    template<typename T>
    inline cxx_mpz pow(T const & x, int e)
    requires std::is_integral_v<T>
    {
        /* I'm lazy */
        return converter_from_mpz<T>(pow(cxx_mpz(x), e));
    }
    /* }}} */
    /* {{{ simple wrappers around std::isnan, + cxx_mpfr extensions */
    /* see comment about ldexpf and friends
    inline bool isnan(float x) { return std::isnan(x); }
    inline bool isnan(double x) { return std::isnan(x); }
    inline bool isnan(long double x) { return std::isnan(x); }
    */
    using std::isnan;
    template<typename T>
    inline bool isnan(std::complex<T> const & x) {
        return std::isnan(x.real()) || std::isnan(x.imag());
    }
#ifdef HAVE_MPFR
    inline bool isnan(cxx_mpfr const & x) { return mpfr_nan_p(mpfr_srcptr(x)); }
#endif
#ifdef HAVE_MPC
    inline bool isnan(cxx_mpc const & x) {
        return mpfr_nan_p(mpc_realref(x)) || mpfr_nan_p(mpc_imagref(x));
    }
#endif
    /* }}} */
    /* {{{ simple wrappers around std::isinf, + cxx_mpfr extensions */
    /* see comment about ldexpf and friends
    inline bool isinf(float x) { return std::isinf(x); }
    inline bool isinf(double x) { return std::isinf(x); }
    inline bool isinf(long double x) { return std::isinf(x); }
    */
    using std::isinf;
    template<typename T>
    inline bool isinf(std::complex<T> const & x) {
        return std::isinf(x.real()) || std::isinf(x.imag());
    }
#ifdef HAVE_MPFR
    inline bool isinf(cxx_mpfr const & x) { return mpfr_inf_p(mpfr_srcptr(x)); }
#endif
#ifdef HAVE_MPC
    inline bool isinf(cxx_mpc const & x) {
        return mpfr_inf_p(mpc_realref(x)) || mpfr_inf_p(mpc_imagref(x));
    }
#endif
    /* }}} */
    /* {{{ simple wrappers around std::ldexp, + cxx_mpfr extensions */

    /* we used to have std::{ldexpf,ldexp,ldexpl} here, but see the drama
     * with https://gcc.gnu.org/bugzilla/show_bug.cgi?id=79700 ; so I
     * think that simply importing std::ldexp _should_ be enough
    inline float ldexp(float x, int e) { return std::ldexp(x, e); }
    inline double ldexp(double x, int e) { return std::ldexp(x, e); }
    inline long double ldexp(long double x, int e) { return std::ldexp(x, e); }
     */
    using std::ldexp;
#ifdef HAVE_MPFR
    inline cxx_mpfr ldexp(cxx_mpfr const & x, int e) {
        cxx_mpfr res = x;
        if (e > 0)
            mpfr_mul_2exp(res, x, e, MPFR_RNDN);
        if (e < 0)
            mpfr_div_2exp(res, x, -e, MPFR_RNDN);
        return res;
    }
#endif
    /* }}} */
    /* {{{ simple wrappers around std::frexp, + cxx_mpfr extensions */
    /* see above
    inline float frexp(float x, int * e) { return std::frexp(x, e); }
    inline double frexp(double x, int * e) { return std::frexp(x, e); }
    inline long double frexp(long double x, int * e) { return std::frexp(x, e); }
    */
    using std::frexp;
#ifdef HAVE_MPFR
    inline cxx_mpfr frexp(cxx_mpfr const & x, int * e) {
        mpfr_exp_t ee;
        cxx_mpfr res = x;
        mpfr_frexp(&ee, res, x, MPFR_RNDN);
        *e = ee;
        return res;
    }
#endif
    /* }}} */
    /* {{{ simple wrappers around std::abs, + cxx_mpfr extensions */
    /* see above
    inline float abs(float x) { return std::abs(x); }
    inline double abs(double x) { return std::abs(x); }
    inline long double abs(long double x) { return std::abs(x); }
    */
    using std::abs;
    /* abs is also defined on std::complex */
#ifdef HAVE_MPFR
    inline cxx_mpfr abs(cxx_mpfr const & x) {
        cxx_mpfr res = x;
        mpfr_abs(res, res, MPFR_RNDN);
        return res;
    }
#endif
#ifdef HAVE_MPC
    inline cxx_mpfr abs(cxx_mpc const & x) {
        cxx_mpfr res;
        mpfr_set_prec(res, x.prec());
        mpc_abs(res, x, MPFR_RNDN);
        return res;
    }
#endif
    /* }}} */
    /* {{{ simple wrappers around std::log and std::exp */
    using std::log;
    using std::exp;
#ifdef HAVE_MPFR
    inline cxx_mpfr exp(cxx_mpfr const & x) {
        cxx_mpfr res;
        mpfr_set_prec(res, x.prec());
        mpfr_exp(res, x, MPFR_RNDN);
        return res;
    }
    inline cxx_mpfr log(cxx_mpfr const & x) {
        cxx_mpfr res;
        mpfr_set_prec(res, x.prec());
        mpfr_log(res, x, MPFR_RNDN);
        return res;
    }
#endif
#ifdef HAVE_MPC
    inline cxx_mpc exp(cxx_mpc const & x) {
        cxx_mpc res;
        mpc_set_prec(res, x.prec());
        mpc_exp(res, x, MPC_RNDNN);
        return res;
    }
    inline cxx_mpc log(cxx_mpc const & x) {
        cxx_mpc res;
        mpc_set_prec(res, x.prec());
        mpc_log(res, x, MPC_RNDNN);
        return res;
    }
#endif
    /* }}} */
    /* {{{ fma/fms operation -- we have to use a function in order to have a
     * unified interface
     */
    /* see above
    inline float fma(float x, float y, float z) { return std::fmaf(x, y, z); }
    inline double fma(double x, double y, double z) { return std::fma(x, y, z); }
    inline long double fma(long double x, long double y, long double z) { return std::fmal(x, y, z); }
    */
    using std::fma;
    inline float fms(float x, float y, float z) { return std::fmaf(x, y, -z); }
    inline double fms(double x, double y, double z) { return std::fma(x, y, -z); }
    inline long double fms(long double x, long double y, long double z) { return std::fma(x, y, -z); }
    inline cxx_mpz fma(cxx_mpz const & x, cxx_mpz const & y, cxx_mpz const & z) {
        cxx_mpz res = z;
        mpz_addmul(res, x, y);
        return res;
    }
    inline cxx_mpz fms(cxx_mpz const & x, cxx_mpz const & y, cxx_mpz const & z) {
        cxx_mpz res = -z;
        mpz_addmul(res, x, y);
        return res;
    }
#ifdef HAVE_MPFR
    /* the working precision is the precision of x */
    inline cxx_mpfr fma(cxx_mpfr const & x, cxx_mpfr const & y, cxx_mpfr const & z) {
        cxx_mpfr res;
        mpfr_set_prec(res, x.prec());
        mpfr_fma(res, x, y, z, MPFR_RNDN);
        return res;
    }
    inline cxx_mpfr fms(cxx_mpfr const & x, cxx_mpfr const & y, cxx_mpfr const & z) {
        cxx_mpfr res;
        mpfr_set_prec(res, x.prec());
        mpfr_fms(res, x, y, z, MPFR_RNDN);
        return res;
    }
#endif
#ifdef HAVE_MPC
    /* the working precision is the precision of x */
    inline cxx_mpc fma(cxx_mpc const & x, cxx_mpc const & y, cxx_mpc const & z) {
        cxx_mpc res;
        mpc_set_prec(res, x.prec());
        mpc_fma(res, x, y, z, MPC_RNDNN);
        return res;
    }
    inline cxx_mpc fms(cxx_mpc const & x, cxx_mpc const & y, cxx_mpc const & z) {
        cxx_mpc res;
        mpc_set_prec(res, x.prec());
        cxx_mpc nz = -z;
        mpc_fms(res, x, y, nz, MPC_RNDNN);
        return res;
    }
#endif
    template<typename T>
    inline std::complex<T> fma(std::complex<T> const & x, std::complex<T> const & y, std::complex<T> const & z) {
        /* res = x * y + z
         * re(res) = re(x) * re(y) - im(x) * im(y) + re(z)
         * im(res) = re(x) * im(y) + im(x) * re(y) + im(z)
         */
        using cado_math_aux::fma;
        using cado_math_aux::fms;
        auto const [ xr, xi ] = std::make_pair<T, T>(x.real(), x.imag());
        auto const [ yr, yi ] = std::make_pair<T, T>(y.real(), y.imag());
        auto const [ zr, zi ] = std::make_pair<T, T>(z.real(), z.imag());
        return {
            fms(xr, yr, fms(xi, yi, zr)),
            fma(xr, yi, fma(xi, yr, zi))
        };
    }
    template<typename T>
    inline std::complex<T> fms(std::complex<T> const & x, std::complex<T> const & y, std::complex<T> const & z) {
        /* res = x * y - z
         * re(res) = re(x) * re(y) - im(x) * im(y) - re(z)
         * im(res) = re(x) * im(y) + im(x) * re(y) - im(z)
         */
        using cado_math_aux::fma;
        using cado_math_aux::fms;
        auto const [ xr, xi ] = std::make_pair<T, T>(x.real(), x.imag());
        auto const [ yr, yi ] = std::make_pair<T, T>(y.real(), y.imag());
        auto const [ zr, zi ] = std::make_pair<T, T>(z.real(), z.imag());
        return {
            fms(xr, yr, fma(xi, yi, zr)),
            fma(xr, yi, fms(xi, yr, zi))
        };
    }
    /* }}} */
    // {{{ divexact
    template<typename T>
        T divexact(T a, T b)
        requires std::is_integral_v<T>
    {
        return a / b;
    }
    template<typename T>
        cxx_mpz divexact(cxx_mpz const & a, T b)
        requires cado::converts_via<T, uint64_t>
    {
        return mpz_divexact_uint64(a, b);
    }
    template<typename T>
        cxx_mpz divexact(cxx_mpz const & a, T b)
        requires cado::converts_via<T, int64_t>
    {
        cxx_mpz res;
        if (b > 0)
            mpz_divexact_uint64(res, a, b);
        else {
            mpz_divexact_uint64(res, a, -b);
            mpz_neg(res, res);
        }
        return res;
    }
    inline cxx_mpz divexact(cxx_mpz const & a, cxx_mpz const & b)
    {
        cxx_mpz res;
        mpz_divexact(res, a, b);
        return res;
    }
    // }}}

    /* This compares two floating point values with a relative error margin
     */

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
        const volatile double d0 = 1.79769313486231570815e+308; // 0x1.fffffffffffffp+1023;
        const volatile double d1 = 1.9958403095347196e+292;     // 0x1.fffffffffffffp+970;
        volatile long double ld = d0;
        ld += d1;
        int e;
        std::frexp(ld, &e);
        return e != 1025;
    }

}       /* namespace cado_math_aux */


#endif	/* CADO_UTILS_CADO_MATH_AUX_HPP */
