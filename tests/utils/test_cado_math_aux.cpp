#include "cado.h"       // IWYU pragma: keep

#include <cstdlib>
#include <cmath>

#include <limits>
#include <type_traits>

#include <gmp.h>
#include "fmt/base.h"

#include "number_context.hpp"
#include "cxx_mpz.hpp"
#ifdef HAVE_MPFR
#include "cxx_mpfr.hpp"
#endif
#include "cado_math_aux.hpp"
#include "cado_compile_time_hacks.hpp"
#include "cado_addsubmul.hpp"
#include "gmp_aux.h"
#include "macros.h"

static void test_cado_math_aux_trivial()
{
    using cado_math_aux::multiply_by_poweroftwo;
    static_assert(multiply_by_poweroftwo<2>(42) == 168);
    ASSERT_ALWAYS(multiply_by_poweroftwo<2>(42) == 168);
}

/* most of the mpfr-related testing here could/should be done within
 * test_cxx_mpfr, really.
 *
 * Manu of the cado_math_aux trampolines that reach the C++ standard library
 * receive testing from other test files.
 */

static void test_cado_math_aux_pow()
{
    using cado_math_aux::pow;
    using cado_math_aux::isnan;
    using cado_math_aux::isinf;

#ifdef HAVE_MPFR
    {
        cxx_mpfr x = 2;
        cxx_mpfr y = 1/cxx_mpfr(2);
        y = pow(x, y);
        cxx_mpfr z = pow(y, 2);
        fmt::print("{} {} {}\n", x, y, z);
        ASSERT_ALWAYS(z > cxx_mpfr(1.99999));
        ASSERT_ALWAYS(z < cxx_mpfr(2.00001));
        ASSERT_ALWAYS(!isnan(z));
        ASSERT_ALWAYS(!isinf(z));
    }
#endif
}

static void test_cado_math_aux_ldexp_frexp()
{
    using cado_math_aux::ldexp;
    using cado_math_aux::frexp;
    using cado_math_aux::pow;

#ifdef HAVE_MPFR
    {
        auto x = pow(ldexp(cxx_mpfr(.615699059), -32), cxx_mpfr(-1)>>3);
        fmt::print("{}\n", x);
        ASSERT_ALWAYS(x > cxx_mpfr(16.99999));
        ASSERT_ALWAYS(x < cxx_mpfr(17.00001));
        int e;
        auto y = frexp(x + 1, &e);
        fmt::print("{} {}\n", y, e);
        ASSERT_ALWAYS(e == 5);
        ASSERT_ALWAYS(y > 0.562500);
        ASSERT_ALWAYS(y < 0.562501);
    }
#endif
    long double const long_pi = 3.141592653589793238462643383279502884L;
    int long_pi_exponent;
    long double const long_pi_mantissa = frexp(long_pi, &long_pi_exponent);
    long double const down = .78539816339744830962L;
    long double const up = .78539816339744830963L;
    ASSERT_ALWAYS(long_pi_mantissa >= down);
    ASSERT_ALWAYS(long_pi_mantissa <= up);
    ASSERT_ALWAYS(long_pi_exponent == 2);
}

static void test_cado_math_aux_abs()
{
    using cado_math_aux::abs;

#ifdef HAVE_MPFR
    {
        cxx_mpfr pi, sqrt2;
        mpfr_const_pi(pi, MPFR_RNDN);
        mpfr_sqrt_ui(sqrt2, 2, MPFR_RNDN);
        cxx_mpfr const z = abs(pi - 3*sqrt2);
        ASSERT_ALWAYS(z > 1.101048);
        ASSERT_ALWAYS(z < 1.101049);
    }
#endif
}

template<typename T>
struct working_precision : public std::integral_constant<int, std::numeric_limits<T>::digits> {};

#ifdef HAVE_MPFR
template<>
struct working_precision<cxx_mpfr> : public std::integral_constant<int, 200> {};
#endif


template<typename T>
static void test_one_fma()
{
    using cado_math_aux::ldexp;
    using cado_math_aux::fma;
    using cado_math_aux::fms;
    using cado_math_aux::addmul;
    using cado_math_aux::submul;
    constexpr int bits = working_precision<T>::value;

    /* It's very important that we use the constant one with the required
     * final precision, because we won't have automatic type promotion
     * for cxx_mpfr's.
     */

    const cado::number_context<T> tr(bits);

    T const one = tr(1);

    constexpr int b = bits - 1;
    T const epsilon = ldexp(one, -b);
    T const rd = ldexp(one - epsilon, b);
    T const ru = ldexp(one + epsilon, b);
    T const t = rd * ru;
    T const bigone = ldexp(one, 2*b);
    ASSERT_ALWAYS(t == bigone);

    fmt::print("{} {:a} {:a} {:a}\n", typeid(T).name(), rd, ru, one);

    {
        /* XXX Attention here. fma *must really* be cado_math_aux::fma, which
         * jumps to std::fma, and **NOT** the bare ::fma, which is defined by
         * cmath and only addresses double precision!!
         */
        T const s = fma(rd, ru, -bigone);
        fmt::print("{} fma {:a}\n", typeid(T).name(), -s);
        ASSERT_ALWAYS(s == -one);
    }
    {
        T const s = fms(rd, ru, bigone);
        fmt::print("{} fms {:a}\n", typeid(T).name(), -s);
        ASSERT_ALWAYS(s == -one);
    }
    {
        T s = -bigone;
        addmul(s, rd, ru);
        fmt::print("{} addmul {:a}\n", typeid(T).name(), -s);
        ASSERT_ALWAYS(s == -one);
        ASSERT_ALWAYS(s + one == tr(0));
    }
    {
        T s = bigone;
        submul(s, rd, ru);
        fmt::print("{} submul {:a}\n", typeid(T).name(), s);
        ASSERT_ALWAYS(s == one);
    }
}

template<>
void test_one_fma<cxx_mpz>()
{
    using cado_math_aux::fma;
    using cado_math_aux::fms;
    using cado_math_aux::addmul;
    using cado_math_aux::submul;
    cxx_gmp_randstate state;
    cxx_mpz a, b, c;
    mpz_rrandomb(a, state, 4 * GMP_LIMB_BITS);
    mpz_rrandomb(b, state, 4 * GMP_LIMB_BITS);
    mpz_rrandomb(c, state, 4 * GMP_LIMB_BITS);
    fmt::print("cxx_mpz {} {} {}\n", a, b, c);
    {
        cxx_mpz z = a;
        addmul(z, b, c);
        fmt::print("cxx_mpz addmul {}\n", z);
        cxx_mpz const zz = fma(b, c, a);
        fmt::print("cxx_mpz fma {}\n", zz);
        ASSERT_ALWAYS(z == zz);
    }
    {
        cxx_mpz z = a;
        submul(z, b, c);
        fmt::print("cxx_mpz submul {}\n", z);
        cxx_mpz const zz = fms(b, c, a);
        fmt::print("cxx_mpz fms {}\n", zz);
        ASSERT_ALWAYS(z + zz == 0);
    }
}


static void test_cado_math_aux_fma()
{
    const bool do_ld = !cado_math_aux::valgrind_long_double_hopeless();

    using cado_math_aux::fma;
    using cado_math_aux::fms;
    using cado_math_aux::addmul;
    using cado_math_aux::submul;

    test_one_fma<float>();
    test_one_fma<double>();

    if (do_ld)
        test_one_fma<long double>();

#ifdef HAVE_MPFR
    test_one_fma<cxx_mpfr>();
#endif
    test_one_fma<cxx_mpz>();
}


int main(int argc, char const * argv[])
{
    cxx_gmp_randstate state;

    if (argc > 1) {
        unsigned long seed;
        seed = strtoul(argv[1], nullptr, 0);
        gmp_randseed_ui(state, seed);
    }

    test_cado_math_aux_trivial();

    test_cado_math_aux_pow();
    test_cado_math_aux_ldexp_frexp();
    test_cado_math_aux_abs();
    test_cado_math_aux_fma();

    return 0;
}
