#include "cado.h"

#include <cfenv>
#include <cmath>
#include <cstdlib>

#include <limits>

#include <algorithm>
#if __cplusplus >= 202302L && defined(__STDCPP_FLOAT128_T__)
#include <stdfloat>
#endif

#include "fmt/core.h"
#include <gmp.h>

#include "cado_math_aux.hpp"
#include "cxx_mpz.hpp"

using namespace cado_math_aux;

/* This way of forming examples is awkward because we need cooperation
 * from the floating point rounding modes
 */
template <typename T> struct example_gen_1 {
    static constexpr char const * name = "gen_1";
    static void gen(T & zd, cxx_mpz & z, cxx_gmp_randstate & rstate)
    {
        cxx_mpz a[3];
        zd = 0;
        mpz_set_si(z, 0);
        temporary_round_mode dummy(FE_TOWARDZERO);

        constexpr int M = std::numeric_limits<double>::digits;

        for (int i = 0; i < 3; i++) {
            mpz_rrandomb(a[i], rstate, M);
            /* The += here must round toward zero if we want to remain
             * consistent.
             */
            zd += std::ldexp((T)mpz_get_d(a[i]), M * i);
            mpz_mul_2exp(a[i], a[i], M * i);
            mpz_add(z, z, a[i]);
        }

        /* Shuffle size and exponent range */
        if (gmp_urandomb_ui(rstate, 1)) {
            zd = -zd;
            mpz_neg(z, z);
        }
        int e = gmp_urandomm_ui(rstate, 4 * M);
        mpz_tdiv_q_2exp(z, z, e);
        zd = std::ldexp(zd, -e);
    }
};

/* This one, after all, is more direct. It does rely on mpz_to<> doing
 * the right thing, though.
 */
template <typename T> struct example_gen_2 {
    static constexpr char const * name = "gen_2";
    static void gen(T & zd, cxx_mpz & z, cxx_gmp_randstate & rstate)
    {
        constexpr int M = std::numeric_limits<double>::digits;
        constexpr int E = std::numeric_limits<double>::max_exponent;
        constexpr int ET = std::numeric_limits<T>::max_exponent;

        mpz_rrandomb(z, rstate, std::min(2 * E, ET));

        zd = mpz_to<T>(z);

        if (gmp_urandomb_ui(rstate, 1)) {
            zd = -zd;
            mpz_neg(z, z);
        }
        int e = gmp_urandomm_ui(rstate, 4 * M);
        mpz_tdiv_q_2exp(z, z, e);
        zd = std::ldexp(zd, -e);
    }
};

template <typename T, template <typename> class G>
static void dotest(cxx_gmp_randstate & rstate)
{
    cxx_mpz z, zc;

    constexpr int M = std::numeric_limits<double>::digits;
    constexpr int MT = std::numeric_limits<T>::digits;

    fmt::print("M={} MT={} gen={}\n", M, MT, G<T>::name);

    for (int i = 0; i < 100; i++) {
        T zd = 0;

        /* form a wide enough integer, together with its closest
         * representative as type T (rounding towards zero, so that we 
         * have 0 <= abs(z)-abs(zd) < ulp(zd)
         */
        G<T>::gen(zd, z, rstate);

        int b = mpz_sizeinbase(z, 2) - (mpz_cmp_ui(z, 0) == 0);

        /* Now z is a random integer, and zd should be its most accurate
         * representative as type T. Check that mpz_from confirms this.
         */

        zc = mpz_from(zd);

        cxx_mpz diff;
        mpz_sub(diff, z, zc);

        int db = mpz_sizeinbase(diff, 2) - (mpz_cmp_ui(diff, 0) == 0);

        // fmt::print(stdout, "Info: b={} db={} z={} zd={}\n", b, db, z, zd);

        if (db > std::max(b - MT, 0)) {
            /* shit happens, but we want to investigate only then.
             */
            fmt::print(stderr, "Error: b={} db={} z={} zc={}\n", b, db, z, zc);
            abort();
        }
    }
}

int main()
{
    cxx_gmp_randstate rstate;

    dotest<double, example_gen_1>(rstate);
    dotest<long double, example_gen_1>(rstate);
    dotest<double, example_gen_2>(rstate);
    dotest<long double, example_gen_2>(rstate);

    // this one can't work because even though gcc boasts a __float128
    // data type, it can't do ldexp on it (why?). It does seem that g++
    // in c++23 mode has __STDCPP_FLOAT128_T__ *AND*
    // _GLIBCXX_HAVE_FLOAT128_MATH, meaning that an ldexp function is
    // available in that case

    /* This float128 (raw hex bytes in memory) yields a wrong result of
     * mpz_from().
     *
     * input
     * 0x4077fffe0000007fffffffc0003fffff
     * 2.65841542675122890537650317297503411e+36
     *
     * expected
     * 0x1fffe0000007fffffffc0003fffff00
     * 2658415426751228905376503172975034112
     *
     * got
     * 2658415426751228905448699303330611456)
     * 0x1fffe000000800001003dffc0008100
     */
#if __cplusplus >= 202302L && defined(__STDCPP_FLOAT128_T__)
    dotest<std::float128_t, example_gen_1>(rstate);
    dotest<std::float128_t, example_gen_2>(rstate);
#endif
}
