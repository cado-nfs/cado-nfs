#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <utility>
#include <complex>
#include <map>

#include <gmp.h>
#include <mpfr.h>
#include <mpc.h>
#include "fmt/base.h"

#include "gmp_aux.h"
#include "cxx_mpc.hpp"
#include "cado_math_aux.hpp"

static void test_basic_arithmetic(cxx_gmp_randstate & state, mpfr_prec_t prec)
{
    std::map<unsigned long, cxx_mpc> v;

    cxx_mpc m;

    mpc_set_prec(m, prec);
    mpc_rootofunity(m, 6, 1, MPC_RNDNN);
    mpfr_frac(mpc_realref(m), mpc_realref(m), MPFR_RNDN);
    mpfr_frac(mpc_imagref(m), mpc_imagref(m), MPFR_RNDN);

    /* create 2000 floating point numbers, all between 1/2 and 1 (thus on
     * average 0.75), place them in a map, and sum all of them
     */
    for(int i = 0 ; i < 2000 ; i++) {
        m *= 1 + gmp_urandomm_ui(state, 255);
        m /= 1 + gmp_urandomm_ui(state, 255);
        cxx_mpc u;
        mpc_set_prec(u, prec);
        mpc_rootofunity(u, 2 + gmp_urandomm_ui(state, 14), gmp_urandomm_ui(state, 16), MPC_RNDNN);
        m *= u;
        mpfr_set_exp(mpc_realref(m), 0);
        mpfr_set_exp(mpc_imagref(m), 0);
        // fmt::print("{}\n", m);
        unsigned int t = gmp_urandomm_ui(state, 100);
        if (v.find(t) == v.end()) {
            mpc_set_prec(v[t], prec);
            mpc_set_ui(v[t], 0, MPC_RNDNN);
        }
        v[t] += m;
    }

    mpc_set_ui(m, 0, MPC_RNDNN);
    for(auto const & x : v) {
        m += x.second;
        fmt::print("[{}] {} {}\n", x.first, x.second, m);
    }
}

static void test_auxiliary_functions(mpfr_prec_t prec)
{
    cxx_mpc m;
    mpc_set_prec(m, prec);
    mpc_rootofunity(m, 5, 1, MPC_RNDNN);
    ASSERT_ALWAYS(m.prec() == prec);
    /* we have operator->, even though we don't really plan to use it */
    ASSERT_ALWAYS(m->re->_mpfr_prec == m.prec());
    ASSERT_ALWAYS(mpc_srcptr(m)->im->_mpfr_prec == m.prec());
}

/* the copy is perfectly intentional here! */
// NOLINTNEXTLINE(performance-unnecessary-value-param)
static void test_simple_copy(cxx_mpc a, cxx_mpc const & b)
{
    ASSERT_ALWAYS(a == b);
}

static void test_copies_and_assignment(mpfr_prec_t prec)
{
    cxx_mpc m, y;
    mpc_set_prec(m, prec);
    mpc_set_prec(y, 2 * prec);
    mpc_rootofunity(m, 5, 1, MPC_RNDNN);

    y = m;
    ASSERT_ALWAYS(y.prec() == prec);

    /* these two deal with real numbers only... */
    mpc_set_prec(y, 2 * prec);
    y = INT64_C(-42);
    ASSERT_ALWAYS(y.prec() == 63);
    ASSERT_ALWAYS(cxx_mpc(INT64_C(-1234)).prec() == 63);

    mpc_set_prec(y, 2 * prec);
    y = UINT64_C(42);
    ASSERT_ALWAYS(y.prec() == 64);
    ASSERT_ALWAYS(cxx_mpc(UINT64_C(1234)).prec() == 64);

    y = cado_math_aux::pow(m, 2);
    ASSERT_ALWAYS(y.prec() == prec);

    test_simple_copy(m, m);

    y = cxx_mpc(static_cast<mpc_srcptr>(m));
    ASSERT_ALWAYS(y == m);

    {
        double d = 42;
        d = std::pow(d, -2);
        m = d;
        ASSERT_ALWAYS(cxx_mpc(d) == m);
    }
    {
        long double d = 42;
        d = std::pow(d, -2);
        m = d;
        ASSERT_ALWAYS(cxx_mpc(d) == m);
    }
}

static void test_operators(mpfr_prec_t prec)
{
    cxx_mpc m, y;
    mpc_set_prec(m, prec);
    mpc_rootofunity(m, 5, 1, MPC_RNDNN);

    ASSERT_ALWAYS(m > 42);
    ASSERT_ALWAYS(-m < 42);

    y = m * m;

    ASSERT_ALWAYS(y != m);

    /* Now here's messy stuff. Imaginary parts win, so _this_ should hold */
    ASSERT_ALWAYS(y < m);
    ASSERT_ALWAYS(y <= m);
    ASSERT_ALWAYS(m > y);
    ASSERT_ALWAYS(m >= y);
    ASSERT_ALWAYS((m <=> y) > 0);

    /* This is slightly more difficult */
    ASSERT_ALWAYS((m - y) < y);
    ASSERT_ALWAYS(m < (y << 1));
    ASSERT_ALWAYS((m >> 1) < y);

    ASSERT_ALWAYS((y >> 1) + (y >> 1) == y);
    ASSERT_ALWAYS(y + y == (y << 1));

    ASSERT_ALWAYS(((4*m)%y) < 0);
    y >>= 2;
    m >>= 1;
    ASSERT_ALWAYS(y == m * m);
    y <<= 10;
    m <<= 5;
    ASSERT_ALWAYS(y == m * m);
    ASSERT_ALWAYS(bool(m));
    ASSERT_ALWAYS(!bool(m-m));
    ASSERT_ALWAYS(bool(static_cast<cxx_mpc const &>(m)));
    ASSERT_ALWAYS(!bool(static_cast<cxx_mpc const &>(m - m)));
    std::swap(y, m);
    ASSERT_ALWAYS(m == y * y);

    {
        std::complex<double> const z { 1, 1 };
        cxx_mpc zz(z);
        m = 2 - zz;
        zz = cado_math_aux::pow(zz, 4);
        ASSERT_ALWAYS(zz + 4 == 0);
        ASSERT_ALWAYS(zz == cado_math_aux::pow(m, 4));
        zz = std::complex<double>(1, -1);
        // fmt::print("{}\n", zz);
        // fmt::print("{}\n", m);
        ASSERT_ALWAYS(zz == m);
    }
    {
        std::complex<long double> const z { 1, 1 };
        cxx_mpc zz(z);
        m = 2 - zz;
        zz = cado_math_aux::pow(zz, 4);
        ASSERT_ALWAYS(zz + 4 == 0);
        ASSERT_ALWAYS(zz == cado_math_aux::pow(m, 4));
        zz = std::complex<long double>(1, -1);
        ASSERT_ALWAYS(zz == m);
    }
}


// coverity[root_function]
int main(int argc, char const * argv[])
{
    cxx_gmp_randstate state;

    if (argc > 1) {
        unsigned long seed;
        seed = strtoul(argv[1], nullptr, 0);
        gmp_randseed_ui(state, seed);
    }

    mpfr_prec_t const prec = 256;

    test_basic_arithmetic(state, prec);
    test_auxiliary_functions(prec);
    test_copies_and_assignment(prec);
    test_operators(prec);

}
