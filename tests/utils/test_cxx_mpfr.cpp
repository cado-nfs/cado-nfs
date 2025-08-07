#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cmath>

#include <utility>
#include <map>

#include <gmp.h>
#include "fmt/base.h"

#include "gmp_aux.h"
#include "cxx_mpfr.hpp"
#include "mpfr_aux.h"
#include "mpfr_auxx.hpp"
#include "cado_math_aux.hpp"
#include "macros.h"

static void test_basic_arithmetic(cxx_gmp_randstate & state, mpfr_prec_t prec)
{
    std::map<unsigned long, cxx_mpfr> v;

    cxx_mpfr m;
    mpfr_set_prec(m, prec);
    mpfr_const_pi(m, MPFR_RNDN);
    mpfr_frac(m, m, MPFR_RNDN);

    /* create 2000 floating point numbers, all between 1/2 and 1 (thus on
     * average 0.75), place them in a map, and sum all of them
     */
    for(int i = 0 ; i < 2000 ; i++) {
        m *= 1 + gmp_urandomm_ui(state, 255);
        m /= 1 + gmp_urandomm_ui(state, 255);

        mpfr_set_exp(m, 0);
        // fmt::print("{}\n", m);
        unsigned int t = gmp_urandomm_ui(state, 100);
        if (v.find(t) == v.end()) {
            mpfr_set_prec(v[t], prec);
            mpfr_set_zero(v[t], 1);
        }
        v[t] += m;
    }

    mpfr_set_zero(m, 1);
    for(auto const & x : v) {
        m += x.second;
        fmt::print("[{}] {} {}\n", x.first, x.second, m);
    }
}

static void test_auxiliary_functions(mpfr_prec_t prec)
{
    cxx_mpfr m;
    mpfr_set_prec(m, prec);
    mpfr_const_pi(m, MPFR_RNDN);
    ASSERT_ALWAYS(m.prec() == prec);
    ASSERT_ALWAYS(m.sgn() == 1);
    ASSERT_ALWAYS((-m).sgn() == -1);
    /* we have operator->, even though we don't really plan to use it */
    ASSERT_ALWAYS(m->_mpfr_prec == m.prec());
    ASSERT_ALWAYS(static_cast<cxx_mpfr const &>(m)->_mpfr_prec == m.prec());
}

/* the copy is perfectly intentional here! */
// NOLINTNEXTLINE(performance-unnecessary-value-param)
static void test_simple_copy(cxx_mpfr a, cxx_mpfr const & b)
{
    ASSERT_ALWAYS(a == b);
}

template<typename T>
static void init_and_compare(T c)
{
    mpfr_t yr;
    mpfr_auxx::cado_mpfr_init_set(yr, c, MPFR_RNDN);
    cxx_mpfr y;
    mpfr_set_prec(y, mpfr_get_default_prec());
    cado_math_aux::similar_set(y, c);
    ASSERT_ALWAYS(y == cxx_mpfr(yr));
    ASSERT_ALWAYS(mpfr_auxx::cado_mpfr_cmp(yr, c) == 0);

    mpfr_auxx::cado_mpfr_add(yr, yr, c, MPFR_RNDN);
    y <<= 1;
    ASSERT_ALWAYS(y == cxx_mpfr(yr));

    mpfr_auxx::cado_mpfr_addmul(yr, yr, 3, MPFR_RNDN);
    y <<= 2;
    ASSERT_ALWAYS(y == cxx_mpfr(yr));

    y = c;
    mpfr_auxx::cado_mpfr_submul(yr, y, 6, MPFR_RNDN);
    mpfr_auxx::cado_mpfr_sub(yr, yr, c, MPFR_RNDN);
    ASSERT_ALWAYS(c == cxx_mpfr(yr));
    y = 2*c - y;
    ASSERT_ALWAYS(c == y);
    y += 31*c;
    mpfr_auxx::cado_mpfr_mul(y, y, 2, MPFR_RNDN);
    ASSERT_ALWAYS(y == 64*c);

    y = c + y;
    y = 2 * y;
    ASSERT_ALWAYS(y == 130*c);

    mpfr_clear(yr);
}

static void test_copies_and_assignment(mpfr_prec_t prec)
{
    cxx_mpfr m, y;
    mpfr_set_prec(m, prec);
    mpfr_set_prec(y, 2 * prec);
    mpfr_const_pi(m, MPFR_RNDN);

    y = m;
    ASSERT_ALWAYS(y.prec() == prec);

    mpfr_set_prec(y, 2 * prec);
    y = INT64_C(-42);
    ASSERT_ALWAYS(y.prec() == 63);
    ASSERT_ALWAYS(cxx_mpfr(INT64_C(-1234)).prec() == 63);

    /* try to reach as many of the trampoline functions as we can. We
     * have several, and some of them are actually only visible on
     * 32-bit.
     */

    init_and_compare<int64_t>(-42);
    init_and_compare<long>(-42L);
    init_and_compare<uint64_t>(1234);
    init_and_compare<unsigned long>(1234UL);


    mpfr_set_prec(y, 2 * prec);
    y = UINT64_C(42);
    ASSERT_ALWAYS(y.prec() == 64);
    ASSERT_ALWAYS(y == cxx_mpfr(UINT64_C(42)));
    ASSERT_ALWAYS(cxx_mpfr(UINT64_C(1234)).prec() == 64);

    y = cado_math_aux::pow(m, 2);
    ASSERT_ALWAYS(y.prec() == prec);

    test_simple_copy(m, m);

    y = static_cast<mpfr_srcptr>(m);
    ASSERT_ALWAYS(y == m);

    {
        double d = 42;
        d = std::pow(d, -2);
        m = d;
        ASSERT_ALWAYS(cxx_mpfr(d) == m);
    }
    {
        long double d = 42;
        d = std::pow(d, -2);
        m = d;
        ASSERT_ALWAYS(cxx_mpfr(d) == m);
    }
}

static void test_operators(mpfr_prec_t prec)
{
    cxx_mpfr m, y;
    mpfr_set_prec(m, prec);
    mpfr_const_pi(m, MPFR_RNDN);

    y = m * m;

    ASSERT_ALWAYS(y != m);
    ASSERT_ALWAYS(y > m);
    ASSERT_ALWAYS(y >= m);
    ASSERT_ALWAYS(m < y);
    ASSERT_ALWAYS(m <= y);
    static_assert((1 <=> 2) < 0);
    ASSERT_ALWAYS((m <=> y) < 0);
    ASSERT_ALWAYS((m <=> (mpfr_srcptr) y) < 0);
    ASSERT_ALWAYS(((mpfr_srcptr) m <=> y) < 0);
    ASSERT_ALWAYS((y <=> 10) < 0);
    ASSERT_ALWAYS((9 <=> y) < 0);
    ASSERT_ALWAYS((y - m) > m);
    ASSERT_ALWAYS(y < (m << 2));
    ASSERT_ALWAYS((y >> 1) > m);
    ASSERT_ALWAYS((y >> 1) + (y >> 1) == y);
    ASSERT_ALWAYS(y + y == (y << 1));
    /* pi^2-3*pi < 0.5 */
    ASSERT_ALWAYS((y % m) * 2 < 1);
    y >>= 2;
    m >>= 1;
    ASSERT_ALWAYS(y == m * m);
    y <<= 10;
    m <<= 5;
    ASSERT_ALWAYS(y == m * m);
    ASSERT_ALWAYS(bool(m));
    ASSERT_ALWAYS(!bool(m-m));
    ASSERT_ALWAYS(bool(static_cast<cxx_mpfr const &>(m)));
    ASSERT_ALWAYS(!bool(static_cast<cxx_mpfr const &>(m - m)));
    std::swap(y, m);
    ASSERT_ALWAYS(m == y * y);
}

static void test_mpfr_aux(mpfr_prec_t prec)
{
    /* mpfr_aux.h has several functions that aren't reached by the test
     * code above, let's try to catch them.
     */
    cxx_mpfr m, y;
    mpfr_set_prec(m, prec);

    {
        mpfr_const_pi(m, MPFR_RNDN);
        y = m;
        mpfr_addmul_ui(m, y, 3, MPFR_RNDN);
        ASSERT_ALWAYS((m >> 2) == y);
        mpfr_submul_ui(m, y, 5, MPFR_RNDN);
        ASSERT_ALWAYS((m + y) == 0);
    }

    {
        mpfr_const_pi(m, MPFR_RNDN);
        y = m * m;
        mpfr_remainder_ui(y, y, 2, MPFR_RNDN);
        mpfr_mul_ui(y, y, 100000000, MPFR_RNDN);
        /* pi^2 == 9.86960440108935861869 */
        ASSERT_ALWAYS(y >= -13039560);
        ASSERT_ALWAYS(y < -13039559);
        y = m * m * m;
        mpfr_remainder_si(y, y, -3, MPFR_RNDN);
        /* pi^3 == 31.00627668029982017479 */
        mpfr_mul_ui(y, y, 100000000, MPFR_RNDN);
        ASSERT_ALWAYS(y >= 100627668);
        ASSERT_ALWAYS(y < 100627669);
    }

    /* What do we want to do with
     * mpfr_{set,cmp,add,sub,mul,addmul,submul,div,remainder}_{int64,uin64}
     * ?
     * These functions are not reached on 64-bit machines.
     */

    {
        mpfr_t x, y;
        mpfr_init_set_uint64(x, 1234, MPFR_RNDN);
        mpfr_init_set_int64(y, -5678, MPFR_RNDN);
        mpfr_add_int64(x, y, -2345, MPFR_RNDN);
        mpfr_add_uint64(y, x, 3456, MPFR_RNDN);
        mpfr_mul_uint64(x, x, 14, MPFR_RNDN);
        mpfr_mul_int64(y, y, -23, MPFR_RNDN);
        mpfr_addmul_int64(y, x, -23, MPFR_RNDN);
        mpfr_addmul_uint64(x, y, 12, MPFR_RNDN);
        mpfr_submul_int64(y, x, -23, MPFR_RNDN);
        mpfr_submul_uint64(x, y, 12, MPFR_RNDN);
        mpfr_div_uint64(x, x, 11, MPFR_RNDN);
        mpfr_remainder_uint64(y, y, 7, MPFR_RNDN);
        mpfr_div_int64(x, x, 11, MPFR_RNDN);
        mpfr_remainder_int64(y, y, 7, MPFR_RNDN);

        // fmt::print("x = {}\n", cxx_mpfr(x));
        // fmt::print("y = {}\n", cxx_mpfr(y));
        ASSERT_ALWAYS(mpfr_cmp_int64(x, -73332627) < 0);
        ASSERT_ALWAYS(mpfr_cmp_int64(x, -73332628) > 0);
        mpfr_set_int64(x, -100, MPFR_RNDN);
        mpfr_sub_int64(x, x, -200, MPFR_RNDN);
        mpfr_add(y, y, x, MPFR_RNDN);
        mpfr_set_uint64(x, 1000, MPFR_RNDN);
        mpfr_sub_uint64(x, x, 200, MPFR_RNDN);
        mpfr_add(y, y, x, MPFR_RNDN);
        ASSERT_ALWAYS(mpfr_cmp_uint64(y, 903) == 0);
        mpfr_clear(x);
        mpfr_clear(y);
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
    test_mpfr_aux(prec);
}
