#include "cado.h" // IWYU pragma: keep

#include <cmath>
#include <cstdlib>

#include <array>
#include <vector>
#include <string>
#include <limits>
#include <type_traits>

#include "fmt/base.h"

#include "cado_math_aux.hpp"
#include "cxx_mpz.hpp"
#include "macros.h"
#include "mpz_poly.h"
#include "polynomial.hpp"
#include "tests_common.h"

/* This verifies that the provided roots are correct to the given
 * accuracy in bits (if positive) or have an accuracy loss wrt the given
 * type that is at most the (negative) given accuracy
 */
template<typename T>
static void 
test_positive_roots(std::string const & poly_str,
        const T bound,
        bool verbose,
        std::vector<T> const & reference,
        int accuracy)
{
    polynomial<T> f(poly_str);

    if (accuracy < 0)
        accuracy += std::numeric_limits<T>::digits;

    if (verbose)
        fmt::print("Testng polynomial {}\n", f);

    if (!reference.empty()) {
        if (bound > 0) {
            const T B = f.bound_positive_roots();
            const T M = *std::max_element(reference.begin(), reference.end());
            ASSERT_ALWAYS(B >= M);
        } else {
            T B = f.bound_positive_roots(true);
            const T M = *std::min_element(reference.begin(), reference.end());
            ASSERT_ALWAYS(B <= M);
        }
    }

    auto pos = f.positive_roots(bound);

    if (pos.size() != reference.size()) {
        fmt::print(stderr,
                "compute_roots produced wrong number of roots {},"
                " reference has {}\n",
                pos.size(), reference.size());
        fmt::print(stderr, "f = {}\n", f);
        fmt::print(stderr, "bound = {}\n", bound);
        abort();
    }
    for(size_t i = 0 ; i < pos.size() ; i++) {
        using cado_math_aux::accurate_bits;
        int a = accurate_bits(pos[i], reference[i]);
        fmt::print("{} vs ref {} : accurate bits: {}\n", pos[i], reference[i], a);
        if (a < accuracy) {
            fmt::print(stderr,
                    "compute_roots produced wrong root {},"
                    " reference has {}\n",
                    pos[i], reference[i]);
            abort();
        }
    }
}

template<typename T>
static void
test_compute_roots(bool verbose)
{
    const int valgrind_penalty = (std::is_same<T, long double>::value && tests_run_under_valgrind()) ? 16 : 0;
    test_positive_roots<T>("1", 3, verbose, {}, -4);

    /* A few roots of 2 */
    test_positive_roots<T>("x-2", 3, verbose,   {2}, -5 - valgrind_penalty);
    test_positive_roots<T>("x^2-2", 3, verbose, {1.4142135623730950488016887242096980786L}, -5 - valgrind_penalty);
    test_positive_roots<T>("x^3-2", 3, verbose, {1.2599210498948731647672106072782283506L}, -5 - valgrind_penalty);
    test_positive_roots<T>("x^4-2", 3, verbose, {1.1892071150027210667174999705604759153L}, -5 - valgrind_penalty);
    test_positive_roots<T>("x^5-2", 3, verbose, {1.1486983549970350067986269467779275894L}, -5 - valgrind_penalty);

    test_positive_roots<T>("(x-1)*(x-2)", 3, verbose, {1, 2}, -8 - valgrind_penalty);
    test_positive_roots<T>("(x-1)*(x-2)*(x-3)", 4, verbose, {1, 2, 3}, -8 - valgrind_penalty);
    test_positive_roots<T>("(x-1)*(x-2)*(x-3)*(x-4)", 5, verbose, {1, 2, 3, 4}, -8 - valgrind_penalty);
    test_positive_roots<T>("(x-1)*(x-2)*(x-3)*(x-4)*(x-5)", 6, verbose, {1, 2, 3, 4, 5}, -8 - valgrind_penalty);

    /* Let f(x+1/x) * x^6 == (x^13-1)/(x-1). Test both positive and
     * negative roots (for negative roots we give a negative bound) */
    test_positive_roots<T>("x^6 + x^5 - 5*x^4 - 4*x^3 + 6*x^2 + 3*x - 1",
            2, verbose, {
            0.24107336051064610669813537490508716455L,
            1.1361294934623116050236151182550332491L,
            1.7709120513064197918007510440301977572L},
            -4 - valgrind_penalty);
    test_positive_roots<T>("x^6 + x^5 - 5*x^4 - 4*x^3 + 6*x^2 + 3*x - 1",
            -2, verbose, {
            -0.70920977408507125193927578520003694863L,
            -1.4970214963422021972692611994027027677L,
            -1.9418836348521040543139645525875784545L,
            },
            -6 - valgrind_penalty);
    /* this is f(x-2). 6 real roots, 0 rational */
    test_positive_roots<T>("x^6 - 11*x^5 + 45*x^4 - 84*x^3 + 70*x^2 - 21*x + 1",
            4, verbose,
            {0.058116365147895945686035447412421545500L,
            0.50297850365779780273073880059729723231L,
            1.2907902259149287480607242147999630514L,
            2.2410733605106461066981353749050871645L,
            3.1361294934623116050236151182550332491L,
            3.7709120513064197918007510440301977572L},
            -8 - valgrind_penalty);

    /* examples below can't be dealt with by the float code because of
     * the limited exponent range
     */
    if (std::is_same<T, float>::value)
        return;

    /* false position needs weighting */
    test_positive_roots<T>(
            "1.1771253911282693e+36"
            " + x * 1.0293658912886811e+39"
            " + x^2 * 1.9888797504712385e+38"
            " - x^3 * 9.7240762665983556e+38"
            " + x^4 * 3.5375219923207381e+37"
            " - x^5 * 1.8078625258002374e+40"
            " - x^6 * 1.6114410711830154e+39",
            1, verbose,
            {0.46942663386512425278156827138724507L},
            -4 - valgrind_penalty);

    /* false position produces b=NaN */
    test_positive_roots<T>(
            "-5.1229871591623088e+251"
            " + x * 4.8231399628079727e+240"
            " - x^2 * 7.5683722678735590e+228"
            " - x^3 * 1.8935837380523070e+224"
            " - x^4 * 3.4853123818766583e+152",
            1e72, verbose,
            {}, -4 - valgrind_penalty);

    /* dichotomy with too few iterations fails */
    test_positive_roots<T>(
            "3396573496846254368813196771328*x^6"
            " + 17192931019341412634837118288585351300728750080*x^5"
            " - 2765086156017059372041966183496747450660061680253272064*x^4"
            " + 4974019329969663881845375223408004510305936664806589566321472066027520*x^2"
            " - 58290863912589135997939905826055669574526326226812326870330700385351723122688*x"
            " + 37718991555021708231785218373859670336563892134301066596876783017325698882362933248",
            9007199254740992.0,
            verbose,
            {
            687391.34048893181184293089912300885635L,
            1.1967277023040306827792365694733244330e7L,
            4.2033871205208627596333114764587438114e7L,
            1.4878252292723770602822450969598401002e8L
            },
            -4 - valgrind_penalty);

    /* false position reaches s==a || s==b very early, need countermeasure */
    test_positive_roots<T>(
            "-4067876448477449548563432969586031412941700056587601428084520911569151852544"
            " + x * 299197339446978982940341333993590148850152155031542808049001636857517178880"
            " - x^2 * 88455602043598760130397057784414636527233402886601450561774258429246832640"
            " + x^4 * 552084484207525010843995104449133201610495729247299431552778240"
            " + x^5 * 386336020186016733435543279553138102808608768"
            " + x^6 * 38629366113825858026209280",
            2199023255552, verbose,
            {400274.63069456928417728406075344817128L},
            -4 - valgrind_penalty);

    /* false position needs many iterations */
    test_positive_roots<T>(
            "416305583514625790805142742552103399071483895746667754265057484511150866432"
            " + x * 202399505763732049099628933992141454064175135567823362216516742332052144128"
            " - x^3 * 2677221347026437285957968988912544408687885411868999680"
            " + x^4 * 21555240319368651153052935288520704"
            " + x^5 * 5123362746908340224",
            1e20, verbose,
            {
            8694859813.2710983618783204718642424277L,
            720776737597677797.82617434751880429252L },
            -4 - valgrind_penalty);
}

template<typename T>
static void test_ctor_and_coeff_access()
{
    polynomial<T> f("42*x^2+17*x-1");
    ASSERT_ALWAYS(f.degree() == 2);
    ASSERT_ALWAYS(f[0] == -1);
    ASSERT_ALWAYS(f[1] == 17);
    ASSERT_ALWAYS(f[2] == 42);
}

template<typename T>
static void test_eval()
{
    polynomial<T> f;

    T w = 0;
    for(int deg = 0 ; deg < std::numeric_limits<T>::digits ; deg++) {
        f[deg] = 1;
        T v = f.eval(2);
        w = 2 * w + f[deg];
        ASSERT_ALWAYS(v == w);
    }
}

template<typename T>
static void test_derivative()
{
    typedef polynomial<T> RX;
    ASSERT_ALWAYS(RX("17+42*x").derivative() == RX("42"));
    ASSERT_ALWAYS(RX("17+42*x+1728*x^2").derivative() == RX("42+3456*x"));
    ASSERT_ALWAYS(RX("17").derivative() == 0);
}

template<typename T>
static void test_reciprocal()
{
    typedef polynomial<T> RX;

    ASSERT_ALWAYS(RX("1+2*x+3*x^2").reciprocal() == RX("3+2*x+x^2"));
    ASSERT_ALWAYS(RX("1+2*x").reciprocal() == RX("2+x"));
}

template<typename T>
static void test_print()
{
    typedef polynomial<T> RX;

    const std::vector<std::array<const char *, 2>> tests {
        { "17", nullptr },
        { "x * 42 + 17", "17+42*x" },
        { "17+42*x+53*x^2", nullptr },
        { "1.7001e+09+53*x^2", nullptr },
        { "17-53.2*x^2+99*x^3", nullptr },
        { "1-x^2+99*x^3", nullptr },
        { "-x+x^2", nullptr },
        { "1-1", "0" },
    };

    for(auto const & t : tests)
        ASSERT_ALWAYS(RX(t[0]).print() == std::string(t[t[1] != nullptr]));
}

template<typename T>
static void test_ctor_mpz_poly()
{
    typedef polynomial<T> RX;
    ASSERT_ALWAYS(RX(cxx_mpz_poly("17*x^2-42-3")) == RX("17*x^2-42-3"));
}

template<typename T>
static void test_resultant()
{
    typedef polynomial<T> RX;
    using namespace cado_math_aux;

    struct test_case {
        std::string f, g;
        int ulps = 0;
        /* default member initializers prevent the struct from being an
         * aggregate until c++14. Oddly enough, it seems to still cause
         * trouble with icpx with c++20.
         * https://stackoverflow.com/questions/39344444/brace-aggregate-initialization-for-structs-with-default-values
         */
    };

    const int valgrind_penalty2 = (std::is_same<T, long double>::value && tests_run_under_valgrind()) ? 8 : 0;

    const std::vector<test_case> test_cases {
        {
            "x^6+13*x^5+13*x^4+9*x^3+7*x+6",
            "128*x^2+128*x+128",
            -8
        },
        {
            "-3-15*x^1-9*x^2+3*x^3-12*x^4-12*x^5-3*x^6-3*x^7-12*x^8-15*x^9+6*x^10",
            "-6-13*x^1+9*x^2+7*x^3-5*x^4-5*x^5+11*x^6+2*x^7",
            -12 - valgrind_penalty2
        },
        {
            "7917871+7917871*x-7916275*x^2-7916275*x^3-7916275*x^4+7917871*x^5+15834944*x^6",
            "128*x^2+128*x+128",
            -1
        },
        {
            "1365*x^6 + 1366*x^5+1368*x^4+1368*x^3+1368*x^2+1366*x+1366",
            "8320*x^2-50560*x-896",
            -12 - valgrind_penalty2
        },
        {
            /* this test case triggers a catastrophic cancellation in the
             * pseudo-division code in double precision, so we fall back
             * on an exact
             * computation instead. Hence in double precision, of course
             * we're accurate within 1 ulp. For long doubles, the
             * situation is not quite the same: we don't have the full
             * cancellation, but the intermediate data only has a few
             * correct bits. Eventually, what happens is that we have a
             * few correct bits, but not many (about 9).
             */
            "1365*x^6 + 1366*x^5+1368*x^4+1368*x^3+1368*x^2+1366*x+1366",
            "15*x^2-43368*x-4753",
            8
        }};

    for(auto const & t : test_cases) {
        T res = RX(t.f).resultant(RX(t.g));
        cxx_mpz val_z;
        cxx_mpz_poly fz(t.f);
        cxx_mpz_poly gz(t.g);
        mpz_poly_resultant(val_z, fz, gz);
        T val = cado_math_aux::mpz_get<T>(val_z);
        int accuracy = t.ulps;
        if (accuracy < 0)
            accuracy += std::numeric_limits<T>::digits;
        int a = accurate_bits(res, val);
        fmt::print("{} vs ref {} : accurate bits: {}\n", res, val, a);
        ASSERT_ALWAYS(a >= accuracy);
    }
}

template<typename T>
static void all_tests()
{
  test_ctor_and_coeff_access<T> ();
  test_eval<T> ();
  test_derivative<T> ();
  test_reciprocal<T> ();
  test_print<T> ();
  test_ctor_mpz_poly<T> ();
  test_resultant<T>();
  test_compute_roots<T>(false);
}


int main()
{
    all_tests<double>();
    all_tests<long double>();
    all_tests<float>();
    return EXIT_SUCCESS;
}
