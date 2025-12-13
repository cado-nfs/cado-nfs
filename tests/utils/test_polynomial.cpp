#include "cado.h" // IWYU pragma: keep

#include <cmath>
#include <cstdlib>

#include <array>
#include <complex>
#include <vector>
#include <string>
#include <limits>
#include <type_traits>
#include <tuple>

#include "fmt/base.h"
#include "fmt/std.h"    // formatter for std::complex

#include "cado_math_aux.hpp"
#include "cxx_mpz.hpp"
#include "macros.h"
#include "mpz_poly.h"
#include "polynomial.hpp"
#include "tests_common.h"
#include "number_context.hpp"
#include "extra_complex_overloads.hpp"

#ifdef HAVE_MPFR
#include "cxx_mpfr.hpp"
#endif
#ifdef HAVE_MPC
#include "cxx_mpc.hpp"
#endif

namespace cado_math_aux {
    template<typename T>
        static int accurate_bits(T reference, T computed)
        requires std::is_floating_point_v<T>
        {
            T c;
            if (reference == 0)
                c = computed;
            else
                c = (computed-reference)/reference;
            c = std::fabs(c);
            return c == 0 ? INT_MAX : -std::ilogb(c);
        }

#ifdef HAVE_MPFR
    template<typename T>
        static int
        accurate_bits(T const & reference, T const & computed)
        requires std::is_same_v<T, cxx_mpfr>
        {
            T c;
            if (reference == 0)
                c = computed;
            else
                c = (computed-reference)/reference;
            mpfr_abs(c, c, MPFR_RNDN);
            if (c == 0) return INT_MAX;
            T mantissa;
            mpfr_exp_t e;
            mpfr_frexp(&e, c, c, MPFR_RNDN);
            return -e;
        }
#endif
} /* namespace cado_math_aux */


/* This verifies that the provided roots are correct to the given
 * accuracy in bits (if positive) or have an accuracy loss wrt the given
 * type that is at most the (negative) given accuracy
 */
template<typename T>
static void
compare_roots(polynomial<T> const & q, std::vector<T> const & roots, std::vector<T> const & reference, int accuracy)
{
    if (roots.size() != reference.size()) {
        fmt::print(stderr,
                "compute_roots produced wrong number of roots {},"
                " reference has {}\n",
                roots.size(), reference.size());
        fmt::print(stderr, "f = {}\n", q);
        abort();
    }
    for(size_t i = 0 ; i < roots.size() ; i++) {
        using cado_math_aux::accurate_bits;
        int a = accurate_bits(roots[i], reference[i]);
        fmt::print("{} vs ref {} : accurate bits: {}\n", roots[i], reference[i], a);
        if (a < accuracy) {
            fmt::print(stderr,
                    "compute_roots produced wrong root {},"
                    " reference has {}\n",
                    roots[i], reference[i]);
            abort();
        }
    }
}

template<typename T>
static void 
test_positive_roots(std::string const & poly_str,
        bool verbose,
        cado::number_context<T> const & tr,
        std::vector<T> const & reference,
        int accuracy)
{
    polynomial<T> f(poly_str);

    if (accuracy < 0)
        accuracy += std::numeric_limits<T>::digits;

    if (verbose)
        fmt::print("Testing polynomial {}\n", f);

    /* test f as well as a few usual transforms */
    compare_roots(f, f.positive_roots(tr), reference, accuracy);
}

template<typename T>
static void
test_all_roots(std::string const & poly_str,
        bool verbose,
        cado::number_context<T> const & tr,
        std::vector<T> const & reference,
        int accuracy)
{
    polynomial<T> f(poly_str);

    if (accuracy < 0)
        accuracy += std::numeric_limits<T>::digits;

    if (verbose)
        fmt::print("Testing polynomial {}\n", f);

    compare_roots(f, f.roots(tr), reference, accuracy);
    {
        auto q = -f;
        compare_roots(q, q.roots(tr), reference, accuracy);
    }

    {
        auto q = f.mirror();
        std::vector<T> roots;
        roots.reserve(reference.size());
        for(auto const & x : reference)
            roots.emplace_back(-x);
        std::ranges::sort(roots);
        compare_roots(q, q.roots(tr), roots, accuracy);
    }

    {
        auto q = f.reciprocal();
        std::vector<T> roots;
        roots.reserve(reference.size());
        for(auto const & x : reference)
            if (x)
                roots.emplace_back(tr(1) / x);
        std::ranges::sort(roots);
        compare_roots(q, q.roots(tr), roots, accuracy);
    }

    {
        auto q = f.inverse_scale(3);
        std::vector<T> roots;
        roots.reserve(reference.size());
        for(auto const & x : reference)
            roots.emplace_back(tr(3) * x);
        std::ranges::sort(roots);
        compare_roots(q, q.roots(tr), roots, accuracy);
    }
}

template<typename T>
static void
test_compute_roots(bool verbose)
    requires cado_math_aux::is_real_v<T>
{
    const cado::number_context<T> tr(128); /* 128 only for cxx_mpfr */
    const int valgrind_penalty = (std::is_same_v<T, long double> && tests_run_under_valgrind()) ? 16 : 0;
    test_positive_roots<T>("1", verbose, tr, {}, -4);

    /* A few roots of 2 */
    test_all_roots<T>("x-2", verbose, tr, {2}, -5 - valgrind_penalty);
    test_positive_roots<T>("x^2-2", verbose, tr, {T(1.4142135623730950488016887242096980786L)}, -5 - valgrind_penalty);
    test_all_roots<T>("x^3-2", verbose, tr, {T(1.2599210498948731647672106072782283506L)}, -5 - valgrind_penalty);
    test_positive_roots<T>("x^4-2", verbose, tr, {T(1.1892071150027210667174999705604759153L)}, -5 - valgrind_penalty);
    test_all_roots<T>("x^5-2", verbose, tr, {T(1.1486983549970350067986269467779275894L)}, -5 - valgrind_penalty);

    test_all_roots<T>("(x-1)*(x-2)", verbose, tr, {1, 2}, -8 - valgrind_penalty);
    test_all_roots<T>("(x-1)*(x-2)*(x-3)", verbose, tr, {1, 2, 3}, -8 - valgrind_penalty);
    test_all_roots<T>("(x-1)*(x-2)*(x-3)*(x-4)", verbose, tr, {1, 2, 3, 4}, -8 - valgrind_penalty);
    test_all_roots<T>("(x-1)*(x-2)*(x-3)*(x-4)*(x-5)", verbose, tr, {1, 2, 3, 4, 5}, -8 - valgrind_penalty);

    test_all_roots<T>("(x-1)*(x-2)*(x-3)*(x-4)*(x-5)", verbose, tr, {1, 2, 3, 4, 5}, -8 - valgrind_penalty);

    /* Let f(x+1/x) * x^6 == (x^13-1)/(x-1). Test both positive and
     * negative roots */
    test_positive_roots<T>("x^6 + x^5 - 5*x^4 - 4*x^3 + 6*x^2 + 3*x - 1",
            verbose, tr, {
            T(0.24107336051064610669813537490508716455L),
            T(1.1361294934623116050236151182550332491L),
            T(1.7709120513064197918007510440301977572L),
            },
            -6 - valgrind_penalty);
    test_all_roots<T>("x^6 + x^5 - 5*x^4 - 4*x^3 + 6*x^2 + 3*x - 1",
            verbose, tr, {
            T(-1.9418836348521040543139645525875784545L),
            T(-1.4970214963422021972692611994027027677L),
            T(-0.70920977408507125193927578520003694863L),
            T(0.24107336051064610669813537490508716455L),
            T(1.1361294934623116050236151182550332491L),
            T(1.7709120513064197918007510440301977572L),
            },
            -6 - valgrind_penalty);
    /* this is f(x-2). 6 real roots, 0 rational */
    test_all_roots<T>("x^6 - 11*x^5 + 45*x^4 - 84*x^3 + 70*x^2 - 21*x + 1",
            verbose, tr,
            {T(0.058116365147895945686035447412421545500L),
            T(0.50297850365779780273073880059729723231L),
            T(1.2907902259149287480607242147999630514L),
            T(2.2410733605106461066981353749050871645L),
            T(3.1361294934623116050236151182550332491L),
            T(3.7709120513064197918007510440301977572L)},
            -8 - valgrind_penalty);

    test_all_roots<T>("4*x^3-1818*x", verbose, tr,
            {tr("-21.319005605327843245456394231442124020"),
             tr(0),
             tr("21.319005605327843245456394231442124020")},
             -8 - valgrind_penalty);

    /* examples below can't be dealt with by the float code because of
     * the limited exponent range
     */
    if (std::is_same_v<T, float>)
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
            verbose, tr,
            {T(0.46942663386512425278156827138724507L)},
            -4 - valgrind_penalty);

    /* false position produces b=NaN */
    test_positive_roots<T>(
            "-5.1229871591623088e+251"
            " + x * 4.8231399628079727e+240"
            " - x^2 * 7.5683722678735590e+228"
            " - x^3 * 1.8935837380523070e+224"
            " - x^4 * 3.4853123818766583e+152",
            verbose, tr,
            {}, -4 - valgrind_penalty);

    /* dichotomy with too few iterations fails */
    test_positive_roots<T>(
            "3396573496846254368813196771328*x^6"
            " + 17192931019341412634837118288585351300728750080*x^5"
            " - 2765086156017059372041966183496747450660061680253272064*x^4"
            " + 4974019329969663881845375223408004510305936664806589566321472066027520*x^2"
            " - 58290863912589135997939905826055669574526326226812326870330700385351723122688*x"
            " + 37718991555021708231785218373859670336563892134301066596876783017325698882362933248",
            verbose,
            tr,
            {
            T(687391.34048893181184293089912300885635L),
            T(1.1967277023040306827792365694733244330e7L),
            T(4.2033871205208627596333114764587438114e7L),
            T(1.4878252292723770602822450969598401002e8L)
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
            verbose, tr,
            {T(400274.63069456928417728406075344817128L)},
            -4 - valgrind_penalty);

    /* false position needs many iterations */
    test_positive_roots<T>(
            "416305583514625790805142742552103399071483895746667754265057484511150866432"
            " + x * 202399505763732049099628933992141454064175135567823362216516742332052144128"
            " - x^3 * 2677221347026437285957968988912544408687885411868999680"
            " + x^4 * 21555240319368651153052935288520704"
            " + x^5 * 5123362746908340224",
            verbose, tr,
            {
            T(8694859813.2710983618783204718642424277L),
            T(720776737597677797.82617434751880429252L) },
            -4 - valgrind_penalty);

    {
        polynomial<T> f("x^2-4");

        if (verbose)
            fmt::print("Testing polynomial {}\n", f);

        compare_roots(f, f.roots(1, tr), {}, 10);
        compare_roots(f, f.roots(3, tr), {-2, 2}, 10);
        compare_roots(f, f.positive_roots(1, tr), {}, 10);
        compare_roots(f, f.positive_roots(3, tr), {2}, 10);
    }
    {
        polynomial<T> f("(x^2-4)*x");

        if (verbose)
            fmt::print("Testing polynomial {}\n", f);

        compare_roots(f, f.roots(1, tr), {0}, 10);
        compare_roots(f, f.roots(3, tr), {-2,0, 2}, 10);
        compare_roots(f, f.positive_roots(1, tr), {}, 10);
        compare_roots(f, f.positive_roots(3, tr), {2}, 10);
    }
}

template<typename T>
static void
test_compute_roots(bool)
    requires cado_math_aux::is_complex_v<T>
{
    for(auto const & s : {"x^3-7*x^2+x+1",
            "(-0.335235471146870 - 0.980233680124672i)*x^5 + (-0.363312420177951 - 0.674311068694115i)*x^4 + (0.555635171051534 + 0.422342866615889i)*x^3 + (-0.119795855286804 + 0.223351587171732i)*x^2 + (0.894110547459500 - 0.368334794919981i)*x - 0.0709390404647527 - 0.447457524669756i",
            "320*x^3-1299080*x^2-218558114*x+99155309961",
            })
    {
        const cado::number_context<T> tr(512); /* 512 only for cxx_mpfr */
        polynomial<T> P(s, tr);
        {
            auto const lo = P.lower_bound_complex_roots();
            auto const hi = P.upper_bound_complex_roots();
            fmt::print("{} -> roots in [{}, {}] (surface={})\n",
                    P, lo, hi, hi*hi-lo*lo);
        }

        {
            auto [ mean, lo, hi ] = P.annulus_complex_roots();
            fmt::print("{} -> roots in [{}, {}, {}] (surface={})\n", P,
                    mean, lo, hi, hi*hi-lo*lo);
        }

        {
            auto zz = P.roots(P.ctx());
            for(auto const & z : zz)
                fmt::print("{}\n", z);
        }
    }
}





template<typename T>
static void
test_compute_roots(bool)
    requires std::is_same_v<T, cxx_mpz>
{
    /* this one is skipped */
}

template<typename T>
static void test_ctor_and_coeff_access()
{
    polynomial<T> f("42*x^2+17*x-1");
    ASSERT_ALWAYS(f.degree() == 2);
    ASSERT_ALWAYS(f[0] == -1);
    ASSERT_ALWAYS(f[1] == 17);
    ASSERT_ALWAYS(f[2] == 42);

    const polynomial<T> g({1, 2, 3, 4});
    ASSERT_ALWAYS(g.degree() == 3);
    ASSERT_ALWAYS(g(1) == 10);
}

template<typename T>
static void test_eval()
{
    ASSERT_ALWAYS(polynomial<T>()(0) == 0);
    ASSERT_ALWAYS(polynomial<int>()(T(0)) == T(0));
    /* long doubles under valgrind basically don't work. They're just
     * doubles in disguise
     */
    if (std::is_same_v<T, long double> && tests_run_under_valgrind())
        return;

    ASSERT_ALWAYS(polynomial<T>()(T()) == T());

    {
        polynomial<T> f;
        T w = 0;

        ASSERT_ALWAYS(f(w) == w);

        T two_n = 1;
        constexpr int D = std::numeric_limits<T>::digits;
        for(int n = 0 ; n < D ; n++) {
            /* f is 1 + 2 * x + ... + (n+1)*x^n */
            /* f(2) is g'(2) with g = 1 + x + ... + x^(n+1) =
             * (x^(n+2)-1)/(x-1), so g'(2) = (n+2)2^(n+1)-2^(n+2)+1 =
             * n * 2^(n+1) + 1. So we must have n*2^(n+1) < 2^digits in
             * order to avoid overflow, iow n >> (digits - (n+1)) == 0.
             * We must make sure that this right shift is less than the
             * type width, though. */

            const int shift = D - (n+1);
            if (shift <= std::numeric_limits<int>::digits)
                if (n >> (D - (n + 1)))
                    break;
            f[n] = n + 1;
            T v = f(2);
            w += two_n * T(f[n]);
            two_n *= 2;
            ASSERT_ALWAYS(v == w);
        }
    }

    /* also test homogeneous evaluation at (2,3) */
    {
        polynomial<T> f;
        T w = 0;
        T two_n = 1;
        constexpr int D = std::numeric_limits<T>::digits;
        for(int n = 0 ; n + 2 < D / 1.58497 ; n++) {
            /* following the same reasoning as above, f(2,3) is
             * g'(2/3)*3^n ==
             * ((n+2)(2/3)^(n+1)*(-1/3)-(2/3)^(n+2)+1)*9*3^n
             * == 3^(n+2)-(n+4)*2^(n+1)
             *
             * Since this value is odd, we only need to make sure that it
             * doesn't exceed 2^D, which boils down to forcing 3^(n+2) <
             * 2^D, so n+2 < D / (log(3)/log(2)), so for example n + 2 <
             * D/1.58497 should suffice.
             */
            f[n] = n + 1;
            T v = f(2, 3);
            T vx = f(T(2), T(3));
            w = 3 * w + two_n * f[n];
            two_n *= 2;
            ASSERT_ALWAYS(v == w);
            ASSERT_ALWAYS(vx == w);
        }
    }


    if constexpr (cado_math_aux::is_real_v<T>) {
        /* the precision 128 here is of course only ever considered by
         * cxx_mpfr. number_context<double> and friends proudly ignore
         * it.  */
        const cado::number_context<T> tr(128);
        const polynomial<int> A("(-x)^2-2");
        T z = tr("1.4142135623730950488016887242096980785696718754");

        auto a = A(z);
        fmt::print("sqrt(2)^2-2 == {:a}\n", a);

        for(auto const & r : A.roots(tr))
            fmt::print(" -> computed root {}\n", r);

        const polynomial<T> B(A, tr);
        auto b = B.eval_safe(z);
        fmt::print("sqrt(2)^2-2 == {:a} (safe eval)\n", b);

        int e;
        auto x = cado_math_aux::frexp(a-b, &e);

        fmt::print("difference = {} * 2^{}\n", x, e);

        /* not 100% sure this looks right in terms of precision. I have a
         * feeling that something odd's going on. Note that 0.5 can be
         * exactly represented by a legit value for x, so something like
         * an _exact_ 2^-21 can be 0.5 * 2^-20
         */
#ifdef HAVE_MPFR
        if constexpr (std::is_same_v<T, cxx_mpfr>)
            ASSERT_ALWAYS(x == 0 || e <= 3 - tr.prec);
        else
#endif
            ASSERT_ALWAYS(x == 0 || e <= 3 - std::numeric_limits<T>::digits);
    }
}

template<typename T>
static void test_derivative()
{
    using RX = polynomial<T>;
    ASSERT_ALWAYS(RX("17+42*x").derivative() == RX("42"));
    ASSERT_ALWAYS(RX("17+42*x+1728*x^2").derivative() == RX("42+3456*x"));
    ASSERT_ALWAYS(RX("17").derivative() == 0);
}

template<typename T>
static void test_reciprocal()
{
    using RX = polynomial<T>;

    ASSERT_ALWAYS(RX("1+2*x+3*x^2").reciprocal() == RX("3+2*x+x^2"));
    ASSERT_ALWAYS(RX("1+2*x").reciprocal() == RX("2+x"));
}

template<typename T>
static void test_print()
{
    using RX = polynomial<T>;

    /* the printing of these polynomial can vary depending on the
     * underlying type, so we'll just check that printing+parsing is the
     * identity map */
    const std::vector<const char *> tests {
        "17",
        "x * 42 + 17",
        "17+42*x+53*x^2",
        "1-x^2+99*x^3",
        "-x+x^2",
        "1-1",
        "1.7001e+12+53*x^2",
        "17-53.2*x^2+99*x^3",
    };
    for(auto const & t : tests) {
        const RX f(t);
        const RX g(fmt::format("{}", f));
        ASSERT_ALWAYS(f == g);
    }

}

template<>
void test_print<cxx_mpz>()
{
    using RX = polynomial<cxx_mpz>;

    const std::vector<std::array<const char *, 2>> tests {
        { "17", nullptr },
        { "x * 42 + 17", "17+42*x" },
        { "17+42*x+53*x^2", nullptr },
        { "1-x^2+99*x^3", nullptr },
        { "-x+x^2", nullptr },
        { "1-1", "0" },
    };

    for(auto const & t : tests) {
        const auto printed = RX(t[0]).print();
        const std::string ref(t[t[1] != nullptr]);
        ASSERT_ALWAYS(printed == ref);
    }
}

#if 0
template<>
void test_print<cxx_mpfr>()
{
    using RX = polynomial<cxx_mpfr>;

    const std::vector<const char *> tests {
        "17",
        "x * 42 + 17",
        "17+42*x+53*x^2",
        "1.7001e+09+53*x^2",
        "17-53.2*x^2+99*x^3",
        "1-x^2+99*x^3",
        "-x+x^2",
        "1-1",
    };

    /* just print */
    for(auto const & t : tests) {
        fmt::print("{}\n", RX(t));
    }
}
#endif

template<typename T>
static void test_ctor_mpz_poly()
{
    using RX = polynomial<T>;
    ASSERT_ALWAYS(RX(cxx_mpz_poly("17*x^2-42-3")) == RX("17*x^2-42-3"));
}

template<typename T>
static void test_resultant()
{
    /* only our integer resultant code makes any sense */
}


template<>
void test_resultant<cxx_mpz>()
{
    using T = cxx_mpz;
    using RX = polynomial<T>;

    const std::vector<std::tuple<std::string, std::string, T>> test_cases {
        {
            "x^6+13*x^5+13*x^4+9*x^3+7*x+6",
            "128*x^2+128*x+128",
            "162727720910848"_mpz
        },
        {
            "-3-15*x^1-9*x^2+3*x^3-12*x^4-12*x^5-3*x^6-3*x^7-12*x^8-15*x^9+6*x^10",
            "-6-13*x^1+9*x^2+7*x^3-5*x^4-5*x^5+11*x^6+2*x^7",
            "-61519394185549843500"_mpz
        },
        {
            "7917871+7917871*x-7916275*x^2-7916275*x^3-7916275*x^4+7917871*x^5+15834944*x^6",
            "128*x^2+128*x+128",
            "1102790158070603587092742144"_mpz
        },
        {
            "1365*x^6 + 1366*x^5+1368*x^4+1368*x^3+1368*x^2+1366*x+1366",
            "8320*x^2-50560*x-896",
            "37263864605996575174727132124282880"_mpz
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
            "11186466747618860118741892404376785"_mpz
        }};

    for(auto const & [f, g, ref] : test_cases) {
        T res = RX(f).resultant(RX(g));
        ASSERT_ALWAYS(res == ref);
    }
}

template<typename T>
static void test_some_arithmetic()
{
    const cado::number_context<T> tr(128); /* 128 only for cxx_mpfr */
    const polynomial<T> f("(x-1) * (x-2) * (x-3)", tr);
    const polynomial<T> g("(x-4) * (x-5) * (x-6)", tr);
    auto fg = f*g;
    ASSERT_ALWAYS(fg[0] == 720 && fg.degree() == 6);
    ASSERT_ALWAYS((f * polynomial<T>()) == 0);

    /* we don't divide with integer polynomials, since that would result
     * in euclidean division, and of course change the evaluation result.
     */
    if constexpr (cado_math_aux::is_real_v<T>)
        fg /= 4;

    for(const int i : {1,2,3,4,5,6})
        ASSERT_ALWAYS(fg(i) == 0);

    fg.addmul(g, polynomial<T>("x-1", tr));
    for(const int i : {1,4,5,6})
        ASSERT_ALWAYS(fg(i) == 0);
    fg.submul(g, polynomial<T>("x-37", tr));
    for(const int i : {4,5,6})
        ASSERT_ALWAYS(fg(i) == 0);
    fg.addmul(g, tr("1234"));
    for(const int i : {4,5,6})
        ASSERT_ALWAYS(fg(i) == 0);
    fg.submul(g, tr("-123"));
    for(const int i : {4,5,6})
        ASSERT_ALWAYS(fg(i) == 0);

    /* some gymnastics. Decompose a polynomial, and then parse the
     * decomposition to form a new polynomial.
     */
    auto q = f;
    for( ; q.degree() >= 0 ; ) {
        const int e = 10;
        auto [ nq, r ] = q.div_qr_xminusr(e);
        {
            auto s = q.div_q_xminusr(e);
            ASSERT_ALWAYS(s == nq);
        }
        {
            auto [ wq, wr ] = q.div_qr(polynomial<T> {-tr(e), tr(1)});
            ASSERT_ALWAYS(nq == wq);
            ASSERT_ALWAYS(wr.degree() <= 0);
            ASSERT_ALWAYS(r == wr[0]);
        }
        {
            auto wq = q.div_q(polynomial<T> {-tr(e), tr(1)});
            ASSERT_ALWAYS(nq == wq);
        }
        {
            auto wr = q.div_r(polynomial<T> {-tr(e), tr(1)});
            ASSERT_ALWAYS(wr.degree() <= 0);
            ASSERT_ALWAYS(r == wr[0]);
        }
        {
            auto s = fmt::format("({}) * (x-{}) + {}", nq, e, r);
            fmt::print("{} == {}\n", q, s);
            ASSERT_ALWAYS(q == polynomial<T>(s, q.ctx()));
        }
        {
            auto fs = q.shift(e);
            auto var = fmt::format("(x+{})", e);
            auto s = fmt::format("{}", q.named(var));
            fmt::print("{} == {}\n", s, fs);
            ASSERT_ALWAYS(fs == polynomial<T>(s, q.ctx()));
        }
        q = nq;
    }

    {
        /* some parsing */
        ASSERT_ALWAYS(polynomial<T>("-(x+1)^0+(x+1)^2-(-x)*-x-2*x") == 0);
    }

    {
        const polynomial<int> A { 1, 2 };
        const polynomial<T> B(A);
        ASSERT_ALWAYS(B.derivative() == 2);
        const cado::number_context<T> tr(128); /* 128 only for cxx_mpfr */
        polynomial<T> C(A, tr);
        ASSERT_ALWAYS(C.derivative() == 2);
        C.set_zero();
        ASSERT_ALWAYS(C == 0);
        ASSERT_ALWAYS(C.pow(12) == 0);
        ASSERT_ALWAYS(C.pow(0) == 1);

        const auto A7 = A.pow(7);
        ASSERT_ALWAYS(A7.degree() == 7);
        const polynomial<int> K(2187);
        ASSERT_ALWAYS(K == A7(1));

        auto [ q, r ] = C.div_qr_xminusr(tr(1));
        ASSERT_ALWAYS(q == 0);
        ASSERT_ALWAYS(r == 0);
    }
}

/* things from polynomial.hpp what are not yet covered (I think)
 *
 * findroot_dichotomy and the code path that leads to it (corner case in
 * findroot_falseposition). test_init_norms might trigger this, though.
 *
 * corner cases in pseudo_division, resultant, and inverse_scale
 *
 * and throwing branches.
 */

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
    test_some_arithmetic<T>();
}


int main()
{
    fmt::print("===== tests with T = double =====\n");
    all_tests<double>();
    fmt::print("===== tests with T = long double =====\n");
    all_tests<long double>();
    fmt::print("===== tests with T = float =====\n");
    all_tests<float>();
    fmt::print("===== tests with T = cxx_mpz =====\n");
    all_tests<cxx_mpz>();
#ifdef HAVE_MPFR
    fmt::print("===== tests with T = cxx_mpfr =====\n");
    all_tests<cxx_mpfr>();
#endif
    fmt::print("===== tests with T = std::complex<double> =====\n");
    all_tests<std::complex<double>>();
#ifdef HAVE_MPC
    fmt::print("===== tests with T = cxx_mpc =====\n");
    all_tests<cxx_mpc>();
#endif
    return EXIT_SUCCESS;
}
