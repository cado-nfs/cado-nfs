#include "cado.h" // IWYU pragma: keep

#include <cstdlib>
#include <cmath>

#include <iostream>
#include <sstream>
#include <algorithm>
#include <map>
#include <vector>

#include <gmp.h>
#include "fmt/base.h"

#include "cxx_mpz.hpp"
#include "tests_common.h"
#include "macros.h"
#include "mpz_poly.h"
#include "mpz_poly_bivariate.hpp"
#include "arithmetic_reductions.hpp"
#include "timing.h"

static void test_mpz_poly_bivariate_trivialities(unsigned long iter)
{
    using T = cxx_mpz_poly_bivariate;

    for(unsigned long i = 0 ; i < iter ; i++) {
        T a, ta, t;
        int const d = gmp_urandomm_ui(state, 256);
        a.set_yi(d);
        T::transpose(ta, a);
        t.set_xi(d);
        ASSERT_ALWAYS(a.degree_x() == 0);
        ASSERT_ALWAYS(a.degree_y() == d);
        ASSERT_ALWAYS(ta.degree_x() == d);
        ASSERT_ALWAYS(ta.degree_y() == 0);
        ASSERT_ALWAYS(t == ta);

        T u, ua;
        u.swap(t);
        ta.swap(ua);
        ASSERT_ALWAYS(u == ua);
    }

    for(unsigned long i = 0 ; i < iter ; i++) {
        cxx_mpz_poly a;
        cxx_mpz_poly_bivariate Ax, Ay, t;
        int const dx = gmp_urandomm_ui(state, 32);
        int const nbits = 1 + gmp_urandomm_ui(state, 255);
        mpz_poly_set_randomb(a, dx, state, nbits,
                MPZ_POLY_UNSIGNED_COEFFICIENTS |
                MPZ_POLY_RRANDOM |
                MPZ_POLY_DEGREE_EXACT);
        Ax = T::lifted_x(a);
        Ay = T::lifted_y(a);
        T::transpose(t, Ax);
        ASSERT_ALWAYS(t == Ay);

        /* default behaviour is to lift x in Z[x] to x in Z[x,y] */
        t = a;
        ASSERT_ALWAYS(t == Ax);
    }
}

static void test_mpz_poly_bivariate_parsing()
{
    using T = cxx_mpz_poly_bivariate;

    char const * example1 = "(11+17*x+42*x^2)+(12+18*x+43*x^2)*y^2";
    char const * example2 =
        "(11+17*x+42*x^2)+(12+18*x+43*x^2)*y^2+(12+18*x+x^5)*y^5";
    char const * example3 = "(11+17*x+42*x^2)-(12+18*x+43*x^2)*y^2";
    {
        T f, g;

        std::istringstream("(11+17*x+42*x^2)+(12+18*x+43*x^2)*y^2") >>
            f.named("x", "y");

        std::ostringstream os;
        os << f << "\n";

        std::istringstream(os.str()) >> g;
        ASSERT_ALWAYS(f == g);
    }

    {
        T f, g;

        std::istringstream(example2) >> f.named("x", "y");

        ASSERT_ALWAYS(f.degree_y() == 5);
        ASSERT_ALWAYS(f.degree_x() == 5);
        ASSERT_ALWAYS(f.degree_in_yi(1) == -1);
        ASSERT_ALWAYS(f.degree_in_yi(2) == 2);
        ASSERT_ALWAYS(f.degree_in_xi(2) == 2);
        ASSERT_ALWAYS(mpz_cmp_ui(mpz_poly_coeff_const(f[2], 1), 18) == 0);

        T::transpose(g, f);

        // read exactly the same string, but transposed.
        std::istringstream(example2) >> f.named("y", "x");

        ASSERT_ALWAYS(f == g);
    }

    {
        T f, g;
        cxx_mpz_poly F, A, B;

        std::istringstream(example3) >> f.named("x", "y");

        ASSERT_ALWAYS(f.degree_y() == 2);
        ASSERT_ALWAYS(f.degree_x() == 2);
        ASSERT_ALWAYS(mpz_cmp_ui(mpz_poly_coeff_const(f[0], 2), 42) == 0);

        std::istringstream("11+17*x+42*x^2") >> F;
        std::istringstream("1+x+x^2") >> A;
        /* f is really (F+A) * (1-y^2) - A */
        mpz_poly_add(F, F, A);

        /* multiply by 1-y^2 in two steps */
        std::istringstream("1-y") >> B.named("y");
        T::mul(g, T::lifted_x(F), T::lifted_y(B));

        /* use literal translation */
        std::istringstream("1+var") >> B.named("var");
        T::mul(g, g, T::lifted_y(B));
        T::sub(g, g, T::lifted_x(A));

        ASSERT_ALWAYS(f == g);

        /* evaluation at y=1 should be -A */
        T::eval_fy(B, f, cxx_mpz_poly(1));
        mpz_poly_add(B, B, A);
        ASSERT_ALWAYS(B == 0);

        /* evaluation at x=-1 should be 36-37*y^2 */
        T::eval_fx(B, f, cxx_mpz(-1));
        cxx_mpz z;
        mpz_poly_add_ui(B, B, 1);
        mpz_poly_content(z, B);
        ASSERT_ALWAYS(z == 37);
        mpz_poly_eval_ui(z, B, 1);
        ASSERT_ALWAYS(z == 0);
    }

    {
        T f, g;
        std::istringstream(example3) >> f.named("x", "y");

        T::pow(f, f, 3);
        std::ostringstream os;
        os << f.named("x", "y") << "\n";
        std::istringstream(os.str()) >> g;
        bool ok = (f == g);
        ASSERT_ALWAYS(ok);

        std::istringstream(
            "-79507*x^6*y^6 + 232974*x^6*y^4 - 227556*x^6*y^2 + 74088*x^6 - "
            "99846*x^5*y^6 + 289347*x^5*y^4 - 279468*x^5*y^2 + 89964*x^5 - "
            "108360*x^4*y^6 + 310821*x^4*y^4 - 297093*x^4*y^2 + 94626*x^4 - "
            "61560*x^3*y^6 + 174672*x^3*y^4 - 165156*x^3*y^2 + 52037*x^3 - "
            "30240*x^2*y^6 + 84924*x^2*y^4 - 79473*x^2*y^2 + 24783*x^2 - "
            "7776*x*y^6 + 21600*x*y^4 - 19998*x*y^2 + 6171*x - 1728*y^6 + "
            "4752*y^4 "
            "- 4356*y^2 + 1331") >>
            g.named("x", "y");
        ok = (f == g);
        ASSERT_ALWAYS(ok);
    }

    {
        T f;
        cxx_mpz_poly A, B;

        std::istringstream(example1) >> f.named("x", "y");

        T::eval_fy(A, f, 42);
        std::istringstream("75894*x^2 + 31769*x + 21179") >> B;
        ASSERT_ALWAYS(A == B);

        T::eval_fy(A, f, 0);
        std::istringstream("42*x^2 + 17*x + 11") >> B;
        ASSERT_ALWAYS(A == B);

        T::eval_fx(A, f, cxx_mpz(42));
        std::istringstream("76620*y^2 + 74813") >> B.named("y");
        ASSERT_ALWAYS(A == B);

        T::eval_fx(A, f, cxx_mpz(0));
        std::istringstream("12*y^2 + 11") >> B.named("y");
        ASSERT_ALWAYS(A == B);

        f = 0;
        T::eval_fy(A, f, 42);
        ASSERT_ALWAYS(A == 0);

        T::eval_fx(A, f, cxx_mpz(1728));
        ASSERT_ALWAYS(A == 0);
    }
}

static void test_mpz_poly_bivariate_basic_arithmetic(unsigned long iter)
{
    fmt::print("{}\n", __func__);
    using T = cxx_mpz_poly_bivariate;
    for(unsigned long i = 0 ; i < iter ; i++) {
        T a, b, a_minus_b, b_minus_a, t;
        /* make sure we don't build anything larger than 2048 bits */
        int const nbits = 1 + gmp_urandomm_ui(state, 255);
        int const dx = gmp_urandomm_ui(state, (mp_limb_t) sqrt(2048./nbits));
        int const dy = gmp_urandomm_ui(state, (mp_limb_t) sqrt(2048./nbits));
        T::set_rrandomb(a, dx, dy, nbits, state);
        T::set_rrandomb(b, dx, dy, nbits, state);
        T::sub(a_minus_b, a, b);
        T::sub(b_minus_a, b, a);
        T::neg(t, a_minus_b);
        ASSERT_ALWAYS(t == b_minus_a);

        T::add(t, a_minus_b, b_minus_a);
        ASSERT_ALWAYS(t == 0);
    }
    for(unsigned long i = 0 ; i < iter ; i++) {
        T a, b, ab;
        /* make sure we don't build anything larger than 1024 bits */
        int const nbits = 1 + gmp_urandomm_ui(state, 63);
        int const dx = gmp_urandomm_ui(state, (mp_limb_t) sqrt(1024./nbits));
        int const dy = gmp_urandomm_ui(state, (mp_limb_t) sqrt(1024./nbits));
        T::set_rrandomb(a, dx, dy, nbits, state);
        T::set_rrandomb(b, dx, dy, nbits, state);
        T::mul(ab, a, b);

        {
            T x;
            T bx,a_bx,ab_x;
            T::set_rrandomb(x, dx, dy, nbits, state);
            T::mul(bx, b, x);
            T::mul(a_bx, a, bx);
            T::mul(ab_x, ab, x);
            ASSERT_ALWAYS(ab_x == a_bx);
        }

        {
            cxx_mpz_poly x;
            T bx,a_bx,ab_x;
            mpz_poly_set_randomb(x, dx, state, nbits,
                    MPZ_POLY_UNSIGNED_COEFFICIENTS |
                    MPZ_POLY_RRANDOM |
                    MPZ_POLY_DEGREE_EXACT);
            T::mul(bx, b, x);
            T::mul(a_bx, a, bx);
            T::mul(ab_x, ab, x);
            ASSERT_ALWAYS(ab_x == a_bx);
        }

        {
            cxx_mpz x;
            T bx,a_bx,ab_x;
            mpz_rrandomb(x, state, nbits);
            T::mul(bx, b, x);
            T::mul(a_bx, a, bx);
            T::mul(ab_x, ab, x);
            ASSERT_ALWAYS(ab_x == a_bx);
        }
    }

    /* test consistency of derivative with the coefficient-wise derivative
     * of the transpose */
    for(unsigned long i = 0 ; i < iter ; i++) {
        T a, dya, ta, dxta, tdxta;
        /* make sure we don't build anything larger than 1024 bits */
        int const nbits = 1 + gmp_urandomm_ui(state, 63);
        int const dx = gmp_urandomm_ui(state, (mp_limb_t) sqrt(1024./nbits));
        int const dy = gmp_urandomm_ui(state, (mp_limb_t) sqrt(1024./nbits));
        T::set_rrandomb(a, dx, dy, nbits, state);
        T::derivative_y(dya, a);
        T::transpose(ta, a);
        dxta = 0;
        T y;
        y.set_yi(1);
        for(int i = ta.degree_y() ; i >= 0 ; i--) {
            cxx_mpz_poly c;
            mpz_poly_derivative(c, ta[i]);
            T::mul(dxta, dxta, y);
            T::add(dxta, dxta, T::lifted_x(c));
        }
        T::transpose(tdxta, dxta);
        ASSERT_ALWAYS(tdxta == dya);
    }
}

template<typename T> using void_t = void;
template <class T, class = void>
struct reduction_is_coefficient_wise {
    static const bool value = true;
};
template<class T>
struct reduction_is_coefficient_wise<T,void_t<decltype(T::fy)>> {
    static const bool value = false;
};

template<typename reducer>
static int test_mpz_poly_bivariate_reduction_operator(reducer const & R)
{
    using T = cxx_mpz_poly_bivariate;
    /* just test commutativity of multiplication and reduction */
    T a, b;
    /* make sure we don't build anything larger than 2048 bits */
    int const nbits = 1 + gmp_urandomm_ui(state, 255);
    int const dx = gmp_urandomm_ui(state, (mp_limb_t) sqrt(2048./nbits));
    int const dy = gmp_urandomm_ui(state, (mp_limb_t) sqrt(2048./nbits));
    T::set_rrandomb(a, dx, dy, nbits, state);
    T::set_rrandomb(b, dx, dy, nbits, state);

    if (reduction_is_coefficient_wise<reducer>::value) {
        /* we test that (a-(a mod R)) is in the ideal R. However, we do
         * so by testing this on evey coefficient in y^i, which won't be
         * the correct thing to do if we're also reducing modulo a
         * polynomial in y
         */
        T r, t;
        R(r, a);
        T::sub(t, a, r);
        for(auto const & c : (T const &) t) {
            cxx_mpz_poly h;
            R(h, c);
            ASSERT_ALWAYS(h == 0);
        }

    }
    {
        T c, t;

        T::mul(c, a, b);
        R(b, b);
        R(a, a);
        T::mul(t, a, b);
        R(t, t);
        R(c, c);
        if (t != c) {
            fmt::print(std::cerr, "a = {}\n", a);
            fmt::print(std::cerr, "b = {}\n", b);
            fmt::print(std::cerr, "t = {}\n", t);
            fmt::print(std::cerr, "c = {}\n", c);
            std::cerr << R.print() << std::endl;
            return 0;
        }
        ASSERT_ALWAYS(t == c);
    }
    return 1;
}

static void test_mpz_poly_bivariate_reduction_functions(unsigned long iter)
{
    using cado::arithmetic_reductions::noop;
    using cado::arithmetic_reductions::mod_p;
    using cado::arithmetic_reductions::mod_fx;
    using cado::arithmetic_reductions::mod_q;
    using cado::arithmetic_reductions::mod_fy_mod_q;

    fmt::print("{}\n", __func__);
    using T = cxx_mpz_poly_bivariate;
    for(unsigned long i = 0 ; i < iter ; i++) {
        test_mpz_poly_bivariate_reduction_operator(noop{});
    }
    /* goal: never reduce mod something larger than 512 bits.
     */
    for(unsigned long i = 0 ; i < iter ; i++) {
        cxx_mpz p;
        const int nbits = 2 + gmp_urandomm_ui(state, 254);
        for( ; p == 0 || !mpz_probab_prime_p(p, 10) ; ) {
            mpz_urandomb(p, state, 4 + nbits / 2);
        }
        const mod_p Rp {p};
        ASSERT_ALWAYS(test_mpz_poly_bivariate_reduction_operator(Rp));

        /* size of fx is dfx * nbits (leading 1 does not count) */
        cxx_mpz_poly fx;
        const int dfx = 1 + std::min((mp_limb_t) (512 / nbits), gmp_urandomm_ui(state, 15));
        mpz_poly_set_randomb(fx, dfx, state, nbits,
                MPZ_POLY_UNSIGNED_COEFFICIENTS |
                MPZ_POLY_RRANDOM |
                MPZ_POLY_MONIC |
                MPZ_POLY_DEGREE_EXACT);
        const mod_fx Rfx{fx};
        ASSERT_ALWAYS(test_mpz_poly_bivariate_reduction_operator(Rfx));

        mod_q Rpfx {Rp, fx};
        ASSERT_ALWAYS(test_mpz_poly_bivariate_reduction_operator(Rpfx));

        /* size of fy is dy * (dx + 1) * nbits */
        {
            T F;
            int const nbits = 4 + gmp_urandomm_ui(state, 32);
            int const dx = 1 + gmp_urandomm_ui(state, (mp_limb_t) sqrt(512.0/nbits));
            int const dy = 2 + gmp_urandomm_ui(state, (mp_limb_t) sqrt(512.0/nbits));
            T::set_rrandomb(F, dx, dy, nbits, state);
            F.setcoeff(dy, 1);
            const mod_fy_mod_q Rpfxfy{Rpfx, F};
            ASSERT_ALWAYS(test_mpz_poly_bivariate_reduction_operator(Rpfxfy));
        }
    }
}

/* resultant is also a good test of evaluation, so we'll leave it to that
 * test */
static void test_mpz_poly_bivariate_resultant(unsigned long iter)
{
    for (unsigned long i = 0; i < iter; i++) {
        int const d = 2 + (int) gmp_urandomm_ui(state, 5);

        cxx_mpz_poly_bivariate f, g;

        for (; f == g;) {
            /* We need two Cab curves, so that there's a single place at
             * infinity, which implies that both resultants have the same
             * degree. */
            cxx_mpz_poly_bivariate::set_rrandomb_cab(f, d - 1, d, 3, state);
            cxx_mpz_poly_bivariate::set_rrandomb_cab(g, d - 1, d, 3, state);
        }

        cxx_mpz_poly rx, ry;

        cxx_mpz_poly_bivariate::resultant_x(rx, f, g);
        cxx_mpz_poly_bivariate::resultant_y(ry, f, g);

        int const ok = (rx != 0) && (ry != 0) &&
                       (mpz_poly_degree(rx) == mpz_poly_degree(ry));
        if (!ok) {
            std::cerr << "// bug in resultant" << "\n";
            std::cerr << "f = " << f << "\n";
            std::cerr << "g = " << g << "\n";
            std::cerr << "rx = " << rx << "\n";
            std::cerr << "ry = " << ry << "\n";
        }
        ASSERT_ALWAYS(ok);
    }

    /* This example is here so that the interpolation must be careful in
     * the choice of interpolation points, because otherwise we get
     * nonsense
     */
    {
        cxx_mpz_poly_bivariate f, g;
        cxx_mpz_poly tmp, rx, ry;
        int ok;

        std::istringstream("(x-1)*y+1") >> f.named("x", "y");
        std::istringstream("x*y^3+x^2+1") >> g.named("x", "y");

        cxx_mpz_poly_bivariate::resultant_y(ry, f, g);
        std::istringstream("x^5 - 3*x^4 + 4*x^3 - 4*x^2 + 2*x - 1") >>
            tmp.named("x");
        ok = (tmp == ry);
        if (!ok) {
            std::cerr << "// bug in resultant" << "\n";
            std::cerr << "f = " << f << "\n";
            std::cerr << "g = " << g << "\n";
            std::cerr << "ry = " << ry << "\n";
            std::cerr << "expected = " << tmp << "\n";
        }
        ASSERT_ALWAYS(ok);

        cxx_mpz_poly_bivariate::resultant_x(rx, f, g);
        std::istringstream("y^5 - y^4 + 2*y^2 - 2*y + 1") >> tmp.named("y");
        ok = (tmp == rx);
        if (!ok) {
            std::cerr << "// bug in resultant" << "\n";
            std::cerr << "f = " << f << "\n";
            std::cerr << "g = " << g << "\n";
            std::cerr << "rx = " << rx << "\n";
            std::cerr << "expected = " << tmp << "\n";
        }
        ASSERT_ALWAYS(ok);
    }

    {
        cxx_mpz_poly_bivariate f, g;
        cxx_mpz_poly tmp, rx, ry;
        int ok;

        std::istringstream("(11+17*x+42*x^2)+(12+18*x+42*x^2)*y^2") >>
            f.named("x", "y");
        std::istringstream("(10+15*x+100*x^2)+(20+41*x+17*x^2)*y^2") >>
            g.named("x", "y");

        cxx_mpz_poly_bivariate::resultant_y(ry, f, g);
        std::istringstream("12152196*x^8 + 2921268*x^7 + 1332913*x^6 - "
                           "2865824*x^5 - 1030822*x^4 "
                           "- 226892*x^3 + 152561*x^2 + 86200*x + 10000") >>
            tmp.named("x");
        ok = (tmp == ry);
        if (!ok) {
            std::cerr << "// bug in resultant" << "\n";
            std::cerr << "f = " << f << "\n";
            std::cerr << "g = " << g << "\n";
            std::cerr << "ry = " << ry << "\n";
            std::cerr << "expected = " << tmp << "\n";
        }
        ASSERT_ALWAYS(ok);

        cxx_mpz_poly_bivariate::resultant_x(rx, f, g);
        std::istringstream(
            "591408*y^8 + 30348*y^6 - 967958*y^4 + 52635*y^2 + 467750") >>
            tmp.named("y");
        ok = (tmp == rx);
        if (!ok) {
            std::cerr << "// bug in resultant" << "\n";
            std::cerr << "f = " << f << "\n";
            std::cerr << "g = " << g << "\n";
            std::cerr << "rx = " << rx << "\n";
            std::cerr << "expected = " << tmp << "\n";
        }
        ASSERT_ALWAYS(ok);
    }

    {
        cxx_mpz_poly_bivariate f, g;
        cxx_mpz_poly rx, ry;
        int ok;

        std::istringstream(
            "(467750+11*x+38960*x^2+17*x^3+42*x^4+44919*x^6+622660*x^8) "
            "+(12+18*x+43*x^2)*y^2 "
            "+(100+11*x+46*x^2+17*x^3+42*x^4+43*x^6+622660*x^8)*y^4 "
            "+(12+11*x+18*x^2+17*x^3+42*x^4+43*x^6+622660*x^8)*y^5") >>
            f.named("x", "y");

        std::istringstream(
            "(10+15*x+100*x^2) "
            "+(10+15*x+100*x^2+17*x^3+42*x^4+43*x^6+622660*x^8)*y "
            "+(20+41*x+17*x^2)*y^2 "
            "+(20+41*x+17*x^2+17*x^3+42*x^4+43*x^6+622660*x^8)*y^4") >>
            g.named("x", "y");

        cxx_mpz_poly_bivariate::resultant_y(ry, f, g);
        ok = (mpz_poly_degree(ry) == 72);
        ok = ok &&
             ry.coeff(42) ==
                 "1494774748820383618990632550811343575455342810669737"_mpz;
        if (!ok) {
            std::cerr << "// bug in resultant" << "\n";
            std::cerr << "f = " << f << "\n";
            std::cerr << "g = " << g << "\n";
            std::cerr << "ry = " << ry << "\n";
        }
        ASSERT_ALWAYS(ok);

        cxx_mpz_poly_bivariate::resultant_x(rx, f, g);
        ok = (mpz_poly_degree(rx) == 72);
        ok =
            ok &&
            rx.coeff(42) ==
                "392307271360601863563100893262849685478724955452018428897254854297277064613707763200"_mpz;
        if (!ok) {
            std::cerr << "// bug in resultant" << "\n";
            std::cerr << "f = " << f << "\n";
            std::cerr << "g = " << g << "\n";
            std::cerr << "rx = " << rx << "\n";
        }
        ASSERT_ALWAYS(ok);
    }
}

/* in reality, most of the complex stuff here should go in a reduction
 * object for the univariate mpz_poly type instead
 */
static void test_mpz_poly_bivariate_frobenius(unsigned long iter)
{
    using cado::arithmetic_reductions::mod_q;

    fmt::print("{}\n", __func__);
    for(unsigned long i = 0 ; i < iter ; i++) {
        cxx_mpz p;
        for( ; p == 0 || !mpz_probab_prime_p(p, 10) ; ) {
            mpz_urandomb(p, state, 2 + gmp_urandomm_ui(state, 8));
        }
        int const d = 1 + gmp_urandomm_ui(state, 5);
        cxx_mpz_poly f;
        for( ; f == 0 || !mpz_poly_is_irreducible(f, p) ; ) {
            mpz_poly_set_randomm(f, d, state, p,
                    MPZ_POLY_DEGREE_EXACT |
                    MPZ_POLY_MONIC |
                    MPZ_POLY_URANDOM);
        }

        const mod_q R { f, p };

        /* Tests on simple mpz_poly's -- this should probably go
         * elsewhere, really. */

        /* make sure that Frobenius and inverse Frobenius are linear */
        for(auto order: { 1, -1 }) {
            cxx_mpz_poly a, b, c, t;
            mpz_poly_set_randomm(a, d, state, p,
                    MPZ_POLY_DEGREE_UPPER_BOUND |
                    MPZ_POLY_URANDOM);
            mpz_poly_set_randomm(b, d, state, p,
                    MPZ_POLY_DEGREE_UPPER_BOUND |
                    MPZ_POLY_URANDOM);
            mpz_poly_add(c, a, b);
            R.frobenius(c, c, order);
            R.frobenius(a, a, order);
            R.frobenius(b, b, order);
            mpz_poly_add(t, a, b);
            R(t, t);
            if (t != c) {
                fmt::print(std::cerr, "a = {}\n", a);
                fmt::print(std::cerr, "b = {}\n", b);
                fmt::print(std::cerr, "t = {}\n", t);
                fmt::print(std::cerr, "c = {}\n", c);
                std::cerr << R.print() << "\n";
            }
            ASSERT_ALWAYS(c == t);
        }

        /* compatibility with scalar multiplication */
        for(auto order: { 1, -1 }) {
            cxx_mpz_poly a, c, t;
            cxx_mpz lambda;
            mpz_poly_set_randomm(a, d-1, state, p,
                    MPZ_POLY_DEGREE_UPPER_BOUND |
                    MPZ_POLY_URANDOM);
            mpz_urandomm(lambda, state, p);
            mpz_poly_mul_mpz(c, a, lambda);
            mpz_poly_mod_mpz(c, c, p, nullptr);
            R.frobenius(c, c, order);
            R.frobenius(a, a, order);
            mpz_poly_mul_mpz(t, a, lambda);
            R(t, t);
            ASSERT_ALWAYS(c == t);
        }

        /* make sure that both operators have order d */
        for(auto order: { 1, -1 }) {
            cxx_mpz_poly a, t;
            mpz_poly_set_randomm(a, d-1, state, p,
                    MPZ_POLY_DEGREE_UPPER_BOUND |
                    MPZ_POLY_URANDOM);
            R(t, a);
            for(int i = 0 ; i < d ; i++) {
                R.frobenius(a, a, order);
            }
            if (t != a) {
                fmt::print(std::cerr, "a = {}\n", a);
                fmt::print(std::cerr, "t = {}\n", t);
                std::cerr << R.print() << "\n";
            }
            ASSERT_ALWAYS(a == t);
        }

        /* Also do some tests on bivariate polynomials. There's not much
         * interesting to test, really. */
        {
            cxx_mpz_poly_bivariate a, b, t;
            cxx_mpz_poly_bivariate::set_urandomm (a, d-1, d, p, state, false, false);
            R.frobenius(b, a, 1);
            R.frobenius(t, b, -1);
            ASSERT_ALWAYS(a == t);
            R.frobenius(t, t, 1);
            ASSERT_ALWAYS(b == t);
        }
    }
}

static void test_mpz_poly_bivariate_factoring(unsigned long iter)
{
    using cado::arithmetic_reductions::mod_q;

    fmt::print("{}\n", __func__);
    for(unsigned long i = 0 ; i < iter ; i++) {
        double t0;

        t0 = -wct_seconds();
        cxx_mpz p;
        for( ; p == 0 || !mpz_probab_prime_p(p, 2) ; ) {
            mpz_urandomb(p, state, 2 + gmp_urandomm_ui(state, 6));
        }
        const int d = 2 + gmp_urandomm_ui(state, 4);
        cxx_mpz_poly f;
        for( ; f == 0 || !mpz_poly_is_irreducible(f, p) ; ) {
            mpz_poly_set_randomm(f, d, state, p,
                    MPZ_POLY_URANDOM |
                    MPZ_POLY_DEGREE_EXACT |
                    MPZ_POLY_MONIC);
        }

        const mod_q R { f, p };

        t0 += wct_seconds();

        fmt::print("{} {} {}\n", i, p, f.degree());
        fmt::print("  setup {:.3f}\n", t0);

        t0 = -wct_seconds();
        /* test square free factorization */
        {
            cxx_mpz_poly_bivariate a, b, t;
            
            /* a is chosen monic */
            cxx_mpz_poly_bivariate::set_urandomm (a, d-1, 10 + gmp_urandomm_ui(state, 32), p, state, false, true);

            auto fl = cxx_mpz_poly_bivariate::factor_sqf(a, R);

            /* all factors must be coprime */
            for(unsigned int i = 0 ; i < fl.size() ; i++) {
                for(unsigned int j = i + 1 ; j < fl.size() ; j++) {
                    cxx_mpz_poly_bivariate::gcd(t, fl[i].first, fl[j].first, R);
                    ASSERT_ALWAYS(t == 1);
                }
            }
            
            /* their product must match a */
            b = 1;
            for(unsigned int i = 0 ; i < fl.size() ; i++) {
                cxx_mpz_poly_bivariate::pow(t, fl[i].first, i);
                cxx_mpz_poly_bivariate::mul(b, b, t);
                cxx_mpz_poly_bivariate::mod(b, b, R);
            }
            ASSERT_ALWAYS(a == b);
        }
        t0 += wct_seconds();
        fmt::print("  sqf {:.3f}\n", t0);

        /* test DDF
         * We'll fabricate irreducible polynomials of degrees between 1
         * and 5, multiply them, and see what comes out of it.
         */
        {
            t0 = -wct_seconds();
            int nfac = 1 + gmp_urandomm_ui(state, 4);
            std::vector<cxx_mpz_poly_bivariate> known;
            std::map<int, int> expected;
            cxx_mpz_poly_bivariate A = 1;
            int tests = 0;
            for(int i = 0 ; i < nfac ; i++) {
                int const dg = 1 + gmp_urandomm_ui(state, 6);
                cxx_mpz_poly_bivariate g;
                do {
                    /* do not impose exactly degree d-1 in x */
                    cxx_mpz_poly_bivariate::set_urandomm (g, d-1, dg, p, state, false, true);
                    tests++;
                } while (!cxx_mpz_poly_bivariate::is_irreducible(g, R));
                known.push_back(g);
                expected[dg]++;
                cxx_mpz_poly_bivariate::mul(A, A, g);
                R(A, A);
            }
            /* sort the known polynomials, while we're at it. */
            std::ranges::sort(known);
            t0 += wct_seconds();
            fmt::print("  ddf setup({}, {}): {:.3f}\n", nfac, tests, t0);

            {
            t0 = -wct_seconds();
            std::map<int, int> got;
            auto fl0 = cxx_mpz_poly_bivariate::factor_sqf(A, R);
            t0 += wct_seconds();
            fmt::print("  ddf compute({}): {:.3f}\n", A.degree(), t0);
            t0 = -wct_seconds();
            for(size_t i = 1 ; i < fl0.size() ; i++) {
                auto fl1 = cxx_mpz_poly_bivariate::factor_ddf(fl0[i].first, R);
                for(size_t j = 1 ; j < fl1.size() ; j++) {
                    /* fl1[j].first is a product of irreducible factors
                     * of degree exactly j, which we're getting i times
                     */
                    if (fl1[j].first.degree() == 0)
                        continue;
                    got[j] += (fl1[j].first.degree() / j) * i;
                }
            }
            ASSERT_ALWAYS(expected == got);
            t0 += wct_seconds();
            fmt::print("  ddf compute({}): {:.3f}\n", A.degree(), t0);
            }
        }
    }
}

int main(int argc, char const * argv[])
{
    unsigned long iter = 500;
    tests_common_cmdline(&argc, &argv, PARSE_SEED | PARSE_ITER);
    tests_common_get_iter(&iter);
    test_mpz_poly_bivariate_trivialities(iter);
    test_mpz_poly_bivariate_parsing();
    test_mpz_poly_bivariate_basic_arithmetic(iter);
    test_mpz_poly_bivariate_reduction_functions(iter);
    test_mpz_poly_bivariate_resultant(iter);
    test_mpz_poly_bivariate_frobenius(iter);
    test_mpz_poly_bivariate_factoring(iter);
    tests_common_clear ();
    exit (EXIT_SUCCESS);
}
