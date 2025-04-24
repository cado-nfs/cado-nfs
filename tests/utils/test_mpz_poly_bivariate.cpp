#include "cado.h" // IWYU pragma: keep
#include <cstdlib>

#include <iostream>
#include <sstream>

#include <gmp.h>

#include "cxx_mpz.hpp"
#include "macros.h"
#include "mpz_poly.h"
#include "mpz_poly_bivariate.hpp"
#include "tests_common.h"

static void test_mpz_poly_bivariate_trivialities()
{
    typedef cxx_mpz_poly_bivariate T;

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

        T::pow_ui(f, f, 3);
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

int main(int argc, char const * argv[])
{
    unsigned long iter = 500;
    tests_common_cmdline(&argc, &argv, PARSE_SEED | PARSE_ITER);
    tests_common_get_iter(&iter);
    test_mpz_poly_bivariate_trivialities();
    test_mpz_poly_bivariate_resultant(iter);
    tests_common_clear();
    return EXIT_SUCCESS;
}
