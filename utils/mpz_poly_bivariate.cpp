#include "cado.h" // IWYU pragma: keep

#include <cstddef>
#include <cctype>

#include <sstream>
#include <istream>
#include <ostream>
#include <string>
#include <algorithm>
#include <vector>
#include <ios>

#include <gmp.h>

#include "cado_expression_parser.hpp"
#include "mpz_poly.h"
#include "cxx_mpz.hpp"
#include "mpz_poly_bivariate.hpp"
#include "runtime_numeric_cast.hpp"
#include "macros.h"
#include "named_proxy.hpp"

/* Polynomial arithmetic */
void cxx_mpz_poly_bivariate::neg(cxx_mpz_poly_bivariate & f,
                                 cxx_mpz_poly_bivariate const & g)
{
    for (size_t i = 0; i < g.size(); i++) {
        mpz_poly_neg(((super &)f)[i], ((super &)f)[i]);
    }
}

void cxx_mpz_poly_bivariate::add(cxx_mpz_poly_bivariate & f,
                                 cxx_mpz_poly_bivariate const & g,
                                 cxx_mpz_poly_bivariate const & h)
{
    const size_t sg = g.size();
    const size_t sh = h.size();
    const size_t sf = std::max(sg, sh);
    f.reserve(sf);
    size_t i = 0;
    if (&f != &g && &f != &h)
        f.clear();
    f.insert(f.end(), sf - f.size(), cxx_mpz_poly()); /* now sf == f.size() */
    for (; i < sg && i < sh; i++)
        mpz_poly_add(((super &)f)[i], g[i], h[i]);
    for (; i < sg; i++)
        mpz_poly_set(((super &)f)[i], g[i]);
    for (; i < sh; i++)
        mpz_poly_set(((super &)f)[i], h[i]);
    f.cleandeg((int)sf - 1);
}

void cxx_mpz_poly_bivariate::sub(cxx_mpz_poly_bivariate & f,
                                 cxx_mpz_poly_bivariate const & g,
                                 cxx_mpz_poly_bivariate const & h)
{
    const size_t sg = g.size();
    const size_t sh = h.size();
    const size_t sf = std::max(sg, sh);
    f.reserve(sf);
    size_t i = 0;
    if (&f != &g && &f != &h)
        f.clear();
    f.insert(f.end(), sf - f.size(), cxx_mpz_poly()); /* now sf == f.size() */
    for (; i < sg && i < sh; i++)
        mpz_poly_sub(((super &)f)[i], g[i], h[i]);
    for (; i < sg; i++)
        mpz_poly_set(((super &)f)[i], g[i]);
    for (; i < sh; i++)
        mpz_poly_neg(((super &)f)[i], h[i]);
    f.cleandeg((int)sf - 1);
}

/*
void cxx_mpz_poly_bivariate::add(cxx_mpz_poly_bivariate & f,
cxx_mpz_poly_bivariate const & g, cxx_mpz_poly_bivariate::lifted_x const & a)
{
    f = g;
    if (f == 0) {
        f = a;
}

*/

void cxx_mpz_poly_bivariate::mul(cxx_mpz_poly_bivariate & f,
                                 cxx_mpz_poly_bivariate const & g,
                                 cxx_mpz_poly_bivariate const & h)
{
    if (g.degree() == -1 || h.degree() == -1) {
        f.clear();
        return;
    }
    if (&f == &g || &f == &h) {
        cxx_mpz_poly_bivariate ff;
        mul(ff, g, h);
        f.swap(ff);
        return;
    }
    f.clear();
    // both f and g are non-zero
    cxx_mpz_poly tmp;
    f.assign(g.degree() + h.degree() + 1, tmp);
    for (int i = 0; i <= g.degree(); i++) {
        for (int j = 0; j <= h.degree(); j++) {
            mpz_poly_mul(tmp, g[i], h[j]);
            mpz_poly_add(((super &)f)[i + j], f[i + j], tmp);
        }
    }
    f.cleandeg(g.degree() + h.degree());
}

/* B = A^n */
void cxx_mpz_poly_bivariate::pow_ui(cxx_mpz_poly_bivariate & B,
                                    cxx_mpz_poly_bivariate const & A,
                                    unsigned long n) /*{{{*/
{
    if (n == 0) {
        B = 1;
        return;
    }
    if (A.degree() < 0) {
        B = 0;
        return;
    }
    if (&B == &A) {
        cxx_mpz_poly_bivariate C;
        pow_ui(C, A, n);
        B.swap(C);
        return;
    }

    unsigned long k = ((~0UL) >> 1) + 1;
    for (; k > n; k >>= 1)
        ;
    B = A;
    for (; k >>= 1;) {
        mul(B, B, B);
        if (n & k)
            mul(B, B, A);
    }
} /*}}}*/

struct mpz_poly_bivariate_parser_traits {
    std::string x, y;
    /*
    mpz_poly_bivariate_parser_traits(std::string const & x,
                                     std::string const & y)
        : x(x)
        , y(y)
    {
    }
    */

    struct unexpected_literal : public cado_expression_parser_details::parse_error {
        std::string msg;
        explicit unexpected_literal(std::string const & v)
            : msg(std::string("unexpected literal " + v))
        {
        }
        char const * what() const noexcept override { return msg.c_str(); }
    };

    static constexpr int const accept_literals = 2;
    typedef cxx_mpz_poly_bivariate type;
    typedef cxx_mpz number_type;

    static void add(cxx_mpz_poly_bivariate & c, cxx_mpz_poly_bivariate const & a,
             cxx_mpz_poly_bivariate const & b)
    {
        type::add(c, a, b);
    }
    static void sub(cxx_mpz_poly_bivariate & c, cxx_mpz_poly_bivariate const & a,
             cxx_mpz_poly_bivariate const & b)
    {
        type::sub(c, a, b);
    }
    static void mul(cxx_mpz_poly_bivariate & c, cxx_mpz_poly_bivariate const & a,
             cxx_mpz_poly_bivariate const & b)
    {
        type::mul(c, a, b);
    }
    static void pow_ui(cxx_mpz_poly_bivariate & c, cxx_mpz_poly_bivariate const & a,
                unsigned long e)
    {
        type::pow_ui(c, a, e);
    }
    static void neg(cxx_mpz_poly_bivariate & a, cxx_mpz_poly_bivariate & b)
    {
        type::neg(a, b);
    }
    static void swap(cxx_mpz_poly_bivariate & a, cxx_mpz_poly_bivariate & b)
    {
        a.swap(b);
    }
    static void set(cxx_mpz_poly_bivariate & a, cxx_mpz const & z)
    {
        a = z;
    }
    /* TODO do something for variable names */
    void set_literal_power(cxx_mpz_poly_bivariate & a, std::string const & v,
                           unsigned long e) const
    {
        if (v == x) {
            a.set_xi(e);
        } else if (v == y) {
            a.set_yi(e);
        } else {
            throw unexpected_literal(v);
        }
    }
};

std::istream &
cado::operator>>(std::istream & in,
           named_proxy<cxx_mpz_poly_bivariate &> F)
{
    std::string const & x = F.x();
    std::string const & y = F.y();
    cxx_mpz_poly_bivariate & f = F.c;

    f = 0;

    std::string line;
    for (;; in.get()) {
        const int c = in.peek();
        if (in.eof() || !isspace(c))
            break;
    }
    if (!getline(in, line))
        return in;
    std::istringstream is(line);

    using poly_parser = cado_expression_parser<mpz_poly_bivariate_parser_traits>
        ;
    poly_parser P(x, y);

    P.tokenize(is);

    try {
        f = P.parse();
    } catch (cado_expression_parser_details::parse_error const & p) {
        in.setstate(std::ios_base::failbit);
        return in;
    }

    return in;
}

std::ostream & cado::operator<<(
    std::ostream & o,
    named_proxy<cxx_mpz_poly_bivariate const &> const & F)
{
    std::string const & x = F.x();
    std::string const & y = F.y();
    cxx_mpz_poly_bivariate const & f = F.c;
    if (f.degree() == -1) {
        return o << 0;
    } else if (f.degree() == 0) {
        return o << f[0].print_poly(x);
    }
    std::ostringstream os;
    for (int i = 0; i <= f.degree(); i++) {
        if (f[i]->deg < 0)
            continue;
        if (f[i] == 1) {
            if (!os.str().empty())
                os << "+";
            if (i == 0)
                os << "1";
        } else if (f[i] == -1) {
            os << "-";
            if (i == 0)
                os << "1";
        } else {
            if (!mpz_poly_is_monomial_multiple(f[i])) {
                /* need parentheses anyway */
                if (!os.str().empty())
                    os << "+";
                os << "(" << f[i].print_poly(x) << ")";
            } else {
                /* if it's a negative multiple of x^i, don't print "+-".
                 */
                if (!os.str().empty() && mpz_sgn(mpz_poly_lc(f[i])) > 0)
                    os << "+";
                os << f[i].print_poly(x);
            }
            if (i > 0)
                os << "*";
        }
        if (i > 0)
            os << y;
        if (i > 1)
            os << "^" << i;
    }
    return o << os.str();
}

void cxx_mpz_poly_bivariate::mod_mpz(cxx_mpz_poly_bivariate & a,
                                     cxx_mpz_poly_bivariate const & b,
                                     mpz_srcptr p)
{
    a = b;
    mpz_ptr invmptr = nullptr;
    cxx_mpz invm;
    if (a.degree() >= 4) {
        barrett_precompute_inverse(invm, p);
        invmptr = invm;
    }
    for (auto & c: a)
        mpz_poly_mod_mpz(c, c, p, invmptr);
}

void cxx_mpz_poly_bivariate::mod_fx(cxx_mpz_poly_bivariate & a,
                                    cxx_mpz_poly_bivariate const & b,
                                    mpz_poly_srcptr fx)
{
    a = b;
    ASSERT_ALWAYS(mpz_poly_is_monic(fx));
    for (auto & c: a)
        mpz_poly_div_r(c, c, fx);
}

void cxx_mpz_poly_bivariate::div_qr(cxx_mpz_poly_bivariate & q,
                                    cxx_mpz_poly_bivariate & r,
                                    cxx_mpz_poly_bivariate const & f,
                                    cxx_mpz_poly_bivariate const & g)
{
    ASSERT_ALWAYS(g.degree() >= 0 && g.lc() == 1);
    ASSERT_ALWAYS(&q != &f && &q != &r && &q != &g);
    ASSERT_ALWAYS(&r != &g);
    /* r == f is allowed */

    const int df = f.degree();
    const int dg = g.degree();
    const int dq = df - dg;

    if (df < dg) /* f is already reduced mod g */
    {
        q = 0;
        r = f;
        return;
    }

    /* now df >= dg */
    ASSERT_FOR_STATIC_ANALYZER(0 <= dg && dg <= df);

    r = f;
    q.assign(dq + 1, 0);
    cxx_mpz_poly tmp;


    for (int k = df - dg; k >= 0; k--) {
        ((super &)q)[k] = r[k + dg];
        for (int j = dg + k - 1; j >= k; j--) {
            mpz_poly_mul(tmp, q[k], g[j - k]);
            mpz_poly_sub(((super &)r)[j], r[j], tmp);
        }
    }
    r.cleandeg(dg - 1);
}

void cxx_mpz_poly_bivariate::mod_fy(cxx_mpz_poly_bivariate & a,
                                    cxx_mpz_poly_bivariate const & b,
                                    cxx_mpz_poly_bivariate const & fy)
{
    cxx_mpz_poly_bivariate q;
    div_qr(q, a, b, fy);
}

void cxx_mpz_poly_bivariate::mul(cxx_mpz_poly_bivariate & a,
                                 cxx_mpz_poly_bivariate const & b,
                                 mpz_poly_srcptr m)
{
    a = b;
    for (auto & c: a)
        mpz_poly_mul(c.x, c.x, m);
}

void cxx_mpz_poly_bivariate::mul(cxx_mpz_poly_bivariate & a,
                                 cxx_mpz_poly_bivariate const & b, mpz_srcptr m)
{
    a = b;
    for (auto & c: a)
        mpz_poly_mul_mpz(c, c, m);
}

void cxx_mpz_poly_bivariate::eval_fy(cxx_mpz_poly & a, self const & f,
                                     cxx_mpz_poly const & e)
{
    ASSERT_ALWAYS(f.degree() >= -1);
    if (f.degree() == -1) {
        a = 0;
        return;
    }
    if (&a == &e) {
        cxx_mpz_poly aa;
        eval_fy(aa, f, e);
        a = aa;
        return;
    }
    a = f.lc();
    for (int i = f.degree() - 1; i >= 0; i--) {
        mpz_poly_mul(a, a, e);
        mpz_poly_add(a, a, f[i]);
    }
}

void cxx_mpz_poly_bivariate::eval_fx(cxx_mpz_poly & a, self const & f,
                                     mpz_srcptr e)
{
    mpz_poly_realloc(a, f.degree() + 1);
    for (size_t i = 0; i < f.size(); i++)
        mpz_poly_eval(mpz_poly_coeff(a, runtime_numeric_cast<int>(i)), f[i], e);
    mpz_poly_cleandeg(a, f.degree());
}

void cxx_mpz_poly_bivariate::transpose(self & a, self && b)
{
    if (&a == &b) {
        self aa = a;
        transpose(aa, b);
        a = aa;
        return;
    }
    const int dx = b.degree_x();
    const int dy = b.degree_y();
    a.assign(dx + 1, 0);
    for (auto & c: a) {
        mpz_poly_realloc(c, dy + 1);
        for (int j = 0; j <= dy; j++)
            mpz_poly_setcoeff_ui(c, j, 0);
    }
    for (int j = 0; j <= dy; j++) {
        for (int i = 0; i <= mpz_poly_degree(((super const &)b)[j]); i++) {
            mpz_swap(mpz_poly_coeff(((super &)a)[i], j),
                     mpz_poly_coeff(((super &)b)[j], i));
        }
    }
    for (int i = 0; i <= dx; i++) {
        mpz_poly_cleandeg(((super &)a)[i], dy);
    }
    a.cleandeg(dx);
    b.clear();
}
void cxx_mpz_poly_bivariate::transpose(self & a, self const & b)
{
    self bb = b;
    transpose(a, (self &&)bb);
}

static constexpr int nth_evaluation_point(int i)
{
    return (i & 1) ? ((i+1)/2) : (-(i/2));
}

void cxx_mpz_poly_bivariate::resultant_y(cxx_mpz_poly & resultant,
                                         self const & f, self const & g)
{
    const size_t nb_eval =
        f.degree_y() * g.degree_x() + f.degree_x() * g.degree_y() + 1;

    std::vector<cxx_mpz> points;
    std::vector<cxx_mpz> evaluations;

    cxx_mpz_poly ef, eg;
    cxx_mpz z;

    for (int i = 0; points.size() < nb_eval; i++) {
        // 0, 1, -1, 2, -2, ...
        const int w = nth_evaluation_point(i);
        eval_fx(ef, f, cxx_mpz(w));
        eval_fx(eg, g, cxx_mpz(w));

        if (ef.degree() < f.degree_y())
            continue;
        if (eg.degree() < g.degree_y())
            continue;

        mpz_poly_resultant(z, ef, eg);
        points.emplace_back(w);
        evaluations.push_back(z);
    }

    const int ok = mpz_poly_interpolate(resultant, points, evaluations);
    if (!ok)
        resultant = 0;
}

void cxx_mpz_poly_bivariate::resultant_x(cxx_mpz_poly & resultant,
                                         self const & f, self const & g)
{
    const size_t nb_eval =
        f.degree_y() * g.degree_x() + f.degree_x() * g.degree_y() + 1;

    std::vector<cxx_mpz> points;
    std::vector<cxx_mpz> evaluations;

    cxx_mpz_poly ef, eg;
    cxx_mpz z;

    for (int i = 0; points.size() < nb_eval; i++) {
        const int w = nth_evaluation_point(i);
        /* simple enough, really. Note that we could also choose
         * polynomials in x as interpolation points, that would work.
         */
        eval_fy(ef, f, cxx_mpz_poly(w));
        eval_fy(eg, g, cxx_mpz_poly(w));

        if (ef.degree() < f.degree_x())
            continue;
        if (eg.degree() < g.degree_x())
            continue;

        mpz_poly_resultant(z, ef, eg);
        points.emplace_back(w);
        evaluations.push_back(z);
    }

    const int ok = mpz_poly_interpolate(resultant, points, evaluations);
    if (!ok)
        resultant = 0;
}

void cxx_mpz_poly_bivariate::set_rrandomb(self & f, int dx, int dy, int bits,
                                          gmp_randstate_ptr rstate)
{
    f.assign(dy + 1, 0);
    for (int i = 0; i <= dy; i++) {
        mpz_poly_set_rrandomb(((super &)f)[i], dx, rstate, bits);
    }
    f.cleandeg(dy);
}

void cxx_mpz_poly_bivariate::set_rrandomb_cab(self & f, int dx, int dy,
                                              int bits,
                                              gmp_randstate_ptr rstate)
{
    // dx and dy must be coprime, but we don't check it.
    // as far as the "Cab" notation goes, dy is a and dx is b.
    set_rrandomb(f, dx, dy, bits, rstate);
    // X^i Y^j appears only if i * dy + j * dx < dx * dy
    ((super &)f)[dy] = 1;
    for (int j = 0; j < dy; j++) {
        // we may have X^i only if i * dy < (dx * dy  - j * dx)
        const int d = (dx * (dy - j) - 1) / dy;
        mpz_poly_cleandeg(((super &)f)[j], d);
        if (j == 0)
            mpz_poly_setcoeff_ui(((super &)f)[j], dx, 1);
    }
}
