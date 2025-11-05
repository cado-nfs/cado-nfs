#include "cado.h" // IWYU pragma: keep

#include <cctype>
#include <cstddef>

#include <algorithm>
#include <ios>
#include <istream>
#include <ostream>
#include <sstream>
#include <string>
#include <utility>

#include <gmp.h>
#include "fmt/base.h"

#include "cado_expression_parser.hpp"
#include "cxx_mpz.hpp"
#include "mpz_poly.h"
#include "arithmetic_reductions.hpp"
#include "mpz_poly_bivariate.hpp"
#include "runtime_numeric_cast.hpp"
#include "macros.h"
#include "mpz_mat.h"
#include "named_proxy.hpp"
#include "timing.h"

/* Polynomial arithmetic.
 *
 * additions / negations ; we don't put any reduction function un the
 * interface because they're always done just as well externally
 *
 * It's the same story for multiplication.
 */
void cxx_mpz_poly_bivariate::neg(cxx_mpz_poly_bivariate & f,
                                 cxx_mpz_poly_bivariate const & g)
{
    f = g;
    for (size_t i = 0; i < g.size(); i++) {
        mpz_poly_neg(f.v()[i], f[i]);
    }
}

void cxx_mpz_poly_bivariate::add(cxx_mpz_poly_bivariate & f,
                                 cxx_mpz_poly_bivariate const & g,
                                 cxx_mpz_poly_bivariate const & h)
{
    size_t const sg = g.size();
    size_t const sh = h.size();
    size_t const sf = std::max(sg, sh);
    f.reserve(sf);
    size_t i = 0;
    if (&f != &g && &f != &h)
        f.clear();
    f.insert(f.end(), sf - f.size(), cxx_mpz_poly()); /* now sf == f.size() */
    for (; i < sg && i < sh; i++)
        mpz_poly_add(f.v()[i], g[i], h[i]);
    for (; i < sg; i++)
        mpz_poly_set(f.v()[i], g[i]);
    for (; i < sh; i++)
        mpz_poly_set(f.v()[i], h[i]);
    f.cleandeg((int)sf - 1);
}

void cxx_mpz_poly_bivariate::sub(cxx_mpz_poly_bivariate & f,
                                 cxx_mpz_poly_bivariate const & g,
                                 cxx_mpz_poly_bivariate const & h)
{
    size_t const sg = g.size();
    size_t const sh = h.size();
    size_t const sf = std::max(sg, sh);
    f.reserve(sf);
    size_t i = 0;
    if (&f != &g && &f != &h)
        f.clear();
    f.insert(f.end(), sf - f.size(), cxx_mpz_poly()); /* now sf == f.size() */
    for (; i < sg && i < sh; i++)
        mpz_poly_sub(f.v()[i], g[i], h[i]);
    for (; i < sg; i++)
        mpz_poly_set(f.v()[i], g[i]);
    for (; i < sh; i++)
        mpz_poly_neg(f.v()[i], h[i]);
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
            mpz_poly_add(f.v()[i + j], f[i + j], tmp);
        }
    }
    f.cleandeg(g.degree() + h.degree());
}

void cxx_mpz_poly_bivariate::mul(cxx_mpz_poly_bivariate & a,
                                 cxx_mpz_poly_bivariate const & b,
                                 mpz_poly_srcptr m)
{
    a = b;
    for (auto & c: a.v())
        mpz_poly_mul(c.x, c.x, m);
}

void cxx_mpz_poly_bivariate::mul(cxx_mpz_poly_bivariate & a,
                                 cxx_mpz_poly_bivariate const & b,
                                 mpz_srcptr m) /*{{{*/
{
    a = b;
    for (auto & c: a.v())
        mpz_poly_mul_mpz(c, c, m);
}
/*}}}*/

/* random generation */
void cxx_mpz_poly_bivariate::set_rrandomb(self & f, int dx, int dy, int bits,
                                          gmp_randstate_ptr rstate)
{
    f.assign(dy + 1, 0);
    for (int i = 0; i <= dy; i++) {
        mpz_poly_set_randomb(f.v()[i], dx, rstate, bits,
                mpz_poly_random_flags::MPZ_POLY_RRANDOM |
                mpz_poly_random_flags::MPZ_POLY_DEGREE_EXACT
                );
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
    f.v()[dy] = 1;
    for (int j = 0; j < dy; j++) {
        // we may have X^i only if i * dy < (dx * dy  - j * dx)
        int const d = (dx * (dy - j) - 1) / dy;
        mpz_poly_cleandeg(f.v()[j], d);
        if (j == 0)
            mpz_poly_setcoeff_ui(f.v()[j], dx, 1);
    }
}

/*
 * exact: degree in x is exact (degree in y always is)
 * monic: monic in y (we never enforce monic in x)
 */
void cxx_mpz_poly_bivariate::set_urandomm(self & f, int dx, int dy,
                                          mpz_srcptr p,
                                          gmp_randstate_ptr rstate,
                                          bool exact,
                                          bool monic) /*{{{*/
{
    f.assign(dy + 1, 0);
    for(auto & c : f.v()) {
        mpz_poly_set_randomm(c, dx, rstate, p,
                MPZ_POLY_URANDOM |
                (exact ? MPZ_POLY_DEGREE_EXACT : 0)
                );
    }
    if (monic)
        f.v()[dy] = 1;
    f.cleandeg(dy);
}
/*}}}*/

void cxx_mpz_poly_bivariate::set_urandomm_cab(self & f, int dx, int dy,
                                              mpz_srcptr p,
                                              gmp_randstate_ptr rstate) /*{{{*/
{
    // dx and dy must be coprime, but we don't check it.
    // as far as the "Cab" notation goes, dy is a and dx is b.
    set_urandomm(f, dx, dy, p, rstate, true, true);
    // X^i Y^j appears only if i * dy + j * dx < dx * dy
    f.v()[dy] = 1;
    for (int j = 0; j < dy; j++) {
        // we may have X^i only if i * dy < (dx * dy  - j * dx)
        int const d = (dx * (dy - j) - 1) / dy;
        mpz_poly_cleandeg(f.v()[j], d);
        if (j == 0)
            mpz_poly_setcoeff_ui(f.v()[j], dx, 1);
    }
}
/*}}}*/

/* transpose */
void cxx_mpz_poly_bivariate::transpose(self & a, self && b)
{
    if (&a == &b) {
        transpose(a, self(b));
        return;
    }
    int const dx = b.degree_x();
    int const dy = b.degree_y();
    a.assign(dx + 1, 0);
    for (auto & c: a.v()) {
        mpz_poly_realloc(c, dy + 1);
        for (int j = 0; j <= dy; j++)
            mpz_poly_setcoeff_ui(c, j, 0);
    }
    for (int j = 0; j <= dy; j++) {
        for (int i = 0; i <= mpz_poly_degree(b[j]); i++) {
            mpz_swap(mpz_poly_coeff(a.v()[i], j),
                     mpz_poly_coeff(b.v()[j], i));
        }
    }
    for (int i = 0; i <= dx; i++) {
        mpz_poly_cleandeg(a.v()[i], dy);
    }
    a.cleandeg(dx);
    b.clear();
}

void cxx_mpz_poly_bivariate::transpose(self & a, self const & b)
{
    transpose(a, self(b));
}

/* parsing */
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

    struct unexpected_literal
        : public cado_expression_parser_details::parse_error {
        std::string msg;
        explicit unexpected_literal(std::string const & v)
            : msg(std::string("unexpected literal " + v))
        {
        }
        char const * what() const noexcept override { return msg.c_str(); }
    };

    static constexpr int const accept_literals = 2;
    using type = cxx_mpz_poly_bivariate;
    using number_type = cxx_mpz;

    static void add(cxx_mpz_poly_bivariate & c,
                    cxx_mpz_poly_bivariate const & a,
                    cxx_mpz_poly_bivariate const & b)
    {
        type::add(c, a, b);
    }
    static void sub(cxx_mpz_poly_bivariate & c,
                    cxx_mpz_poly_bivariate const & a,
                    cxx_mpz_poly_bivariate const & b)
    {
        type::sub(c, a, b);
    }
    static void mul(cxx_mpz_poly_bivariate & c,
                    cxx_mpz_poly_bivariate const & a,
                    cxx_mpz_poly_bivariate const & b)
    {
        type::mul(c, a, b);
    }
    static void pow(cxx_mpz_poly_bivariate & c,
                       cxx_mpz_poly_bivariate const & a, unsigned long e)
    {
        type::pow(c, a, e);
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
        int const c = in.peek();
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

std::string cxx_mpz_poly_bivariate::print_poly(std::string const& x, std::string const & y) const
{
    std::ostringstream os;
    os << named(x, y);
    return os.str();
}

/* this proxy is used for pretty printing in gdb */
std::string cxx_mpz_poly_bivariate::print_poly() const
{
    return print_poly("x", "y");
}

/* derivative_y. sort of one of a kind, really. We use it for square free
 * factorization. */
void cxx_mpz_poly_bivariate::derivative_y(self & B, self const & A)
{
    int degree_of_B = (A.degree_y() < 0) ? -1 : A.degree() - 1;
    B = A;
    /* take coefficients in increasing order, so that B==A is allowed */
    for (int i = 1; i <= A.degree(); i++)
        mpz_poly_mul_si(B.v()[i - 1], B[i], i);
    B.cleandeg(degree_of_B);
}

/* evaluation */
void cxx_mpz_poly_bivariate::eval_fy(cxx_mpz_poly & a, self const & f,
                                     cxx_mpz_poly const & e)
{
    ASSERT_ALWAYS(f.degree() >= -1);
    if (f.degree() == -1) {
        a = 0;
        return;
    }
    if (&a == &e) {
        eval_fy(a, f, cxx_mpz_poly(e));
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

/* resultant */
static constexpr int nth_evaluation_point(int i)
{
    return (i & 1) ? ((i + 1) / 2) : (-(i / 2));
}

void cxx_mpz_poly_bivariate::resultant_y(cxx_mpz_poly & resultant,
                                         self const & f, self const & g)
{
    size_t const nb_eval =
        f.degree_y() * g.degree_x() + f.degree_x() * g.degree_y() + 1;

    std::vector<cxx_mpz> points;
    std::vector<cxx_mpz> evaluations;

    cxx_mpz_poly ef, eg;
    cxx_mpz z;

    for (int i = 0; points.size() < nb_eval; i++) {
        // 0, 1, -1, 2, -2, ...
        int const w = nth_evaluation_point(i);
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

    int const ok = mpz_poly_interpolate(resultant, points, evaluations);
    if (!ok)
        resultant = 0;
}

void cxx_mpz_poly_bivariate::resultant_x(cxx_mpz_poly & resultant,
                                         self const & f, self const & g)
{
    size_t const nb_eval =
        f.degree_y() * g.degree_x() + f.degree_x() * g.degree_y() + 1;

    std::vector<cxx_mpz> points;
    std::vector<cxx_mpz> evaluations;

    cxx_mpz_poly ef, eg;
    cxx_mpz z;

    for (int i = 0; points.size() < nb_eval; i++) {
        int const w = nth_evaluation_point(i);
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

    int const ok = mpz_poly_interpolate(resultant, points, evaluations);
    if (!ok)
        resultant = 0;
}

/* basic reduction functions */
void cxx_mpz_poly_bivariate::mod_mpz(cxx_mpz_poly_bivariate & a,
                                     cxx_mpz_poly_bivariate const & b,
                                     mpz_srcptr p) /*{{{*/
{
    if (b.degree() >= 4) {
        cxx_mpz invp;
        barrett_precompute_inverse(invp, p);
        mod_mpz(a, b, p, invp);
    } else {
        mod_mpz(a, b, p, NULL);
    }
}
/*}}}*/
void cxx_mpz_poly_bivariate::mod_mpz(cxx_mpz_poly_bivariate & a,
                                     cxx_mpz_poly_bivariate const & b,
                                     mpz_srcptr p, mpz_srcptr invp) /*{{{*/
{
    a = b;
    for (auto & c: a.v())
        mpz_poly_mod_mpz(c, c, p, invp);
}
/*}}}*/
void cxx_mpz_poly_bivariate::mod_fx(cxx_mpz_poly_bivariate & a,
                                    cxx_mpz_poly_bivariate const & b,
                                    mpz_poly_srcptr fx) /*{{{*/
{
    a = b;
    ASSERT_ALWAYS(mpz_poly_is_monic(fx));
    for (auto & c: a.v())
        mpz_poly_div_r(c, c, fx);
    a.cleandeg(a.degree());
}

/*}}}*/
void cxx_mpz_poly_bivariate::mod_fy(cxx_mpz_poly_bivariate & a,
                                    cxx_mpz_poly_bivariate const & b,
                                    cxx_mpz_poly_bivariate const & fy) /*{{{*/
{
    cxx_mpz_poly_bivariate q;
    div_qr(q, a, b, fy);
}
/*}}}*/

/* division */
template <typename T>
void cxx_mpz_poly_bivariate::div_qr(cxx_mpz_poly_bivariate & q,
                                    cxx_mpz_poly_bivariate & r,
                                    cxx_mpz_poly_bivariate const & f,
                                    cxx_mpz_poly_bivariate const & g,
                                    T const & reducer) /*{{{*/
{
    ASSERT_ALWAYS(g.degree() >= 0 && g.lc() == 1);
    ASSERT_ALWAYS(&q != &r && &q != &g);
    ASSERT_ALWAYS(&r != &g);
    /* r == f or q == f are both allowed */

    int const df = f.degree();
    int const dg = g.degree();
    int const dq = df - dg;

    if (df < dg) /* f is already reduced mod g */
    {
        r = f;
        q = 0;
        return;
    }

    /* now df >= dg */

    r = f;
    q.assign(dq + 1, 0);
    cxx_mpz_poly tmp;

    for (int k ; (k = r.degree() - dg) >= 0 ; ) {
        q.v()[k] = r[k + dg];
        for (int j = dg + k - 1; j >= k; j--) {
            mpz_poly_mul(tmp, q[k], g[j - k]);
            mpz_poly_sub(r.v()[j], r[j], tmp);
        }
        r.cleandeg(k + dg - 1);
        reducer(r, r);
    }
}
/*}}}*/

void cxx_mpz_poly_bivariate::div_qr(self & q, self & r, self const & f, self const & g)
{ /*{{{*/
    div_qr(q, r, f, g, cado::arithmetic_reductions::noop {});
} /*}}}*/
void cxx_mpz_poly_bivariate::div_qr_mod_mpz(self & q, self & r, self const & f,
        self const & g, mpz_srcptr p)
{ /*{{{*/
    div_qr(q, r, f, g, cado::arithmetic_reductions::mod_mpz {p});
} /*}}}*/
void cxx_mpz_poly_bivariate::div_qr_mod_fx_mod_p(self & q, self & r, self const & f,
        self const & g, mpz_poly_srcptr fx,
        mpz_srcptr p)
{ /*{{{*/
    div_qr(q, r, f, g, cado::arithmetic_reductions::mod_fx_mod_p {fx, p});
} /*}}}*/
void cxx_mpz_poly_bivariate::div_q(self & q, self const & f, self const & g)
{ /*{{{*/
    div_q(q, f, g, cado::arithmetic_reductions::noop {});
} /*}}}*/
void cxx_mpz_poly_bivariate::div_q_mod_mpz(self & q, self const & f, self const & g,
        mpz_srcptr p)
{ /*{{{*/
    div_q(q, f, g, cado::arithmetic_reductions::mod_mpz {p});
} /*}}}*/
void cxx_mpz_poly_bivariate::div_q_mod_fx_mod_p(self & q, self const & f, self const & g,
        mpz_poly_srcptr fx, mpz_srcptr p)
{ /*{{{*/
    div_q(q, f, g, cado::arithmetic_reductions::mod_fx_mod_p {fx, p});
} /*}}}*/
void cxx_mpz_poly_bivariate::div_r(self & r, self const & f, self const & g)
{ /*{{{*/
    div_r(r, f, g, cado::arithmetic_reductions::noop {});
} /*}}}*/
void cxx_mpz_poly_bivariate::div_r_mod_mpz(self & r, self const & f, self const & g,
        mpz_srcptr p)
{ /*{{{*/
    div_r(r, f, g, cado::arithmetic_reductions::mod_mpz {p});
} /*}}}*/
void cxx_mpz_poly_bivariate::div_r_mod_fx_mod_p(self & r, self const & f, self const & g,
        mpz_poly_srcptr fx, mpz_srcptr p)
{ /*{{{*/
    div_r(r, f, g, cado::arithmetic_reductions::mod_fx_mod_p {fx, p});
} /*}}}*/
void cxx_mpz_poly_bivariate::divexact(self & q, self const & f, self const & g)
{ /*{{{*/
    divexact(q, f, g, cado::arithmetic_reductions::noop {});
} /*}}}*/
void cxx_mpz_poly_bivariate::divexact_mod_mpz(self & q, self const & f, self const & g,
        mpz_srcptr p)
{ /*{{{*/
    divexact(q, f, g, cado::arithmetic_reductions::mod_mpz {p});
} /*}}}*/
void cxx_mpz_poly_bivariate::divexact_mod_fx_mod_p(self & q, self const & f,
        self const & g, mpz_poly_srcptr fx,
        mpz_srcptr p)
{ /*{{{*/
    divexact(q, f, g, cado::arithmetic_reductions::mod_fx_mod_p {fx, p});
} /*}}}*/

/* powering. same pattern.  */
void cxx_mpz_poly_bivariate::pow(cxx_mpz_poly_bivariate & B,
        cxx_mpz_poly_bivariate const & A,
        unsigned long n) /*{{{*/
{
    pow(B, A, n, cado::arithmetic_reductions::noop {});
}
/*}}}*/
void cxx_mpz_poly_bivariate::pow_mod_mpz(cxx_mpz_poly_bivariate & B,
        cxx_mpz_poly_bivariate const & A,
        unsigned long n, mpz_srcptr p) /*{{{*/
{
    pow(B, A, n, cado::arithmetic_reductions::mod_mpz {p});
}
/*}}}*/
void cxx_mpz_poly_bivariate::pow_mod_fy_mod_q(self & B, self const & A,
        unsigned long n,
        mpz_poly_srcptr fy,
        mpz_poly_srcptr fx,
        mpz_srcptr p) /*{{{*/
{
    pow(B, A, n, cado::arithmetic_reductions::mod_fy_mod_q {fy, fx, p});
}
/*}}}*/



void cxx_mpz_poly_bivariate::mul_mod_fy_mod_q(self & C, self const & A,
                                                       self const & B,
                                                       mpz_poly_srcptr fy,
                                                       mpz_poly_srcptr fx,
                                                       mpz_srcptr p)
{
    mul(C, A, B);
    cado::arithmetic_reductions::mod_fy_mod_q {fy, fx, p}(C, C);
}

/* gcd. Well it's only about being able to run the Euclidean algorithm,
 * so we're talking about univariate polynomials over finite fields in
 * disguise. So we need a reducer, and that reducer has to be both mod a
 * polynomial and mod a prime. No point in defining gcd with more general
 * reducers, that would not make much sense.
 */
void cxx_mpz_poly_bivariate::gcd(self & C, self const & A, self const & B,
                                 cado::arithmetic_reductions::mod_q const & Rx) /*{{{*/
{
    self R0 = A;
    self R1 = B;
    self R2;
    for (; R1.degree_y() >= 0;) {
        Rx.make_monic(R1);
        div_r(R2, R0, R1, Rx);
        std::swap(R0, R1);
        Rx(R1, R2);
    }
    C = R0;
}
/*}}}*/

static int sqf_inner(cxx_mpz_poly_bivariate::factor_list & fl, cxx_mpz_poly_bivariate const & f, int stride, cado::arithmetic_reductions::mod_q const & Rx)
{
    int r = 0;

    using self = cxx_mpz_poly_bivariate;

    self g, mi, mi1;
    self t0, t1, T, tmp;

    self::derivative_y(t0, f);
    self::gcd(g, t0, f, Rx);
    self::divexact(mi, f, g, Rx);
    /* mi is f/gcd(f,f') == all repeated prime factors of f whose
     * multiplicity isn't a multiple of the field characteristic.
     */

    T.set_yi(0);
    for (int i = 1; mi.degree_y() > 0; i++) {
        /* note the weird argument ordering */
        self::pow(t0, mi, i, Rx);
        self::divexact(t1, f, t0, Rx);
        /* t1 = all polynomials in mi taken out from f with multiplicity i */
        self::gcd(mi1, t1, mi, Rx);
        /* mi1 = almost like mi, but since factors with multiplicity i
         * are no longer in t1, there's absent from mi1 too. Whence
         * mi/mi1 is exactly the product of factors of multiplicity 1.
         */
        self::divexact(tmp, mi, mi1, Rx);
        if ((size_t) (i * stride) >= fl.size()) {
            /* insert unit polynomials for consistency */
            fl.insert(fl.end(), i*stride + 1 - fl.size(), {1, 1});
        }
        /* Use tmp so that we don't absurdly keep storage within
         * lf->factors */
        fl[i * stride] = { tmp, 1 }; /* multiplicity field in fact unused */
        self::pow(t0, fl[i * stride].first, i, Rx);
        self::mul(T, T, t0);
        Rx(T, T);
        mi.swap(mi1);
        r = i * stride;
    }

    if (fl.empty())
        fl.insert(fl.end(), {1, 1});

    self::divexact(fl[0].first, f, T, Rx);
    fl[0].second = 1;

    return r;
}

/* returns a list of square-free factors of f. In the returned
 * factor_list, item i has multiplicity i. Note that f0 is made monic
 * before any work, and we do not keep track of the leading coefficient
 */
cxx_mpz_poly_bivariate::factor_list cxx_mpz_poly_bivariate::factor_sqf(cxx_mpz_poly_bivariate const & f0, cado::arithmetic_reductions::mod_q const & R)
{
    using self = cxx_mpz_poly_bivariate;

    /* factoring 0 doesn't make sense, really */
    ASSERT(f0.degree() >= 0);

    /* We'll call mpz_poly_factor_sqf_inner, possibly several times if
     * we are in small characteristic.
     */
    self f = f0;
    R.make_monic(f);
    ASSERT_ALWAYS(f.lc() == 1);

    int m = 0;
    // coverity[zero_return]
    int pu = mpz_get_ui(R.p);  // see below

    cxx_mpz_poly_bivariate::factor_list lf;

    for(int stride = 1 ; ; stride *= pu) {
        int r = sqf_inner(lf, f, stride, R);
        if (r > m) m = r;
        cxx_mpz_poly_bivariate repeated = 1;
        repeated.swap(lf[0].first);
        if (repeated.degree() == 0) {
            // if p is LAAAARGE, then of course we'll never have a linear
            // polynomial out of sqf_inner, thus we'll break early here.
            break;
        }
        /* lf[0].first is a p-th power. Take the p-th root, which also
         * means taking the p-th root of each coefficient (mod fx).
         */
        R.frobenius(f, repeated, -1);

    }
    /* Now make sure that all factors in the factor list are non-zero */
    for(auto & c : lf) {
        ASSERT_ALWAYS(c.first.degree_y() >= 0);
    }
    return lf;
}

/* This performs distinct degree factorization */
/* Input polynomial must be squarefree -- otherwise repeated factors
 * probably won't show up in the factor list, or maybe at the wrong place
 * as parasites. */
/* Coefficients of f0 need not be reduced mod p.
 * Coefficients of all polynomials stored in lf are reduced mod p.
 */
cxx_mpz_poly_bivariate::factor_list cxx_mpz_poly_bivariate::factor_ddf(cxx_mpz_poly_bivariate const & f0, cado::arithmetic_reductions::mod_q const & R, bool only_check_irreducible)
{
    /* factoring 0 doesn't make sense, really */
    ASSERT_ALWAYS(f0.degree() >= 0);

    ASSERT_ALWAYS(mpz_size(R.p) == 1);

    // mp_limb_t pp = mpz_get_ui(R.p);

    self f = f0;
    R.make_monic(f);
    ASSERT_ALWAYS(f.lc() == 1);

    self::factor_list lf;

    self g, gmy, y, tmp;

    y.set_yi(1);
    g.set_yi(1);

    double t0 = 0;

    lf.emplace_back(1, 1);

    auto Rf = std::make_unique<cado::arithmetic_reductions::mod_fy_mod_q>(R, f);

    for (int i = 1; i <= f.degree() ; ++i) {
        ASSERT_ALWAYS(i == (int) lf.size());
        if (2 * i > f.degree()) {
            /* Then we know that the remaining f is irreducible.  */
            /* multiplicity field still unused at this point */
            lf.insert(lf.end(), f.degree() - i, {1, 1});
            lf.emplace_back( std::move(f), 1 );
            break;
        }

        /* g <- g^(p^(deg(fx))) mod fx, p */
        t0 -= wct_seconds();
        Rf->frobenius(g, g, 1);
        t0 += wct_seconds();
        fmt::print("g^(pp^{})-g mod (pp, fx, f): {:.3f}\n",
                R.fx.degree(), t0);

        /* subtract y */
        self::sub(gmy, g, y);

        R(gmy, gmy);

        /* lf[i] <- gcd (f, y^(p^(deg(fx)*i))-y) */

        self::gcd(tmp, f, gmy, R);

        lf.emplace_back(tmp, 1);

        if (tmp.degree() == 0)
            continue;

        fmt::print("deg={}: factors degree = {}\n", i, tmp.degree());

        self::divexact(f, f, tmp, R);
        
        /* Rf has a _reference_ to f, and no derived data. It's
         * dangerous, I know, but for the moment it seems fine. So we
         * don't have to reinitialize it (which wouldn't be possible
         * anyway because of embedded refs). */
        // Rf = self::mod_fy_mod_q(R, f);
        Rf = std::make_unique<cado::arithmetic_reductions::mod_fy_mod_q>(R, f);

        /* Note for a mere irreducibility test: the length of the loop in
         * the irreducible case would still be deg(f)/2, and the penalty
         * caused by storing factors can be neglected.
         */
        if (only_check_irreducible && lf[i].first.degree() > 0)
            break;

        if (f.degree() == 0)
            break;
    }
    return lf;
}

