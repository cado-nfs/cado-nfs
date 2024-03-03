#include "cado.h" // IWYU pragma: keep

#include <cctype>
#include <cstddef>

#include <algorithm>
#include <ios>
#include <istream>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include <gmp.h>

#include "cado_expression_parser.hpp"
#include "cxx_mpz.hpp"
#include "mpz_poly.h"
#include "mpz_poly_bivariate.hpp"
#include "runtime_numeric_cast.hpp"
#include "macros.h"
#include "mpz_mat.h"
#include "named_proxy.hpp"

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
        mpz_poly_neg(((super &)f)[i], ((super &)f)[i]);
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
    size_t const sg = g.size();
    size_t const sh = h.size();
    size_t const sf = std::max(sg, sh);
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

void cxx_mpz_poly_bivariate::mul(cxx_mpz_poly_bivariate & a,
                                 cxx_mpz_poly_bivariate const & b,
                                 mpz_poly_srcptr m)
{
    a = b;
    for (auto & c: a)
        mpz_poly_mul(c.x, c.x, m);
}

void cxx_mpz_poly_bivariate::mul(cxx_mpz_poly_bivariate & a,
                                 cxx_mpz_poly_bivariate const & b,
                                 mpz_srcptr m) /*{{{*/
{
    a = b;
    for (auto & c: a)
        mpz_poly_mul_mpz(c, c, m);
}
/*}}}*/

/* random generation */
void cxx_mpz_poly_bivariate::set_rrandomb(self & f, int dx, int dy, int bits,
                                          gmp_randstate_ptr rstate)
{
    f.assign(dy + 1, 0);
    for (int i = 0; i <= dy; i++) {
        mpz_poly_set_randomb(((super &)f)[i], dx, rstate, bits,
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
    ((super &)f)[dy] = 1;
    for (int j = 0; j < dy; j++) {
        // we may have X^i only if i * dy < (dx * dy  - j * dx)
        int const d = (dx * (dy - j) - 1) / dy;
        mpz_poly_cleandeg(((super &)f)[j], d);
        if (j == 0)
            mpz_poly_setcoeff_ui(((super &)f)[j], dx, 1);
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
    for (int i = 0; i <= dy; i++) {
        mpz_poly_set_randomm(((super &)f)[i], dx, rstate, p,
                MPZ_POLY_URANDOM |
                (exact ? MPZ_POLY_DEGREE_EXACT : 0)
                );
    }
    if (monic) {
        ((super &)f)[dy] = 1;
    }
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
    ((super &)f)[dy] = 1;
    for (int j = 0; j < dy; j++) {
        // we may have X^i only if i * dy < (dx * dy  - j * dx)
        int d = (dx * (dy - j) - 1) / dy;
        mpz_poly_cleandeg(((super &)f)[j], d);
        if (j == 0)
            mpz_poly_setcoeff_ui(((super &)f)[j], dx, 1);
    }
}
/*}}}*/

/* transpose */
void cxx_mpz_poly_bivariate::transpose(self & a, self && b)
{
    if (&a == &b) {
        self bb = b;
        transpose(a, bb);
        return;
    }
    int const dx = b.degree_x();
    int const dy = b.degree_y();
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
    typedef cxx_mpz_poly_bivariate type;
    typedef cxx_mpz number_type;

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
    static void pow_ui(cxx_mpz_poly_bivariate & c,
                       cxx_mpz_poly_bivariate const & a, unsigned long e)
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

/* derivative_y. sort of one of a kind, really. We use it for square free
 * factorization. */
void cxx_mpz_poly_bivariate::derivative_y(self & B, self const & A)
{
    int degree_of_B = (A.degree_y() < 0) ? -1 : A.degree() - 1;
    B = A;
    /* take coefficients in increasing order, so that B==A is allowed */
    for (int i = 1; i <= A.degree(); i++)
        mpz_poly_mul_si(((super &)B)[i - 1], B[i], i);
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
        cxx_mpz_poly ee = e;
        eval_fy(a, f, ee);
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
    for (auto & c: a)
        mpz_poly_mod_mpz(c, c, p, invp);
}
/*}}}*/
void cxx_mpz_poly_bivariate::mod_fx(cxx_mpz_poly_bivariate & a,
                                    cxx_mpz_poly_bivariate const & b,
                                    mpz_poly_srcptr fx) /*{{{*/
{
    a = b;
    ASSERT_ALWAYS(mpz_poly_is_monic(fx));
    for (auto & c: a)
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

/* reduction operator inside function objects */
void cxx_mpz_poly_bivariate::reducer_mod_mpz::operator()(
    self & b, self const & a) const /*{{{*/
{
    mpz_srcptr inv = invp.get() ? (mpz_srcptr)*invp : nullptr;
    b = a;
    for (auto & c: b)
        mpz_poly_mod_mpz(c, c, p, inv);
    b.cleandeg(b.degree());
}
/*}}}*/
void cxx_mpz_poly_bivariate::reducer_mod_fx::operator()(
    self & B, self const & A) const /*{{{*/
{
    self::mod_fx(B, A, fx);
}
/*}}}*/
void cxx_mpz_poly_bivariate::reducer_mod_fx_mod_mpz::operator()(
    self & B, self const & A) const /*{{{*/
{
    // assume that we are are in the quadratic range, and that the
    // following conditions hold:
    //  - deg(fx)=dx ; fx is monic. Coefficients are reduced mod p
    //  - p is np words
    //  - coeffs in y^i of A are polynomials of degree k*dx with l*np words
    //  - a Barrett preinverse of p is known (so that reducing l*np words
    //    mod p costs (l-1) mulhi(np,np) + 1 mullo(np,np) -- we don't
    //    implement mulhi and mullo, so we count (l-1)*2*np^2).
    //
    // Furthermore, we assume that dx>=2, k>=1, l>=1
    //
    // In this case, the costs and new values for s,k,l are:
#if 0
    sage: QP.<k,l,dx,np> = QQ[]
    sage: def mod_p(k,l):
    sage:     return (k*dx*(l-1)*2*np^2, k, 1)
    sage: def mod_fx(k, l):
    sage:     return ((k-1)*dx^2*l*np^2, 1, l+1)
#endif
    // We can compare different strategies.
#if 0
    sage: def chain(*functions):
    sage:     A=(k,l); cost=0
    sage:     for f in functions:
    sage:         (c,*A) = f(*A); cost += c
    sage:     return cost / (dx * np^2)
    sage: cost0 = chain(mod_fx, mod_p)
    sage: cost1 = chain(mod_p, mod_fx, mod_p)
    sage: assert cost1 - cost0 == - (dx - 2) * (l - 1) * (k - 1)
#endif
    // So cost1 wins over cost0
    reducer_mod_mpz::operator()(B, A);
    self::mod_fx(B, B, fx);
    reducer_mod_mpz::operator()(B, B);
}
/*}}}*/
void cxx_mpz_poly_bivariate::reducer_mod_mpz::operator()(
    cxx_mpz_poly & b, cxx_mpz_poly const & a) const /*{{{*/
{
    mpz_srcptr inv = invp.get() ? (mpz_srcptr)*invp : nullptr;
    mpz_poly_mod_mpz(b, a, p, inv);
}
/*}}}*/
void cxx_mpz_poly_bivariate::reducer_mod_fx::operator()(
    cxx_mpz_poly & B, cxx_mpz_poly const & A) const /*{{{*/
{
    mpz_poly_div_r(B, A, fx);
}
/*}}}*/
void cxx_mpz_poly_bivariate::reducer_mod_fx_mod_mpz::operator()(
    cxx_mpz_poly & B, cxx_mpz_poly const & A) const /*{{{*/
{
    reducer_mod_mpz::operator()(B, A);
    mpz_poly_div_r(B, B, fx);
    reducer_mod_mpz::operator()(B, B);
}
/*}}}*/
void cxx_mpz_poly_bivariate::reducer_mod_fy_mod_fx_mod_mpz::operator()(
    self & B, self const & A) const /*{{{*/
{
    // assume that we are are in the quadratic range, and that the
    // following conditions hold:
    //  - deg_y(fy) = dy ; fy is monic ; coefficients are reduced mod fx
    //  - deg(fx)=dx ; fx is monic. Coefficients are reduced mod p
    //  - p is np words
    //  - deg_y(A) = s*dy
    //  - coeffs in y^i of A are polynomials of degree k*dx with l*np words
    //  - a Barrett preinverse of p is known (so that reducing l*np words
    //    mod p costs (l-1) mulhi(np,np) + 1 mullo(np,np) -- we don't
    //    implement mulhi and mullo, so we count (l-1)*2*np^2).
    // Furthermore, we assume that dx>=2, dy>=2, s>=1, k>=1, l>=1
    //
    // In this case, the costs and new values for s,k,l are:
#if 0
    sage: QP.<s,k,l,dx,dy,np> = QQ[]
    sage: def mod_p(s,k,l):
    sage:     return (s*dy*k*dx*(l-1)*2*np^2, s, k, 1)
    sage: def mod_fx(s, k, l):
    sage:     return (s*dy*(k-1)*dx^2*l*np^2, s, 1, l+1)
    sage: def mod_fy(s, k, l):
    sage:     return ((s-1)*dy^2*k*dx^2*l*np^2, 1, k+1, l+1)
#endif
    // We can compare different strategies.
#if 0
    sage: def chain(*functions):
    sage:     A=(s,k,l); cost=0
    sage:     for f in functions:
    sage:         (c,*A) = f(*A); cost += c
    sage:     return cost / (dx * dy * np^2)
    sage: cost0 = chain(mod_fy, mod_fx, mod_p)
    sage: cost1 = chain(mod_fy, mod_p, mod_fx, mod_p)
    sage: cost3 = chain(mod_p, mod_fy, mod_p, mod_fx, mod_p)
    sage: assert ((cost3-cost1)) == (l-1)*((s-1)*k*(2-dx*dy) - 2)
#endif
    // So cost1 wins over cost0 if dx >= 2, and cost3 always wins over cost1
#if 0
    sage: cost4 = chain(mod_fx, mod_p, mod_fy, mod_p, mod_fx, mod_p)
    sage: cost5 = chain(mod_p, mod_fx, mod_p, mod_fy, mod_p, mod_fx, mod_p)
    sage: assert cost5-cost3 == (k-1)*(s-1)*(dx+2-dx*dy)-(k-2)*2*s
#endif
    // if k>=2, cost5 always wins. If k=1 it doesn't
    reducer_mod_fx_mod_mpz::operator()(B, A);
    self::mod_fy(B, B, fy);
    reducer_mod_fx_mod_mpz::operator()(B, B);
}
/*}}}*/

std::string cxx_mpz_poly_bivariate::reducer_noop::print() const
{
    return "dummy no-reduction";
}

std::string cxx_mpz_poly_bivariate::reducer_mod_mpz::print() const
{
    std::ostringstream os;
    os << "reduction mod " << p;
    return os.str();
}

std::string cxx_mpz_poly_bivariate::reducer_mod_fx::print() const
{
    std::ostringstream os;
    os << "reduction mod " << fx.named("x");
    return os.str();
}

std::string cxx_mpz_poly_bivariate::reducer_mod_fx_mod_mpz::print() const
{
    std::ostringstream os;
    os << "reduction mod " << fx.named("x") << ", mod " << p;
    return os.str();
}

std::string cxx_mpz_poly_bivariate::reducer_mod_fy_mod_fx_mod_mpz::print() const
{
    std::ostringstream os;
    os << "reduction mod " << fy.named("x", "y") << ", mod " << fx.named("x")
       << ", mod " << p;
    return os.str();
}

void cxx_mpz_poly_bivariate::reducer_mod_fx_mod_mpz::compute_frobenius_matrices() const
{
    unsigned int d = fx->deg;
    ASSERT_ALWAYS(fx->deg > 0);
    if (frobenius_matrix->m) {
        ASSERT_ALWAYS(frobenius_matrix->m == d);
        ASSERT_ALWAYS(frobenius_matrix->n == d);
        ASSERT_ALWAYS(inverse_frobenius_matrix->m == d);
        ASSERT_ALWAYS(inverse_frobenius_matrix->n == d);
        return;
    }
    mpz_mat_realloc(frobenius_matrix, d, d);
    for(unsigned int i = 0 ; i < d ; i++) {
        cxx_mpz_poly R;
        mpz_poly_set_xi(R, i);
        mpz_poly_pow_mod_f_mod_mpz(R, R, fx, p, p);
        ASSERT_ALWAYS(R->deg >= 0);
        ASSERT_ALWAYS((unsigned int) R->deg < d);
        for(unsigned int j = 0 ; j <= (unsigned int) R->deg ; j++)
            mpz_set(mpz_mat_entry(frobenius_matrix, i, j), R->_coeff[j]);
    }
    mpz_mat_inv_mod_mpz(inverse_frobenius_matrix, frobenius_matrix, p);
}

void cxx_mpz_poly_bivariate::reducer_mod_fx_mod_mpz::frobenius(
        cxx_mpz_poly & B,
        cxx_mpz_poly const & A,
        int order) const
{
    if (A == 0) {
        B = A;
        return;
    }

    /* It's just a linear operation, all we care about is gettint the
     * correct (precomputed) matrix */

    mpz_mat_ptr mat;
    compute_frobenius_matrices();
    if (order == 1) {
        mat = frobenius_matrix;
    } else if (order == -1) {
        mat = inverse_frobenius_matrix;
    } else {
        /* for the moment the others are unimplemented. Of course we could at
         * least reduce the order mod the degree of fx (thereby asserting
         * that fx is irreducible), but we have no strong use case for such
         * an interface at the moment. */
        ASSERT_ALWAYS(0);
    }

    /* there might be situations where the input isn't reduced.
     */
    (*this)(B, A);

    cxx_mpz_poly h;
    /* It's weird, but we don't have row times matrix in mpz_mat...  */
    for(int j = 0 ; j < fx->deg ; j++) {
        cxx_mpz s;
        ASSERT_ALWAYS(B->deg < fx->deg);
        for(int i = 0 ; i <= B->deg ; i++) {
            mpz_addmul(s, B->_coeff[i], mpz_mat_entry_const(mat, i, j));
        }
        mpz_mod(s, s, p);
        mpz_poly_setcoeff(h, j, s);
    }
    B = h;
}

void cxx_mpz_poly_bivariate::reducer_mod_fx_mod_mpz::frobenius(
        cxx_mpz_poly_bivariate & B,
        cxx_mpz_poly_bivariate const & A,
        int order) const
{
    ASSERT_ALWAYS(mpz_size(p) == 1);
    mp_limb_t pp = mpz_get_ui(p);

    /* for the moment the others are unimplemented */
    ASSERT_ALWAYS(order == 1 || order == -1);

    B = A;
    int n = B.degree_y();

    if (n == -1)
        return;

    if (order == -1) {
        /* scale down powers for inverse Frobenius */
        ASSERT_ALWAYS(n % pp == 0);
        n = n / pp;
        for(int k = 0 ; k <= n ; k++) {
            if (k < n) {
                for(mp_limb_t l = 1 ; l < pp ; l++)
                    ASSERT_ALWAYS(B[k * pp + l] == 0);
            }
            mpz_poly_swap(((super&)B)[k], ((super&)B)[k*pp]);
        }
        B.erase(B.begin() + n + 1, B.end());
    }

    for(int k = 0 ; k <= n ; k++)
        frobenius(((super&)B)[k], ((super&)B)[k], order);

    if (order == 1) {
        /* scale up powers for forward Frobenius */

        B.insert(B.end(), n * (pp - 1), cxx_mpz_poly());
        for(int k = n ; k >= 0 ; k--)
            mpz_poly_swap(((super&)B)[k * pp], ((super&)B)[k]);
        B.cleandeg(n * pp);
        ASSERT_ALWAYS(B.degree_y() >= 0);
        ASSERT_ALWAYS((unsigned int) B.degree_y() == n * pp);
    }
}

/* powering. Both the generic template, and the instantiations that go
 * through the reduction operator */

template <typename T> struct pow_ui_inner { /*{{{*/
    T const & reducer;
    pow_ui_inner(T const & reducer)
        : reducer(reducer)
    {
    }

  private:
    void wrapped(cxx_mpz_poly_bivariate & B, cxx_mpz_poly_bivariate const & A,
                 unsigned long n) const
    {
        ASSERT_ALWAYS(&B != &A);
        unsigned long k = ((~0UL) >> 1) + 1;
        for (; k > n; k >>= 1)
            ;
        B = A;
        for (; k >>= 1;) {
            cxx_mpz_poly_bivariate::mul(B, B, B);
            if (n & k)
                cxx_mpz_poly_bivariate::mul(B, B, A);
            reducer(B, B);
        }
    }

  public:
    void operator()(cxx_mpz_poly_bivariate & B,
                    cxx_mpz_poly_bivariate const & A, unsigned long n) const
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
            wrapped(C, A, n);
            B.swap(C);
        } else {
            wrapped(B, A, n);
        }
    }
}; /*}}}*/

template <typename T>
inline void cxx_mpz_poly_bivariate::pow_ui(cxx_mpz_poly_bivariate & B,
                                           cxx_mpz_poly_bivariate const & A,
                                           unsigned long n,
                                           T const & reducer) /*{{{*/
{
    pow_ui_inner {reducer}(B, A, n);
}
/*}}}*/

template void
cxx_mpz_poly_bivariate::pow_ui<cxx_mpz_poly_bivariate::reducer_noop>(
    cxx_mpz_poly_bivariate & B, cxx_mpz_poly_bivariate const &, unsigned long,
    reducer_noop const &);

/* division */
template <typename T>
void cxx_mpz_poly_bivariate::div_qr(cxx_mpz_poly_bivariate & q,
                                    cxx_mpz_poly_bivariate & r,
                                    cxx_mpz_poly_bivariate const & f,
                                    cxx_mpz_poly_bivariate const & g,
                                    T const & reducer) /*{{{*/
{
    ASSERT_ALWAYS(g.degree() >= 0 && g.lc() == 1);
    ASSERT_ALWAYS(&q != &f && &q != &r && &q != &g);
    ASSERT_ALWAYS(&r != &g);
    /* r == f is allowed */

    int j, df = f.degree(), dg = g.degree(), dq = df - dg;

    if (df < dg) /* f is already reduced mod g */
    {
        q = 0;
        r = f;
        return;
    }

    /* now df >= dg */

    r = f;
    q.assign(dq + 1, 0);
    cxx_mpz_poly tmp;

    for (int k ; (k = r.degree() - dg) >= 0 ; ) {
        ((super &)q)[k] = r[k + dg];
        for (j = dg + k - 1; j >= k; j--) {
            mpz_poly_mul(tmp, q[k], g[j - k]);
            mpz_poly_sub(((super &)r)[j], r[j], tmp);
        }
        r.cleandeg(k + dg - 1);
        reducer(r, r);
    }
}
/*}}}*/

void cxx_mpz_poly_bivariate::mul_mod_fx_mod_fy_mod_mpz(self & C, self const & A,
                                                       self const & B,
                                                       mpz_poly_srcptr fx,
                                                       mpz_poly_srcptr fy,
                                                       mpz_srcptr p)
{
    mul(C, A, B);
    reducer_mod_fy_mod_fx_mod_mpz {fy, fx, p}(C, C);
}

/* gcd. Well it's only about being able to run the Euclidean algorithm,
 * so we're talking about univariate polynomials over finite fields in
 * disguise. So we need a reducer, and that reducer has to be both mod a
 * polynomial and mod a prime. No point in defining gcd with more general
 * reducers, that would not make much sense.
 */
void cxx_mpz_poly_bivariate::gcd(self & C, self const & A, self const & B,
                                 reducer_mod_fx_mod_mpz const & Rx) /*{{{*/
{
    self R0 = A;
    self R1 = B;
    self R2;
    for (; R1.degree_y() >= 0;) {
        Rx.make_monic(R1);
        mod_fy(R2, R0, R1);
        std::swap(R0, R1);
        Rx(R1, R2);
    }
    C = R0;
}
/*}}}*/

bool cxx_mpz_poly_bivariate::reducer_mod_fx_mod_mpz::make_monic(cxx_mpz_poly_bivariate & f) const
{
    if (f.lc() == 1)
        return true;
    cxx_mpz_poly u,v, d;
    mpz_poly_xgcd_mpz(d, f.lc(), fx, u, v, p);
    cxx_mpz_poly_bivariate::mul(f, f, u);
    (*this)(f, f);
    return d == 1;
}

int sqf_inner(cxx_mpz_poly_bivariate::factor_list & fl,
              cxx_mpz_poly_bivariate const & f, int stride,
              cxx_mpz_poly_bivariate::reducer_mod_fx_mod_mpz const & Rx)

{
    int r = 0;

    typedef cxx_mpz_poly_bivariate self;

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
        self::pow_ui(t0, mi, i, Rx);
        self::divexact(t1, f, t0, Rx);
        /* t1 = all polynomials in mi taken out from f with multiplicity i */
        self::gcd(mi1, t1, mi, Rx);
        /* mi1 = almost like mi, but since factors with multiplicity i
         * are no longer in t1, there's absent from mi1 too. Whence
         * mi/mi1 is exactly the product of factors of multiplicity 1.
         */
        if ((size_t)(i * stride) >= fl.size()) {
            fl.insert(fl.end(), i * stride + 1 - fl.size(), {});
        }
        /* Use tmp so that we don't absurdly keep storage within
         * lf->factors */
        self::divexact(tmp, mi, mi1, Rx);
        fl[i * stride] = {tmp, 1}; /* multiplicity field in fact unused */
        self::pow_ui(t0, fl[i * stride].first, i, Rx);
        self::mul(T, T, t0);
        Rx(T, T);
        mi.swap(mi1);
        r = i * stride;
    }

    self::divexact(fl[0].first, f, T, Rx);
    fl[0].second = 1;

    return r;
}
