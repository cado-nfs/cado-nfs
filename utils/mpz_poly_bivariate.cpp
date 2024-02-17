#include "cado.h"
#include "mpz_poly_bivariate.hpp"
#include "cado_expression_parser.hpp"
#include <sstream>
#include <string>

/* Polynomial arithmetic */
void cxx_mpz_poly_bivariate::neg(cxx_mpz_poly_bivariate & f, cxx_mpz_poly_bivariate const & g)
{
    for(size_t i = 0 ; i < g.size() ; i++) {
        mpz_poly_neg(((super&)f)[i], ((super&)f)[i]);
    }
}

void cxx_mpz_poly_bivariate::add(cxx_mpz_poly_bivariate & f, cxx_mpz_poly_bivariate const & g, cxx_mpz_poly_bivariate const & h)
{
    size_t sg = g.size();
    size_t sh = h.size();
    size_t sf = std::max(sg, sh);
    f.reserve(sf);
    size_t i = 0;
    if (&f != &g && &f != &h)
        f.clear();
    f.insert(f.end(), sf - f.size(), cxx_mpz_poly());   /* now sf == f.size() */
    for( ; i < sg && i < sh ; i++)
        mpz_poly_add(((super&)f)[i], g[i], h[i]);
    for( ; i < sg ; i++)
        mpz_poly_set(((super&)f)[i], g[i]);
    for( ; i < sh ; i++)
        mpz_poly_set(((super&)f)[i], h[i]);
    f.cleandeg((int) sf - 1);
}

void cxx_mpz_poly_bivariate::sub(cxx_mpz_poly_bivariate & f, cxx_mpz_poly_bivariate const & g, cxx_mpz_poly_bivariate const & h)
{
    size_t sg = g.size();
    size_t sh = h.size();
    size_t sf = std::max(sg, sh);
    f.reserve(sf);
    size_t i = 0;
    if (&f != &g && &f != &h)
        f.clear();
    f.insert(f.end(), sf - f.size(), cxx_mpz_poly());   /* now sf == f.size() */
    for( ; i < sg && i < sh ; i++)
        mpz_poly_sub(((super&)f)[i], g[i], h[i]);
    for( ; i < sg ; i++)
        mpz_poly_set(((super&)f)[i], g[i]);
    for( ; i < sh ; i++)
        mpz_poly_neg(((super&)f)[i], h[i]);
    f.cleandeg((int) sf - 1);
}

/*
void cxx_mpz_poly_bivariate::add(cxx_mpz_poly_bivariate & f, cxx_mpz_poly_bivariate const & g, cxx_mpz_poly_bivariate::lifted_x const & a)
{
    f = g;
    if (f == 0) {
        f = a;
}

*/

void cxx_mpz_poly_bivariate::mul(cxx_mpz_poly_bivariate & f, cxx_mpz_poly_bivariate const & g, cxx_mpz_poly_bivariate const & h)
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
    for(int i = 0 ; i <= g.degree() ; i++) {
        for(int j = 0 ; j <= h.degree() ; j++) {
            mpz_poly_mul(tmp, g[i], h[j]);
            mpz_poly_add(((super&)f)[i+j], f[i+j], tmp);
        }
    }
    f.cleandeg(g.degree() + h.degree());
}

/* B = A^n */
void cxx_mpz_poly_bivariate::pow_ui(cxx_mpz_poly_bivariate & B, cxx_mpz_poly_bivariate const & A, unsigned long n)/*{{{*/
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

    unsigned long k = ((~0UL)>>1) + 1;
    for( ; k > n ; k >>= 1);
    B = A;
    for( ; k >>= 1 ; ) {
        mul(B, B, B);
        if (n & k)
            mul(B, B, A);
    }
}/*}}}*/


struct mpz_poly_bivariate_parser_traits {
    std::string x, y;
    mpz_poly_bivariate_parser_traits(std::string const & x, std::string const & y)
        : x(x), y(y)
    {}

    struct parse_error : public std::exception {};
    struct unexpected_literal : public parse_error {
        std::string msg;
        unexpected_literal(std::string const & v)
            : msg(std::string("unexpected literal " + v))
        {}
        const char *what() const noexcept override {
            return msg.c_str();
        }
    };

    static constexpr const int accept_literals = 2;
    typedef cxx_mpz_poly_bivariate type;
    void add(cxx_mpz_poly_bivariate & c, cxx_mpz_poly_bivariate const & a, cxx_mpz_poly_bivariate const & b) {
        type::add(c, a, b);
    }
    void sub(cxx_mpz_poly_bivariate & c, cxx_mpz_poly_bivariate const & a, cxx_mpz_poly_bivariate const & b) {
        type::sub(c, a, b);
    }
    void mul(cxx_mpz_poly_bivariate & c, cxx_mpz_poly_bivariate const & a, cxx_mpz_poly_bivariate const & b) {
        type::mul(c, a, b);
    }
    void pow_ui(cxx_mpz_poly_bivariate & c, cxx_mpz_poly_bivariate const & a, unsigned long e) {
        type::pow_ui(c, a, e);
    }
    void swap(cxx_mpz_poly_bivariate & a, cxx_mpz_poly_bivariate & b) {
        a.swap(b);
    }
    void set_mpz(cxx_mpz_poly_bivariate & a, cxx_mpz const & z) {
        a = z;
    }
    /* TODO do something for variable names */
    void set_literal_power(cxx_mpz_poly_bivariate & a, std::string const & v, unsigned long e) {
        if (v == x) {
            a.set_xi(e);
        } else if (v == y) {
            a.set_yi(e);
        } else {
            throw unexpected_literal(v);
        }
    }
};

std::istream& operator>>(std::istream& in, cxx_mpz_poly_bivariate::named_proxy<cxx_mpz_poly_bivariate &> F)
{
    std::string const & x = F.x;
    std::string const & y = F.y;
    cxx_mpz_poly_bivariate & f = F.c;

    f = 0;

    std::string line;
    for(;;in.get()) {
        int c = in.peek();
        if (in.eof() || !isspace(c)) break;
    }
    if (!getline(in, line)) return in;
    std::istringstream is(line);

    typedef cado_expression_parser<mpz_poly_bivariate_parser_traits> poly_parser;
    poly_parser P(x, y);

    P.tokenize(is);

    try {
        f = P.parse();
    } catch (poly_parser::parse_error const & p) {
        in.setstate(std::ios_base::failbit);
        return in;
    }

    return in;
}

std::ostream& operator<<(std::ostream& o, cxx_mpz_poly_bivariate::named_proxy<cxx_mpz_poly_bivariate const &> F) {
    std::string const & x = F.x;
    std::string const & y = F.y;
    cxx_mpz_poly_bivariate const & f = F.c;
    if (f.degree() == -1) {
        return o << 0;
    } else if (f.degree() == 0) {
        return o << f[0].print_poly(x);
    }
    std::ostringstream os;
    for(int i = 0 ; i <= f.degree() ; i++) {
        if (f[i]->deg < 0) continue;
        if (f[i] == 1) {
            if (os.str().size())
                os << "+";
            if (i == 0) os << "1";
        } else if (f[i] == -1) {
            os << "-";
            if (i == 0) os << "1";
        } else {
            if (!mpz_poly_is_monomial_multiple(f[i])) {
                /* need parentheses anyway */
                if (os.str().size())
                    os << "+";
                os << "(" << f[i].print_poly(x) << ")";
            } else {
                /* if it's a negative multiple of x^i, don't print "+-".
                */
                if (os.str().size() && mpz_sgn(mpz_poly_lc(f[i])) > 0)
                    os << "+";
                os << f[i].print_poly(x);
            }
            if (i > 0) os << "*";
        }
        if (i > 0) os << y;
        if (i > 1) os << "^" << i;
    }
    return o << os.str();
}

void cxx_mpz_poly_bivariate::mod_mpz(cxx_mpz_poly_bivariate & a, cxx_mpz_poly_bivariate const & b, mpz_srcptr p)
{
    a = b;
    mpz_ptr invmptr = NULL;
    cxx_mpz invm;
    if (a.degree() >= 4) {
        barrett_precompute_inverse(invm, p);
        invmptr = invm;
    }
    for(auto & c : a)
        mpz_poly_mod_mpz(c, c, p, invmptr);
}

void cxx_mpz_poly_bivariate::mod_fx(cxx_mpz_poly_bivariate & a, cxx_mpz_poly_bivariate const & b, mpz_poly_srcptr fx)
{
    a = b;
    ASSERT_ALWAYS(mpz_poly_is_monic(fx));
    for(auto & c : a)
        mpz_poly_div_r(c, c, fx);
}

void cxx_mpz_poly_bivariate::div_qr(cxx_mpz_poly_bivariate & q, cxx_mpz_poly_bivariate & r, cxx_mpz_poly_bivariate const & f, cxx_mpz_poly_bivariate const & g)
{
    ASSERT_ALWAYS(g.degree() >= 0 && g.lc() == 1);
    ASSERT_ALWAYS(&q != &f && &q != &r && &q != &g);
    ASSERT_ALWAYS(&r != &g);
    /* r == f is allowed */

    int k, j, df = f.degree(), dg = g.degree(), dq = df - dg;

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

    for (k = df-dg ; k >=0 ; k--) {
        ((super&)q)[k] = r[k+dg];
        for (j = dg+k-1 ; j >= k ; j--) {
            mpz_poly_mul(tmp, q[k], g[j-k]);
            mpz_poly_sub(((super&)r)[j], r[j], tmp);
        }
    }
    r.cleandeg(dg - 1);
}

void cxx_mpz_poly_bivariate::mod_fy(cxx_mpz_poly_bivariate & a, cxx_mpz_poly_bivariate const & b, cxx_mpz_poly_bivariate const & fy)
{
    cxx_mpz_poly_bivariate q;
    div_qr(q, a, b, fy);
}

void cxx_mpz_poly_bivariate::mul(cxx_mpz_poly_bivariate & a, cxx_mpz_poly_bivariate const & b, mpz_poly_srcptr m)
{
    a = b;
    for(auto & c : a)
        mpz_poly_mul(c.x, c.x, m);
}

void cxx_mpz_poly_bivariate::mul(cxx_mpz_poly_bivariate & a, cxx_mpz_poly_bivariate const & b, mpz_srcptr m)
{
    a = b;
    for(auto & c : a)
        mpz_poly_mul_mpz(c, c, m);
}

void cxx_mpz_poly_bivariate::eval_fy(cxx_mpz_poly & a, self const & f, cxx_mpz_poly const & e)
{
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
    for(int i = f.degree() - 1 ; i >= 0 ; i--) {
        mpz_poly_mul(a, a, e);
        mpz_poly_add(a, a, f[i]);
    }
}

void cxx_mpz_poly_bivariate::eval_fx(cxx_mpz_poly & a, self const & f, mpz_srcptr e)
{
    mpz_poly_realloc(a, f.degree() + 1);
    for(size_t i = 0 ; i < f.size() ; i++)
        mpz_poly_eval(a->coeff[i], f[i], e);
    mpz_poly_cleandeg(a, f.degree());
}

void cxx_mpz_poly_bivariate::transpose(self & a, self && b)
{
    if (&a == &b) {
        self bb = b;
        transpose(a, bb);
        return;
    }
    int dx = b.degree_x();
    int dy = b.degree_y();
    a.assign(dx + 1, 0);
    for(auto & c : a) {
        mpz_poly_realloc(c, dy + 1);
        for(int j = 0 ; j <= dy ; j++)
            mpz_set_ui(c->coeff[j], 0);
    }
    for(int j = 0 ; j <= dy ; j++) {
        for(int i = 0 ; i <= mpz_poly_degree(((super const &)b)[j]) ; i++) {
            mpz_swap(((super&)a)[i]->coeff[j], ((super&)b)[j]->coeff[i]);
        }
    }
    for(int i = 0 ; i <= dx ; i++) {
        mpz_poly_cleandeg(((super&)a)[i], dy);
    }
    a.cleandeg(dx);
    b.clear();
}
void cxx_mpz_poly_bivariate::transpose(self & a, self const & b)
{
    self bb = b;
    transpose(a, (self &&) bb);
}


void cxx_mpz_poly_bivariate::resultant_y(
        cxx_mpz_poly & resultant,
        self const & f,
        self const & g)
{
    size_t nb_eval = f.degree_y() * g.degree_x() + f.degree_x() * g.degree_y() + 1;

    std::vector<cxx_mpz> points;
    std::vector<cxx_mpz> evaluations;

    cxx_mpz_poly ef, eg;
    cxx_mpz z;

    for(int i = 0 ; points.size() < nb_eval ; i++) {
        int w = i ? ((i & 1) ? (i + 1) / 2 : -(i / 2)) : 0;
        eval_fx(ef, f, cxx_mpz(w));
        eval_fx(eg, g, cxx_mpz(w));

        if (ef.degree() < f.degree_y()) continue;
        if (eg.degree() < g.degree_y()) continue;

        mpz_poly_resultant(z, ef, eg);
        points.push_back(w);
        evaluations.push_back(z);
    }

    int ok = mpz_poly_interpolate(resultant, points, evaluations);
    if (!ok)
        resultant = 0;
}

void cxx_mpz_poly_bivariate::resultant_x(
        cxx_mpz_poly & resultant,
        self const & f,
        self const & g)
{
    size_t nb_eval = f.degree_y() * g.degree_x() + f.degree_x() * g.degree_y() + 1;

    std::vector<cxx_mpz> points;
    std::vector<cxx_mpz> evaluations;

    cxx_mpz_poly ef, eg;
    cxx_mpz z;

    for(int i = 0 ; points.size() < nb_eval ; i++) {
        int w = i ? ((i & 1) ? (i + 1) / 2 : -(i / 2)) : 0;
        /* simple enough, really. Note that we could also choose
         * polynomials in x as interpolation points, that would work.
         */
        eval_fy(ef, f, cxx_mpz_poly(w));
        eval_fy(eg, g, cxx_mpz_poly(w));

        if (ef.degree() < f.degree_x()) continue;
        if (eg.degree() < g.degree_x()) continue;

        mpz_poly_resultant(z, ef, eg);
        points.push_back(w);
        evaluations.push_back(z);
    }

    int ok = mpz_poly_interpolate(resultant, points, evaluations);
    if (!ok)
        resultant = 0;
}


void cxx_mpz_poly_bivariate::set_rrandomb(self & f, int dx, int dy, int bits, gmp_randstate_ptr rstate)
{
    f.assign(dy + 1, 0);
    for(int i = 0 ; i <= dy ; i++) {
        mpz_poly_set_rrandomb(((super&)f)[i], dx, bits, rstate);
    }
    f.cleandeg(dy);
}

void cxx_mpz_poly_bivariate::set_rrandomb_cab(self & f, int dx, int dy, int bits, gmp_randstate_ptr rstate)
{
    // dx and dy must be coprime, but we don't check it.
    // as far as the "Cab" notation goes, dy is a and dx is b.
    set_rrandomb(f, dx, dy, bits, rstate);
    // X^i Y^j appears only if i * dy + j * dx < dx * dy
    ((super&)f)[dy] = 1;
    for(int j = 0 ; j < dy ; j++) {
        // we may have X^i only if i * dy < (dx * dy  - j * dx)
        int d = (dx * (dy - j) - 1) / dy;
        mpz_poly_cleandeg(((super&)f)[j], d);
        if (j == 0) {
            mpz_set_ui(((super&)f)[j]->coeff[dx], 1);
            mpz_poly_cleandeg(((super&)f)[j], dx);
        }
    }
}

