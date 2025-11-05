#include "cado.h"       // IWYU pragma: keep

#include <string>
#include <sstream>
#include <tuple>
#include <utility>

#include <gmp.h>

#include "mpz_poly.h"
#include "mpz_mat.h"
#include "cxx_mpz.hpp"
#include "arithmetic_reductions.hpp"
#include "mpz_poly_bivariate.hpp"
#include "macros.h"
#include "matrix.hpp"

namespace cado::arithmetic_reductions {

/* reduction operator inside function objects */
void mod_mpz::operator()(
    cxx_mpz_poly_bivariate & b, cxx_mpz_poly_bivariate const & a) const /*{{{*/
{
    mpz_srcptr inv = invp_lazy();
    b = a;
    for (auto & c: b.v())
        mpz_poly_mod_mpz(c, c, p, inv);
    b.cleandeg(b.degree());
}
/*}}}*/
void mod_fx::operator()(
    cxx_mpz_poly_bivariate & B, cxx_mpz_poly_bivariate const & A) const /*{{{*/
{
    cxx_mpz_poly_bivariate::mod_fx(B, A, fx);
}
/*}}}*/
void mod_fx_mod_p::operator()(
    cxx_mpz_poly_bivariate & B, cxx_mpz_poly_bivariate const & A) const /*{{{*/
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
    mod_mpz::operator()(B,A);
#if 0
    cxx_mpz_poly_bivariate::mod_fx(B, B, fx);
    mod_mpz::operator()(B,B);
#else
    cxx_mpz_poly q;
    for(auto & c : B.v())
        mpz_poly_div_qr_mod_mpz(q, c, c, fx, p);
    B.cleandeg(B.degree());
#endif
}
/*}}}*/

void mod_mpz::operator()(
    cxx_mpz_poly & b, cxx_mpz_poly const & a) const /*{{{*/
{
    mpz_srcptr inv = invp_lazy();
    mpz_poly_mod_mpz(b, a, p, inv);
}
/*}}}*/
void mod_mpz::operator()(
    cado::matrix<cxx_mpz> & b, cado::matrix<cxx_mpz> const & a) const /*{{{*/
{
    if (&b != &a) {
        b.resize(a.m, a.n);
    }
    for(unsigned int i = 0 ; i < a.m * a.n ; i++)
        mpz_mod(b.coeffs[i], a.coeffs[i], p);
}
/*}}}*/
void mod_mpz::operator()(
    cado::matrix<cxx_mpz_poly> & b, cado::matrix<cxx_mpz_poly> const & a) const /*{{{*/
{
    if (&b != &a) {
        b.resize(a.m, a.n);
    }
    for(unsigned int i = 0 ; i < a.m * a.n ; i++) {
        (*this)(b.coeffs[i], a.coeffs[i]);
    }
}
/*}}}*/
void mod_fx::operator()(
    cxx_mpz_poly & B, cxx_mpz_poly const & A) const /*{{{*/
{
    mpz_poly_div_r(B, A, fx);
}
/*}}}*/
void mod_fx_mod_p::operator()(
    cxx_mpz_poly & B, cxx_mpz_poly const & A) const /*{{{*/
{
    mod_mpz::operator()(B,A);
    // mpz_poly_div_r_mod_mpz_clobber(B, fx, p);
    cxx_mpz_poly Q;
    mpz_poly_div_qr_mod_mpz(Q, B, B, fx, p);
    mod_mpz::operator()(B,B);
}
/*}}}*/
void mod_fx_mod_p::operator()(
    cado::matrix<cxx_mpz_poly> & b, cado::matrix<cxx_mpz_poly> const & a) const /*{{{*/
{
    if (&b != &a) {
        b.resize(a.m, a.n);
    }
    for(unsigned int i = 0 ; i < a.m * a.n ; i++)
        (*this)(b.coeffs[i], a.coeffs[i]);
}
/*}}}*/

void mod_fy_mod_q::operator()(
    cxx_mpz_poly_bivariate & B, cxx_mpz_poly_bivariate const & A) const /*{{{*/
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
    mod_q const & R0 = *this;
    R0(B,A);
    // cxx_mpz_poly_bivariate::mod_fy(B, B, fy);
    cxx_mpz_poly_bivariate::div_r(B, B, fy, R0);
    R0(B,B);
}
/*}}}*/

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
std::string noop::print() const
{
    return "dummy no-reduction";
}

std::string mod_mpz::print() const
{
    std::ostringstream os;
    os << "reduction mod " << p;
    return os.str();
}

std::string mod_fx::print() const
{
    std::ostringstream os;
    os << "reduction mod " << fx.named("x");
    return os.str();
}

std::string mod_fx_mod_p::print() const
{
    std::ostringstream os;
    os << "reduction mod " << fx.named("x") << ", mod " << p;
    return os.str();
}

std::string mod_fy_mod_q::print() const
{
    std::ostringstream os;
    os << "reduction mod " << fy.named("x", "y") << ", mod " << fx.named("x")
       << ", mod " << p;
    return os.str();
}

bool mod_fx_mod_p::make_monic(cxx_mpz_poly_bivariate & f) const
{
    if (f.lc() == 1)
        return true;
    cxx_mpz_poly u, v, d;
    mpz_poly_xgcd_mpz(d, f.lc(), fx, u, v, p);
    cxx_mpz_poly_bivariate::mul(f, f, u);
    (*this)(f, f);
    return d == 1;
}

/* f being a polynomial over Z/pZ of degree d, compute the matrix of the
 * p-th powering linear map modulo f. It's a d*d matrix.
 *
 * Cost: one evaluation to the p-th power, followed by d-1 products, so
 * in total O(log(p)+d) products of polynomials mod p, each of which costs
 * M(d) operations mod p. ---> O(d + log(p)) * M(d)
 */
cado::matrix<cxx_mpz> mod_q::compute_frobenius_matrix() const
{
    int const d = fx->deg;
    ASSERT_ALWAYS(d >= 0);
    cado::matrix<cxx_mpz> M(d, d);
    if (d == 0)
        return M;
    cxx_mpz_poly xp;
    mpz_poly_set_xi(xp, 1);
    mpz_poly_pow_mod_f_mod_mpz(xp, xp, fx, p, p);

    cxx_mpz_poly R = 1;
    for(int i = 0 ; i < d ; i++) {
        ASSERT_ALWAYS(R->deg >= 0);
        ASSERT_ALWAYS(R->deg < d);
        for(int j = 0 ; j <= R->deg ; j++)
            M(i, j) = R->_coeff[j];
        if (i + 1 < d)
            mpz_poly_mul_mod_f_mod_mpz(R, R, xp, fx, p, nullptr, invp_lazy());
    }
    return M;
}

/* f is a polynomial of degree n-1 over Fq, where Fq is represented by
 * polynomials over Z/pZ modulo fx of degree d (thus q=p^d). The p-th
 * powering map (call it F) is a skew linear map on polynomials, since
 * for a\in Fq and u a polynomial in Fq[x] mod f, we have have F(a*u) =
 * a^p F(u). Nevertheless, we're interested in how the the basis
 * (1,x,...,x^(n-1)) is transformed by F. This will be a matrix, but it
 * must *NOT* be construed as the matrix of anything super useful.
 * However it will be of interest for the computation of the frobenius
 * map on polynomials of Fq[x].
 *
 * A special case is when f is defined over Fp. Then the matrix that we
 * compute here is a matrix over Fp, which actually matches the Frobenius
 * matrix of Fp[x]/f(x).
 *
 * The cost of the computation is one p-th power mod f, followed by n-1
 * products. In total, it means O(log(p) + n) products of polynomials mod
 * f, so O(log(p) + n) * M(n) * M(d)
 *
 */
cado::matrix<cxx_mpz_poly> mod_fy_mod_q::compute_small_frobenius_matrix() const
{
    int const n = fy.degree_y();
    ASSERT_ALWAYS(n >= 0);
    cado::matrix<cxx_mpz_poly> M(n, n);
    if (n == 0)
        return M;

    if (fy.degree_x() <= 0) {
        cxx_mpz_poly_bivariate g;
        cxx_mpz_poly_bivariate::transpose(g, cxx_mpz_poly_bivariate(fy));
        auto M0 = mod_q { g[0], p }.compute_frobenius_matrix();
        for(int i = 0 ; i < n ; i++)
            for(int j = 0 ; j < n ; j++)
                M(i, j) = M0(i, j);
    } else {
        auto yp = cxx_mpz_poly_bivariate::yi(1);
        cxx_mpz_poly_bivariate::pow(yp, yp, p, *this);
        cxx_mpz_poly_bivariate R = 1;

        for(int i = 0 ; i < n ; i++) {
            for(int j = 0 ; j <= R.degree_y() ; j++)
                M(i, j) = R[j];
            if (i + 1 < n) {
                cxx_mpz_poly_bivariate::mul(R, R, yp);
                (*this)(R, R);
            }
        }
    }
    return M;
}

/* Here we compute the matrix of the Fq-linear map Q that takes a
 * polynomial in Fq[x] mod f to its q-th power. f has degree n, and q is
 * p^d.
 *
 * We assume that the following data has been computed already.
 *
 * - F (compute_small_frobenius_matrix) costs O(n + log(p)) * M(n) * M(d)
 *   F is represented by M_1 (n*n matrix with coefficients in Fq)
 * - sigma (mod_q::compute_frobenius_matrix) costs O(d + log(p)) * M(d)
 *   sigma is represented by M_0 (d*d matrix with coefficients in Fp)
 *
 * We want to compute the matrix M that represents Q (n*n matrix with
 * coefficients in Fq).  We have Q = F^d, however because F is just the
 * p-th power map, which is a _skew_ linear map, we do not have M =
 * M_1^d.
 *
 * Instead, look at F over an F_p basis. Since F(a*u) = a^p*F(u), a
 * matrix representation of F in dimension n*d can be written as
 * (I⊗M_0) * M_1, where each coefficient of M_1 is implicitly expanded to 
 * the d*d matrix that represents the multiplication by it. Now the
 * matrix of F is ((I⊗M_0) * M_1)^d.
 *
 * Skew-linearity can also be expressed as M_1 * (I⊗M_0) = (I⊗M_0) * M_1^\sigma.
 * In other words, when we restrict to what happens for one coefficient
 * of Fq if x and c are in Fq, Mc is the multiplication-by-c matrix, and
 * elements of Fq are identified to row vectors in Fp^d, we have
 * x*M_c*M_0 = sigma(x*M_c) = sigma(x*c) = sigma(x)*sigma(c) =
 * x*M_0*M_{c^\sigma}.
 *
 * A consequence is that the matrix of F is
 * (I⊗M_0)^d * M_1^{\sigma^{d-1}} * ... * M_1^\sigma * M_1
 *
 * and (I⊗M_0)^d is actually the identity matrix. Therefore, we merely
 * have to compute the "matrix Galois norm" above.
 *
 * We do this by divide-and-conquer. To compute X^(\sigma^i) for a matrix
 * X, we need the d*d matrix that represents \sigma^i.
 *
 * For k a power of two, the algorithm A that returns
 * M_1^{\sigma^{k-1}} * ... * M_1^\sigma * M_1 as well as the matrix of
 * \sigma^k goes recursively as:
 *
 * def A(k):
 *     if k==1:
 *         return M_1, M_0
 *     elif k % 2 == 0:  # always true if the input k is a power of two
 *         X, H = A(k//2)
 *         Y = apply H to all n^2 entries of X
 *         return Y*X, H^2
 * The odd case generalization is simply a third branch in the above code:
 *     elif k % 2 == 1:
 *         X, H = A(k-1)
 *         Y = apply M_0 to all n^2 entries of X
 *         return Y * M_1, H * M_0
 *
 * The cost to compute the matrix of F is thus
 *  - O(log(d)) products of n*n matrices over Fq, 
 *  - O(n^2*log(d)) applications of powers of sigma, which are
 *    vector-times-matrix products in dimension d (thus similar to M(d))
 *  - O(log(d)) multiplications of d*d matrices.
 * So in total the first cost very likely dominates and we should have
 * something like O(log(d) * max(n^3, d) * M(d)) operations over Fp.
 */

std::tuple<cado::matrix<cxx_mpz_poly>, cado::matrix<cxx_mpz>> mod_fy_mod_q::matrix_galois_norm(int k) const
{
    unsigned int const n = fy.degree();
    if (k == 1) {
        return { small_frobenius_matrix(), mod_q::frobenius_matrix() };
    } else if ((k % 2) == 0) {
        auto [ X, H ] = matrix_galois_norm(k/2);
        cado::matrix<cxx_mpz_poly> Y(n, n);
        for(unsigned int i = 0 ; i < n ; i++)
            for(unsigned int j = 0 ; j < n ; j++)
                mod_q::matrix_apply(Y(i, j), X(i, j), H);
        return { Y * X, H * H };
    } else {
        auto [ X, H ] = matrix_galois_norm(k - 1);
        cado::matrix<cxx_mpz_poly> Y(n, n);
        for(unsigned int i = 0 ; i < n ; i++)
            for(unsigned int j = 0 ; j < n ; j++)
                mod_q::frobenius(Y(i, j), X(i, j), 1);
        return { Y * small_frobenius_matrix(), H * mod_q::frobenius_matrix() };
    }
}

cado::matrix<cxx_mpz_poly> mod_fy_mod_q::compute_frobenius_matrix() const
{
    auto [ Y, H ] = matrix_galois_norm(mod_q::fx.degree());
    return Y;
}

cado::matrix<cxx_mpz> mod_q::compute_inverse_frobenius_matrix() const
{
    cxx_mpz_mat Zi;
    mpz_mat_inv_mod_mpz(Zi, cxx_mpz_mat(frobenius_matrix()), p);
    return cado::matrix<cxx_mpz>(std::move(Zi));
}

void mod_q::matrix_apply(
        cxx_mpz_poly & B,
        cxx_mpz_poly const & A,
        cado::matrix<cxx_mpz> const & mat) const
{
    /* there might be situations where the input isn't reduced.
     */
    (*this)(B, A);

    cxx_mpz_poly h;
    /* It's weird, but we don't have row times matrix in mpz_mat...  */
    for(int j = 0 ; j < fx->deg ; j++) {
        cxx_mpz s;
        ASSERT_ALWAYS(B->deg < fx->deg);
        for(int i = 0 ; i <= B->deg ; i++) {
            mpz_addmul(s, B->_coeff[i], mat(i, j));
        }
        mpz_mod(s, s, p);
        mpz_poly_setcoeff(h, j, s);
    }
    B = h;
}

void mod_q::frobenius(
        cxx_mpz_poly & B,
        cxx_mpz_poly const & A,
        int order) const
{
    if (A == 0) {
        B = A;
        return;
    }

    /* It's just a linear operation, all we care about is getting the
     * correct (precomputed) matrix */

    if (order == 1) {
        matrix_apply(B, A, frobenius_matrix());
    } else if (order == -1) {
        matrix_apply(B, A, inverse_frobenius_matrix());
    }
}

/*
 * Given A(x) = \sum a_i y^i with a_i in Fp[x]/f(x),
 * return the polynomial B^p = \sum (a_i^p) y^(pi) (if order==1). With
 * order=-1, do the reverse operation.
 *
 * This is used in square-free factorization, with order==-1
 *
 * This must not be confused with the frobenius map modulo fx *and* fy.
 */
void mod_q::frobenius(
        cxx_mpz_poly_bivariate & B,
        cxx_mpz_poly_bivariate const & A,
        int order) const
{
    ASSERT_ALWAYS(mpz_size(p) == 1);
    const mp_limb_t pp = mpz_get_ui(p);

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
            mpz_poly_swap(B.v()[k], B.v()[k*pp]);
        }
        B.erase(B.begin() + n + 1, B.end());
    }

    for(auto & c : B.v())
        frobenius(c, c, order);

    if (order == 1) {
        /* scale up powers for forward Frobenius */

        B.insert(B.end(), n * (pp - 1), cxx_mpz_poly());
        for(int k = n ; k >= 0 ; k--)
            mpz_poly_swap(B.v()[k * pp], B.v()[k]);
        B.cleandeg(n * pp);
        ASSERT_ALWAYS(B.degree_y() >= 0);
        ASSERT_ALWAYS((unsigned int) B.degree_y() == n * pp);
    }
}

void mod_fy_mod_q::frobenius(
        cxx_mpz_poly_bivariate & B,
        cxx_mpz_poly_bivariate const & A,
        int order) const
{
    ASSERT_ALWAYS(mpz_size(p) == 1);

    /* for the moment the others are unimplemented */
    ASSERT_ALWAYS(order == 1 || order == -1);

    /* there might be situations where the input isn't reduced.
     */
    (*this)(B, A);

    int n = B.degree_y();

    if (n == -1)
        return;

    /* for the moment we don't have code to compute the inverse frobenius
     * matrix */
    ASSERT_ALWAYS(order == 1);

    cado::matrix<cxx_mpz_poly> const & mat = frobenius_matrix();

    cxx_mpz_poly_bivariate h;
    /* It's weird, but we don't have row times matrix in mpz_mat...  */
    for(int j = 0 ; j < fy.degree() ; j++) {
        cxx_mpz_poly s;
        ASSERT_ALWAYS(B.degree() < fy.degree());
        for(int i = 0 ; i <= B.degree() ; i++) {
            s += B[i] * mat(i, j);
        }
        (*this)(s, s);
        h.setcoeff(j, s);
    }
    B = h;
}

} /* namespace cado::arithmetic_reductions */

