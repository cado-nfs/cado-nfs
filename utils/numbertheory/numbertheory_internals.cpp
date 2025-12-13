#include "cado.h" // IWYU pragma: keep
/* 
 * Authors: Joshua Peignier and Emmanuel Thomé
 */
#include <climits>
#include <cstdio>

#include <algorithm>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <gmp.h>

#include "cxx_mpz.hpp"
#include "gmp_aux.h"
#include "macros.h"
#include "mpz_mat.h"
#include "mpz_mat_accessors.h"
#include "mpz_poly.h"
#include "numbertheory_internals.hpp"
#include "numbertheory_fwd_types.hpp"
#include "runtime_numeric_cast.hpp"


/*{{{ commonly used wrappers around HNF functions */

static cxx_mpz_mat join_HNF(cxx_mpz_mat const& K, cxx_mpz const& p)//{{{
{
    cxx_mpz_mat J(K->n, K->n, p);
    cxx_mpz_mat I;

    mpz_mat_vertical_join(I, J, K);
    mpz_mat_hermite_form(I);
    mpz_mat_submat_swap(I, 0, 0, J, 0, 0, K->n, K->n);
    return J;
}
//}}}

#if 0
static cxx_mpz_mat join_HNF_rev(cxx_mpz_mat const& K, cxx_mpz const& p)//{{{
{
    // Builds the block matrix containing p*identity in the top, and K in
    // the bottom Then computes its HNF and stores it in I
    cxx_mpz_mat J(K->n, K->n, p);
    cxx_mpz_mat T0;
    cxx_mpz_mat I;

    mpz_mat_vertical_join(I, J, K);
    mpz_mat_hermite_form_rev(I, T0);
    mpz_mat_submat_swap(I, 0, 0, J, 0, 0, K->n, K->n);
    return J;
}//}}}
#endif
/*}}}*/

/* {{{ multiplication_table_of_order */
/* The correct interface is
number_field_order::number_field_order(class number_field & K, cxx_mpq_mat const & mat)
 */
cxx_mpz_mat numbertheory_internals::multiplication_table_of_order(cxx_mpq_mat const& O,
                                  cxx_mpz_poly const& g)
{
    /* Let O be an order, with basis written with respect to the
     * polynomial basis defined by g (of degree denoted below by n). O is
     * thus an n*n matrix.
     *
     * This function computes the integer matrix of size n*n^2 such that
     * the n coordinates at position (i,j*n) to (i,j*n+n-1) are the
     * coordinates of the i-th times the j-th generator of O, expressed
     * as combinations of the generators of O.
     */
    unsigned int const n = g->deg;
    ASSERT_ALWAYS(O->m == n);
    ASSERT_ALWAYS(O->n == n);
    cxx_mpq_mat R;
    mpq_mat_inv(R, O);
    cxx_mpz_mat M(n, n * n);

    for(unsigned int i = 0; i < n; i++) {
        cxx_mpz_poly w;
        cxx_mpz dw;
        mpq_mat_row_to_poly(w, dw, O, i);
        cxx_mpq_mat T(n, n);
        for (unsigned int j = 0; j < n; j++) {
            cxx_mpz_poly c;
            cxx_mpz dc;
            mpq_mat_row_to_poly(c, dc, O, j);
            mpz_poly_mul_mpz(c, c, mpz_poly_lc(g));
            mpz_poly_mul_mod_f(c, c, w, g);
            mpz_mul(dc, dc, dw);
            mpz_mul(dc, dc, mpz_poly_lc(g));
            mpq_poly_to_mat_row(T, j, c, dc);
        }
        mpq_mat_mul(T, T, R);
        cxx_mpz_mat Tz;
        int const rc = mpq_mat_numden(Tz, nullptr, T);
        if (!rc) {
            throw element_not_integral();
        }
        for (unsigned int j = 0; j < n; j++) {
            mpz_mat_submat_swap(M, i, j*n, Tz, j, 0, 1, n);
        }
    }
    return M;
}/*}}}*/

/* {{{ multiplication_table_of_ideal*/
static cxx_mpz_mat multiplication_table_of_ideal(cxx_mpz_mat const& M,
				  cxx_mpz_mat const& I)
{
    /* Let O be an order of a degree n number field. Let M be the n*n^2
     * multiplication matrix of O.
     *
     * Let I be an ideal of O, given as an n*n matrix with entries
     * expressed as coordinate vectors with respect to the basis O.
     *
     * This function computes the integer matrix of size n*n^2 such that
     * the n coordinates at position (i,j*n) to (i,j*n+n-1) are the
     * coordinates of the i-th generator of O times the j-th generator of
     * I, expressed as combinations of the generators of I.
     */
    unsigned int const n = M->m;
    ASSERT_ALWAYS(M->n == n * n);
    ASSERT_ALWAYS(I->m == n);
    ASSERT_ALWAYS(I->n == n);

    cxx_mpq_mat R;
    mpq_mat_inv(R, cxx_mpq_mat(I));
    cxx_mpz_mat MI(n, n * n);
    for(unsigned int i = 0 ; i < n ; i++) {
        cxx_mpz_mat Tz(n, n);
        for (unsigned int j = 0; j < n; j++) {
            mpz_mat_submat_set(Tz, j, 0, M, i, j*n, 1, n);
        }
        /* We have the matrix of multiplication by w (some generator of
         * O). Row i is w*wi.
         *
         * We need to compute I*T*I^-1 in order to have that represent
         * how multiplication by w affects the generators of I.
         */
        mpz_mat_mul(Tz, I, Tz);
        cxx_mpq_mat T(Tz);
        mpq_mat_mul(T, T, R);
        int const rc = mpq_mat_numden(Tz, nullptr, T);
        ASSERT_ALWAYS(rc);
        for (unsigned int j = 0; j < n; j++) {
            mpz_mat_submat_swap(MI, i, j*n, Tz, j, 0, 1, n);
        }
    }
    return MI;
}/*}}}*/

/*{{{ multiply_elements_in_order */
/* The correct interface is number_field_order_element::operator* */
cxx_mpz_mat numbertheory_internals::multiply_elements_in_order(cxx_mpz_mat const& M, cxx_mpz_mat const& E, cxx_mpz_mat const& F)
{
    /* Let O be an order in a degree n number field.
     *
     * Given the n*n^2 matrix M of the multiplication within that order,
     * given matrices E and F of size k*n with coordinates of k elements
     * of O, return the matrix of size k*n with i-th row giving
     * coordinates of e_i times f_i in O.
     */
    unsigned int const n = M->m;
    unsigned int const k = E->m;
    ASSERT_ALWAYS(M->n == n * n);
    ASSERT_ALWAYS(E->n == n);
    ASSERT_ALWAYS(F->n == n);
    ASSERT_ALWAYS(F->m == k);

    cxx_mpz_mat EM;
    mpz_mat_mul(EM, E, M);
    cxx_mpz_mat R(k, n);
    for(unsigned int ell = 0 ; ell < k ; ell++) {
        for (unsigned int j = 0; j < n; j++) {
            mpz_ptr Rlj = mpz_mat_entry(R, ell, j);
            for(unsigned int i = 0 ; i < n ; i++) {
                mpz_addmul(Rlj,
                        mpz_mat_entry_const(F, ell, i),
                        mpz_mat_entry_const(EM, ell, i*n+j));
            }
        }
    }
    return R;
}/*}}}*/

//{{{ frobenius_matrix
static cxx_mpz_mat frobenius_matrix(cxx_mpz_mat const& M, cxx_mpz const& p)
{
    // frobenius_matrix ; utility function for computing the p-radical
    // Takes a matrix B containing the generators (w_0, ... w_{n-1}) of an
    // order, expressed as polynomials in the root of the polynomial g,
    // return the matrix U containing ((w_0)^p, ..., (w_{n-1})^p), and
    // expressed as a linear transformation within the order B.

    /* This version uses the multiplication table of the order (and does all
     * arithmetic mod p), and binary powering. */
    unsigned int const n = M->m;
    ASSERT_ALWAYS(M->n == n * n);
    cxx_mpz_mat Mp = M;
    mpz_mat_mod_mpz(Mp, Mp, p);

    const cxx_mpz_mat E(n, n, 1);
    cxx_mpz_mat F(E);

    int k = mpz_sizeinbase(p, 2) - 1;
    for( ; k-- ; ) {
        F = numbertheory_internals::multiply_elements_in_order(M, F, F);
        if (mpz_tstbit(p, k))
            F = numbertheory_internals::multiply_elements_in_order(M, E, F);
        mpz_mat_mod_mpz(F, F, p);
    }
    return F;
}//}}}

// {{{ cxx_mpz_mat p_radical_of_order

/* the correct interface is number_field_order::p_radical */

// Stores in I the p-radical of the order whose multiplication matrix is
// given by M. I is expressed with respect to the basis of the order.
cxx_mpz_mat numbertheory_internals::p_radical_of_order(cxx_mpz_mat const& M, cxx_mpz const& p)
{
    unsigned int const n = M->m;
    ASSERT_ALWAYS(M->n == n * n);

    cxx_mpz_mat K;

    // Now building the matrix U, containing all generators to the power
    // of p Generators are polynomials, stored in the matrix B
    cxx_mpz_mat T = frobenius_matrix(M, p);

    /* Which power of p ? */
    int k = 1;
    for (cxx_mpz pk = p; mpz_cmp_ui(pk, n) < 0; k++)
        mpz_mul(pk, pk, p);
    mpz_mat_pow_ui_mod_mpz(T, T, k, p);

    // Storing in K a basis of Ker((z -> (z^%p mod g mod p))^k)
    mpz_mat_kernel_mod_mpz(K, T, p);

    // Getting generators of the p radical from Ker(X) by computing HNF
    // of the vertical block matrix (p*Id, K);
    return join_HNF(K, p);
}// }}}


cxx_mpq_mat numbertheory_internals::p_maximal_order(cxx_mpq_mat const & B, cxx_mpz_poly const& g, cxx_mpz const& p)// {{{
{
    /* Given the basis of an order that is expressed with respect to the
     * polynomial basis of the monic polynomial g, return the basis of a
     * p-maximal order with respect to this same basis.
     */

    ASSERT_ALWAYS(mpz_poly_is_monic(g));

    cxx_mpq_mat D = B;

    for(;;) {
        cxx_mpz_mat const M = multiplication_table_of_order(D, g);

        // Getting the p-radical of the order generated by D in I.
        cxx_mpz_mat const I = p_radical_of_order(M, p);

        // Building the (n,n^2) matrix containing the integers mod p
        // associated to all generators of O
        cxx_mpz_mat M2 = multiplication_table_of_ideal(M, I);
        mpz_mat_mod_mpz(M2, M2, p);
        //printf("M is :\n"); mpz_mat_fprint(stdout, M); printf("\n");
        
        // Computing Ker(M)
        cxx_mpz_mat K;	        // Ker(M)
        mpz_mat_kernel_mod_mpz(K, M2, p);
        // printf("Ker(M) is :\n"); mpz_mat_fprint(stdout, K_M); printf("\n");

        // Getting generators of p*O' by computing HNF of the vertical block matrix (p*Id, K_M);
        cxx_mpz_mat const J = join_HNF(K, p);

        cxx_mpq_mat new_D;

        // Converting in the basis which is used to express elements of D
        mpq_mat_mul(new_D, cxx_mpq_mat(J), D);

        // Dividing by p
        mpq_mat_div_mpz(new_D, new_D, p);

        if (D == new_D)
            return D;

        std::swap(D, new_D);
    }
}// }}}

//{{{ mpz_mat_minpoly_mod_ui
static cxx_mpz_poly mpz_mat_minpoly_mod_mpz(cxx_mpz_mat M, cxx_mpz const& p)
{
    ASSERT_ALWAYS(M->m == M->n);
    unsigned long const n = M->n;
    cxx_mpz_mat B(n+1, n*n);
    cxx_mpz_mat T(n,n,1);
    for(unsigned int k = 0 ; k < n + 1 ; k++) {
        for(unsigned int i = 0 ; i < n ; i++) {
            mpz_mat_submat_set(B, n - k, i * n, T, i, 0, 1, n);
        }
        mpz_mat_mul_mod_mpz(T, T, M, p);
    }
    cxx_mpz_mat K;
    mpz_mat_kernel_mod_mpz(K, B, p);
    mpz_mat_gauss_backend_mod_mpz(K, nullptr, p);

    /* TODO: write down exactly what we need in terms of ordering from
     * mpz_mat_gauss_backend_mod_ui. I take it that we're happy if it
     * computes a RREF, so that the last row is the one with the largest
     * number of column coordinates killed, corresponding to highest
     * degree coefficients absent ?
     */
    cxx_mpz_poly f;
    mpz_mat_row_to_poly_rev(f,K,K->m-1);

    return f;
}
/*}}}*/

/*{{{ matrix_of_multmap */
static cxx_mpz_mat matrix_of_multmap(
        cxx_mpz_mat const& M,
        cxx_mpz_mat const& J, 
        cxx_mpz_mat const& c,
        cxx_mpz const& p)
{
    /* Let O be an order, represented here simply by its multiplication
     * matrix. Let J be an m-dimensional subalgebra of O/pO. Let c be an
     * element of J (a linear combination of the rows of J).
     * This returns the m times m matrix which expresses multiplication by c
     * within J.
     */
    /* We write the coordinates in O/pO of all products of the for c*j_k,
     * j_k being one of the generators of J, and denote this matrix by
     * CJ.  Now we look for a matrix M such that M*J = CJ. To do this, we
     * need to do gaussian reduction on transpose(J): Find a matrix K
     * such that K*transpose(J) is verticaljoin(identity_matrix(m),
     * zero_matrix(n-m,m)). Then, taking K' as the first m rows of K, we
     * will have M = CJ * transpose(K').
     */
    unsigned int const m = J->m;
    unsigned int const n = M->m;
    ASSERT_ALWAYS(M->n == n*n);
    ASSERT_ALWAYS(J->n == n);
    ASSERT_ALWAYS(c->n == n);
    ASSERT_ALWAYS(c->m == 1);
    cxx_mpz_mat Mc;
    cxx_mpz_mat Jt;
    cxx_mpz_mat const K;
    mpz_mat_transpose(Jt, J);
    cxx_mpz_mat T;
    mpz_mat_gauss_backend_mod_mpz(Jt, T, p);
    cxx_mpz_mat Kt(m, n);
    mpz_mat_submat_swap(Kt, 0, 0, T, 0, 0, m, n);
    mpz_mat_transpose(Kt, Kt);  /* Kt is now n times m */
    cxx_mpz_mat CJ(m, n);
    for(unsigned int k = 0 ; k < m ; k++) {
        cxx_mpz_mat jk(1,n);
        mpz_mat_submat_set(jk,0,0,J,k,0,1,n);
        cxx_mpz_mat w = numbertheory_internals::multiply_elements_in_order(M, c, jk);
        mpz_mat_submat_swap(CJ,k,0,w,0,0,1,n);
    }
    mpz_mat_mul_mod_mpz(Mc, CJ, Kt, p);
    return Mc;
}
/*}}}*/

/*{{{ template <typename T> void append_move(vector<T> &a, vector<T> &b) */
template <typename T> static void append_move(std::vector<T> &a, std::vector<T> &b)
{
    a.reserve(a.size() + b.size());
    size_t const na = a.size();
    size_t const nb = b.size();
    if (&a == &b) {
        for(size_t i = 0 ; i < nb ; i++) {
            a.push_back(b[i]);
        }
    } else {
        a.insert(a.end(), nb, T());
        for(size_t i = 0 ; i < nb ; i++) {
            std::swap(a[na + i], b[i]);
        }
    }
}
/*}}}*/

struct ideal_comparator {
    using Im_t = std::pair<cxx_mpz_mat,int>;
    bool operator()(Im_t const& a, Im_t const& b) const {
        int const r = mpz_mat_cmp(a.first, b.first);
        if (r) return r < 0;
        return a.second < b.second;     /* should never happen */
    }
};

// {{{ factorization_of_prime
static std::vector<std::pair<cxx_mpz_mat, int> > factorization_of_prime_inner(
        cxx_mpq_mat const & B,
        cxx_mpz_mat const & M,
        cxx_mpz const& p,
        cxx_mpz_mat const& Ip,
        cxx_mpz_mat const& I,
        cxx_mpz_mat const& J,
        gmp_randstate_t state)
{
    unsigned int const m = J->m;
    unsigned int const n = J->n;
    ASSERT_ALWAYS(I->m == n);
    ASSERT_ALWAYS(I->n == n);
    /* J represents an m-dimensional subalgebra of O/pO */

    /* Pick a random element of J */
    cxx_mpz_mat c(1, m);
    for(unsigned int i = 0 ; i < m ; i++) {
        mpz_urandomm(mpz_mat_entry(c,0,i), state, p);
    }
    mpz_mat_mul_mod_mpz(c, c, J, p);

    cxx_mpz_mat Mc = matrix_of_multmap(M, J, c, p);

    /* Now we would like to find the minimal polynomial of Mc, which is
     * simply the matrix of multiplication by c. For this, we write an
     * (m+1) times m^2 matrix, and compute a left kernel.
     */
    cxx_mpz_poly Pc = mpz_mat_minpoly_mod_mpz(Mc, p);

    std::vector<std::pair<cxx_mpz_poly, int> > facP = mpz_poly_factor(Pc, p, state);

    std::vector<std::pair<cxx_mpz_mat, int> > ideals;

    std::vector<cxx_mpz_mat> characteristic_subspaces;

    for(auto const & [ f, e ] : facP) {
        /* We need the basis of the kernel of f(Mc) */
        cxx_mpz_mat E;
        mpz_poly_eval_mpz_mat_mod_mpz(E, f, Mc, p);
        mpz_mat_pow_ui_mod_mpz(E, E, e, p);
        mpz_mat_kernel_mod_mpz(E, E, p);
        /* This line is just to be exactly in line with what magma says
         */
        mpz_mat_gauss_backend_mod_mpz(E, nullptr, p);
        mpz_mat_mul_mod_mpz(E, E, J, p);
        characteristic_subspaces.push_back(E);
    }

    for(unsigned int i = 0 ; i < facP.size() ; i++) {
        auto const & [ f, e ] = facP[i];
        cxx_mpz_mat const& Ci(characteristic_subspaces[i]);
        /* We need to find an ideal which is smaller than Ip, and whose
         * generators are p*the generators of Ci, as well as the
         * unmodified generators of I and of the other characteristic
         * subspaces.
         */
        cxx_mpz_mat Ix(n + J->m - Ci->m, n);
        mpz_mat_submat_set(Ix,0,0,I,0,0,n,n);
        for(unsigned int r = n, k = 0; k < facP.size() ; k++) {
            if (k == i) continue;
            cxx_mpz_mat const& Ck(characteristic_subspaces[k]);
            unsigned int const mk = Ck->m;
            mpz_mat_submat_set(Ix, r, 0, Ck, 0, 0, mk, n);
            r += mk;
        }
        cxx_mpz_mat Ihead(n, n);
        if (Ci->m == runtime_numeric_cast<unsigned int>(e * f->deg)) {
            mpz_mat_vertical_join(Ix, Ix, Ip);
            mpz_mat_hermite_form_rev(Ix);
            mpz_mat_submat_swap(Ihead,0,0,Ix,0,0,n,n);
            ideals.emplace_back(Ihead, e);
        } else {
            mpz_mat_hermite_form_rev(Ix);
            mpz_mat_submat_swap(Ihead,0,0,Ix,0,0,n,n);
            std::vector<std::pair<cxx_mpz_mat, int> > more_ideals;
            more_ideals = factorization_of_prime_inner(B,M,p,Ip,Ihead,Ci,state);
            append_move(ideals, more_ideals);
        }
    }
    std::ranges::sort(ideals, ideal_comparator());
    return ideals;
}

/* Let B be an order, with basis written with respect to the
 * polynomial basis defined by g (of degree denoted below by n). B is
 * thus an n*n matrix.
 *
 * This returns the factorization of the ideal generated by p in the
 * order B.
 *
 * As it turns out, this function is only ever called with B as returned
 * by the p_maximal_order function. No guarantees as to how it would
 * behave is B is not a b-maximal order. Maybe your chair will transform
 * into a coconut.
 *
 * the correct interface is number_field_order::factor
 */
std::vector<std::pair<cxx_mpz_mat, int> > numbertheory_internals::factorization_of_prime(
        cxx_mpq_mat const & B, cxx_mpz_poly const& g,
        cxx_mpz const& p,
        gmp_randstate_ptr state)
{
    int const n = g->deg;
    cxx_mpz_mat const M = multiplication_table_of_order(B, g);
    cxx_mpz_mat const Ip = p_radical_of_order(M, p);
    return factorization_of_prime_inner(B, M, p, Ip,
            cxx_mpz_mat(n, n, p),
            cxx_mpz_mat(n, n, 1), state);
}
//}}}


// {{{ valuation_helper_for_ideal
//
// compute an uniformizing element a for the prime ideal I of the order O,
// where I is above the rational prime p, and O is p-maximal and given by
// means of its multiplication matrix M.
//
// The uniformizing element is such that (a/p)*I is in O, yet a is not in
// p*O.
//
// a is returned as a 1 times n matrix, and the coefficients are a
// polynomial representation of a with respect to the basis of the order
// that M is a multiplication matrix of.
cxx_mpz_mat numbertheory_internals::valuation_helper_for_ideal(cxx_mpz_mat const& M, cxx_mpz_mat const& I, cxx_mpz const& p)
{
    unsigned int const n = M->m;
    ASSERT_ALWAYS(M->n == n * n);
    ASSERT_ALWAYS(I->m == n);
    ASSERT_ALWAYS(I->n == n);

    /* We begin by something which is very similar to
     * multiplication_table_of_ideal. It's a bit like computing I*M,
     * except that we want it in a different order.
     */
    cxx_mpz_mat MI(n, n * n);
    for(unsigned int i = 0 ; i < n ; i++) {
        cxx_mpz_mat Tz(n, n);
        for (unsigned int j = 0; j < n; j++) {
            mpz_mat_submat_set(Tz, j, 0, M, i, j*n, 1, n);
        }
        mpz_mat_mul_mod_mpz(Tz, I, Tz, p);
        for (unsigned int j = 0; j < n; j++) {
            mpz_mat_submat_swap(MI, i, j*n, Tz, j, 0, 1, n);
        }
    }

    cxx_mpz_mat ker;
    mpz_mat_kernel_mod_mpz(ker, MI, p);
    mpz_mat_hermite_form(ker);

    cxx_mpz_mat res(1, n);
    mpz_mat_submat_swap(res, 0, 0, ker, 0, 0, 1, n);
    return res;
}// }}}

// {{{ generate_ideal -- create an ideal from a set of generators
// generators (gens) are here given as elements of the field, not
// elements of the order. In case the ideal is fractional, its
// denominator is also returned. For an integral ideal, the denominator
// is always 1.
std::pair<cxx_mpz_mat, cxx_mpz> numbertheory_internals::generate_ideal(cxx_mpq_mat const& O, cxx_mpz_mat const& M, cxx_mpq_mat const& gens)
{
    unsigned int const n = M->m;
    ASSERT_ALWAYS(M->n == n * n);
    ASSERT_ALWAYS(O->m == n);
    ASSERT_ALWAYS(O->n == n);
    ASSERT_ALWAYS(gens->n == n);
    cxx_mpq_mat R;
    mpq_mat_inv(R, O);
    cxx_mpq_mat I0q;
    mpq_mat_mul(I0q, gens, R);
    cxx_mpz_mat I0;
    cxx_mpz denom;
    mpq_mat_numden(I0, denom, I0q);
    /* Now create the full list of elements of O which are generated by
     * the products I_i * O_j */
    cxx_mpz_mat products(gens->m * O->m, n);
    for(unsigned int i = 0 ;  i < I0->m ; i++) {
        cxx_mpz_mat a(1,n);
        mpz_mat_submat_set(a,0,0,I0,i,0,1,n);
        mpz_mat_mul(a,a,M);
        for(unsigned int j = 0 ; j < n ; j++){
            mpz_mat_submat_swap(a,0,j*n,products,i*n+j,0,1,n);
        }
    }
    /* And put this in HNF */
    mpz_mat_hermite_form_rev(products);
    cxx_mpz_mat I(n,n);
    mpz_mat_submat_swap(I,0,0,products,0,0,n,n);
    return { I, denom };
}//}}}

int numbertheory_internals::prime_ideal_inertia_degree(cxx_mpz_mat const& I)/*{{{*/
{
    unsigned int const n = I->m;
    ASSERT_ALWAYS(I->n == n);
    int f = 0;
    for(unsigned int i = 0 ; i < n ; i++) {
        f += mpz_cmp_ui(mpz_mat_entry_const(I, i, i), 1) != 0;
    }
    return f;
}
/*}}}*/
int numbertheory_internals::valuation_of_ideal_at_prime_ideal(cxx_mpz_mat const& M, cxx_mpz_mat const& I, cxx_mpz_mat const& a, cxx_mpz const& p)/*{{{*/
{
    /* M is the multiplication table of the order. I is an ideal. We want
     * to compute the fkp-valuation of I, where fkp is an ideal above p
     * with a helper element given by a
     */
    unsigned int const n = M->m;
    ASSERT_ALWAYS(M->n == n * n);
    ASSERT_ALWAYS(I->m == n);
    ASSERT_ALWAYS(I->n == n);
    ASSERT_ALWAYS(a->m == 1);
    ASSERT_ALWAYS(a->n == n);
    
    if (mpz_mat_is_zero(I))
        return INT_MAX;

    cxx_mpz_mat Ia(I);
    int val = 0;
    for(;;val++) {
        for(unsigned int i = 0 ; i < n ; i++) {
            cxx_mpz_mat b(1,n);
            mpz_mat_submat_set(b,0,0,Ia,i,0,1,n);
            b = numbertheory_internals::multiply_elements_in_order(M, a, b);
            mpz_mat_submat_swap(b,0,0,Ia,i,0,1,n);
        }
        if (mpz_mat_p_valuation(Ia, p) < 1)
            return val;
        mpz_mat_divexact_mpz(Ia, Ia, p);
    }
}
/*}}}*/
int numbertheory_internals::valuation_of_ideal_at_prime_ideal(cxx_mpz_mat const& M, std::pair<cxx_mpz_mat,cxx_mpz> const& Id, cxx_mpz_mat const& a, int e, cxx_mpz const& p)/*{{{*/
{
    /* M is the multiplication table of the order. Id is a pair (ideal,
     * denominator of ideal). We want to compute the fkp-valuation of I,
     * where fkp is an ideal above p with a helper element given by a. e
     * is the fkp-valuation of p (that is, the ramification index of
     * fkp).
     *
     * INT_MAX is returned if Id.first is the zero ideal, but please
     * don't try to use it.
     */
    int const v = numbertheory_internals::valuation_of_ideal_at_prime_ideal(M, Id.first, a, p);
    if (v == INT_MAX) return v;

    int const w = mpz_p_valuation(Id.second, p);
    return v-w*e;
}
/*}}}*/

struct hypercube_walk {/*{{{*/
    struct iterator {
        int B;
        std::vector<int> v;
        std::vector<int> speed;
        std::pair<int, int> last;
        iterator& operator++() {
            for(unsigned int j = 0 ; j < v.size() ; j++) {
                int const s = v[j] + speed[j]; 
                if (s >= 0 && s <= B) {
                    v[j] = s; 
                    last = std::make_pair(j, speed[j]);
                    return *this;
                }
                speed[j]=-speed[j];
            }
            last = std::make_pair(-1, -1);
            return *this;
        }
        bool operator!=(iterator const& x) const { return last != x.last; }
        std::vector<int> & operator*() { return v; }
    };
    int n,B;
    hypercube_walk(int n, int B) : n(n), B(B) {}
    iterator begin() const {
        iterator z;
        z.B = B; 
        z.v.assign(n, 0);
        z.speed.assign(n, 1);
        return z;
    }
    iterator middle(std::vector<int> const& x) const {
        ASSERT_ALWAYS(x.size() == (unsigned int) n);
        iterator z;
        z.B = B; 
        z.v = x;
        z.speed.assign(n, 1);
        return z;
    }
    iterator end() const {
        iterator z;
        z.B = B; 
        z.v.assign(n,0);
        z.speed.assign(n, 1);
        z.last = std::make_pair(-1,-1);
        return z;
    }
};/*}}}*/
 
/* {{{ prime_ideal_two_element */
/* I must be in HNF */
std::pair<cxx_mpz, cxx_mpz_mat> numbertheory_internals::prime_ideal_two_element(cxx_mpq_mat const& O, cxx_mpz_poly const& f, cxx_mpz_mat const& M, cxx_mpz_mat const& I)
{
    cxx_mpz p;
    unsigned int const n = M->m;
    ASSERT_ALWAYS(M->n == n * n);
    ASSERT_ALWAYS(I->m == n);
    ASSERT_ALWAYS(I->n == n);
    int inertia = 0;
    for(unsigned int i = 0 ; i < n ; i++) {
        mpz_srcptr di = mpz_mat_entry_const(I, i, i);
        if (mpz_cmp_ui(di, 1) != 0) {
            mpz_set(p, di);
            inertia++;
        }
    }
    ASSERT_ALWAYS(mpz_size(p) != 0);

    /* We use a homebrew algorithm, loosely inspired on Cohen's. We
     * prefer an enumeration which favors small generators. Based on
     * [Cohen93, lemma 4.7.9], we're simply going to search for an
     * element whose norm has p-valuation exactly f.
     *
     * Note though that we're only achieving small element size in terms
     * of coordinates on the order.
     */
    /* There's a question about which generators we should combine in
     * order to find our winning generator. Cohen takes "the generators
     * of p*O, and the generators of I". In fact, we can simply take the
     * generators of I as a Z-basis: that generates the right set. Then,
     * we'd better work on some ordering for these. Elements in the HNF
     * form of I with a p on the diagonal have in fact all other elements
     * on that row equal to zero. That is because since I contains p*O,
     * the subtraction of this generating row and the p-th multiple of
     * the corresponding generator has to be zero, or I would not be in
     * HNF. Therefore, these elements are of secondary importance in
     * generating I (e.g.: alone, they can't).
     */
    unsigned int const m = n; /* number of generators we take */
    cxx_mpz_mat OI(m,n);
    for(unsigned int i = 0, r = 0, s = 0 ; i < n ; i++) {
        mpz_srcptr di = mpz_mat_entry_const(I, i, i);
        if (mpz_cmp_ui(di, 1) != 0) {
            mpz_mat_submat_set(OI,n - inertia + s++,0,I,i,0,1,n);
        } else {
            mpz_mat_submat_set(OI,r++,0,I,i,0,1,n);
        }
    }

    cxx_mpq_mat OIOq;
    mpq_mat_mul(OIOq, cxx_mpq_mat(OI), O); 
    for(int B = 1 ; ; B++) {
        hypercube_walk const H(m, B);
        std::vector<int> const H0(m,0);
        cxx_mpq_mat gen(1, n);
        // H0[n]=1;
        // mpq_mat_submat_set(gen, 0, 0, OIOq, n, 0, 1, n);
        cxx_mpq_mat temp(1, n);
        for(hypercube_walk::iterator it = H.middle(H0) ; it != H.end() ; ++it) {
            int const j = it.last.first;
            int const s = it.last.second;
            /* add (s=1) or subtract (s=-1) the j-th generator to gen */
            mpq_mat_submat_swap(temp,0,0,OIOq,j,0,1,n);
            if (s==1)
                mpq_mat_add(gen, gen, temp);
            else if (s==-1)
                mpq_mat_sub(gen, gen, temp);
            mpq_mat_submat_swap(temp,0,0,OIOq,j,0,1,n);

            cxx_mpz_poly pgen;
            cxx_mpz dgen;
            cxx_mpz res;
            mpq_mat_row_to_poly(pgen, dgen, gen, 0);
            if (pgen->deg < 0) continue;
            mpz_poly_resultant(res, pgen, f);
            /* We want the absolute norm to have p-valuation equal to the
             * inertia degree.
             * absolute norm is galois norm. galois norm is product of
             * all conjugate of pgen(alpha). Resultant(f,pgen) is
             * lc(f)^deg(pgen) times the galois norm.
             */
            int const v = mpz_p_valuation(res, p) - pgen->deg * mpz_p_valuation(mpz_poly_lc(f), p) - f->deg * mpz_p_valuation(dgen, p);
            ASSERT_ALWAYS(v >= inertia);
            if (v == inertia) {
                cxx_mpz_mat lambda(1, m);
                std::vector<int> & v(*it);
                for(unsigned int i = 0 ; i < m ; i++) {
                    mpz_set_si(mpz_mat_entry(lambda,0,i),v[i]);
                }
                mpz_mat_mul(lambda, lambda, OI);
                return std::make_pair(p, lambda);
            }
        }
    }
}
// }}}

std::string numbertheory_internals::write_element_as_polynomial(cxx_mpq_mat const& theta_q, std::string const& var)
{
    ASSERT_ALWAYS(theta_q->m == 1);
    cxx_mpz theta_denom;
    cxx_mpz_poly theta;
    mpq_mat_row_to_poly(theta, theta_denom, theta_q, 0);

    /* first write the numerator as a string */
    std::string num = theta.print_poly(var);
    if (mpz_cmp_ui(theta_denom, 1) == 0) {
        return num;
    } else {
        std::ostringstream os2;
        os2 << "(" << num << ")/" << theta_denom;
        return os2.str();
    }
}

std::string numbertheory_internals::write_order_element_as_vector(cxx_mpz_mat const& z)
{
    ASSERT_ALWAYS(z->m == 1);
    std::ostringstream s;
    s << "[";
    for(unsigned int i = 0 ; i < z->n ; i++) {
        if (i) s << ", ";
        s << z(0,i);
    }
    s << "]";
    return s.str();
}

std::vector<cxx_mpz>
numbertheory_internals::write_element_as_list_of_integers(cxx_mpq_mat const& theta_q)
{
    ASSERT_ALWAYS(theta_q->m == 1);
    cxx_mpz theta_denom;
    cxx_mpz_mat theta;
    mpq_mat_numden(theta, theta_denom, theta_q);
    std::vector<cxx_mpz> res;
    res.push_back(theta_denom);
    res.reserve(theta->n + 1);
    for(unsigned int i = 0 ; i < theta->n ; i++)
        res.emplace_back(mpz_mat_entry(theta, 0, i));
    return res;
}
