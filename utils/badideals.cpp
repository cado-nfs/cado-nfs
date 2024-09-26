/*
 * Authors: Joshua Peignier and Emmanuel Thomé
 */
#include "cado.h" // IWYU pragma: keep
// IWYU pragma: no_include <memory>
// iwyu wants it for allocator_traits<>::value_type, which seems weird
#include <ostream>
#include <iomanip>      // std::dec // IWYU pragma: keep
#include <sstream>
#include <utility>

#include <cstddef> // size_t
#include <gmp.h>
#include <fmt/format.h>

#include "gmp_aux.h"  // mpz_p_valuation
#include "mpz_mat.h"
#include "numbertheory.hpp"
#include "badideals.hpp"
#include "rootfinder.h"
#include "getprime.h"  // for getprime_mt, prime_info_clear, prime_info_init
#include "macros.h" // ASSERT_ALWAYS // IWYU pragma: keep
#include "misc.h"

using namespace std;


    std::ostream& badideal::print_dot_badideals_file(std::ostream & o, int side) const {/*{{{*/
        o << p
          << "," << r
          << ":" << side
          << ": " << nbad << std::endl;
        return o;
    }/*}}}*/
    std::ostream& badideal::print_dot_badidealinfo_file(std::ostream& o, int side) const {/*{{{*/
        o << comments;
        for(unsigned int j = 0 ; j < branches.size() ; j++) {
            badideal::branch const& br(branches[j]);
            o << p << " " << br.k << " " << br.r << " " << side;
            for(unsigned int k = 0 ; k < br.v.size() ; k++) {
                o << " " << br.v[k];
            }
            o << std::endl;
        }
        return o;
    }/*}}}*/

std::ostream& badideal::operator<<(std::ostream& os) const
{
    std::ios_base::fmtflags ff = os.flags();
    os << comments;
    // p and r are printed according to the current format flags, but the
    // number of ideals and the branches are always in decimal.
    os << p
        << " " << r
        << dec
        << "  " << nbad
        << "  " << branches.size()
        << std::endl;
    for(unsigned int j = 0 ; j < branches.size() ; j++) {
        badideal::branch const& br(branches[j]);
        os << p << " " << br.k << " " << br.r;
        ASSERT_ALWAYS(br.v.size() == (size_t) nbad);
        for(unsigned int k = 0 ; k < br.v.size() ; k++) {
            os << " " << br.v[k];
        }
        os << std::endl;
    }
    os.flags(ff);
    return os;
}

std::istream& badideal::operator>>(std::istream& is)
{
    badideal c(is);
    std::swap(*this, c);
    return is;
}

badideal::badideal(std::istream& is)
{
    for(std::string s; std::ws(is).peek() == '#' ; getline(is, s) ) ;
    size_t nbranches;
    // coverity[tainted_argument]
    is >> p >> r >> nbad >> nbranches;
    if (!is) return;
    for(unsigned int j = 0 ; j < nbranches ; j++) {
        badideal::branch br;
        is >> p >> br.k >> br.r;
        ASSERT_ALWAYS_OR_THROW(is, std::invalid_argument);
        br.v.assign(nbad, 0);
        for(unsigned int k = 0 ; k < br.v.size() ; k++) {
            is >> br.v[k];
            ASSERT_ALWAYS_OR_THROW(is, std::invalid_argument);
        }
        branches.emplace_back(std::move(br));
    }
    return;
}

struct all_valuations_above_p {/*{{{*/
    cxx_mpz_poly f;
    cxx_mpz p;
private:
    cxx_mpq_mat O;
    cxx_mpz_mat M;
    vector<pair<cxx_mpz_mat, int> > F;
    vector<int> inertia;
    vector<int> ramification; // valuation of the ideal in the
                              // factorization of the underlying prime
                              // ideal.
    pair<cxx_mpz_mat, cxx_mpz> jjinv;
    vector<cxx_mpz_mat> helpers;
    vector<int> val_base;

public:
    all_valuations_above_p(cxx_mpz_poly const& f, cxx_mpz const& p, gmp_randstate_t state) : f(f), p(p) {/*{{{*/
        O = p_maximal_order(f, p);
        M = multiplication_table_of_order(O, f);
        F = factorization_of_prime(O, f, p, state);
        for(unsigned int k = 0 ; k < F.size() ; k++) {
            cxx_mpz_mat const& fkp(F[k].first);
            inertia.push_back(prime_ideal_inertia_degree(fkp));
            ramification.push_back(F[k].second);
            helpers.push_back(valuation_helper_for_ideal(M, fkp, p));
        }
        cxx_mpq_mat jjinv_gen(2, f->deg);
        mpq_set_ui(mpq_mat_entry(jjinv_gen,0,0),1,1);
        mpq_set_ui(mpq_mat_entry(jjinv_gen,1,1),1,1);
        jjinv = ::generate_ideal(O,M,jjinv_gen);

        val_base.assign(f->deg, 0);
        val_base = (*this)(jjinv);
    }/*}}}*/
    vector<int> operator()(pair<cxx_mpz_mat, cxx_mpz> const& Id) const {/*{{{*/
        int w = mpz_p_valuation(Id.second, p);
        vector<int> res;
        for(unsigned int k = 0 ; k < F.size() ; k++) {
            cxx_mpz_mat const& a(helpers[k]);
            int v = valuation_of_ideal_at_prime_ideal(M, Id.first, a, p);
            int e = F[k].second;
            res.push_back(v - w * e - val_base[k]);
        }
        return res;
    }/*}}}*/
    std::string sagemath_string(int k, int side) {
        cxx_mpz_mat const& fkp(F[k].first);
        pair<cxx_mpz, cxx_mpz_mat> two = prime_ideal_two_element(O, f, M, fkp);
        /* Write the uniformizer as a polynomial with respect to the
         * polynomial basis defined by f */
        cxx_mpq_mat theta_q;
        {
            mpq_mat_set_mpz_mat(theta_q, two.second);
            mpq_mat_mul(theta_q, theta_q, O);
        }

        /* That's only for debugging, so it's not terribly important.
         * But we may have a preference towards giving ideal info in
         * a more concise way, based on alpha_hat for instance */
        ostringstream alpha;
        alpha << "alpha" << side;
        string uniformizer = write_element_as_polynomial(theta_q, alpha.str());

        return fmt::format("OK{}.ideal({}, {})", side, two.first, uniformizer);
    }
    void print_info(ostream& o, int k, cxx_mpz const& r MAYBE_UNUSED, int side) const {/*{{{*/
        cxx_mpz_mat const& fkp(F[k].first);
        pair<cxx_mpz, cxx_mpz_mat> two = prime_ideal_two_element(O, f, M, fkp);
        /* Write the uniformizer as a polynomial with respect to the
         * polynomial basis defined by f */
        cxx_mpq_mat theta_q;
        {
            mpq_mat_set_mpz_mat(theta_q, two.second);
            mpq_mat_mul(theta_q, theta_q, O);
        }

        /* That's only for debugging, so it's not terribly important.
         * But we may have a preference towards giving ideal info in
         * a more concise way, based on alpha_hat for instance */
        ostringstream alpha;
        alpha << "alpha" << side;
        string uniformizer = write_element_as_polynomial(theta_q, alpha.str());

        /* This prints magma code. */
        int e = F[k].second;
        o << "# I" << k
            << ":=ideal<O" << side << "|" << two.first << "," << uniformizer << ">;"
            << " // f=" << prime_ideal_inertia_degree(fkp)
            << " e="<< e
            << endl;
        o << "# I_" << two.first << "_" << r << "_" << side << "_" << k
            << " " << two.first
            << " " << r
            << " " << side
            << " " << theta_q
            << " // " << prime_ideal_inertia_degree(fkp)
            << " " << e
            << endl;
    }/*}}}*/
    pair<cxx_mpz_mat, cxx_mpz> generate_ideal(cxx_mpq_mat const& gens) const {/*{{{*/
        return ::generate_ideal(O, M, gens);
    }/*}}}*/
    pair<cxx_mpz_mat, cxx_mpz> generate_ideal(cxx_mpz_mat const& gens) const {/*{{{*/
        return ::generate_ideal(O, M, cxx_mpq_mat(gens));
    }/*}}}*/
    /* create ideal I=<p^k,p^k*alpha,v*alpha-u> and decompose I*J */
    vector<int> operator()(int k, cxx_mpz const& r) const {/*{{{*/
        cxx_mpz pk;
        mpz_pow_ui(pk, p, k);
        cxx_mpz_mat Igens(3, f->deg);
        mpz_set(mpz_mat_entry(Igens,0,0),pk);
        if (mpz_cmp(r, pk) < 0) {
            mpz_neg(mpz_mat_entry(Igens,1,0),r);
            mpz_set_ui(mpz_mat_entry(Igens,1,1),1);
        } else {
            mpz_set_si(mpz_mat_entry(Igens,1,0),-1);
            mpz_sub(mpz_mat_entry(Igens,1,1), r, pk);
        }
        /* hell, do I _really_ need p*alpha here ??? */
        mpz_set(mpz_mat_entry(Igens,2,1),pk);
        pair<cxx_mpz_mat, cxx_mpz> I = generate_ideal(Igens);
        return (*this)(I);
    }/*}}}*/
    vector<int> multiply_inertia(vector<int> const& v) const {/*{{{*/
        ASSERT_ALWAYS(v.size() == inertia.size());
        vector<int> res(v.size(),0);
        for(unsigned int i = 0 ; i < v.size() ; i++) {
            res[i] = v[i] * inertia[i];
        }
        return res;
    }/*}}}*/
};/*}}}*/

vector<cxx_mpz> lift_p1_elements(cxx_mpz const& p, int k, cxx_mpz const& x)/*{{{*/
{
    /* Given x which represents an element of P^1(Z/p^kZ), return all the p
     * lifts of x in P^1(Z/p^(k+1)Z), all following the same representation
     * convention (detailed somewhat above)
     */
    /* Note how we have complexity O(p) here ! */
    ASSERT_ALWAYS(k >= 1);
    vector<cxx_mpz> res;
    cxx_mpz pk, pkp1;
    mpz_pow_ui(pk, p, k);
    mpz_mul(pkp1, pk, p);
    cxx_mpz xx = x;
    if (mpz_cmp(xx, pk) >= 0) {
        mpz_sub(xx, xx, pk);
        mpz_add(xx, xx, pkp1);
    }
    for(unsigned int w = 0 ; mpz_cmp_ui(p, w) > 0 ; w++) {
        res.push_back(xx);
        mpz_add(xx, xx, pk);
    }
    return res;
}/*}}}*/

vector<badideal::branch> lift_root(all_valuations_above_p const& A, int k0, cxx_mpz const& Q, vector<int> v)/*{{{*/
{
    vector<badideal::branch> dead_branches_reports;
    vector<pair<cxx_mpz, vector<int> > > live_branches;
    vector<int> live_ideals;
    cxx_mpz const& p(A.p);

    //int k = k0 + 1;

    vector<cxx_mpz> rootlifts = lift_p1_elements(p, k0, Q);

    for(unsigned int i = 0 ; i < rootlifts.size() ; i++) {
        cxx_mpz const& nQ(rootlifts[i]);
        vector<int> newvals = A(k0 + 1, nQ);
        vector<int> alive;
        /* All valuations are expected to be >= the base valuation. For
         * some, it's going to lift higher. Find which ones.  */
        for(unsigned int j = 0 ; j < v.size() ; j++) {
            if (newvals[j] != v[j]) alive.push_back(j);
        }
        if (alive.empty()) {
            badideal::branch br;
            br.k = k0 + 1;
            br.r = nQ;
            br.v = A.multiply_inertia(v);
            dead_branches_reports.push_back(br);
        } else {
            live_branches.push_back(make_pair(nQ, newvals));
        }
        live_ideals.insert(live_ideals.end(), alive.begin(), alive.end());
    }
    if (live_ideals.size() <= 1) {
        /* Then the decision is basically taken */
        vector<int> vv = A.multiply_inertia(v);
        int sumvv = 0;
        for(unsigned int j = 0 ; j < vv.size() ; sumvv+=vv[j++]);
        for(unsigned int j = 0 ; j < live_ideals.size() ; j++) {
            int jj = live_ideals[j];
            vv[jj] = vv[jj] - sumvv;
        }
        badideal::branch br;
        br.k = k0;
        br.r = Q;
        br.v = vv;
        return vector<badideal::branch>(1, br);
    }

    vector<badideal::branch> res = dead_branches_reports;

    for(unsigned int i = 0 ; i < live_branches.size() ; i++) {
        cxx_mpz const& nQ(live_branches[i].first);
        vector<int> const& nv(live_branches[i].second);

        vector<badideal::branch> add = lift_root(A, k0 + 1, nQ, nv);

        res.insert(res.end(), add.begin(), add.end());
    }
    return res;
}/*}}}*/

vector<cxx_mpz> projective_roots_modp(cxx_mpz_poly const& f, cxx_mpz const& p, gmp_randstate_ptr rstate)/*{{{*/
{
    /* p must be prime */
    vector<cxx_mpz> roots;
    mpz_t * rr = new mpz_t[f->deg];
    for(int i = 0 ; i < f->deg ; i++) mpz_init(rr[i]);

    int d = mpz_poly_roots(rr, f, p, rstate);
    for(int i = 0 ; i < d ; i++) {
        cxx_mpz a;
        mpz_set(a, rr[i]);
        roots.push_back(a);
    }
    if (mpz_divisible_p(f->coeff[f->deg], p)) {
        roots.push_back(p);
    }
    for(int i = 0 ; i < f->deg ; i++) mpz_clear(rr[i]);
    delete[] rr;
    return roots;
}/*}}}*/

vector<badideal> badideals_above_p(cxx_mpz_poly const& f, int side, cxx_mpz const& p, gmp_randstate_t state)/*{{{*/
{
    vector<badideal> badideals;

    all_valuations_above_p A(f, p, state);

    vector<cxx_mpz> roots = projective_roots_modp(f, p, state);

    for(unsigned int i = 0 ; i < roots.size() ; i++) {
        /* first try to decompose <p,(v*alpha-u)>*J */
        vector<int> vals = A(1, roots[i]);

        vector<int> nonzero;
        for(unsigned int k = 0 ; k < vals.size() ; k++) {
            if (vals[k]) nonzero.push_back(k);
        }
        if (nonzero.size() == 1)
            continue;

        badideal b(p,roots[i], nonzero.size());

        vector<badideal::branch> lifts = lift_root(A, 1, roots[i], vals);

        ostringstream cmt;
        cmt << "# p=" << p << ", r=" << roots[i] << " : " << nonzero.size() << " ideals among " << vals.size() << " are bad\n";
        for(unsigned int j = 0 ; j < nonzero.size() ; j++) {
            A.print_info(cmt, nonzero[j], roots[i], side);
            b.sagemath_string.push_back(A.sagemath_string(nonzero[j], side));
        }
        cmt << "# " << lifts.size() << " branch"
            << (lifts.size() == 1 ? "" : "es") << " found\n";
        b.comments = cmt.str();

        ASSERT_ALWAYS(b.sagemath_string.size() == (size_t) b.nbad);

        /* compres all branches so that we keep only the valuations in
         * the nonzero indirection table */
        for(unsigned int j = 0 ; j < lifts.size() ; j++) {
            badideal::branch & br(lifts[j]);
            vector<int> w;
            for(unsigned int k = 0 ; k < nonzero.size() ; k++) {
                w.push_back(br.v[nonzero[k]]);

            }
            swap(br.v, w);
            b.branches.push_back(br);
        }
        badideals.push_back(b);
    }
    return badideals;
}/*}}}*/

vector<badideal> badideals_for_polynomial(cxx_mpz_poly const& f, int side, gmp_randstate_t state)/*{{{*/
{
    vector<badideal> badideals;

    if (f->deg == 1) return badideals;

    ASSERT_ALWAYS(f->deg > 1);

    cxx_mpz disc;
    mpz_poly_discriminant(disc, f);
    mpz_mul(disc, disc, f->coeff[f->deg]);

    /* We're not urged to use ecm here */
    vector<pair<cxx_mpz,int> > small_primes = trial_division(disc, 10000000, disc);


    typedef vector<pair<cxx_mpz,int> >::const_iterator vzci_t;

    for(vzci_t it = small_primes.begin() ; it != small_primes.end() ; it++) {
        vector<badideal> tmp = badideals_above_p(f, side, it->first, state);
        badideals.insert(badideals.end(), tmp.begin(), tmp.end());
    }

    return badideals;
}/*}}}*/

vector<badideal> badideals_for_polynomial(cxx_mpz_poly const& f, int side)
{
    gmp_randstate_t state;
    gmp_randinit_default(state);
    auto x = badideals_for_polynomial(f, side, state);
    gmp_randclear(state);
    return x;
}

vector<badideal> badideals_above_p(cxx_mpz_poly const& f, int side, cxx_mpz const & p)
{
    gmp_randstate_t state;
    gmp_randinit_default(state);
    auto x = badideals_above_p(f, side, p, state);
    gmp_randclear(state);
    return x;
}

cxx_mpz badideal::r_from_rk(cxx_mpz const & p, int k, cxx_mpz const & rk)
{
    cxx_mpz pk, r;
    mpz_pow_ui(pk, p, k);
    if (mpz_cmp(rk, pk) < 0) {
        mpz_mod(r, rk, p);
        return r;
    } else {
        /* then u=rk-p^k is divisible by p, and represents (1:u). Reduced
         * mod p, this means (1:0), which is encoded by p */
        return p;
    }
}

std::string generic_sagemath_string(cxx_mpz_poly const & f, int side, cxx_mpz const & p, cxx_mpz const & r, gmp_randstate_t state)
{
    /* This will crash for non-prime ideals, **on purpose** */
    auto A = all_valuations_above_p(f, p, state);
    auto v = A(1, r);
    int k = -1;
    for(unsigned x = 0 ; x < v.size() ; x++) {
        if (v[x] == 0)
            continue;
        if (k != -1)
            throw std::runtime_error("ideal is not prime");
        k = x;
    }
    if (k == -1)
        throw std::runtime_error("valuations of ideal not found");
    return A.sagemath_string(k, side);
}

