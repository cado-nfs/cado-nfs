/*
 * Authors: Joshua Peignier and Emmanuel Thomé
 */
#include "cado.h" // IWYU pragma: keep

#include <cstddef>

#include <iostream>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <gmp.h>
#include "fmt/format.h"

#include "badideals.hpp"
#include "cxx_mpz.hpp"
#include "gmp_aux.h"
#include "macros.h"
#include "misc.h"
#include "mpz_poly.h"
#include "numbertheory/all_valuations_above_p.hpp"
#include "rootfinder.h"
#include "utils_cxx.hpp"

using namespace std;

std::ostream& badideal::print_dot_badideals_file(std::ostream & o, int side) const {/*{{{*/
    o << p
        << "," << r
        << ":" << side
        << ": " << nbad << "\n";
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
        o << "\n";
    }
    return o;
}/*}}}*/

std::ostream& badideal::operator<<(std::ostream& os) const
{
    std::ios_base::fmtflags const ff = os.flags();
    os << comments;
    // p and r are printed according to the current format flags, but the
    // number of ideals and the branches are always in decimal.
    os << p
        << " " << r
        << dec
        << "  " << nbad
        << "  " << branches.size()
        << "\n";
    for(auto const & br: branches) {
        os << p << " " << br.k << " " << br.r;
        ASSERT_ALWAYS(br.v.size() == (size_t) nbad);
        for(auto v : br.v) os << " " << v;
        os << "\n";
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
    branches.reserve(nbranches);
    for(unsigned int j = 0 ; j < nbranches ; j++) {
        badideal::branch br;
        is >> p >> br.k >> br.r;
        ASSERT_ALWAYS_OR_THROW(is, std::invalid_argument);
        br.v.assign(nbad, 0);
        for(auto & v : br.v) {
            is >> v;
            ASSERT_ALWAYS_OR_THROW(is, std::invalid_argument);
        }
        branches.emplace_back(std::move(br));
    }
}

static vector<cxx_mpz> lift_p1_elements(cxx_mpz const& p, int k, cxx_mpz const& x)/*{{{*/
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

static vector<badideal::branch> lift_root(numbertheory_internals::all_valuations_above_p const& A, int k0, cxx_mpz const& Q, vector<int> v)/*{{{*/
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
            live_branches.emplace_back(nQ, newvals);
        }
        live_ideals.insert(live_ideals.end(), alive.begin(), alive.end());
    }
    if (live_ideals.size() <= 1) {
        /* Then the decision is basically taken */
        vector<int> vv = A.multiply_inertia(v);
        int sumvv = 0;
        for(unsigned int j = 0 ; j < vv.size() ; sumvv+=vv[j++]);
        for(int const jj : live_ideals) {
            vv[jj] -= sumvv;
        }
        badideal::branch br;
        br.k = k0;
        br.r = Q;
        br.v = vv;
        return vector<badideal::branch>(1, br);
    }

    vector<badideal::branch> res = dead_branches_reports;

    for(auto const & br : live_branches) {
        cxx_mpz const& nQ(br.first);
        vector<int> const& nv(br.second);

        vector<badideal::branch> add = lift_root(A, k0 + 1, nQ, nv);

        res.insert(res.end(), add.begin(), add.end());
    }
    return res;
}/*}}}*/

static vector<cxx_mpz> projective_roots_modp(cxx_mpz_poly const& f, cxx_mpz const& p, gmp_randstate_ptr rstate)/*{{{*/
{
    /* p must be prime */
    vector<cxx_mpz> roots;
    mpz_t * rr = new mpz_t[f->deg];
    for(int i = 0 ; i < f->deg ; i++) mpz_init(rr[i]);

    int const d = mpz_poly_roots(rr, f, p, rstate);
    for(int i = 0 ; i < d ; i++) {
        cxx_mpz a;
        mpz_set(a, rr[i]);
        roots.push_back(a);
    }
    if (mpz_divisible_p(mpz_poly_lc(f), p)) {
        roots.push_back(p);
    }
    for(int i = 0 ; i < f->deg ; i++) mpz_clear(rr[i]);
    delete[] rr;
    return roots;
}/*}}}*/

vector<badideal> badideals_above_p(cxx_mpz_poly const& f, int side, cxx_mpz const& p, cxx_gmp_randstate & state)/*{{{*/
{
    vector<badideal> badideals;

    numbertheory_internals::all_valuations_above_p A(f, p, state);

    A.bless_side(side);

    vector<cxx_mpz> roots = projective_roots_modp(f, p, state);

    for(auto const & r : roots) {
        /* first try to decompose <p,(v*alpha-u)>*J */
        vector<int> vals = A(1, r);

        vector<int> nonzero;
        for(unsigned int k = 0 ; k < vals.size() ; k++) {
            if (vals[k]) nonzero.push_back(k);
        }
        if (nonzero.size() == 1)
            continue;

        badideal b(p,r, nonzero.size());

        vector<badideal::branch> lifts = lift_root(A, 1, r, vals);

        ostringstream cmt;
        cmt << "# p=" << p << ", r=" << r << " : " << nonzero.size() << " ideals among " << vals.size() << " are bad\n";
        for(auto nz : nonzero) {
            A.print_info(cmt, nz, r, side);
            b.sagemath_string.push_back(A.sagemath_string(nz, side));
            b.machine_description.push_back(A.machine_description(nz));
        }
        cmt << "# " << lifts.size() << " branch"
            << (lifts.size() == 1 ? "" : "es") << " found\n";
        b.comments = cmt.str();

        ASSERT_ALWAYS(b.sagemath_string.size() == (size_t) b.nbad);
        ASSERT_ALWAYS(b.machine_description.size() == (size_t) b.nbad);

        /* compres all branches so that we keep only the valuations in
         * the nonzero indirection table */
        for(auto & br : lifts) {
            vector<int> w;
            w.reserve(nonzero.size());
            for(auto nz : nonzero)
                w.push_back(br.v[nz]);
            swap(br.v, w);
            b.branches.push_back(br);
        }
        badideals.push_back(b);
    }
    return badideals;
}/*}}}*/

vector<badideal> badideals_for_polynomial(cxx_mpz_poly const& f, int side, cxx_gmp_randstate & state)/*{{{*/
{
    vector<badideal> badideals;

    if (f->deg == 1) return badideals;

    ASSERT_ALWAYS(f->deg > 1);

    cxx_mpz disc;
    mpz_poly_discriminant(disc, f);
    mpz_mul(disc, disc, mpz_poly_lc(f));

    /* We're not urged to use ecm here */
    vector<pair<cxx_mpz,int> > const small_primes = trial_division(disc, 10000000, disc);

    for(auto const & pe : small_primes) {
        vector<badideal> tmp = badideals_above_p(f, side, pe.first, state);
        badideals.insert(badideals.end(), tmp.begin(), tmp.end());
    }

    return badideals;
}/*}}}*/

vector<badideal> badideals_for_polynomial(cxx_mpz_poly const& f, int side)
{
    cxx_gmp_randstate state;
    auto x = badideals_for_polynomial(f, side, state);
    return x;
}

vector<badideal> badideals_above_p(cxx_mpz_poly const& f, int side, cxx_mpz const & p)
{
    cxx_gmp_randstate state;
    auto x = badideals_above_p(f, side, p, state);
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

std::string generic_sagemath_string(cxx_mpz_poly const & f, int side, cxx_mpz const & p, cxx_mpz const & r)
{
    cxx_gmp_randstate state;
    /* This will crash for non-prime ideals, **on purpose** */
    auto A = numbertheory_internals::all_valuations_above_p(f, p, state);
    A.bless_side(side);
    auto v = A(1, r);
    int k = -1;
    for(unsigned x = 0 ; x < v.size() ; x++) {
        if (v[x] == 0)
            continue;
        int inertia = A.get_inertia_degree(x);
        if (inertia != 1) {
            std::cerr << fmt::format(
                        "# note: seemingly innocuous prime ideal ({},{}) on side {} has non-trivial residue class degree {}\n",
                    p, r, side, inertia);
        }
        if (k != -1)
            throw std::runtime_error("ideal is not prime");
        k = x;
    }
    if (k == -1)
        throw std::runtime_error("valuations of ideal not found");
    return A.sagemath_string(k, side);
}

std::string generic_machine_description(cxx_mpz_poly const & f, int, cxx_mpz const & p, cxx_mpz const & r)
{
    cxx_gmp_randstate state;
    /* This will crash for non-prime ideals, **on purpose** */
    auto A = numbertheory_internals::all_valuations_above_p(f, p, state);
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
    return A.machine_description(k);
}

int get_inertia_of_prime_ideal(cxx_mpz_poly const & f, cxx_mpz const & p, cxx_mpz const & r)
{
    cxx_gmp_randstate state;
    auto A = numbertheory_internals::all_valuations_above_p(f, p, state);
    auto v = A(1, r);
    for(unsigned x = 0 ; x < v.size() ; x++) {
        if (v[x] == 0)
            continue;
        return A.get_inertia_degree(x);
    }
    return 1;
}
