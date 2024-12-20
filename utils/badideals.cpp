/*
 * Authors: Joshua Peignier and Emmanuel Thomé
 */
#include "cado.h" // IWYU pragma: keep
// IWYU pragma: no_include <memory>
// iwyu wants it for allocator_traits<>::value_type, which seems weird
#include <ostream>
#include <iomanip>      // std::dec // IWYU pragma: keep
#include <iostream>
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
#include "numbertheory/all_valuations_above_p.hpp"

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
    std::ios_base::fmtflags const ff = os.flags();
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
            int const jj = live_ideals[j];
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
            b.machine_description.push_back(A.machine_description(nonzero[j]));
        }
        cmt << "# " << lifts.size() << " branch"
            << (lifts.size() == 1 ? "" : "es") << " found\n";
        b.comments = cmt.str();

        ASSERT_ALWAYS(b.sagemath_string.size() == (size_t) b.nbad);
        ASSERT_ALWAYS(b.machine_description.size() == (size_t) b.nbad);

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
