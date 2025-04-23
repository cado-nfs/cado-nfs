#include "cado.h" // IWYU pragma: keep

#include <cstdint>
#include <climits>
#include <cstdio>
#include <cstdlib>

#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <gmp.h>
#include "fmt/base.h"
#include "fmt/format.h"

#include "cxx_mpz.hpp"
#include "getprime.h"
#include "gmp_aux.h"
#include "macros.h"
#include "misc.h"
#include "mpz_mat.h"
#include "mpz_poly.h"
#include "numbertheory.hpp"
#include "params.h"
#include "portability.h"
#include "timing.h"

using namespace std;

/* we'd like to get rid of this! */
using namespace numbertheory_internals;

static const char ** original_argv;

static void decl_usage(cxx_param_list & pl)/*{{{*/
{
    param_list_decl_usage(pl, "test", "which test to run");
    param_list_decl_usage(pl, "prime", "prime");
    param_list_decl_usage(pl, "poly", "polynomial (string)");
    param_list_decl_usage(pl, "out", "output file");
    param_list_decl_usage(pl, "batch", "batch input file with test vectors and expected results");
    param_list_decl_usage(pl, "seed", "seed used for random picks");
    param_list_decl_usage(pl, "elements", "ideal generators (separated by ;)");
}/*}}}*/

static void usage(cxx_param_list & pl, char const ** argv, const char * msg = nullptr)/*{{{*/
{
    param_list_print_usage(pl, argv[0], stderr);
    if (msg) {
        fprintf(stderr, "%s\n", msg);
    }
    exit(EXIT_FAILURE);
}/*}}}*/

static int do_p_maximal_order(cxx_param_list & pl) /*{{{*/
{
    cxx_mpz p;
    if (!param_list_parse_mpz(pl, "prime", p)) usage(pl, original_argv, "missing prime argument");

    std::string polystr;
    if (!param_list_parse(pl, "poly", polystr)) usage(pl, original_argv, "missing poly argument");
    cxx_mpz_poly const f(polystr);

    cxx_mpq_mat M = p_maximal_order(f, p);
    cxx_mpz D;
    cxx_mpz_mat A;
    mpq_mat_numden(A, D, M);
    std::cout << "1/" << D << "*\n" << A << "\n";

    return 1;
}
/*}}}*/

static bool sl_equivalent_matrices(cxx_mpq_mat const& M, cxx_mpq_mat const& A, cxx_mpz const& p)/*{{{*/
{
    /* This is over SL_n(Z_p) */
    if (M->m != A->m) return false;
    if (M->n != A->n) return false;
    cxx_mpq_mat Mi;
    mpq_mat_inv(Mi, M);
    cxx_mpq_mat AMi;
    mpq_mat_mul(AMi, A, Mi);
    /* check that the p-valuation is zero */
    for(unsigned int i = 0 ; i < AMi->m ; i++) {
        for(unsigned int j = 0 ; j < AMi->n ; j++) {
            mpq_srcptr mij = mpq_mat_entry_const(AMi, i, j);
            if (mpz_divisible_p(mpq_denref(mij), p)) return false;
        }
    }
    return true;
}/*}}}*/

#if 0
bool sl_equivalent_matrices(cxx_mpq_mat const& M, cxx_mpq_mat const& A)/*{{{*/
{
    /* unimplemented for the moment, since unneeded.  We would compute
     * A*M^-1, and see whether we have a denominator. */
    if (M->m != A->m) return false;
    if (M->n != A->n) return false;
    cxx_mpq_mat Mi;
    mpq_mat_inv(Mi, M);
    cxx_mpq_mat AMi;
    mpq_mat_mul(AMi, A, Mi);
    for(unsigned int i = 0 ; i < AMi->m ; i++) {
        for(unsigned int j = 0 ; j < AMi->n ; j++) {
            mpq_srcptr mij = mpq_mat_entry_const(AMi, i, j);
            if (mpz_cmp_ui(mpq_denref(mij), 1) != 0) return false;
        }
    }
    return true;
}/*}}}*/
#endif

static cxx_mpq_mat batch_read_order_basis(istream & in, unsigned int n)/*{{{*/
{
    cxx_mpq_mat O;
    string keyword;
    if (!(in >> keyword) || keyword != "order")
        throw invalid_argument("Parse error");
    cxx_mpz_mat A(n, n);
    cxx_mpz d;
    if (!(in >> d))
        throw invalid_argument("Parse error");
    for(unsigned int i = 0 ; i < A->m ; i++) {
        for(unsigned int j = 0 ; j < A->n ; j++) {
            if (!(in >> mpz_mat_entry(A, i, j)))
                throw invalid_argument("Parse error");
        }
    }
    mpq_mat_set_mpz_mat_denom(O, A, d);
    return O;
}/*}}}*/

static vector<pair<cxx_mpz_mat, int> > batch_read_prime_factorization(istream & in, unsigned int n, cxx_mpz const& p, cxx_mpq_mat const& O, cxx_mpz_mat const& M)/*{{{*/
{
    vector<pair<cxx_mpz_mat, int> > ideals;
    string keyword;
    if (!(in >> keyword) || keyword != "ideals")
        throw invalid_argument("Parse error");
    unsigned int nideals;
    // coverity[tainted_argument]
    if (!(in >> nideals))
        throw invalid_argument("Parse error");
    for(unsigned int k = 0 ; k < nideals ; k++) {
        cxx_mpq_mat A(2, n);
        cxx_mpz den;
        int e;
        mpq_set_z(mpq_mat_entry(A,0,0),p);
        if (!(in >> den))
            throw invalid_argument("Parse error");
        for(unsigned int j = 0 ; j < A->n ; j++) {
            cxx_mpz num;
            if (!(in >> num))
                throw invalid_argument("Parse error");
            mpq_ptr aa = mpq_mat_entry(A, 1, j);
            mpz_set(mpq_numref(aa), num);
            mpz_set(mpq_denref(aa), den);
            mpq_canonicalize(aa);
        }
        if (!(in >> e))
            throw invalid_argument("Parse error");
        pair<cxx_mpz_mat, cxx_mpz> Id = generate_ideal(O,M,A);
        ASSERT_ALWAYS(mpz_cmp_ui(Id.second, 1) == 0);
        ideals.emplace_back(Id.first,e);
    }
    return ideals;
}/*}}}*/


static int do_p_maximal_order_batch(cxx_param_list & pl) /*{{{*/
{
    const char * tmp;

    if ((tmp = param_list_lookup_string(pl, "batch")) == nullptr)
        usage(pl, original_argv, "missing batch argument");

    ifstream is(tmp);
    string s;
    int nok = 0;
    int nfail = 0;

    for(int test = 0; getline(is, s, '\n') ; ) {
        if (s.empty()) continue;
        if (s[0] == '#') continue;
        const string exc = string("Parse error on input ") + s;

        cxx_mpz_poly f;
        if (!mpz_poly_set_from_expression(f, s))
            usage(pl, original_argv, "cannot parse polynomial");


        if (!(getline(is, s, '\n')))
            throw invalid_argument(exc);
        cxx_mpz const d;
        cxx_mpz_mat const A(f->deg, f->deg);
        istringstream is1(s);
        cxx_mpz p;
        if (!(is1 >> p))
            throw invalid_argument(exc);

        cxx_mpq_mat const O = batch_read_order_basis(is1, f->deg);
        cxx_mpq_mat const my_O = p_maximal_order(f, p);

        bool const ok = sl_equivalent_matrices(O, my_O, p);

        cout << ok_NOK(ok) << " test " << test
            << " (degree " << f->deg << ", p=" << p << ")"
            << "\n";
        test++;
        nok += ok;
        nfail += !ok;
    }
    cout << nok << " tests passed";
    if (nfail) cout << ", " << nfail << " TESTS FAILED";
    cout << "\n";
    return nfail == 0;
}
/*}}}*/

static int do_factorization_of_prime(cxx_param_list & pl) /*{{{*/
{
    cxx_mpz p;
    if (!param_list_parse_mpz(pl, "prime", p)) usage(pl, original_argv, "missing prime argument");

    std::string polystr;
    if (!param_list_parse(pl, "poly", polystr)) usage(pl, original_argv, "missing poly argument");
    cxx_mpz_poly const f(polystr);

    number_field K(f);
    K.bless("K", "alpha");
    number_field_order O = K.p_maximal_order(p);
    O.bless(fmt::format("O{}", p));

    cxx_gmp_randstate state;
    unsigned long seed = 0;
    if (param_list_parse_ulong(pl, "seed", &seed)) {
        gmp_randseed_ui(state, seed);
    }
    fmt::print("{}.<{}lpha>={:S}\n", K.name, K.varname, K);
    fmt::print("{}={:S}\n", O.name, O);
    auto F = O.factor(p, state);

    // number_field_fractional_ideal I = O.fractional_ideal({K(11)});

    fmt::print("assert {}.fractional_ideal({}) == prod([\n", O.name, p);
    for(auto const & fe : F) {
        fmt::print("\t{}", fe.first);
        if (fe.second > 1)
            fmt::print("^{}", fe.second);
        fmt::print(",\n");
    }
    fmt::print("])\n");
    return 1;
}
/*}}}*/

static int do_factorization_of_prime_batch(cxx_param_list & pl) /*{{{*/
{
    const char * tmp;

    if ((tmp = param_list_lookup_string(pl, "batch")) == nullptr)
        usage(pl, original_argv, "missing batch argument");

    ifstream is(tmp);
    string s;
    int nok = 0;
    int nfail = 0;
    for(int test = 0; getline(is, s, '\n') ; ) {
        if (s.empty()) continue;
        if (s[0] == '#') continue;

        const string exc = string("Parse error on input ") + s;

        cxx_mpz_poly f;
        if (!mpz_poly_set_from_expression(f, s))
            usage(pl, original_argv, "cannot parse polynomial");

        if (!(getline(is, s, '\n')))
            throw invalid_argument(exc);

        istringstream is1(s);

        cxx_mpz p;
        if (!(is1 >> p)) throw invalid_argument(exc);

        cxx_mpq_mat const O = batch_read_order_basis(is1, f->deg);
        cxx_mpz_mat const M = multiplication_table_of_order(O, f);

        vector<pair<cxx_mpz_mat, int> > ideals = batch_read_prime_factorization(is1, f->deg, p, O, M);

        cxx_mpq_mat const my_O = p_maximal_order(f, p);

        bool ok = sl_equivalent_matrices(O, my_O, p);

        gmp_randstate_t state;
        gmp_randinit_default(state);
        unsigned long seed = 0;
        if (param_list_parse_ulong(pl, "seed", &seed)) {
            gmp_randseed_ui(state, seed);
        }

        /* What if we simply ask our code to compute the factorization of
         * p with respect to the order basis which was chosen by magma ?
         * */
        vector<pair<cxx_mpz_mat, int> > my_ideals = factorization_of_prime(O, f, p, state);
        gmp_randclear(state);

        // sort magma ideals. Ours are sorted already.
        sort(ideals.begin(), ideals.end(), ideal_comparator());
        ok = ok && (ideals.size() == my_ideals.size());
        for(unsigned int k = 0 ; ok && k < ideals.size() ; k++) {
            ok = (ideals[k] == my_ideals[k]);
        }
        cout << ok_NOK(ok) << " test " << test
            << " (degree " << f->deg << ", p=" << p << ")"
            << "\n";
        test++;
        nok += ok;
        nfail += !ok;
    }
    cout << nok << " tests passed";
    if (nfail) cout << ", " << nfail << " TESTS FAILED";
    cout << "\n";
    return nfail == 0;
}
/*}}}*/

static int do_valuations_of_ideal(cxx_param_list & pl) /*{{{*/
{
    cxx_mpz p;
    if (!param_list_parse_mpz(pl, "prime", p)) usage(pl, original_argv, "missing prime argument");

    std::string polystr;
    if (!param_list_parse(pl, "poly", polystr)) usage(pl, original_argv, "missing poly argument");
    cxx_mpz_poly f(polystr);


    /* Now read the element description */
    vector<cxx_mpz_poly> elements; 
    {
        std::string tmp;
        if (!param_list_parse(pl, "elements", tmp))
            usage(pl, original_argv, "missing ideal generators");
        for( ; !tmp.empty() ; ) {
            auto nq = tmp.find(';');
            auto desc = tmp;
            if (nq != std::string::npos) {
                desc = tmp.substr(0, nq);
                tmp = tmp.substr(nq + 1);
            } else {
                tmp.clear();
            }
            istringstream is(desc);
            cxx_mpz_poly x;
            if (!(is >> x)) usage(pl, original_argv, "cannot parse ideal generators");
            elements.push_back(x);
        }
    }


    cxx_mpq_mat O = p_maximal_order(f, p);
    cxx_mpz_mat const M = multiplication_table_of_order(O, f);

    /* We need to reduce our generating elements modulo f, at the expense
     * of creating denominators all over the place.
     *
     * Note that we are *not* changing conventions at all.
     *
     * We start with something in the basis [alpha^i]. We rewrite it in
     * the basis [alpha_hat^i]. We reduce modulo f_hat, which is monic.
     * And then we rewrite that back to the basis [alpha^i].
     */
    pair<cxx_mpz_mat, cxx_mpz> Id;
    {
        cxx_mpq_mat generators(elements.size(), f->deg);
        cxx_mpz_poly fh;
        mpz_poly_to_monic(fh, f);
        for(unsigned int i = 0 ; i < elements.size() ; i++) {
            cxx_mpz_poly& e(elements[i]);
            /* e is e0+e1*x+e2*x^2+...
             * f(alpha) is zero.
             * let y = lc(f)*alpha
             * g is monic, and g(y) = 0.
             * e(alpha) is also e'(lc*alpha), with
             * lc(f)^deg(e)*e'(y) = e0*lc^deg_e + e1*lc^(deg_e-1)*x + ...
             * e' can be reduced mod fh.
             * now we need to compute ([y^i]e')*lc^i / lc^deg_e.
             */
            cxx_mpz denom, c;
            mpz_set_ui(denom, 1);
            mpz_set_ui(c, 1);
            for(int k = e->deg ; k-- ; ) {
                mpz_mul(denom, denom, mpz_poly_lc(f));
                mpz_mul(mpz_poly_coeff(e, k), mpz_poly_coeff_const(e, k), denom);
            }
            mpz_poly_div_r(e, e, fh);
            for(int k = 0 ; k <= e->deg ; k++) {
                mpz_mul(mpz_poly_coeff(e, k), mpz_poly_coeff_const(e, k), c);
                mpz_mul(c, c, mpz_poly_lc(f));
                mpq_ptr gik = mpq_mat_entry(generators, i, k);
                mpz_set(mpq_numref(gik), mpz_poly_coeff_const(e, k));
                mpz_set(mpq_denref(gik), denom);
                mpq_canonicalize(gik);
            }
        }
        Id = generate_ideal(O, M, generators);
    }

    gmp_randstate_t state;
    gmp_randinit_default(state);
    unsigned long seed = 0;
    if (param_list_parse_ulong(pl, "seed", &seed)) {
        gmp_randseed_ui(state, seed);
    }
    vector<pair<cxx_mpz_mat, int> > F = factorization_of_prime(O, f, p, state);
    gmp_randclear(state);

    for(unsigned int k = 0 ; k < F.size() ; k++) {
        cxx_mpz_mat const& fkp(F[k].first);
        cxx_mpz_mat const a = valuation_helper_for_ideal(M, fkp, p);
        pair<cxx_mpz, cxx_mpz_mat> two = prime_ideal_two_element(O, f, M, fkp);

        cxx_mpq_mat theta_q;
        {
            mpq_mat_set_mpz_mat(theta_q, two.second);
            mpq_mat_mul(theta_q, theta_q, O);
        }
        string const uniformizer = write_element_as_polynomial(theta_q, "alpha");

        int const e = F[k].second;
        int const v = valuation_of_ideal_at_prime_ideal(M, Id, a, e, p);

        cout << "# (p=" << p
            << ", k=" << k
            << ", f="<< prime_ideal_inertia_degree(fkp)
            << ", e="<< e
            << "; "
            << "ideal<O|" << two.first << "," << uniformizer << ">"
            << ")"
            << "^" << v
            << ";" << "\n";
    }
    return 1;
}
/*}}}*/

static int do_valuations_of_ideal_batch(cxx_param_list & pl) /*{{{*/
{
    const char * tmp;

    if ((tmp = param_list_lookup_string(pl, "batch")) == nullptr)
        usage(pl, original_argv, "missing batch argument");

    ifstream is(tmp);
    string s;
    int nok = 0;
    int nfail = 0;

    for(int test = 0; getline(is, s, '\n') ; ) {
        if (s.empty()) continue;
        if (s[0] == '#') continue;

        const string exc = string("Parse error on input ") + s;

        cxx_mpz_poly f;
        if (!mpz_poly_set_from_expression(f, s))
            usage(pl, original_argv, "cannot parse polynomial");

        if (!(getline(is, s, '\n')))
            throw invalid_argument(exc);

        istringstream is1(s);

        cxx_mpz p;
        if (!(is1 >> p)) throw invalid_argument(exc);

        cxx_mpq_mat const O = batch_read_order_basis(is1, f->deg);
        cxx_mpz_mat const M = multiplication_table_of_order(O, f);

        vector<pair<cxx_mpz_mat, int> > ideals = batch_read_prime_factorization(is1, f->deg, p, O, M);

        cxx_mpq_mat const my_O = p_maximal_order(f, p);

        bool ok = sl_equivalent_matrices(O, my_O, p);

        gmp_randstate_t state;
        gmp_randinit_default(state);
        unsigned long seed = 0;
        if (param_list_parse_ulong(pl, "seed", &seed)) {
            gmp_randseed_ui(state, seed);
        }

        /* What if we simply ask our code to compute the factorization of
         * p with respect to the order basis which was chosen by magma ?
         * */
        vector<pair<cxx_mpz_mat, int> > my_ideals = factorization_of_prime(O, f, p, state);
        gmp_randclear(state);

        /* compute matching table */
        ok = ok && (ideals.size() == my_ideals.size());
        vector<unsigned int> magma_to_mine(ideals.size());
        vector<unsigned int> mine_to_magma(my_ideals.size());
        for(unsigned int k = 0 ; ok && k < ideals.size() ; k++) {
            bool found = false;
            for(unsigned int ell = 0 ; ell < my_ideals.size() ; ell++) {
                if (ideals[k] == my_ideals[ell]) {
                    magma_to_mine[k] = ell;
                    mine_to_magma[ell] = k;
                    found = true;
                    break;
                }
            }
            ok = found;
        }

        /* now read the list of composites */
        for( ; ok ; ) {
            string keyword;
            if (!(is1 >> keyword) || keyword != "composite") throw invalid_argument(exc);
            int ngens;
            // coverity[-taint_source]
            if (!(is1 >> ngens)) throw invalid_argument(exc);
            if (ngens < 0) throw invalid_argument(exc);
            if (ngens == 0) break;
            cxx_mpz_mat gens(ngens, f->deg);
            for(unsigned int i = 0 ; i < gens->m ; i++) {
                for(unsigned int j = 0 ; j < gens->n ; j++) {
                    if (!(is1 >> mpz_mat_entry(gens, i, j)))
                        throw invalid_argument(exc);
                }
            }
            pair<cxx_mpz_mat, cxx_mpz> const Id = generate_ideal(O, M, cxx_mpq_mat(gens));
            vector<int> my_vals;
            for(auto const & Ie : my_ideals) {
                cxx_mpz_mat const& fkp(Ie.first);
                int const e = Ie.second;
                cxx_mpz_mat const a = valuation_helper_for_ideal(M, fkp, p);
                int const v = valuation_of_ideal_at_prime_ideal(M, Id, a, e, p);
                my_vals.push_back(v);
            }
            if (!(is1 >> keyword) || keyword != "valuations") throw invalid_argument(exc);
            for(unsigned int k = 0 ; ok && k < ideals.size() ; k++) {
                std::string s;
                is1 >> s;
                int v;
                if (s == "Infinity") {
                    v = INT_MAX;
                } else {
                    std::istringstream is2(s);
                    if (!(is2 >> v)) throw invalid_argument(exc);
                }
                ok = (v == my_vals[magma_to_mine[k]]);
            }
        }


        cout << ok_NOK(ok) << " test " << test
            << " (degree " << f->deg << ", p=" << p << ")"
            << "\n";
        test++;
        nok += ok;
        nfail += !ok;
    }
    cout << nok << " tests passed";
    if (nfail) cout << ", " << nfail << " TESTS FAILED";
    cout << "\n";
    return nfail == 0;
}
/*}}}*/

static int do_maximal_order(cxx_param_list & pl)
{
    cxx_gmp_randstate state;
    unsigned long seed = 0;
    if (param_list_parse_ulong(pl, "seed", &seed)) {
        gmp_randseed_ui(state, seed);
    }

    unsigned long bound;
    if (!param_list_parse(pl, "bound", bound)) {
        usage(pl, original_argv, "missing bound argument for maximal_order");
    }

    std::string polystr;
    if (!param_list_parse(pl, "poly", polystr)) usage(pl, original_argv, "missing poly argument");
    const cxx_mpz_poly f(polystr);

    number_field K(f);
    K.bless("K", "alpha");
    fmt::print("print(\"working with {}\")\n", K);

    const number_field_order OK = K.maximal_order(bound);
    // OK.bless("OK");

    for(auto const & e : OK.basis())
        fmt::print("{}\n", e);

    return 0;
}

static int do_number_theory_object_interface(cxx_param_list & pl)
{
    cxx_gmp_randstate state;
    unsigned long seed = 0;
    if (param_list_parse_ulong(pl, "seed", &seed)) {
        gmp_randseed_ui(state, seed);
    }

    std::string polystr;
    if (!param_list_parse(pl, "poly", polystr)) usage(pl, original_argv, "missing poly argument");
    cxx_mpz_poly f(polystr);


    number_field K(f);

    K.bless("K", "alpha");

    fmt::print("print(\"working with {}\")\n", K);
    fmt::print("ZP.<x> = ZZ[]\n");
    fmt::print("{}.<{}> = NumberField({})\n", K.name, K.varname, f);

    auto alpha = K.gen();
    auto x = K(1);

    const int N = 2 * f.degree();

    for(int i = 0 ; i < N ; i++) {
        fmt::print("assert alpha^{} == {}\n", i, x);
        x = x * alpha;
    }

    fmt::print("assert (alpha^{}).trace() == {}\n", N, x.trace());

    prime_info pi;
    prime_info_init(pi);
    for(unsigned long pp ; (pp = getprime_mt(pi)) < 12 ; pp++) {
        fmt::print("print(\"tests mod p={}\")\n", pp);
        cxx_mpz p = pp;
        number_field_order Op = K.p_maximal_order(p);
        Op.bless(fmt::format("O{}", p));
        fmt::print("print(\"computed {}\")\n", Op);
        fmt::print("{} = {:S}\n", Op.name, Op);
        fmt::print("assert {}.is_maximal({})\n", Op.name, p);

        auto fac_p = Op.factor(p);
        fmt::print("assert {}.fractional_ideal({}) == prod([\n", Op.name, p);
        for(auto const & Ie : fac_p) {
            number_field_prime_ideal const & I(Ie.first);
            int e = Ie.second;
            if (e == 1) {
                fmt::print(" {},\n", I);
            } else {
                fmt::print(" ({})^{},\n", I, e);
            }
        }
        fmt::print("])\n");

        for(int i = 0 ; i < 4 ; i++) {
            cxx_mpz_poly phi;
            int64_t a = gmp_urandomm_ui(state, 1000);
            uint64_t b = gmp_urandomm_ui(state, 1000);
            fmt::print("print(\"computing valuations above {} for a={} b={}\")\n", p, a, b);
            mpz_poly_set_ab(phi, a, b);
            number_field_element const z = K(phi);
            number_field_fractional_ideal I = Op.fractional_ideal({z});

            fmt::print("I_{0}_{1} = {2}.fractional_ideal([{0}-{1}*{3}])\n",
                    a, b, Op.name, K.varname);
            fmt::print("assert I_{}_{} == {}\n", a, b, I);
            for(auto const & Ie : fac_p) {
                fmt::print("assert I_{}_{}.valuation({}) == {}\n",
                        a, b, Ie.first, 
                        I.valuation(Ie.first));
                fmt::print("print(\"{}\")\n", I.valuation(Ie.first));
            }
        }
    }
    prime_info_clear(pi);

    return 1;
}

static int do_linear_algebra_timings(cxx_param_list & pl)/*{{{*/
{
    unsigned int m = 8;
    unsigned int n = 5;
    param_list_parse_uint(pl, "m", &m);
    param_list_parse_uint(pl, "n", &n);

    gmp_randstate_t state;
    gmp_randinit_default(state);
    unsigned long seed = 0;
    if (param_list_parse_ulong(pl, "seed", &seed)) {
        gmp_randseed_ui(state, seed);
    }

    cxx_mpq_mat M(m,n);
    cxx_mpq_mat T(m,m);
    cxx_mpz_mat Mz(m,n);
    cxx_mpz_mat Tz(m,m);
    cxx_mpz p;

    mpz_set_ui(p, 19);
    param_list_parse_mpz(pl, "prime", p);

    {
        printf("\n\nCase 0.1\n\n");
        mpq_mat_urandomm(M, state, p);
        mpq_mat_fprint(stdout, M);
        printf("\n");
        mpq_mat_gauss_backend(M, T);
        mpq_mat_fprint(stdout, M);
        printf("\n");
        mpq_mat_fprint(stdout, T);
        printf("\n");
    }

    {
        printf("\n\nCase 0.2\n\n");
        mpz_mat_urandomm(Mz, state, p);
        mpz_mat_fprint(stdout, Mz);
        printf("\n");
        mpz_mat_gauss_backend_mod_mpz(Mz, Tz, p);
        mpz_mat_fprint(stdout, Mz);
        printf("\n");
        mpz_mat_fprint(stdout, Tz);
        printf("\n");
    }

    {
        printf("\n\nCase 1\n\n");
        mpz_mat_realloc(Mz, m, n);
        mpz_mat_urandomm(Mz, state, p);
        mpz_mat_fprint(stdout, Mz); printf("\n");
        double t = seconds();
        mpz_mat_hermite_form(Mz, Tz);
        t = seconds()-t;
        mpz_mat_fprint(stdout, Mz); printf("\n");
        mpz_mat_fprint(stdout, Tz); printf("\n");

        printf("%1.4f\n", t);
    }

    gmp_randclear(state);
    return 1;
}/*}}}*/

// coverity[root_function]
int main(int argc, char const * argv[])
    /*{{{ */
{
    cxx_param_list pl;

    param_list_configure_alias(pl, "prime", "p");

    decl_usage(pl);

    original_argv = argv;

    argv++,argc--;

    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        fprintf (stderr, "Unknown option: %s\n", argv[0]);
        usage(pl, original_argv, "unexpected argument");
    }

    std::string tmp;
    if (!param_list_parse(pl, "test", tmp))
        usage(pl, original_argv, "missing --test argument");

    int rc = 0; /* placate gcc */

    if (tmp == "p-maximal-order") {
        rc = do_p_maximal_order(pl);
    } else if (tmp == "p-maximal-order-batch") {
        rc = do_p_maximal_order_batch(pl);
    } else if (tmp == "factorization-of-prime") {
        rc = do_factorization_of_prime(pl);
    } else if (tmp == "factorization-of-prime-batch") {
        rc = do_factorization_of_prime_batch(pl);
    } else if (tmp == "valuations-of-ideal") {
        rc = do_valuations_of_ideal(pl);
    } else if (tmp == "valuations-of-ideal-batch") {
        rc = do_valuations_of_ideal_batch(pl);
    } else if (tmp == "linear-algebra-timings") {
        rc = do_linear_algebra_timings(pl);
    } else if (tmp == "nt-object-interface") {
        rc = do_number_theory_object_interface(pl);
    } else if (tmp == "maximal-order") {
        rc = do_maximal_order(pl);
    } else {
        usage(pl, original_argv, "unknown test");
    }
    return rc ? EXIT_SUCCESS : EXIT_FAILURE;

}

/*}}}*/
