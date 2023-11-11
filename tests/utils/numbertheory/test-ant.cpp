#include "cado.h" // IWYU pragma: keep
// IWYU pragma: no_include <ext/alloc_traits.h>
#include <cstring>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <cstdio>
#include <cstdlib>
#include <memory>              // for allocator_traits<>::value_type
#include <string>              // for string, getline, operator!=, operator>>
#include <utility>             // for pair, operator==, swap, make_pair
#include <vector>              // for vector
#include <gmp.h>               // for mpz_t, operator>>, gmp_randclear, gmp_...
#include "cxx_mpz.hpp"
#include "gmp_aux.h"           // for mpz_p_valuation
#include "macros.h"
#include "mpz_mat.h"           // for cxx_mpz_mat, cxx_mpq_mat, mpz_mat_fprint
#include "mpz_poly.h"
#include "numbertheory.hpp"
#include "params.h"
#include "portability.h" // strdup // IWYU pragma: keep
#include "timing.h"     // seconds
#include "misc.h"

using namespace std;
static char ** original_argv;

namespace straightforward_poly_io {
istream& operator>>(istream& is, cxx_mpz_poly& f)/*{{{*/
{
    vector<cxx_mpz> v;
    cxx_mpz a;
    for( ; is >> a ; v.push_back(a)) ;
    mpz_poly_realloc(f, v.size());
    for(unsigned int i = 0 ; i < v.size() ; i++) {
        mpz_set(f->coeff[i], v[i]);
    }
    mpz_poly_cleandeg(f, v.size()-1);
    is.clear();
    return is;
}
/*}}}*/

ostream& operator<<(ostream& o, cxx_mpz_poly const& v)/*{{{*/
{
    /* note that we can't cheat and use cxx_mpz here */
    ostream_iterator<mpz_t> it(o, " ");
    if (v->deg>=0) {
        copy(v->coeff, v->coeff + v->deg, it);
        o << v->coeff[v->deg];
    } else {
        o << "0";
    }
    return o;
}/*}}}*/
}

static void decl_usage(param_list_ptr pl)/*{{{*/
{
    param_list_decl_usage(pl, "test", "which test to run");
    param_list_decl_usage(pl, "prime", "prime");
    param_list_decl_usage(pl, "polystr", "polynomial (string)");
    param_list_decl_usage(pl, "polyfile", "polynomial (file)");
    param_list_decl_usage(pl, "out", "output file");
    param_list_decl_usage(pl, "batch", "batch input file with test vectors and expected results");
    param_list_decl_usage(pl, "seed", "seed used for random picks");
    param_list_decl_usage(pl, "elements", "ideal generators (separated by ;)");
}/*}}}*/

void usage(param_list_ptr pl, char ** argv, const char * msg = NULL)/*{{{*/
{
    param_list_print_usage(pl, argv[0], stderr);
    if (msg) {
        fprintf(stderr, "%s\n", msg);
    }
    exit(EXIT_FAILURE);
}/*}}}*/

struct iowrap {/*{{{*/
    streambuf * ibuf, * obuf;
private:
    ifstream ifs;
    ofstream ofs;
    istringstream iss;
public:
    iowrap(param_list_ptr pl) {
        const char * tmp;
        if ((tmp = param_list_lookup_string(pl, "polystr")) != NULL) {
            iss.str(tmp);
            ibuf = iss.rdbuf();
        } else if ((tmp = param_list_lookup_string(pl, "polyfile")) != NULL) {
            if (strcmp(tmp, "-") == 0) {
                ibuf = cin.rdbuf();
            } else {
                ifs.open(tmp);
                ibuf = ifs.rdbuf();
            }
        } else {
            usage(pl, original_argv, "Please provide either --polyfile or --polystr");
        }
        if ((tmp = param_list_lookup_string(pl, "out")) != NULL && strcmp(tmp, "-") != 0) {
            ofs.open(tmp);
            obuf = ofs.rdbuf();
        } else {
            obuf = cout.rdbuf();
        }
    }
};/*}}}*/

int do_p_maximal_order(param_list_ptr pl) /*{{{*/
{
    cxx_mpz_poly f;
    cxx_mpz p;

    if (!param_list_parse_mpz(pl, "prime", p)) usage(pl, original_argv, "missing prime argument");

    iowrap io(pl);
    istream in(io.ibuf);
    ostream out(io.obuf);

    if (!(straightforward_poly_io::operator>>(in, f))) usage(pl, original_argv, "cannot parse polynomial");
    cxx_mpq_mat M = p_maximal_order(f, p);
    cxx_mpz D;
    cxx_mpz_mat A;
    mpq_mat_numden(A, D, M);
    out << "1/" << D << "*\n" << A << endl;

    return 1;
}
/*}}}*/

bool sl_equivalent_matrices(cxx_mpq_mat const& M, cxx_mpq_mat const& A, cxx_mpz const& p)/*{{{*/
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

cxx_mpq_mat batch_read_order_basis(istream & in, unsigned int n)/*{{{*/
{
    invalid_argument exc(string("Parse error"));
    cxx_mpq_mat O;
    string keyword;
    if (!(in >> keyword) || keyword != "order") throw exc;
    cxx_mpz_mat A(n, n);
    cxx_mpz d;
    if (!(in >> d))
        throw exc;
    for(unsigned int i = 0 ; i < A->m ; i++) {
        for(unsigned int j = 0 ; j < A->n ; j++) {
            if (!(in >> mpz_mat_entry(A, i, j)))
                throw exc;
        }
    }
    mpq_mat_set_mpz_mat_denom(O, A, d);
    return O;
}/*}}}*/

vector<pair<cxx_mpz_mat, int> > batch_read_prime_factorization(istream & in, unsigned int n, cxx_mpz const& p, cxx_mpq_mat const& O, cxx_mpz_mat const& M)/*{{{*/
{
    invalid_argument exc(string("Parse error"));
    vector<pair<cxx_mpz_mat, int> > ideals;
    string keyword;
    if (!(in >> keyword) || keyword != "ideals") throw exc;
    unsigned int nideals;
    // coverity[tainted_argument]
    if (!(in >> nideals))
        throw exc;
    for(unsigned int k = 0 ; k < nideals ; k++) {
        cxx_mpq_mat A(2, n);
        cxx_mpz den;
        int e;
        mpq_set_z(mpq_mat_entry(A,0,0),p);
        if (!(in >> den)) throw exc;
        for(unsigned int j = 0 ; j < A->n ; j++) {
            cxx_mpz num;
            if (!(in >> num)) throw exc;
            mpq_ptr aa = mpq_mat_entry(A, 1, j);
            mpz_set(mpq_numref(aa), num);
            mpz_set(mpq_denref(aa), den);
            mpq_canonicalize(aa);
        }
        if (!(in >> e)) throw exc;
        pair<cxx_mpz_mat, cxx_mpz> Id = generate_ideal(O,M,A);
        ASSERT_ALWAYS(mpz_cmp_ui(Id.second, 1) == 0);
        ideals.push_back(make_pair(Id.first,e));
    }
    return ideals;
}/*}}}*/


int do_p_maximal_order_batch(param_list_ptr pl) /*{{{*/
{
    cxx_mpz_poly f;
    const char * tmp;

    if ((tmp = param_list_lookup_string(pl, "batch")) == NULL)
        usage(pl, original_argv, "missing batch argument");

    ifstream is(tmp);
    string s;
    int nok = 0;
    int nfail = 0;

    for(int test = 0; getline(is, s, '\n') ; ) {
        if (s.empty()) continue;
        if (s[0] == '#') continue;
        invalid_argument exc(string("Parse error on input") + s);
        istringstream is0(s);
        if (!(straightforward_poly_io::operator>>(is0, f))) usage(pl, original_argv, "cannot parse polynomial");

        if (!(getline(is, s, '\n')))
            throw exc;
        cxx_mpz d;
        cxx_mpz_mat A(f->deg, f->deg);
        istringstream is1(s);
        cxx_mpz p;
        if (!(is1 >> p))
            throw exc;

        cxx_mpq_mat O = batch_read_order_basis(is1, f->deg);
        cxx_mpq_mat my_O = p_maximal_order(f, p);

        bool ok = sl_equivalent_matrices(O, my_O, p);

        cout << ok_NOK(ok) << " test " << test
            << " (degree " << f->deg << ", p=" << p << ")"
            << endl;
        test++;
        nok += ok;
        nfail += !ok;
    }
    cout << nok << " tests passed";
    if (nfail) cout << ", " << nfail << " TESTS FAILED";
    cout << endl;
    return nfail == 0;
}
/*}}}*/

int do_factorization_of_prime(param_list_ptr pl) /*{{{*/
{
    cxx_mpz_poly f;
    cxx_mpz p;

    if (!param_list_parse_mpz(pl, "prime", p)) usage(pl, original_argv, "missing prime argument");

    iowrap io(pl);
    istream in(io.ibuf);
    ostream out(io.obuf);

    if (!(straightforward_poly_io::operator>>(in, f))) usage(pl, original_argv, "cannot parse polynomial");
    cxx_mpq_mat M = p_maximal_order(f, p);
    gmp_randstate_t state;
    gmp_randinit_default(state);
    unsigned long seed = 0;
    if (param_list_parse_ulong(pl, "seed", &seed)) {
        gmp_randseed_ui(state, seed);
    }
    vector<pair<cxx_mpz_mat, int> > F = factorization_of_prime(M, f, p, state);
    gmp_randclear(state);
    out << F.size() << "\n";
    for(unsigned int k = 0 ; k < F.size() ; k++) {
        out << F[k].first << " " << F[k].second << endl;
    }
    return 1;
}
/*}}}*/

int do_factorization_of_prime_batch(param_list_ptr pl) /*{{{*/
{
    cxx_mpz_poly f;
    const char * tmp;

    if ((tmp = param_list_lookup_string(pl, "batch")) == NULL)
        usage(pl, original_argv, "missing batch argument");

    ifstream is(tmp);
    string s;
    int nok = 0;
    int nfail = 0;
    for(int test = 0; getline(is, s, '\n') ; ) {
        if (s.empty()) continue;
        if (s[0] == '#') continue;

        invalid_argument exc(string("Parse error on input") + s);

        istringstream is0(s);
        if (!(straightforward_poly_io::operator>>(is0, f))) usage(pl, original_argv, "cannot parse polynomial");

        if (!(getline(is, s, '\n')))
            throw exc;

        istringstream is1(s);

        cxx_mpz p;
        if (!(is1 >> p)) throw exc;

        string keyword;
        cxx_mpq_mat O = batch_read_order_basis(is1, f->deg);
        cxx_mpz_mat M = multiplication_table_of_order(O, f);

        vector<pair<cxx_mpz_mat, int> > ideals = batch_read_prime_factorization(is1, f->deg, p, O, M);

        cxx_mpq_mat my_O = p_maximal_order(f, p);

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
            << endl;
        test++;
        nok += ok;
        nfail += !ok;
    }
    cout << nok << " tests passed";
    if (nfail) cout << ", " << nfail << " TESTS FAILED";
    cout << endl;
    return nfail == 0;
}
/*}}}*/

int do_valuations_of_ideal(param_list_ptr pl) /*{{{*/
{
    cxx_mpz_poly f;
    cxx_mpz p;

    if (!param_list_parse_mpz(pl, "prime", p)) usage(pl, original_argv, "missing prime argument");

    iowrap io(pl);
    istream in(io.ibuf);
    ostream out(io.obuf);

    /* Now read the element description */
    vector<cxx_mpz_poly> elements; 
    {
        const char * tmp;
        if ((tmp = param_list_lookup_string(pl, "elements")) == NULL)
            usage(pl, original_argv, "missing ideal generators");
        char * desc = strdup(tmp);
        for(char * q = desc, * nq ; q ; q = nq) {
            nq = strchr(desc, ';');
            if (nq) *nq++='\0';
            istringstream is(q);
            cxx_mpz_poly x;
            if (!(is >> x)) usage(pl, original_argv, "cannot parse ideal generators");
            elements.push_back(x);
        }
    }


    if (!(straightforward_poly_io::operator>>(in, f))) usage(pl, original_argv, "cannot parse polynomial");
    cxx_mpq_mat O = p_maximal_order(f, p);
    cxx_mpz_mat M = multiplication_table_of_order(O, f);

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
                mpz_mul(denom, denom, f->coeff[f->deg]);
                mpz_mul(e->coeff[k], e->coeff[k], denom);
            }
            mpz_poly_div_r_z(e, e, fh);
            for(int k = 0 ; k <= e->deg ; k++) {
                mpz_mul(e->coeff[k], e->coeff[k], c);
                mpz_mul(c, c, f->coeff[f->deg]);
                mpq_ptr gik = mpq_mat_entry(generators, i, k);
                mpz_set(mpq_numref(gik), e->coeff[k]);
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
        cxx_mpz_mat a = valuation_helper_for_ideal(M, fkp, p);
        pair<cxx_mpz, cxx_mpz_mat> two = prime_ideal_two_element(O, f, M, fkp);

        cxx_mpq_mat theta_q;
        {
            mpq_mat_set_mpz_mat(theta_q, two.second);
            mpq_mat_mul(theta_q, theta_q, O);
        }
        string uniformizer = write_element_as_polynomial(theta_q, "alpha");

        int e = F[k].second;
        int v = valuation_of_ideal_at_prime_ideal(M, Id, a, e, p);

        cout << "# (p=" << p
            << ", k=" << k
            << ", f="<< prime_ideal_inertia_degree(fkp)
            << ", e="<< e
            << "; "
            << "ideal<O|" << two.first << "," << uniformizer << ">"
            << ")"
            << "^" << v
            << ";" << endl;
    }
    return 1;
}
/*}}}*/

int do_valuations_of_ideal_batch(param_list_ptr pl) /*{{{*/
{
    cxx_mpz_poly f;
    const char * tmp;

    if ((tmp = param_list_lookup_string(pl, "batch")) == NULL)
        usage(pl, original_argv, "missing batch argument");

    ifstream is(tmp);
    string s;
    int nok = 0;
    int nfail = 0;

    for(int test = 0; getline(is, s, '\n') ; ) {
        if (s.empty()) continue;
        if (s[0] == '#') continue;

        invalid_argument exc(string("Parse error on input") + s);

        istringstream is0(s);
        if (!(straightforward_poly_io::operator>>(is0, f))) usage(pl, original_argv, "cannot parse polynomial");

        if (!(getline(is, s, '\n')))
            throw exc;

        istringstream is1(s);

        cxx_mpz p;
        if (!(is1 >> p)) throw exc;

        string keyword;
        cxx_mpq_mat O = batch_read_order_basis(is1, f->deg);
        cxx_mpz_mat M = multiplication_table_of_order(O, f);

        vector<pair<cxx_mpz_mat, int> > ideals = batch_read_prime_factorization(is1, f->deg, p, O, M);

        cxx_mpq_mat my_O = p_maximal_order(f, p);

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
            if (!(is1 >> keyword) || keyword != "composite") throw exc;
            int ngens;
            // coverity[-taint_source]
            if (!(is1 >> ngens)) throw exc;
            if (ngens < 0) throw exc;
            if (ngens == 0) break;
            cxx_mpz_mat gens(ngens, f->deg);
            for(unsigned int i = 0 ; i < gens->m ; i++) {
                for(unsigned int j = 0 ; j < gens->n ; j++) {
                    if (!(is1 >> mpz_mat_entry(gens, i, j)))
                        throw exc;
                }
            }
            pair<cxx_mpz_mat, cxx_mpz> Id = generate_ideal(O, M, cxx_mpq_mat(gens));
            vector<int> my_vals;
            for(unsigned int ell = 0 ; ell < my_ideals.size() ; ell++) {
                cxx_mpz_mat const& fkp(my_ideals[ell].first);
                int e = my_ideals[ell].second;
                cxx_mpz_mat a = valuation_helper_for_ideal(M, fkp, p);
                int v = valuation_of_ideal_at_prime_ideal(M, Id, a, e, p);
                my_vals.push_back(v);
            }
            if (!(is1 >> keyword) || keyword != "valuations") throw exc;
            for(unsigned int k = 0 ; ok && k < ideals.size() ; k++) {
                std::string s;
                is1 >> s;
                int v;
                if (s == "Infinity") {
                    v = INT_MAX;
                } else {
                    std::istringstream is2(s);
                    if (!(is2 >> v)) throw exc;
                }
                ok = (v == my_vals[magma_to_mine[k]]);
            }
        }


        cout << ok_NOK(ok) << " test " << test
            << " (degree " << f->deg << ", p=" << p << ")"
            << endl;
        test++;
        nok += ok;
        nfail += !ok;
    }
    cout << nok << " tests passed";
    if (nfail) cout << ", " << nfail << " TESTS FAILED";
    cout << endl;
    return nfail == 0;
}
/*}}}*/

int do_linear_algebra_timings(param_list_ptr pl)/*{{{*/
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

    if (0) {
        printf("\n\nCas 0.1\n\n");
        mpq_mat_urandomm(M, state, p);
        mpq_mat_fprint(stdout, M);
        printf("\n");
        mpq_mat_gauss_backend(M, T);
        mpq_mat_fprint(stdout, M);
        printf("\n");
        mpq_mat_fprint(stdout, T);
        printf("\n");
    }

    if (0) {
        printf("\n\nCas 0.2\n\n");
        mpz_mat_urandomm(Mz, state, p);
        mpz_mat_fprint(stdout, Mz);
        printf("\n");
        mpz_mat_gauss_backend_mod_mpz(Mz, Tz, p);
        mpz_mat_fprint(stdout, Mz);
        printf("\n");
        mpz_mat_fprint(stdout, Tz);
        printf("\n");
    }

    if (1) {
        printf("\n\nCas 1\n\n");
        mpz_mat_realloc(Mz, m, n);
        mpz_mat_urandomm(Mz, state, p);
        mpz_mat_fprint(stdout, Mz); printf("\n");
        double t = seconds();
        mpz_mat_hnf_backend(Mz, Tz);
        t = seconds()-t;
        mpz_mat_fprint(stdout, Mz); printf("\n");
        mpz_mat_fprint(stdout, Tz); printf("\n");

        printf("%1.4f\n", t);
    }

    gmp_randclear(state);
    return 1;
}/*}}}*/

// coverity[root_function]
int main(int argc, char *argv[]) /*{{{ */
{
    cxx_param_list pl;

    param_list_configure_alias(pl, "prime", "p");
    param_list_configure_alias(pl, "polystr", "f");
    param_list_configure_alias(pl, "out", "o");

    decl_usage(pl);

    original_argv = argv;

    argv++,argc--;

    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        fprintf (stderr, "Unknown option: %s\n", argv[0]);
        usage(pl, original_argv, "unexpected argument");
    }

    const char * tmp = param_list_lookup_string(pl, "test");
    if (!tmp) usage(pl, original_argv, "missing --test argument");

    int rc = 0; /* placate gcc */

    if (strcmp(tmp, "p-maximal-order") == 0) {
        rc = do_p_maximal_order(pl);
    } else if (strcmp(tmp, "p-maximal-order-batch") == 0) {
        rc = do_p_maximal_order_batch(pl);
    } else if (strcmp(tmp, "factorization-of-prime") == 0) {
        rc = do_factorization_of_prime(pl);
    } else if (strcmp(tmp, "factorization-of-prime-batch") == 0) {
        rc = do_factorization_of_prime_batch(pl);
    } else if (strcmp(tmp, "valuations-of-ideal") == 0) {
        rc = do_valuations_of_ideal(pl);
    } else if (strcmp(tmp, "valuations-of-ideal-batch") == 0) {
        rc = do_valuations_of_ideal_batch(pl);
    } else if (strcmp(tmp, "linear-algebra-timings") == 0) {
        rc = do_linear_algebra_timings(pl);
    } else {
        usage(pl, original_argv, "unknown test");
    }
    return rc ? EXIT_SUCCESS : EXIT_FAILURE;

}

/*}}}*/
