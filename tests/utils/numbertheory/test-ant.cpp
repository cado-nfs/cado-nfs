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

static const char ** original_argv;

static void decl_usage(cxx_param_list & pl)/*{{{*/
{
    param_list_decl_usage(pl, "test", "which test to run");
    param_list_decl_usage(pl, "prime", "prime");
    param_list_decl_usage(pl, "poly", "polynomial (std::string)");
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
    number_field const K(f);
    number_field_order O = K.p_maximal_order(p);

    fmt::print("{}\n", O);

    return 1;
}
/*}}}*/

static cxx_mpq_mat batch_read_order_basis(std::istream & in, unsigned int n)/*{{{*/
{
    cxx_mpq_mat O;
    std::string keyword;
    if (!(in >> keyword) || keyword != "order")
        throw std::invalid_argument("Parse error");
    cxx_mpz_mat A(n, n);
    cxx_mpz d;
    if (!(in >> d))
        throw std::invalid_argument("Parse error");
    for(unsigned int i = 0 ; i < A->m ; i++) {
        for(unsigned int j = 0 ; j < A->n ; j++) {
            if (!(in >> mpz_mat_entry(A, i, j)))
                throw std::invalid_argument("Parse error");
        }
    }
    mpq_mat_set_mpz_mat_denom(O, A, d);
    return O;
}/*}}}*/

static std::vector<number_field_prime_ideal> batch_read_prime_factorization(std::istream & in, number_field_order const & O, cxx_mpz const & p) /*{{{*/
{
    number_field const & K(O.number_field());

    std::vector<number_field_prime_ideal> ideals;
    std::string keyword;
    if (!(in >> keyword) || keyword != "ideals")
        throw std::invalid_argument("Parse error");
    unsigned int nideals;
    // coverity[tainted_argument]
    if (!(in >> nideals))
        throw std::invalid_argument("Parse error");
    for(unsigned int k = 0 ; k < nideals ; k++) {
        std::vector<number_field_element> gens;
        gens.push_back(K(p));
        cxx_mpz den;
        int e;
        if (!(in >> den))
            throw std::invalid_argument("Parse error");
        cxx_mpq_mat A(1, K.degree());
        for(unsigned int j = 0 ; j < A->n ; j++) {
            cxx_mpz num;
            if (!(in >> num))
                throw std::invalid_argument("Parse error");
            mpq_ptr aa = mpq_mat_entry(A, 0, j);
            mpz_set(mpq_numref(aa), num);
            mpz_set(mpq_denref(aa), den);
            mpq_canonicalize(aa);
        }
        gens.push_back(K(A));
        if (!(in >> e))
            throw std::invalid_argument("Parse error");
        auto I = O.fractional_ideal(gens);
        ASSERT_ALWAYS(I.is_integral());
        ideals.emplace_back(I, p, e);
    }
    return ideals;
}/*}}}*/


static int do_p_maximal_order_batch(cxx_param_list & pl) /*{{{*/
{
    const char * tmp;

    if ((tmp = param_list_lookup_string(pl, "batch")) == nullptr)
        usage(pl, original_argv, "missing batch argument");

    std::ifstream is(tmp);
    std::string s;
    int nok = 0;
    int nfail = 0;

    for(int test = 0; getline(is, s, '\n') ; ) {
        if (s.empty()) continue;
        if (s[0] == '#') continue;
        const auto exc = std::string("Parse error on input ") + s;

        cxx_mpz_poly f;
        if (!mpz_poly_set_from_expression(f, s))
            usage(pl, original_argv, "cannot parse polynomial");


        if (!(getline(is, s, '\n')))
            throw std::invalid_argument(exc);
        cxx_mpz const d;
        cxx_mpz_mat const A(f->deg, f->deg);
        std::istringstream is1(s);
        cxx_mpz p;
        if (!(is1 >> p))
            throw std::invalid_argument(exc);

        number_field const K(f);
        number_field_order const O(K, batch_read_order_basis(is1, f->deg));
        number_field_order const my_O = K.p_maximal_order(p);
        bool const ok = O.equal_mod(my_O, p);

        std::cout << ok_NOK(ok) << " test " << test
            << " (degree " << f->deg << ", p=" << p << ")"
            << "\n";
        test++;
        nok += ok;
        nfail += !ok;
    }
    std::cout << nok << " tests passed";
    if (nfail) std::cout << ", " << nfail << " TESTS FAILED";
    std::cout << "\n";
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

    std::ifstream is(tmp);
    std::string s;
    int nok = 0;
    int nfail = 0;
    for(int test = 0; getline(is, s, '\n') ; ) {
        if (s.empty()) continue;
        if (s[0] == '#') continue;

        const auto exc = std::string("Parse error on input ") + s;

        cxx_mpz_poly f;
        if (!mpz_poly_set_from_expression(f, s))
            usage(pl, original_argv, "cannot parse polynomial");

        if (!(getline(is, s, '\n')))
            throw std::invalid_argument(exc);

        std::istringstream is1(s);

        cxx_mpz p;
        if (!(is1 >> p)) throw std::invalid_argument(exc);

        number_field const K(f);

        number_field_order const O(K, batch_read_order_basis(is1, f->deg));

        std::vector<number_field_prime_ideal> ideals = batch_read_prime_factorization(is1, O, p);

        number_field_order const my_O = K.p_maximal_order(p);

        bool ok = O.equal_mod(my_O, p);

        cxx_gmp_randstate state;
        unsigned long seed = 0;
        if (param_list_parse_ulong(pl, "seed", &seed)) {
            gmp_randseed_ui(state, seed);
        }

        /* What if we simply ask our code to compute the factorization of
         * p with respect to the order basis which was chosen by magma ?
         * */
        auto my_ideals = O.factor(p, state);

        // sort magma ideals. Ours are sorted already.
        std::sort(ideals.begin(), ideals.end());
        ok = ok && (ideals.size() == my_ideals.size());
        for(unsigned int k = 0 ; ok && k < ideals.size() ; k++) {
            ok = (ideals[k] == my_ideals[k].first);
        }
        std::cout << ok_NOK(ok) << " test " << test
            << " (degree " << f->deg << ", p=" << p << ")"
            << "\n";
        test++;
        nok += ok;
        nfail += !ok;
    }
    std::cout << nok << " tests passed";
    if (nfail) std::cout << ", " << nfail << " TESTS FAILED";
    std::cout << "\n";
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
    std::vector<cxx_mpz_poly> elements; 
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
            std::istringstream is(desc);
            cxx_mpz_poly x;
            if (!(is >> x)) usage(pl, original_argv, "cannot parse ideal generators");
            elements.push_back(x);
        }
    }

    number_field K(f);
    number_field_order O = K.p_maximal_order(p);

    K.bless("alpha");

    std::vector<number_field_element> gens;
    gens.reserve(elements.size());
    for(auto const & e : elements)
        gens.emplace_back(K(e));

    auto I = O.fractional_ideal(gens);

    cxx_gmp_randstate state;
    unsigned long seed = 0;
    if (param_list_parse_ulong(pl, "seed", &seed)) {
        gmp_randseed_ui(state, seed);
    }
    std::vector<std::pair<number_field_prime_ideal, int> > F = O.factor(p, state);

    for(unsigned int k = 0 ; k < F.size() ; k++) {
        auto const& fkp(F[k].first);
        number_field_prime_ideal::two_element two(fkp);
        int const e = F[k].second;
        int const v = fkp.valuation(I);

        fmt::print("# (p={}, k={}, f={}. e={}; ideal<O|{},{}>)^{};\n",
                p, k, fkp.inertia_degree(), e, two.first, K(two.second), v);
    }
    return 1;
}
/*}}}*/

static int do_valuations_of_ideal_batch(cxx_param_list & pl) /*{{{*/
{
    const char * tmp;

    if ((tmp = param_list_lookup_string(pl, "batch")) == nullptr)
        usage(pl, original_argv, "missing batch argument");

    std::ifstream is(tmp);
    std::string s;
    int nok = 0;
    int nfail = 0;

    for(int test = 0; getline(is, s, '\n') ; ) {
        if (s.empty()) continue;
        if (s[0] == '#') continue;

        const std::string exc = std::string("Parse error on input ") + s;

        cxx_mpz_poly f;
        if (!mpz_poly_set_from_expression(f, s))
            usage(pl, original_argv, "cannot parse polynomial");

        if (!(getline(is, s, '\n')))
            throw std::invalid_argument(exc);

        std::istringstream is1(s);

        cxx_mpz p;
        if (!(is1 >> p)) throw std::invalid_argument(exc);

        number_field const K(f);

        number_field_order const O(K, batch_read_order_basis(is1, f->deg));

        auto ideals = batch_read_prime_factorization(is1, O, p);

        number_field_order const my_O = K.p_maximal_order(p);

        bool ok = O.equal_mod(my_O, p);

        cxx_gmp_randstate state;
        unsigned long seed = 0;
        if (param_list_parse_ulong(pl, "seed", &seed)) {
            gmp_randseed_ui(state, seed);
        }

        /* What if we simply ask our code to compute the factorization of
         * p with respect to the order basis which was chosen by magma ?
         * */
        auto my_ideals = O.factor(p, state);

        /* compute matching table */
        ok = ok && (ideals.size() == my_ideals.size());
        std::vector<unsigned int> magma_to_mine(ideals.size());
        std::vector<unsigned int> mine_to_magma(my_ideals.size());
        for(unsigned int k = 0 ; ok && k < ideals.size() ; k++) {
            bool found = false;
            for(unsigned int ell = 0 ; ell < my_ideals.size() ; ell++) {
                if (ideals[k] == my_ideals[ell].first) {
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
            std::string keyword;
            if (!(is1 >> keyword) || keyword != "composite")
                throw std::invalid_argument(exc);
            int ngens;
            // coverity[-taint_source]
            if (!(is1 >> ngens)) throw std::invalid_argument(exc);
            if (ngens < 0) throw std::invalid_argument(exc);
            if (ngens == 0) break;

            std::vector<number_field_element> gens;
            for(int i = 0 ; i < ngens ; i++) {
                /* read n integers */
                cxx_mpq_mat A(1, K.degree());
                for(unsigned int j = 0 ; j < A->n ; j++) {
                    cxx_mpz num;
                    if (!(is1 >> num))
                        throw std::invalid_argument(exc);
                    mpq_ptr aa = mpq_mat_entry(A, 0, j);
                    mpz_set(mpq_numref(aa), num);
                    mpz_set_ui(mpq_denref(aa), 1);
                    // mpz_set(mpq_denref(aa), den);
                    mpq_canonicalize(aa);
                }
                gens.push_back(K(A));
            }
            auto I = O.fractional_ideal(gens);
            std::vector<int> my_vals;
            my_vals.reserve(my_ideals.size());
            for(auto const & Ie : my_ideals)
                my_vals.emplace_back(Ie.first.valuation(I));
            if (!(is1 >> keyword) || keyword != "valuations")
                throw std::invalid_argument(exc);
            for(unsigned int k = 0 ; ok && k < ideals.size() ; k++) {
                std::string s;
                is1 >> s;
                int v;
                if (s == "Infinity") {
                    v = INT_MAX;
                } else {
                    std::istringstream is2(s);
                    if (!(is2 >> v)) throw std::invalid_argument(exc);
                }
                ok = (v == my_vals[magma_to_mine[k]]);
            }
        }


        std::cout << ok_NOK(ok) << " test " << test
            << " (degree " << f->deg << ", p=" << p << ")"
            << "\n";
        test++;
        nok += ok;
        nfail += !ok;
    }
    std::cout << nok << " tests passed";
    if (nfail) std::cout << ", " << nfail << " TESTS FAILED";
    std::cout << "\n";
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
