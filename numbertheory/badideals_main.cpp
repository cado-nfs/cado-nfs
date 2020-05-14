#include "cado.h"       // IWYU pragma: keep
#include <iostream>
#include <sstream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include "badideals.hpp"

using namespace std;

char ** original_argv;
gmp_randstate_t state;

/* This program is intended to replicate exactly the behaviour of the
 * scripts/badideals.mag program of old.
 *
 * Now, we acknowledge several bizarre things in the form of the
 * .badideals and .badidealinfo file.
 *
 * For reference, here is a test case:
 *
 * > clear; DEBUG:=1; load "scripts/badideals.mag";
 * Loading "scripts/badideals.mag"
 * 2,2:0: 2
 * d,5:0: 2
 * 2,1:1: 2
 * 2,2:1: 2
 * 3,1:1: 2
 * # bad ideals for poly[0]=240*x^4 + 4327846*x^3 - 11463289949*x^2 -
 * 48524924823332*x + 99623823815957215
 * 2 2 4 0 1 -1
 * 2 2 6 0 -1 1
 * 13 2 18 0 1 1
 * 13 2 31 0 1 1
 * 13 2 44 0 1 1
 * 13 2 57 0 1 1
 * 13 2 70 0 1 1
 * 13 2 83 0 1 1
 * 13 2 109 0 1 1
 * 13 2 122 0 1 1
 * 13 2 135 0 1 1
 * 13 2 148 0 1 1
 * 13 2 161 0 1 1
 * 13 2 5 0 1 -1
 * 13 2 96 0 -1 1
 * # bad ideals for poly[1]=840*x^5 - 16210*x^4 + 2610983*x^3 -
 * 2560484269*x^2 - 34656366009*x + 253976285986533
 * 2 1 1 1 1 -1
 * 2 2 4 1 -1 1
 * 2 2 6 1 1 -1
 * 3 2 7 1 1 1
 * 3 2 1 1 1 -1
 * 3 2 4 1 -1 1
 *
 *
 * The bizarre points are the following:
 *
 *  - That we have both a .badideals and a .badidealinfo file is very
 *    odd.
 *    Actually maybe this info would be just as well saved in the .roots file,
 *    perhaps.
 *  - The "side" info in the .badidealinfo file is weird.
 *  - The 13 lines above are inelegant. A more concise way would be:
 *        13 2 5 0 1 -1
 *        13 2 96 0 -1 1
 *        13 1 5 0 1 1
 *    But that would be non commutative, so this would be a clear change
 *    in semantics.
 *  - The file format, for now, is not ready to handle larger-degree
 *    prime ideals.
 */ 

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

void badideals_declare_usage(param_list_ptr pl)/*{{{*/
{
    param_list_decl_usage(pl, "badideals", "badideals file");
    param_list_decl_usage(pl, "badidealinfo", "badidealinfo file");
    param_list_decl_usage(pl, "polystr", "polynomial (string)");
    param_list_decl_usage(pl, "poly", "polynomial file");
    param_list_decl_usage(pl, "seed", "random seed");
    param_list_decl_usage(pl, "ell", "ell (for computing default number of maps)");
}/*}}}*/

void usage(param_list_ptr pl, char ** argv, const char * msg = NULL)/*{{{*/
{
    param_list_print_usage(pl, argv[0], stderr);
    if (msg) {
        fprintf(stderr, "%s\n", msg);
    }
    exit(EXIT_FAILURE);
}/*}}}*/

int main(int argc, char * argv[])
{
    param_list pl;
    param_list_init(pl);

    badideals_declare_usage(pl);
    param_list_configure_alias(pl, "polystr", "f");

    original_argv = argv;

    argv++,argc--;
    /* switches, if any. See below */
    /* aliases, if any. See below */

    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        /* Do perhaps some other things on the argument that haven't
         * been eaten at all. Like check whether it is a valid file to
         * source in order to get more options. See
         * param_list_read_stream and param_list_read_file for that. */
        fprintf (stderr, "Unknown option: %s\n", argv[0]);
        usage(pl, original_argv);
    }

    gmp_randinit_default(state);
    unsigned long seed = 1;
    if (param_list_parse_ulong(pl, "seed", &seed)) {
        gmp_randseed_ui(state, seed);
    }

    typedef vector<badideal>::const_iterator vbci_t;

    const char * tmp;
    if ((tmp = param_list_lookup_string(pl, "polystr")) != NULL) {
        int side = 0;
        cxx_mpz_poly f;
        string stmp(tmp);
        for(unsigned int i = 0 ; i < stmp.size() ; i++) {
            if (stmp[i]==',') stmp[i]=' ';
        }
        istringstream is(stmp);
        if (!(straightforward_poly_io::operator>>(is, f)))
            usage(pl, original_argv, "cannot parse polynomial");

        vector<badideal> badideals = badideals_for_polynomial(f, side);
        cout << "--- .badideals data ---\n";
        for(vbci_t it = badideals.begin() ; it != badideals.end() ; it++) {
            badideal const& b(*it);
            b.print_dot_badideals_file(cout, side);
        }

        cout << "--- .badidealinfo data ---\n";
        cout << "# bad ideals for poly"<<side<<"=" << f.print_poly("x") << endl;
        for(vbci_t it = badideals.begin() ; it != badideals.end() ; it++) {
            badideal const& b(*it);
            b.print_dot_badidealinfo_file(cout, side);
        }
    } else if ((tmp = param_list_lookup_string(pl, "poly")) != NULL) {
        cado_poly cpoly;
        cado_poly_init(cpoly);
        cado_poly_read(cpoly, tmp);
        const char * fbname = param_list_lookup_string(pl, "badideals");
        const char * fbiname = param_list_lookup_string(pl, "badidealinfo");
        if (!fbname || !fbiname)
            usage(pl, original_argv, "-poly requires both -badideals and -badidealinfo");

        ofstream fb(fbname);
        ofstream fbi(fbiname);

        for(int side = 0 ; side < cpoly->nb_polys ; side++) {
            cxx_mpz_poly f(cpoly->pols[side]);
            if (f->deg == 1) continue;
            vector<badideal> badideals = badideals_for_polynomial(f, side);
            for(vbci_t it = badideals.begin() ; it != badideals.end() ; it++) {
                badideal const& b(*it);
                b.print_dot_badideals_file(fb, side);
            }

            fbi << "# bad ideals for poly"<<side<<"=" << f.print_poly("x") << endl;
            for(vbci_t it = badideals.begin() ; it != badideals.end() ; it++) {
                badideal const& b(*it);
                b.print_dot_badidealinfo_file(fbi, side);
            }
        }
        cxx_mpz ell;
        if (param_list_parse_mpz(pl, "ell", ell)) {
            for(int side = 0 ; side < cpoly->nb_polys ; side++) {
                sm_side_info sm;
                sm_side_info_init(sm, cpoly->pols[side], ell);
                cout << "# nmaps" << side << " " << sm->nsm << endl;
                sm_side_info_clear(sm);
            }
        }
        cado_poly_clear(cpoly);
    } else {
        usage(pl, original_argv, "-poly or -polystr are mandatory");
    }

    gmp_randclear(state);
    param_list_clear(pl);
}
