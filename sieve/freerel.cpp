/*
 * Program: free relations
 * Original author : F. Morain
 * Purpose: creating free relations in a suitable format
 * Modified / rewritten by C. Bouvier (and others)
 * Multi-thread code by A. Filbois

This file is part of CADO-NFS.

CADO-NFS is free software; you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 2.1 of the License, or (at your option)
any later version.

CADO-NFS is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
details.

You should have received a copy of the GNU Lesser General Public License
along with CADO-NFS; see the file COPYING.  If not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include "cado.h" // IWYU pragma: keep
#include <algorithm>
#include <vector>
#include <memory>        // for unique_ptr, allocator_traits<>::value_type
#include <utility>       // for pair
#include <cstdio>       // fprintf
#include <cstdlib>       // exit
#include <climits>       // CHAR_BIT
#include "cado_poly.h"   // for cxx_cado_poly, cado_poly_s, cado_poly_read
#include "gzip.h"       // ofstream_maybe_compressed
#include "macros.h"      // for ASSERT_ALWAYS
#include "mpz_poly.h"   // mpz_poly_srcptr
#include "omp_proxy.h"
#include "params.h"      // for cxx_param_list, param_list_decl_usage, param...
#include "renumber.hpp" // renumber_t
#include "typedefs.h"    // for index_t, p_r_values_t, PRid, SIZEOF_INDEX, PRpr
#include "verbose.h"    // verbose_interpret_parameters
#include "fmt/core.h"
#include "fmt/format.h"

char * argv0;

static void
usage(cxx_param_list & pl, char* argv0)
{
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
}

struct freerel_data_t : public renumber_t::hook {
    ofstream_maybe_compressed sink;
    unsigned long pmin = 2;
    unsigned long pmax = 0;
    unsigned long nfree = 0;
    bool print_this_p(unsigned long p) const {
        return pmin <= p && p <= pmax;
    }
    freerel_data_t(cxx_param_list & pl, cxx_cado_poly const & cpoly, std::vector<unsigned int> const & lpb);
    void operator()(renumber_t & R, p_r_values_t p, index_t idx, renumber_t::cooked const & C) override;
    static void declare_usage(cxx_param_list & pl) {
        param_list_decl_usage(pl, "out", "output file for free relations");
        param_list_decl_usage(pl, "pmin", "do not create freerel below this bound");
        param_list_decl_usage(pl, "pmax", "do not create freerel beyond this bound");
    }
    static void lookup_parameters(cxx_param_list & pl) {
        param_list_lookup_string(pl, "out");
        param_list_lookup_string(pl, "pmin");
        param_list_lookup_string(pl, "pmax");
    }
    ~freerel_data_t() override { }
};

freerel_data_t::freerel_data_t(cxx_param_list & pl, cxx_cado_poly const & cpoly, std::vector<unsigned int> const & lpb) : sink(param_list_lookup_string(pl, "out"))
{
    param_list_parse_ulong(pl, "pmin", &pmin);
    param_list_parse_ulong(pl, "pmax", &pmax);
    /* if pmax is not equal to 0 (i.e., was not given on the command line),
     * set pmax to the largest integer that is less than or equal to
     * two of the large prime bounds.
     */
    if (!pmax && cpoly->nb_polys > 1) {
        std::vector<unsigned int> lpb_copy = lpb;
        std::sort(lpb_copy.begin(), lpb_copy.end());
        pmax = 1UL << lpb_copy[cpoly->nb_polys-2];
    }
}

void freerel_data_t::operator()(renumber_t & R, p_r_values_t p, index_t idx, renumber_t::cooked const & C)
{
    std::vector<std::pair<int, index_t>> full_sides;

    if (!print_this_p(p)) return;

    for(int side = 0 ; side < R.get_nb_polys() ; ++side) {
        mpz_poly_srcptr f = R.get_poly(side);
        /* Check if p corresponds to free relations.
           For 2 polynomials, this happens when the number of roots modulo p
           is maximal for both polynomials.
           For 3 or more polynomials, let t be the number of polynomials
           such that the number of roots modulo p is maximal. Then the
           number of free relations is t-1 (see Section 4.6 from the PhD
           thesis from Marije Elkenbracht-Huizing, "Factoring integers
           with the Number Field Sieve") */
        if ((int) C.nroots[side] == f->deg)
            full_sides.emplace_back(side, idx);

        idx += C.nroots[side];
    }

    if (full_sides.size() > 1) {
        for(size_t i = 1 ; i < full_sides.size() ; i++) {
            /* print a new free relation */
            sink << fmt::format(FMT_STRING("{:x},0:"), p);
            int side0 = full_sides[i-1].first;
            index_t i0 = full_sides[i-1].second;
            unsigned int n0 = C.nroots[side0];
            int side1 = full_sides[i].first;
            index_t i1 = full_sides[i].second;
            unsigned int n1 = C.nroots[side1];
            bool first = true;
            sink << std::hex;
            for(unsigned int k = 0 ; k < n0 ; k++, first=false) {
                if (!first) sink << ',';
                sink << i0 + k;
            }
            for(unsigned int k = 0 ; k < n1 ; k++, first=false) {
                if (!first) sink << ',';
                sink << i1 + k;
            }
            sink << '\n';
            nfree++;
        }
    }
}

static void
declare_usage(param_list pl)
{
    param_list_decl_usage(pl, "poly", "input polynomial file");
    param_list_decl_usage(pl, "lpb0", "large primes bound on side 0");
    param_list_decl_usage(pl, "lpb1", "large primes bound on side 1");
    param_list_decl_usage(pl,
                          "lpbs",
                          "large primes bounds (comma-separated list) "
                          "(for MNFS)");
    param_list_decl_usage(pl, "t", "number of threads");
    param_list_decl_usage(pl,
                          "dl",
                          "Add ideals for the leading "
                          "coeffs of the polynomials (for DL)");

}

// coverity[root_function]
int
main(int argc, char* argv[])
{
    argv0 = argv[0];
    cxx_param_list pl;
    cxx_cado_poly cpoly;
    int for_dl = 0;

    /* {{{ parse cmdline */
    declare_usage(pl);
    freerel_data_t::declare_usage(pl);
    verbose_decl_usage(pl);
    renumber_t::builder_declare_usage(pl);
    param_list_configure_switch(pl, "dl", &for_dl);

    argv++, argc--;
    if (argc == 0)
        usage(pl, argv0);
    for (; argc;) {
        if (param_list_update_cmdline(pl, &argc, &argv))
            continue;
        FILE* f;
        if ((f = fopen(argv[0], "r")) != NULL) {
            param_list_read_stream(pl, f, 0);
            fclose(f);
            argv++, argc--;
            continue;
        }
        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        usage(pl, argv0);
    }
    /* print command-line arguments */
    verbose_interpret_parameters(pl);
    param_list_print_command_line(stdout, pl);
    fflush(stdout);
    /* }}} */
    /* {{{ interpret cmdline parameters and catch errors */
    const char * polyfilename = param_list_lookup_string(pl, "poly");

    renumber_t::builder_lookup_parameters(pl);

    int nthreads = 0;
    if (param_list_parse_int (pl, "t", &nthreads)) {
        fprintf(stderr, "Warning: the -t argument to freerel is kept for compatibility, but you should rather take it out and let openmp deal with it\n");
        omp_set_num_threads(nthreads);
    }

    freerel_data_t::lookup_parameters(pl);

    param_list_lookup_string(pl, "lpb0");
    param_list_lookup_string(pl, "lpb1");
    param_list_lookup_string(pl, "lpbs");


    if (param_list_warn_unused(pl))
        usage(pl, argv0);

    if (polyfilename == NULL) {
        fprintf(stderr, "Error, missing -poly command line argument\n");
        usage(pl, argv0);
    }
    if (param_list_lookup_string(pl, "renumber") == NULL) {
        fprintf(stderr, "Error, missing -renumber command line argument\n");
        usage(pl, argv0);
    }
    if (!cado_poly_read(cpoly, polyfilename)) {
        fprintf(stderr, "Error reading polynomial file\n");
        exit(EXIT_FAILURE);
    }

    std::vector<unsigned int> lpb(cpoly->nb_polys, 0);

    if (!param_list_parse_uint_args_per_side(pl, "lpb", lpb.data(), cpoly->nb_polys, ARGS_PER_SIDE_DEFAULT_COPY_PREVIOUS)) {
        fprintf(stderr,
                "Error, could not obtain values for the lpb bounds (or not for all polynomials)\n");
        usage(pl, argv0);
    }

    /* }}} */

    std::unique_ptr<freerel_data_t> F;

    if (param_list_lookup_string(pl, "out")) {
        F.reset(new freerel_data_t(pl, cpoly, lpb));

        printf("Considering freerels for %lu <= p <= %lu\n", 
                F->pmin,
                F->pmax);
        fflush(stdout);
    }

    renumber_t renumber_table(cpoly);
    renumber_table.set_lpb(lpb);

    int max_lpb = renumber_table.get_max_lpb();
    if (max_lpb >= (int) sizeof(unsigned long) * CHAR_BIT) {
      fprintf (stderr, "Error, cannot handle lpb >= %zu (max(lpbs) is %d)\n",
                       sizeof(unsigned long) * CHAR_BIT, max_lpb);
      abort();
    }
    unsigned long lpbmax = 1UL << max_lpb;
    printf("Generating renumber table for 2 <= p <= %lu\n", lpbmax);

    /* This reads the options:
     *
     * renumber
     *
     */
    index_t R_max_index = renumber_table.build(pl, for_dl, F.get());

    if (F.get()) {
        /* /!\ Needed by the Python script. /!\ */
        fprintf(stderr, "# Free relations: %lu\n", F->nfree);
    }
    fprintf(stderr, "Renumbering struct: nprimes=%lu\n", (unsigned long) R_max_index);


    /* produce an error when index_t is too small to represent all ideals */
    if ((SIZEOF_INDEX < 8) && renumber_table.get_size() >> (8 * SIZEOF_INDEX)) {
        fprintf(stderr, "Error, please increase SIZEOF_INDEX\n");
        fprintf(stderr, "(see local.sh.example)\n");
        exit(1);
    }

    return 0;
}
