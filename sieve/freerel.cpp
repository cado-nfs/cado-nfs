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

#include <cstdio>
#include <cstdlib>
#include <climits>

#include <algorithm>
#include <vector>
#include <memory>
#include <utility>
#include <ios>

#include "fmt/base.h"
#include "fmt/format.h"

#include "cado_poly.h"
#include "gzip.h"
#include "macros.h"
#include "mpz_poly.h"
#include "omp_proxy.h"
#include "params.h"
#include "renumber.hpp"
#include "typedefs.h"
#include "verbose.h"

static char const * argv0;

static void
usage(cxx_param_list & pl, char const * argv0)
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
        std::ranges::sort(lpb_copy);
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
            sink << fmt::format("{:x},0:", p);
            int const side0 = full_sides[i-1].first;
            index_t const i0 = full_sides[i-1].second;
            unsigned int const n0 = C.nroots[side0];
            int const side1 = full_sides[i].first;
            index_t const i1 = full_sides[i].second;
            unsigned int const n1 = C.nroots[side1];
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
main(int argc, char const * argv[])
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
        if ((f = fopen(argv[0], "r")) != nullptr) {
            param_list_read_stream(pl, f, 0);
            fclose(f);
            argv++, argc--;
            continue;
        }
        fmt::print(stderr, "Unhandled parameter {}\n", argv[0]);
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
        fmt::print(stderr, "Warning: the -t argument to freerel is kept for compatibility, but you should rather take it out and let openmp deal with it\n");
        omp_set_num_threads(nthreads);
    }

    freerel_data_t::lookup_parameters(pl);

    param_list_lookup_string(pl, "lpb0");
    param_list_lookup_string(pl, "lpb1");
    param_list_lookup_string(pl, "lpbs");


    if (param_list_warn_unused(pl))
        usage(pl, argv0);

    if (polyfilename == nullptr) {
        fmt::print(stderr, "Error, missing -poly command line argument\n");
        usage(pl, argv0);
    }
    if (!param_list_lookup_string(pl, "renumber")) {
        fmt::print(stderr, "Error, missing -renumber command line argument\n");
        usage(pl, argv0);
    }
    if (!cado_poly_read(cpoly, polyfilename)) {
        fmt::print(stderr, "Error reading polynomial file\n");
        exit(EXIT_FAILURE);
    }

    std::vector<unsigned int> lpb(cpoly->nb_polys, 0);

    if (!param_list_parse_uint_args_per_side(pl, "lpb", lpb.data(), cpoly->nb_polys, ARGS_PER_SIDE_DEFAULT_COPY_PREVIOUS)) {
        fmt::print(stderr,
                "Error, could not obtain values for the lpb bounds (or not for all polynomials)\n");
        usage(pl, argv0);
    }

    /* }}} */

    std::unique_ptr<freerel_data_t> F;

    if (param_list_lookup_string(pl, "out")) {
        F = std::make_unique<freerel_data_t>(pl, cpoly, lpb);

        fmt::print("Considering freerels for {} <= p <= {}\n", 
                F->pmin,
                F->pmax);
        fflush(stdout);
    }

    renumber_t renumber_table(cpoly);
    renumber_table.set_lpb(lpb);

    auto const max_lpb = renumber_table.get_max_lpb();
    if (max_lpb >= sizeof(unsigned long) * CHAR_BIT) {
      fmt::print (stderr, "Error, cannot handle lpb >= {} (max(lpbs) is {})\n",
                       sizeof(unsigned long) * CHAR_BIT, max_lpb);
      abort();
    }
    unsigned long const lpbmax = 1UL << max_lpb;
    fmt::print("Generating renumber table for 2 <= p <= {}\n", lpbmax);

    /* This reads the options:
     *
     * renumber
     *
     */
    index_t const R_max_index = renumber_table.build(pl, for_dl, F.get());

    if (F) {
        /* /!\ Needed by the Python script. /!\ */
        fmt::print(stderr, "# Free relations: {}\n", F->nfree);
    }
    fmt::print(stderr, "Renumbering struct: nprimes={}\n", (unsigned long) R_max_index);


    /* produce an error when index_t is too small to represent all ideals */
    if ((SIZEOF_INDEX < 8) && renumber_table.size() >> (8 * SIZEOF_INDEX)) {
        fmt::print(stderr, "Error, please increase SIZEOF_INDEX\n");
        fmt::print(stderr, "(see local.sh.example)\n");
        exit(1);
    }

    return 0;
}
