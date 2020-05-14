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

/* Model of this code : One producer -> Many consumers/producers -> one
   consumer.

   1. One thread produces NB_PRIMERS_PER_BUFFER primes per buffer, in
      NB_BUFFERS_PER_THREAD buffers (named bufs) for each working thread. This
      thread uses the function pthread_primes_producer.

      When all the primes are produced, this thread produces one buffer with
      only one "false" prime = FREEREL_END_MARKER for each consumer. It's the
      end marker.

   2. 'nthreads' threads load a primes buffer by a classical one producer/many
      consumer models (2 pseudo semaphores BY pair of producer/consumer).
      Each thread produces 2 buffers: one buffer containing the local
      renumbering table on the form of a ASCII array and of one buffer
      containing the information to retrieve the free relations. Those threads
      use the function pthread_primes_consumer.

      When a thread loads a primes buffer which begins and contains only
      FREEREL_END_MARKER, it produces a free relations buffer which contains
      only FREEREL_END_MARKER and exits.

   3. With a many producers/one consumer model, the principal programs (the
      function generate_renumber_and_freerels) loads the data produced by the
      'nthreads' pthread_primes_consumer, writes sequentially the renumbering
      table, and print the free relations.

      When a roots buffer contains exactly one FREEREL_END_MARKER, the job is
      done.
*/

#include "cado.h"
#include <algorithm>
#include <vector>
#include <list>
#include <fstream>
#include <gmp.h>

#include "utils.h"

struct freerel_data_t {
    FILE * sink = NULL;
    const char * filename = NULL;
    unsigned long pmin = 2;
    unsigned long pmax = 0;
    unsigned long nfree = 0;
    bool print_this_p(unsigned long p) const {
        return sink != NULL && pmin <= p && p <= pmax;
    }
};

struct renumber_data_t {
    std::ostream & os;
    renumber_t & R;
    renumber_data_t(renumber_t & R, std::ostream & os)
        : os(os)
        , R(R)
    {}
};

struct synchronous_data {
    stats_data_t stats;
    freerel_data_t F;
    renumber_data_t R;
    unsigned long R_max; // sigh... *must* be ulong for stats().
    synchronous_data(freerel_data_t F, renumber_t & R, std::ostream& os_r)
        : F(F)
        , R(R, os_r)
        , R_max(R.get_max_index())
    {
        /* will print report at 2^10, 2^11, ... 2^23 computed primes
         * and every 2^23 primes after that */
        stats_init(stats, stdout, &R_max, 23,
                "Processed", "primes", "", "p");
    }
    void progress() {
        if (stats_test_progress(stats))
            stats_print_progress(stats, R_max, 0, 0, 0);
    }
    ~synchronous_data() {
        stats_print_progress(stats, R_max, 0, 0, 1);
    }
};

renumber_t::cooked write_cooked_data_for_prime(renumber_t const & renumber_table, unsigned long p)/*{{{*/
{
    std::vector<std::vector<unsigned long>> all_roots;
    for (unsigned int side = 0; side < renumber_table.get_nb_polys(); side++) {
        std::vector<unsigned long> roots;
        mpz_poly_srcptr f = renumber_table.get_poly(side);

        if (UNLIKELY(p >> renumber_table.get_lpb(side)))
            roots.clear();
        else if (f->deg == 1)
            roots.assign(1, 0);
        else {
            /* Note that roots are sorted */
            roots = mpz_poly_roots(f, p);
        }

        /* Check for a projective root ; append it (so that the list of roots
         * is still sorted)
         */

        if ((int) roots.size() != renumber_table.get_poly_deg(side)
                && mpz_divisible_ui_p(f->coeff[f->deg], p))
            roots.push_back(p);

        /* take off bad ideals from the list, if any. */
        if (p <= renumber_table.get_max_bad_p()) { /* can it be a bad ideal ? */
            for (size_t i = 0; i < roots.size() ; i++) {
                unsigned long r = roots[i];
                if (!renumber_table.is_bad(p, r, side))
                    continue;
                /* bad ideal -> remove this root from the list */
                roots.erase(roots.begin() + i);
                i--;
            }
        }
        all_roots.emplace_back(roots);
    }

    /* Data is written in the temp buffer in a way that is not quite
     * similar to the renumber table, but still close enough.
     */
    return renumber_table.cook(p, all_roots);
}/*}}}*/

void use_cooked_data_for_prime(synchronous_data & S, unsigned long p, renumber_t::cooked const & C)
{
    index_t old_nprimes = S.R_max;

    index_t current_nprimes = S.R_max;
    std::vector<std::pair<int, index_t>> full_sides;

    S.R_max = S.R.R.use_cooked_nostore(S.R_max, C);

    for(unsigned int side = 0 ; side < S.R.R.get_nb_polys() ; ++side) {
        mpz_poly_srcptr f = S.R.R.get_poly(side);
        /* Check if p corresponds to free relations.
           For 2 polynomials, this happens when the number of roots modulo p
           is maximal for both polynomials.
           For 3 or more polynomials, let t be the number of polynomials
           such that the number of roots modulo p is maximal. Then the
           number of free relations is t-1 (see Section 4.6 from the PhD
           thesis from Marije Elkenbracht-Huizing, "Factoring integers
           with the Number Field Sieve") */
        if ((int) C.nroots[side] == f->deg && S.F.print_this_p(p))
            full_sides.emplace_back(side, current_nprimes);

        current_nprimes += C.nroots[side];
    }

    if (full_sides.size() > 1) {
        for(size_t i = 1 ; i < full_sides.size() ; i++) {
            /* print a new free relation */
            fprintf(S.F.sink, "%" PRpr ",0:%" PRid, (p_r_values_t) p, old_nprimes);
            int side0 = full_sides[i-1].first;
            index_t i0 = full_sides[i-1].second;
            unsigned int n0 = C.nroots[side0];
            int side1 = full_sides[i].first;
            index_t i1 = full_sides[i].second;
            unsigned int n1 = C.nroots[side1];
            bool first = true;
            for(unsigned int k = 0 ; k < n0 ; k++, first=false)
                fprintf(S.F.sink, "%s%" PRid, first ? "" : ",", i0 + k);
            for(unsigned int k = 0 ; k < n1 ; k++, first=false)
                fprintf(S.F.sink, "%s%" PRid, first ? "" : ",", i1 + k);
            fputc('\n', S.F.sink);

            S.F.nfree++;
        }
    }
}

struct prime_chunk {
    bool done = false;
    std::vector<unsigned long> primes;
    std::vector<renumber_t::cooked> C;
    void process(renumber_t const & renumber_table) {
        ASSERT_ALWAYS(!done);
        /* change x (list of input primes) into the list of integers that go
         * to the renumber table, and then set "done" to true.
         * This is done asynchronously.
         */
        for(auto p : primes) {
            C.emplace_back(write_cooked_data_for_prime(renumber_table, p));
        }
#pragma omp atomic write
        done = true;
    }
    void postprocess(synchronous_data & S) {
        ASSERT_ALWAYS(!done);
        /* put all entries from x into the renumber table, and also print
         * to freerel_file any free relation encountered. This is done
         * synchronously.
         *
         * (if freerel_file is NULL, store only into the renumber table)
         */
        for(size_t i = 0; i < primes.size() ; i++)
            use_cooked_data_for_prime(S, primes[i], C[i]);
        /* free memory ! */
        primes.clear();
        C.clear();
        S.progress();
    }
    prime_chunk(std::vector<unsigned long> && primes) : primes(primes) {}
};

/* Generate all the free relations and the renumbering table.
 * Return the number of free relations */
static void
generate_renumber_and_freerels(freerel_data_t & freerel_data,
                               renumber_t& renumber_table,
                               std::ostream& os_r)
{
    /* open freerel_file */
    if (freerel_data.filename) {
        freerel_data.sink = fopen_maybe_compressed(freerel_data.filename, "w");
        ASSERT_ALWAYS(freerel_data.sink != NULL);
    }

    synchronous_data sync(freerel_data, renumber_table, os_r);

    unsigned long lpbmax = 1UL << renumber_table.get_max_lpb();

    printf("Generating renumber table for 2 <= p <= %lu\n", lpbmax);
    printf("Considering freerels for %lu <= p <= %lu\n", 
            freerel_data.pmin,
            freerel_data.pmax);
    fflush(stdout);


    prime_info pi;
    prime_info_init(pi);
    unsigned long p = 2;
    constexpr const unsigned int granularity = 1024;
    std::list<prime_chunk> chunks;
    std::list<prime_chunk>::iterator next = chunks.begin();

#pragma omp parallel
    {
#pragma omp single
        {
            for (; p <= lpbmax;) {
                std::vector<unsigned long> pp;
                pp.reserve(granularity);
                for (; p <= lpbmax && pp.size() < granularity;) {
                    pp.push_back(p);
                    p = getprime_mt(pi); /* get next prime */
                }
                chunks.emplace_back(std::move(pp));
#pragma omp task
                chunks.back().process(renumber_table);

                for (; next != chunks.end(); ) {
                    bool ok;
#pragma omp atomic read
                    ok = next->done;
                    if (!ok)
                        break;
                    next->postprocess(sync);
                    auto garbage = next++;
                    chunks.erase(garbage);
                }
            }
        }
    }
    prime_info_clear(pi);

    if (freerel_data.filename)
        fclose_maybe_compressed(freerel_data.sink, freerel_data.filename);
}

static void
declare_usage(param_list pl)
{
    param_list_decl_usage(pl, "poly", "input polynomial file");
    param_list_decl_usage(pl, "renumber", "output file for renumbering table");
    param_list_decl_usage(pl, "out", "output file for free relations");
    param_list_decl_usage(pl, "lpb0", "large primes bound on side 0");
    param_list_decl_usage(pl, "lpb1", "large primes bound on side 1");
    param_list_decl_usage(pl,
                          "lpbs",
                          "large primes bounds (comma-separated list) "
                          "(for MNFS)");
    param_list_decl_usage(pl, "pmin", "do not create freerel below this bound");
    param_list_decl_usage(
      pl, "pmax", "do not create freerel beyond this bound");
    // param_list_decl_usage(pl, "badideals", "file describing bad ideals (for
    // DL)");
    param_list_decl_usage(pl,
                          "lcideals",
                          "Add ideals for the leading "
                          "coeffs of the polynomials (for DL)");
    // param_list_decl_usage(pl, "t", "number of threads");
}

static void
usage(param_list pl, char* argv0)
{
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
}

int
main(int argc, char* argv[])
{
    char* argv0 = argv[0];
    freerel_data_t freerel_data;
    int lcideals = 0;
    int lpb_arg[NB_POLYS_MAX] = { 0 };
    std::vector<unsigned int> lpb;
    cxx_param_list pl;
    cxx_cado_poly cpoly;

    /* {{{ parse cmdline */
    declare_usage(pl);
    param_list_configure_switch(pl, "-lcideals", &lcideals);
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
    const char* polyfilename = param_list_lookup_string(pl, "poly");
    freerel_data.filename = param_list_lookup_string(pl, "out");
    // const char * badidealsfilename = param_list_lookup_string(pl,
    // "badideals");
    const char* renumberfilename = param_list_lookup_string(pl, "renumber");

    int has_lpb01 = 0;
    has_lpb01 += param_list_parse_int(pl, "lpb0", &(lpb_arg[0]));
    has_lpb01 += param_list_parse_int(pl, "lpb1", &(lpb_arg[1]));
    int has_nlpbs =
      param_list_parse_int_list(pl, "lpbs", lpb_arg, NB_POLYS_MAX, ",");
    param_list_parse_ulong(pl, "pmin", &freerel_data.pmin);
    param_list_parse_ulong(pl, "pmax", &freerel_data.pmax);
    // param_list_parse_uint (pl, "t", &nthreads);

    if (param_list_warn_unused(pl))
        usage(pl, argv0);

    if (polyfilename == NULL) {
        fprintf(stderr, "Error, missing -poly command line argument\n");
        usage(pl, argv0);
    }
    if (renumberfilename == NULL) {
        fprintf(stderr, "Error, missing -renumber command line argument\n");
        usage(pl, argv0);
    }
    if (freerel_data.filename == NULL) {
        fprintf(stderr, "Error, missing -out command line argument\n");
        usage(pl, argv0);
    }
    if (!cado_poly_read(cpoly, polyfilename)) {
        fprintf(stderr, "Error reading polynomial file\n");
        exit(EXIT_FAILURE);
    }
    if (has_nlpbs && has_lpb01) {
        fprintf(stderr, "Error, lpb[01] and lpbs are incompatible\n");
        exit(EXIT_FAILURE);
    }

    if (has_nlpbs == 0) /* lpbs were given as -lpb0 and -lpb1 */
    {
        has_nlpbs = 2;
        if (cpoly->nb_polys > 2) /* With more than 2 polys, must use -lpbs. */
        {
            fprintf(stderr, "Error, missing -lpbs command line argument\n");
            usage(pl, argv0);
        }
    }

    if (has_nlpbs != cpoly->nb_polys) {
        fprintf(stderr,
                "Error, the number of values given in -lpbs does not "
                "correspond to the number of polynomials\n");
        usage(pl, argv0);
    }
    lpb.assign(&lpb_arg[0], &lpb_arg[has_nlpbs]);
    for (auto l : lpb) {
        if (l <= 0) {
            fprintf(stderr,
                    "Error, -lpbs command line argument cannot contain "
                    "non-positive values\n");
            usage(pl, argv0);
        }
    }
    /* }}} */

    renumber_t renumber_table(cpoly);
    renumber_table.set_lpb(lpb);
    if (lcideals)
        renumber_table.use_additional_columns_for_dl();
    renumber_table.compute_bad_ideals();

    std::ofstream of_r(renumberfilename);

    renumber_table.write_header(of_r);
    renumber_table.write_bad_ideals(of_r);

    /* if pmax is not equal to 0 (i.e., was not given on the command line),
     * set pmax to the *minimum* of the large prime bounds, since larger primes
     * will never occur on both sides.
     * We generate the renumbering table from the first prime (2) and
     * up to the *maximum* of the large prime bounds.
     */
    if (!freerel_data.pmax)
        freerel_data.pmax = 1UL << renumber_table.get_min_lpb();

    generate_renumber_and_freerels(freerel_data, renumber_table, of_r);

    of_r.close();

    /* /!\ Needed by the Python script. /!\ */
    fprintf(stderr, "# Free relations: %" PRIu64 "\n", freerel_data.nfree);
    fprintf(stderr,
            "Renumbering struct: nprimes=%" PRIu64 "\n",
            renumber_table.get_size());

    /* produce an error when index_t is too small to represent all ideals */
    if ((SIZEOF_INDEX < 8) && renumber_table.get_size() >> (8 * SIZEOF_INDEX)) {
        fprintf(stderr, "Error, please increase SIZEOF_INDEX\n");
        fprintf(stderr, "(see local.sh.example)\n");
        exit(1);
    }

    return 0;
}
