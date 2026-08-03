#include "cado.h" // IWYU pragma: keep

#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#ifdef HAVE_HWLOC
#include <hwloc.h>
#endif

#include <gmp.h>

#include "cado_poly.hpp"
#include "las-cofactor.hpp"
#include "las-info.hpp"
#include "las-side-config.hpp"
#include "las-siever-config.hpp"
#include "las-sieve-shared-data.hpp"
#include "las-special-q-task-collection.hpp"
#include "las-todo-list.hpp"
#include "macros.h"
#include "params.hpp"
#include "relation_cache.hpp"

/* las_info stuff */

/* Note that both sieve_shared_data and las_info have a poly field, which
 * is quite awkward. This should be cleaned up.
 */
void las_info::configure(cxx_param_list & pl)
{
    cxx_cado_poly::declare_usage(pl);
    cxx_cado_poly::configure_aliases(pl);
    cxx_cado_poly::configure_switches(pl);

    siever_config_pool::declare_usage(pl);

    sieve_shared_data::declare_usage(pl);

    las_dlog_base::declare_usage(pl);

    cofactorization_statistics::declare_usage(pl);

    batch_side_config::declare_usage(pl);
    pl.declare_usage("batch", "use batch cofactorization");
    pl.declare_usage("batch-print-survivors", "just print survivors to files with the given basename for an external cofactorization");
    pl.declare_usage("batch-print-survivors-filesize", "write that many survivors per file");
    pl.declare_usage("batch-print-survivors-number-of-printers", "use this number of I/O threads to write survivor files. defaults to 1, and should not be changed except in very unusual cases");
    pl.configure_switch("-batch");

    relation_cache::configure(pl);

    pl.declare_usage("galois", "depending on the specified galois automorphism, sieve only part of the q's");

    pl.declare_usage("dup", "suppress duplicate relations");
    pl.declare_usage("dup-qmin", "lower limit of global q-range for 2-sided duplicate removal");
    pl.declare_usage("dup-qmax", "upper limit of global q-range for 2-sided duplicate removal");
    pl.configure_switch("-dup");


    pl.declare_usage("smallset-purge", "use experimental 'smallset' code in purge_buckets");
    pl.configure_switch("-smallset-purge");


    pl.declare_usage("dumpfile", "Dump entire sieve region to file for debugging.");


    las_todo_list::declare_usage(pl);
    las_todo_list::configure_switches(pl);
}


void las_info::prepare_sieve_shared_data(cxx_param_list & pl)
{
#ifdef HAVE_HWLOC
    shared_structure_cache.clear();
    set_loose_binding();
    for(int k = 0 ; k < number_of_subjobs_total() ; k+= number_of_subjobs_per_memory_binding_zone()) {
        set_subjob_mem_binding(k);
        const auto nn = current_memory_binding();
        shared_structure_cache.emplace(nn, sieve_shared_data(cpoly, pl));
        /* for this one, the default ctor is good enough */
        las_memory_accessor_cache[nn];
    }
    set_loose_binding();
#else
    shared_structure_private = sieve_shared_data(cpoly, pl);
    /* the default-constructed memory accessor is fine */
#endif
}

void las_info::load_factor_base(cxx_param_list & pl)
{
#ifdef HAVE_HWLOC
    set_loose_binding();
    for(int k = 0 ; k < number_of_subjobs_total() ; k+= number_of_subjobs_per_memory_binding_zone()) {
        set_subjob_mem_binding(k);
        /* right, so at this point we would probably need to 
         * compute the factor base just once, and copy it in ram, isn't
         * it ?
         */
        local_cache().load_factor_base(pl, number_of_threads_loose());
    }
    set_loose_binding();
#else
    shared_structure_private.load_factor_base(pl, number_of_threads_loose());
#endif
}

template<sieve_method Algo>
las_info::las_info(cxx_param_list & pl, Algo)
    : galois(pl.lookup_old("galois"))
    , suppress_duplicates(pl.parse<bool>("-dup"))
    , use_smallset_purge(pl.parse<bool>("-smallset-purge"))
    , cpoly(pl)
    , config_pool(pl, cpoly.nsides(), Algo{})
#ifndef HAVE_HWLOC
    , shared_structure_private(cpoly, pl)
#endif
    , dlog_base(cpoly, pl)
    , tree(special_q_task_collection_base::create<Algo>(cpoly, pl))
    , rel_cache(pl)
    , cofac_stats(pl)
      /*{{{*/
{
    int const nsides = cpoly.nsides();

    /* We strive to initialize things in the exact order they're written
     * in the struct */
    // ----- general operational flags {{{


    if (const char * tmp = pl.lookup_old("bkmult")) {
        bk_multiplier = bkmult_specifier(tmp);
    }


    // }}}


    // ----- stuff roughly related to the descent {{{
    descent_helper = nullptr;
    // }}}

    /* {{{ duplicate suppression */
    dupqmin.assign(nsides, ULONG_MAX);
    dupqmax.assign(nsides, ULONG_MAX);
    if (suppress_duplicates) {
        std::vector<unsigned long> v;
        if (pl.parse_per_side("dup-qmin", v, nsides, ULONG_MAX)) {
            dupqmin = v;
        } else {
            fprintf(stderr, "Error: -dup-qmin is mandatory with -dup\n");
            exit(EXIT_FAILURE);
        }
        if (pl.parse_per_side("dup-qmax", v, nsides, ULONG_MAX)) {
            dupqmax = v;
        }
        /* The command-line value 0 also means ULONG_MAX */
        for (auto & x : dupqmin) if (x == 0) x = ULONG_MAX;
    }

    /* }}} */

    // ----- batch mode {{{
    batch = pl.parse<bool>("-batch");

    if (batch) {
        batch_side_config::parse(pl, bsides, nsides);
        ASSERT_ALWAYS(config_pool.default_config_ptr);
        siever_config const & sc0(*config_pool.default_config_ptr);

        /* Set some defaults. I agree that this logic is a little bit
         * quirky.
         */
        for(int side = 0 ; side < nsides ; side++) {
            if (bsides[side].batchlpb == UINT_MAX)
                bsides[side].batchlpb = sc0.sides[side].lpb;
            if (bsides[side].batchmfb == UINT_MAX)
                bsides[side].batchmfb = sc0.sides[side].lpb;
        }

        for(int side = 0 ; side < nsides ; side++) {
            auto const & bS = bsides[side];
            // the product of primes up to B takes \log2(B)-\log\log 2 /
            // \log 2 bits. The added constant is 0.5287.
            if (bS.batchlpb + 0.5287 >= 31 + log2(GMP_LIMB_BITS)) {
                fprintf(stderr, "Gnu MP cannot deal with primes product that large (max 37 bits, asked for batchlpb%d=%d)\n", side, bS.batchlpb);
                abort();
            } else if (bS.batchlpb + 0.5287 >= 34) {
                fprintf(stderr, "Gnu MP's mpz_inp_raw and mpz_out_raw functions are limited to integers of at most 34 bits (asked for batchlpb%d=%d)\n",side,bS.batchlpb);
                abort();
            }
        }
    }


    batch_print_survivors.filename = pl.lookup_old("batch-print-survivors");
    if (batch_print_survivors.filename) {
        batch_print_survivors.counter = 0;
        batch_print_survivors.number_of_printers = 1;
        batch_print_survivors.filesize = 1000000;
        pl.parse("batch-print-survivors-filesize", batch_print_survivors.filesize);
        pl.parse("batch-print-survivors-number-of-printers", batch_print_survivors.number_of_printers);
    }
    // }}} 

    dump_filename = pl.lookup_old("dumpfile");
}/*}}}*/

template las_info::las_info(cxx_param_list & pl, NFS);
template las_info::las_info(cxx_param_list & pl, SIQS);
