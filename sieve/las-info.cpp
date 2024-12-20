#include "cado.h" // IWYU pragma: keep

#include <climits>           // for ULONG_MAX
#include <cmath>             // for log2
#include <cstdio>             // for fprintf, NULL, stderr
#include <cstdlib>           // for abort, exit, EXIT_FAILURE
#include <gmp.h>              // for GMP_LIMB_BITS, gmp_randclear, gmp_randi...
#include "ecm/batch.hpp"          // for cofac_list
#include "las-info.hpp"
#include "las-todo-list.hpp"  // for las_todo_list
#include "macros.h"           // for ASSERT_ALWAYS
#include "params.h"

/* las_info stuff */

void las_info::configure_aliases(cxx_param_list & pl)
{
    cxx_cado_poly::configure_aliases(pl);
}

void las_info::configure_switches(cxx_param_list & pl)
{
    cxx_cado_poly::configure_switches(pl);
    las_todo_list::configure_switches(pl);
    param_list_configure_switch(pl, "-allow-compsq", NULL);
    param_list_configure_switch(pl, "-dup", NULL);
    param_list_configure_switch(pl, "-batch", NULL);
}

void las_info::declare_usage(cxx_param_list & pl)
{
    cxx_cado_poly::declare_usage(pl);
    siever_config_pool::declare_usage(pl);
    sieve_shared_data::declare_usage(pl);
    las_dlog_base::declare_usage(pl);
    cofactorization_statistics::declare_usage(pl);
    batch_side_config::declare_usage(pl);


    param_list_decl_usage(pl, "seed", "Use this seed for random state seeding (currently used only by --random-sample)");

    param_list_decl_usage(pl, "galois", "depending on the specified galois automorphism, sieve only part of the q's");

    /* Note: also declared by las_todo_list ! */
    param_list_decl_usage(pl, "allow-compsq", "allows composite special-q");
    param_list_decl_usage(pl, "qfac-min", "factors of q must be at least that");
    param_list_decl_usage(pl, "qfac-max", "factors of q must be at most that");




    param_list_decl_usage(pl, "dup", "suppress duplicate relations");
    param_list_decl_usage(pl, "dup-qmin", "lower limit of global q-range for 2-sided duplicate removal");
    param_list_decl_usage(pl, "dup-qmax", "upper limit of global q-range for 2-sided duplicate removal");
    param_list_decl_usage(pl, "adjust-strategy", "strategy used to adapt the sieving range to the q-lattice basis (0 = logI constant, J so that boundary is capped; 1 = logI constant, (a,b) plane norm capped; 2 = logI dynamic, skewed basis; 3 = combine 2 and then 0) ; default=0");


    param_list_decl_usage(pl, "batch", "use batch cofactorization");
    param_list_decl_usage(pl, "batch-print-survivors", "just print survivors to files with the given basename for an external cofactorization");
    param_list_decl_usage(pl, "batch-print-survivors-filesize", "write that many survivors per file");
    param_list_decl_usage(pl, "batch-print-survivors-number-of-printers", "use this number of I/O threads to write survivor files. defaults to 1, and should not be changed except in very unusual cases");

    param_list_decl_usage(pl, "relation_cache", "Directory with cache of collected relation for sampling within a known data set. Useful only with --random-sample\n");
    param_list_decl_usage(pl, "dumpfile", "Dump entire sieve region to file for debugging.");
}


void las_info::prepare_sieve_shared_data(cxx_param_list & pl)
{
#ifdef HAVE_HWLOC
    shared_structure_cache.clear();
    set_loose_binding();
    for(int k = 0 ; k < number_of_subjobs_total() ; k+= number_of_subjobs_per_memory_binding_zone()) {
        set_subjob_mem_binding(k);
        /* gcc pre 4.7 does not have map::emplace ; while it's
         * admittedly not a c++11-complete compiler anyway, we still have
         * hope to use such an ancient thing on a few machines (old *bsd,
         * for example).
         *
         * See:
         * https://gcc.gnu.org/bugzilla/show_bug.cgi?id=44436#c30
         */
#if 0
        shared_structure_cache.emplace(
                current_memory_binding(),
                sieve_shared_data(cpoly, pl));
#else
        std::pair<cxx_hwloc_nodeset, sieve_shared_data> v(
                current_memory_binding(),
                sieve_shared_data(cpoly, pl));
        shared_structure_cache.insert(std::move(v));
        /* for this one, the default ctor is good enough */
        las_memory_accessor_cache[current_memory_binding()];
#endif
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
        shared_structure_cache.at(current_memory_binding()).load_factor_base(pl, number_of_threads_loose());
    }
    set_loose_binding();
#else
    shared_structure_private.load_factor_base(pl, number_of_threads_loose());
#endif
}

las_info::las_info(cxx_param_list & pl)
    : cpoly(pl)
    , config_pool(pl, cpoly->nb_polys)
#ifdef HAVE_HWLOC
    , shared_structure_cache()
#else
    , shared_structure_private(cpoly, pl)
#endif
    , dlog_base(cpoly, pl)
    , cofac_stats(pl)
      /*{{{*/
{
    int const nsides = cpoly->nb_polys;

    /* We strive to initialize things in the exact order they're written
     * in the struct */
    // ----- general operational flags {{{
    gmp_randinit_default(rstate);
    unsigned long seed = 0;
    if (param_list_parse_ulong(pl, "seed", &seed))
        gmp_randseed_ui(rstate, seed);

    galois = param_list_lookup_string(pl, "galois");
    suppress_duplicates = param_list_parse_switch(pl, "-dup");

    if (const char * tmp = param_list_lookup_string(pl, "bkmult")) {
        bk_multiplier = bkmult_specifier(tmp);
    }

    param_list_parse_int(pl, "adjust-strategy", &adjust_strategy);

    // }}}


    /* composite special-q ? Note: this block is present both in
     * las-todo-list.cpp and las-info.cpp */
    if ((allow_composite_q = param_list_parse_switch(pl, "-allow-compsq"))) {
        /* defaults are set in the class description */
        param_list_parse_uint64(pl, "qfac-min", &qfac_min);
        param_list_parse_uint64(pl, "qfac-max", &qfac_max);
    }

    param_list_parse(pl, "relation_cache", relation_cache);

    // ----- stuff roughly related to the descent {{{
    descent_helper = NULL;
    // }}}

    /* {{{ duplicate suppression */
    dupqmin.assign(nsides, ULONG_MAX);
    dupqmax.assign(nsides, ULONG_MAX);
    if (suppress_duplicates) {
        if (!param_list_parse_per_side<unsigned long>(pl, "dup-qmin", dupqmin.data(), nsides, ARGS_PER_SIDE_DEFAULT_AS_IS)) {
            fprintf(stderr, "Error: -dup-qmin is mandatory with -dup\n");
            exit(EXIT_FAILURE);
        }
        param_list_parse_per_side<unsigned long>(pl, "dup-qmax", dupqmax.data(), nsides, ARGS_PER_SIDE_DEFAULT_AS_IS);
        /* The command-line value 0 also means ULONG_MAX */
        for (auto & x : dupqmin) if (x == 0) x = ULONG_MAX;
    }

    /* }}} */

    // ----- batch mode {{{
    batch = param_list_parse_switch(pl, "-batch");

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


    batch_print_survivors.filename = param_list_lookup_string(pl, "batch-print-survivors");
    if (batch_print_survivors.filename) {
        batch_print_survivors.counter = 0;
        batch_print_survivors.number_of_printers = 1;
        batch_print_survivors.filesize = 1000000;
        param_list_parse_uint64(pl, "batch-print-survivors-filesize", &batch_print_survivors.filesize);
        param_list_parse_int(pl, "batch-print-survivors-number-of-printers", &batch_print_survivors.number_of_printers);
    }
    // }}} 

    dump_filename = param_list_lookup_string(pl, "dumpfile");
}/*}}}*/


las_info::~las_info()/*{{{*/
{

    // ----- general operational flags {{{
    gmp_randclear(rstate);

    // }}}

}/*}}}*/
